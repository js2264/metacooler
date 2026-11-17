"""Add, replace or drop the metadata columns of an existing map.

``/genomes`` and ``/contigs`` carry, beside the columns metacooler derives,
any number of *metadata* columns: completeness, taxonomy, coverage, ...
``create`` and ``collect`` take them from the ``--genomes`` / ``--contigs``
tables at creation time; this module attaches them afterwards, when the
numbers arrive (CheckM, GTDB-Tk and depth tables usually come after the map).

What it guarantees:

- **Same meaning as at creation.**  The tables are read by the same parsers
  (``parse_genomes_table``, ``collect``'s per-contig metadata reader), columns
  colliding with metacooler's own names get the same ``user_`` prefix, and a
  column is typed by the same rule (:func:`_coerce_numeric_columns`).  So a
  column added here holds what the same table would have given ``create``,
  and a ``--filter`` written against it means the same thing either way.
- **Only metadata moves.**  Rows are matched by name onto the stored order and
  never added, removed or reordered; derived columns (``name``, ``length``,
  ``genome_id``, ``n_contigs``, ``n_bins``, ``qc_*``) are never written; bins,
  pixels and indexes are not read.  The cost is O(n_contigs).
- **Nothing half-written.**  Every table is checked -- names, column clashes,
  ASCII in ``/contigs`` -- before anything is written; a table is rebuilt in a
  scratch group and swapped into place; the map is validated afterwards.

A bare ``.mmetacool`` path annotates every entry of the collection (each entry
carries its own ``/genomes`` and ``/contigs``: ``zoomify`` levels, stacked
replicates); an entry URI annotates that entry alone.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import (
    UNMATCHED_CHOICES,
    Unmatched,
    check_output_path,
    report_unmatched_rows,
    unmatched_rows_error,
    write_column,
)
from metacooler.io.colspec import CONTIG_METADATA_FORMAT_DEFAULT, GENOMES_FORMAT_DEFAULT
from metacooler.qc import is_qc_column
from metacooler.schema import (
    MetaCoolerFormatError,
    format_ascii_violations,
    is_mmetacooler,
    list_metacoolers,
    non_ascii_values,
    validate,
)
from metacooler.uri import parse_uri

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

Missing = Literal["keep", "na"]
TABLES = ("genomes", "contigs")

#: Scratch group, at the map's root, where a table is rebuilt before being
#: swapped into place.  Removed on success; a leftover from an interrupted run
#: is cleared by the next one.
_SCRATCH = ".annotate-scratch"

# How many names or values to quote in a message.
_MAX_QUOTED = 5


@dataclasses.dataclass
class _TableEdit:
    """What annotating does to one table of one map.

    Attributes:
        table: ``"genomes"`` or ``"contigs"``.
        values: New content of every added or replaced column, in stored row
            order.
        added: Columns that did not exist, in table order.
        replaced: Existing columns whose values change.
        unchanged: Existing columns given the values they already hold.
        dropped: Columns removed.
        matched: How many rows of the user table matched a row of the map.
    """

    table: str
    values: dict[str, np.ndarray] = dataclasses.field(default_factory=dict)
    added: list[str] = dataclasses.field(default_factory=list)
    replaced: list[str] = dataclasses.field(default_factory=list)
    unchanged: list[str] = dataclasses.field(default_factory=list)
    dropped: list[str] = dataclasses.field(default_factory=list)
    matched: int = 0

    @property
    def changed(self) -> bool:
        """Whether the table is rewritten at all."""
        return bool(self.added or self.replaced or self.dropped)


@dataclasses.dataclass
class _MapEdit:
    """What annotating does to one map.

    Attributes:
        uri: The map, as ``file::group``.
        group: HDF5 group of the map inside its file.
        tables: The edit of each table that was touched.
        conflicts: ``table.column`` -> description, for existing columns the
            tables would change while ``replace`` is off.
        absent_drops: ``table.column`` names ``drop`` asked for that this map
            does not have.
    """

    uri: str
    group: str
    tables: dict[str, _TableEdit] = dataclasses.field(default_factory=dict)
    conflicts: dict[str, str] = dataclasses.field(default_factory=dict)
    absent_drops: list[str] = dataclasses.field(default_factory=list)

    @property
    def changed(self) -> bool:
        """Whether anything is written to this map."""
        return any(edit.changed for edit in self.tables.values())


def annotate(
    uri: str | Path,
    output: str | Path | None = None,
    *,
    genomes_table: str | Path | None = None,
    contigs_table: str | Path | None = None,
    genomes_table_format: str = GENOMES_FORMAT_DEFAULT,
    contigs_table_format: str = CONTIG_METADATA_FORMAT_DEFAULT,
    columns: Sequence[str] | None = None,
    drop: Sequence[str] | None = None,
    replace: bool = False,
    missing: Missing = "keep",
    unmatched: Unmatched = "warn",
    inplace: bool = False,
    metadata: dict | None = None,
    force: bool = False,
    keep_partial: bool = True,
) -> None:
    """Add, replace or drop metadata columns of ``/genomes`` and ``/contigs``.

    Args:
        uri: A ``.metacool`` path, an ``.mmetacool`` entry URI
            (``file.mmetacool::entries/<name>``), or a bare ``.mmetacool``
            path, which annotates every entry of the collection.
        output: Where to write the annotated copy.  A whole file in gives a
            whole file out (a collection stays a collection); an entry URI
            gives a standalone ``.metacool``.  Required unless *inplace*.
        genomes_table: Genomes metadata TSV.  Its genome-name column (see
            *genomes_table_format*) is matched against ``/genomes/name``; every
            other column is metadata.
        contigs_table: Per-contig metadata TSV, as ``collect`` reads it: its
            contig-name column (see *contigs_table_format*) names a contig as
            the map stores it (or, in a collected map, as in its input cool);
            an optional ``genome_name`` column is checked, never applied;
            every other column is metadata.
        genomes_table_format: Which column of *genomes_table* holds the genome
            names: a 1-based position or a header name, plus ``,noheader``.
        contigs_table_format: Which column of *contigs_table* holds the contig
            names, in the same grammar.
        columns: Take only these metadata columns, as ``genomes.<name>`` /
            ``contigs.<name>`` (named as in the tables).  Default: all of them.
        drop: Remove these columns, as ``genomes.<name>`` / ``contigs.<name>``.
        replace: Allow changing the values of a column the map already has.
            A column given the values it already holds is never a change.
        missing: For a replaced column, what the rows the table does not list
            get: ``"keep"`` (their current value: the table is a patch) or
            ``"na"`` (an empty value: the table is authoritative).  New
            columns always leave unlisted rows empty.
        unmatched: What to do with the rows of a table that name no genome
            (contig) of a map: ``"warn"`` (default) ignores them with a
            WARNING, ``"ignore"`` without one (DEBUG only), ``"error"``
            refuses the tables before anything is written.  A table that
            matches nothing at all is refused whatever this says.
        inplace: Edit *uri* itself instead of writing *output*.
        metadata: Keys merged into each edited map's ``metadata`` attribute; a
            ``commands`` list is appended to the provenance already there.
        force: Replace *output* if it exists (only once the new file
            validates).
        keep_partial: On a failure while the output is being written, leave
            ``<output>.partial`` behind for inspection, as the default does
            (logged at ERROR).  ``False`` removes it instead, whatever stopped
            the run (``KeyboardInterrupt`` included), and logs that at INFO.
            The output path itself is never written by a failed run.  Unused
            with *inplace*, which stages nothing.

    Raises:
        ValueError: On inconsistent arguments, a table that matches nothing (or,
            under ``unmatched="error"``, a row that does), a column the tables
            would change without *replace*, a derived column in *drop*, or a
            contig table contradicting the genome assignment.
        MetaCoolerFormatError: If a ``/contigs`` column would hold non-ASCII
            text (cooler could no longer open the map).
        FileNotFoundError: If the map or a table does not exist.
    """
    from metacooler.api import MetaCooler

    if inplace and output is not None:
        raise ValueError(
            "inplace edits the input, so no output path is expected; pass one "
            "or the other."
        )
    if not inplace and output is None:
        raise ValueError(
            "An output path is required. Pass output=... to write an annotated "
            "copy, or inplace=True to edit the input."
        )
    if missing not in ("keep", "na"):
        raise ValueError(f"missing must be 'keep' or 'na', got {missing!r}.")
    if unmatched not in UNMATCHED_CHOICES:
        raise ValueError(
            f"unmatched must be 'warn', 'ignore' or 'error', got {unmatched!r}."
        )
    wanted = _qualified_names(columns, "columns")
    dropping = _qualified_names(drop, "drop")
    if genomes_table is None and contigs_table is None and not dropping:
        raise ValueError(
            "Nothing to do: give a genomes table, a contigs table, or columns to drop."
        )
    _check_drops(dropping, wanted)
    for table, given in (("genomes", genomes_table), ("contigs", contigs_table)):
        if given is None and wanted.get(table):
            raise ValueError(
                f"columns names {table} column(s) "
                f"({', '.join(wanted[table])}) but no {table} table was given."
            )
    inputs = [uri, genomes_table, contigs_table]
    if output is not None:
        from metacooler.io.create import _refuse_staging_an_input

        check_output_path(output, force, inputs=inputs, defer_removal=True)
        _refuse_staging_an_input(output, inputs)

    # --- Read the tables once; they do not depend on the map ---
    genomes_extra = (
        _read_genomes(genomes_table, genomes_table_format, wanted)
        if genomes_table is not None
        else None
    )
    contig_table = (
        _read_contigs(contigs_table, contigs_table_format, wanted)
        if contigs_table is not None
        else None
    )

    # --- Plan every map before writing any ---
    file_path, group_path = parse_uri(str(uri))
    groups = _map_groups(file_path, group_path)
    edits = []
    # Rows refused under unmatched="error", raised once every map is planned:
    # a table that matches nothing at all gets its own, clearer refusal first.
    unmatched_rows: list[str] = []
    for group in groups:
        map_uri = f"{file_path}::{group}"
        if len(groups) > 1:
            logger.info("Annotating entry %s", group)
        edits.append(
            _plan_map(
                MetaCooler(map_uri),
                group,
                genomes_extra,
                contig_table,
                contigs_table,
                dropping,
                replace=replace,
                missing=missing,
                unmatched=unmatched,
                unmatched_rows=unmatched_rows,
            )
        )
    _check_plans(edits, genomes_table, contigs_table, dropping, unmatched_rows)

    history = {
        "genomes_table": None if genomes_table is None else str(genomes_table),
        "contigs_table": None if contigs_table is None else str(contigs_table),
        "replace": replace,
        "missing": missing,
    }
    if inplace:
        for edit in edits:
            _write_map(file_path, edit.group, edit, metadata, history)
            if edit.changed:
                validate(edit.uri)
        return

    from metacooler.io.create import _staged_output

    assert output is not None  # guaranteed by the argument checks above
    target = Path(output)
    with _staged_output(target, keep_partial=keep_partial, inputs=inputs) as staged:
        if group_path.strip("/") and len(groups) == 1:
            # One entry out of a collection: a standalone map, as reassign does.
            from metacooler.io.reassign import _materialise_copy

            _materialise_copy(MetaCooler(f"{file_path}::{group_path}"), staged)
            moved = [dataclasses.replace(edits[0], group="/")]
        else:
            shutil.copy2(file_path, staged)
            moved = edits
        for edit in moved:
            _write_map(staged, edit.group, edit, metadata, history, always=True)
            validate(f"{staged}::{edit.group}")
    logger.info("Annotated map written to %s", target)


def _qualified_names(names: Sequence[str] | None, what: str) -> dict[str, list[str]]:
    """Split ``table.column`` names by table, refusing an unknown table."""
    out: dict[str, list[str]] = {}
    for name in names or ():
        table, dot, column = str(name).partition(".")
        if not dot or table not in TABLES or not column:
            raise ValueError(
                f"{what} takes columns as 'genomes.<column>' or "
                f"'contigs.<column>', got {name!r}."
            )
        out.setdefault(table, []).append(column)
    return out


def _check_drops(dropping: dict[str, list[str]], wanted: dict[str, list[str]]) -> None:
    """Refuse to drop a derived column, or one also asked for in *columns*."""
    from metacooler.io.create import (
        _CONTIGS_RESERVED_COLUMNS,
        _GENOMES_RESERVED_COLUMNS,
    )

    protected = {
        "genomes": _GENOMES_RESERVED_COLUMNS,
        "contigs": _CONTIGS_RESERVED_COLUMNS,
    }
    for table, names in dropping.items():
        derived = [n for n in names if n in protected[table] or is_qc_column(table, n)]
        if derived:
            raise ValueError(
                f"Cannot drop {', '.join(f'{table}.{n}' for n in derived)}: "
                "metacooler derives these columns. The genome assignment "
                "changes with 'metacooler reassign', and the qc_* columns with "
                "'metacooler qc --store'."
            )
        both = sorted(set(names) & set(wanted.get(table, ())))
        if both:
            raise ValueError(
                f"{', '.join(f'{table}.{n}' for n in both)} is both dropped "
                "and asked for in columns."
            )


def _select_columns(
    frame: pd.DataFrame,
    keys: Iterable[str],
    wanted: dict[str, list[str]],
    table: str,
    path: str | Path,
) -> pd.DataFrame:
    """Keep the key columns of *frame* and the metadata columns *wanted* asks for."""
    keys = list(keys)
    asked = wanted.get(table)
    if not asked:
        return frame
    available = [c for c in frame.columns if c not in keys]
    unknown = [c for c in asked if c not in available]
    if unknown:
        raise ValueError(
            f"The {table} table {str(path)!r} has no column "
            f"{', '.join(repr(c) for c in unknown)}; its metadata columns are: "
            f"{', '.join(map(str, available)) or 'none'}."
        )
    return frame[[*(k for k in keys if k in frame.columns), *asked]]


def _read_genomes(
    path: str | Path, fmt: str, wanted: dict[str, list[str]]
) -> pd.DataFrame:
    """Read a genomes table as ``create`` does: deduplicated, typed, renamed."""
    from metacooler.io.create import (
        _GENOMES_RESERVED_COLUMNS,
        _GENOMES_RESERVED_PREFIXES,
        _coerce_numeric_columns,
        _load_extra_genomes_table,
        _protect_reserved_columns,
    )

    extra = _load_extra_genomes_table(path, fmt)
    assert extra is not None  # path is not None
    extra = _select_columns(extra, ["genome_id"], wanted, "genomes", path)
    extra = _protect_reserved_columns(
        extra,
        reserved=_GENOMES_RESERVED_COLUMNS,
        keys=frozenset({"genome_id"}),
        context="Genomes table",
        reserved_prefixes=_GENOMES_RESERVED_PREFIXES,
    )
    extra = extra.drop_duplicates(subset="genome_id", keep="last")
    _coerce_numeric_columns(extra, [c for c in extra.columns if c != "genome_id"])
    return extra.reset_index(drop=True)


def _read_contigs(
    path: str | Path, fmt: str, wanted: dict[str, list[str]]
) -> tuple[pd.DataFrame, bool]:
    """Read a per-contig metadata table as ``collect`` does."""
    from metacooler.io.collect import _read_contig_metadata

    table, has_genome = _read_contig_metadata(path, fmt)
    table = _select_columns(
        table, ["contig_name", "genome_name"], wanted, "contigs", path
    )
    return table, has_genome


def _map_groups(file_path: Path, group_path: str) -> list[str]:
    """The map groups a URI addresses: one map, or every entry of a collection."""
    if group_path.strip("/"):
        return [group_path]
    if is_mmetacooler(file_path):
        entries = list_metacoolers(file_path)
        if not entries:
            raise ValueError(f"{file_path} is a collection with no entry.")
        return entries
    return ["/"]


def _catalog(
    mc: MetaCooler, contigs: pd.DataFrame, genome_names: np.ndarray
) -> pd.DataFrame:
    """Describe a map's contigs the way ``collect``'s metadata reader expects.

    ``original`` is the name without its ``<genome>__`` prefix when it carries
    one -- as every contig of a ``collect``-ed map does -- so a row may name a
    contig as in its input cool, exactly as it could at collect time.
    """
    from metacooler.io.collect import CONTIG_NAME_SEPARATOR

    names = contigs["name"].astype(str).to_numpy()
    genome = genome_names[contigs["genome_id"].to_numpy(np.int64)].astype(str)
    prefix = np.char.add(genome, CONTIG_NAME_SEPARATOR)
    names_u = names.astype(str)
    prefixed = np.char.startswith(names_u, prefix)
    stripped = np.array(
        [n[len(p) :] for n, p in zip(names_u[prefixed], prefix[prefixed], strict=True)],
        dtype=object,
    )
    original = names.astype(object)
    original[prefixed] = stripped
    return pd.DataFrame(
        {"name": names, "original": original, "genome": genome, "source": mc.uri}
    )


def _plan_map(
    mc: MetaCooler,
    group: str,
    genomes_extra: pd.DataFrame | None,
    contig_table: tuple[pd.DataFrame, bool] | None,
    contigs_path: str | Path | None,
    dropping: dict[str, list[str]],
    *,
    replace: bool,
    missing: Missing,
    unmatched: Unmatched,
    unmatched_rows: list[str],
) -> _MapEdit:
    """Work out, without writing, what annotating does to one map."""
    from metacooler.io.collect import _resolve_contig_metadata
    from metacooler.io.create import _coerce_numeric_columns

    edit = _MapEdit(uri=mc.uri, group=group)
    genomes = mc.genomes
    genome_names = genomes["name"].astype(str).to_numpy()
    supplied: dict[str, set[str]] = {}

    if genomes_extra is not None:
        keys = genomes_extra["genome_id"].astype(str)
        listed = set(keys)
        absent = [g for g in genome_names if g not in listed and g != "__unbinned__"]
        unknown = sorted(listed - set(genome_names))
        _report_genomes(mc, absent, unknown, unmatched, unmatched_rows)
        genome_values = {
            column: pd.Series(genomes_extra[column].to_numpy(), index=keys.to_numpy())
            for column in genomes_extra.columns
            if column != "genome_id"
        }
        supplied["genomes"] = set(genome_values)
        edit.tables["genomes"] = _plan_table(
            "genomes",
            genomes,
            genome_names,
            genome_values,
            matched=len(listed & set(genome_names)),
            replace=replace,
            missing=missing,
            conflicts=edit.conflicts,
        )

    if contig_table is not None:
        contigs = mc.contigs
        table, has_genome = contig_table
        assert contigs_path is not None  # contig_table comes from a path
        frame = _resolve_contig_metadata(
            table,
            has_genome,
            _catalog(mc, contigs, genome_names),
            contigs_path,
            command="annotate",
            unmatched=unmatched,
            deferred=unmatched_rows,
        )
        names = contigs["name"].astype(str).to_numpy()
        values: dict[str, pd.Series] = {}
        if frame is None:
            matched = int(table["contig_name"].isin(set(names)).sum())
        else:
            extra_columns = [
                c for c in frame.columns if c not in ("contig_name", "genome_id")
            ]
            _coerce_numeric_columns(frame, extra_columns)
            index = frame["contig_name"].to_numpy()
            values = {
                column: pd.Series(frame[column].to_numpy(), index=index)
                for column in extra_columns
            }
            matched = len(frame)
        supplied["contigs"] = set(values)
        edit.tables["contigs"] = _plan_table(
            "contigs",
            contigs,
            names,
            values,
            matched=matched,
            replace=replace,
            missing=missing,
            conflicts=edit.conflicts,
        )

    for table, columns in dropping.items():
        clash = sorted(set(columns) & supplied.get(table, set()))
        if clash:
            raise ValueError(
                f"{', '.join(f'{table}.{c}' for c in clash)} is both dropped and "
                f"given by the {table} table."
            )
        stored = mc.genomes if table == "genomes" else mc.contigs
        table_edit = edit.tables.setdefault(table, _TableEdit(table=table))
        for column in columns:
            if column in stored.columns:
                table_edit.dropped.append(column)
            else:
                edit.absent_drops.append(f"{table}.{column}")
    return edit


def _report_genomes(
    mc: MetaCooler,
    absent: Sequence[str],
    unknown: Sequence[str],
    unmatched: Unmatched = "warn",
    unmatched_rows: list[str] | None = None,
) -> None:
    """Report names the genomes table and the map disagree on.

    A genome the table leaves out is normal -- a patch, a QC run over the
    high-quality bins only -- so it is only logged; a table row naming no
    genome of the map is more often a typo or a renamed bin, so it warns,
    unless *unmatched* says otherwise.
    """
    if absent:
        logger.info(
            "%s: %d genome(s) have no row in the genomes table, e.g. %s; they "
            "get no new value from it.",
            mc.uri,
            len(absent),
            _quote(absent),
        )
    if unknown:
        report_unmatched_rows(
            unmatched,
            f"{mc.uri}: {len(unknown)} row(s) of the genomes table name no "
            "genome of the map",
            _quote(unknown),
            deferred=unmatched_rows,
        )


def _plan_table(
    table: str,
    stored: pd.DataFrame,
    row_names: np.ndarray,
    values: dict[str, pd.Series],
    *,
    matched: int,
    replace: bool,
    missing: Missing,
    conflicts: dict[str, str],
) -> _TableEdit:
    """Decide, column by column, what one table gets.

    Args:
        table: ``"genomes"`` or ``"contigs"``.
        stored: The table as the map holds it.
        row_names: Its row names, in stored order.
        values: New values per column, indexed by row name (deduplicated).
        matched: How many table rows matched a row of the map.
        replace: Whether an existing column may change.
        missing: What unlisted rows of a replaced column get.
        conflicts: Filled with the columns that would change without *replace*.

    Returns:
        The table's edit.
    """
    edit = _TableEdit(table=table, matched=matched)
    for column, series in values.items():
        aligned = series.reindex(row_names)
        new = aligned.to_numpy()
        if column not in stored.columns:
            edit.values[column] = new
            edit.added.append(column)
            continue
        old = stored[column].to_numpy()
        if missing == "keep":
            listed = pd.Index(series.index).get_indexer(row_names) >= 0
            new = _combine(old, new, listed)
        if _same_values(old, new):
            edit.unchanged.append(column)
            continue
        if not replace:
            conflicts[f"{table}.{column}"] = _describe_change(old, new, row_names)
            continue
        edit.values[column] = new
        edit.replaced.append(column)
    return edit


def _combine(old: np.ndarray, new: np.ndarray, listed: np.ndarray) -> np.ndarray:
    """New values on the listed rows, old ones elsewhere, typed as ``create`` would."""
    merged = np.where(listed, new.astype(object), old.astype(object))
    try:
        return pd.to_numeric(pd.Series(merged), errors="raise").to_numpy()
    except (ValueError, TypeError):
        return merged


def _is_numeric(values: np.ndarray) -> bool:
    return values.dtype.kind in "iufb"


def _as_text(values: np.ndarray) -> np.ndarray:
    return np.array(["" if pd.isna(v) else str(v) for v in values], dtype=object)


def _same_values(old: np.ndarray, new: np.ndarray) -> bool:
    """Whether a column would hold the same values after the edit."""
    if _is_numeric(old) and _is_numeric(new):
        return bool(
            np.array_equal(
                old.astype(np.float64), new.astype(np.float64), equal_nan=True
            )
        )
    if _is_numeric(old) or _is_numeric(new):
        return False
    return bool(np.array_equal(_as_text(old), _as_text(new)))


def _describe_change(old: np.ndarray, new: np.ndarray, row_names: np.ndarray) -> str:
    """Say how many rows a refused replacement would change, with one example."""
    if _is_numeric(old) and _is_numeric(new):
        differs = ~(
            (old.astype(np.float64) == new.astype(np.float64))
            | (pd.isna(old) & pd.isna(new))
        )
    else:
        differs = _as_text(old) != _as_text(new)
    rows = np.flatnonzero(differs)
    if not rows.size:
        return "its type would change"
    first = int(rows[0])
    return (
        f"{rows.size} row(s) would change, e.g. {str(row_names[first])!r}: "
        f"{_show(old[first])} -> {_show(new[first])}"
    )


def _show(value: object) -> str:
    """Render one stored or table value for a message."""
    if pd.isna(value):
        return "(empty)"
    if isinstance(value, np.generic):
        value = value.item()
    return repr(value)


def _check_plans(
    edits: Sequence[_MapEdit],
    genomes_table: str | Path | None,
    contigs_table: str | Path | None,
    dropping: dict[str, list[str]],
    unmatched_rows: Sequence[str] = (),
) -> None:
    """Refuse what would be wrong in every map, before any is written.

    *unmatched_rows* are the rows refused under ``unmatched="error"``, which
    only come after a table that matches nothing at all.
    """
    for table, path in (("genomes", genomes_table), ("contigs", contigs_table)):
        if path is None:
            continue
        if not any(edit.tables[table].matched for edit in edits):
            where = edits[0].uri if len(edits) == 1 else "any entry of the collection"
            raise ValueError(
                f"The {table} table {str(path)!r} matches no {table[:-1]} of "
                f"{where}: it is the wrong table, or --{table}-format points at "
                "the wrong column."
            )
    if unmatched_rows:
        raise ValueError(unmatched_rows_error(unmatched_rows))
    conflicts = {
        f"{name} ({edit.uri})" if len(edits) > 1 else name: detail
        for edit in edits
        for name, detail in edit.conflicts.items()
    }
    if conflicts:
        detail = "\n".join(f"  {name}: {what}" for name, what in conflicts.items())
        raise ValueError(
            f"The tables would change {len(conflicts)} existing column(s):\n"
            f"{detail}\nPass --replace to change them (and --missing na to empty "
            "the rows the table does not list), or --columns to leave them out."
        )
    asked = {f"{t}.{c}" for t, cs in dropping.items() for c in cs}
    nowhere = sorted(
        name for name in asked if all(name in edit.absent_drops for edit in edits)
    )
    if nowhere:
        raise ValueError(
            f"Cannot drop {', '.join(nowhere)}: no such column in "
            f"{edits[0].uri if len(edits) == 1 else 'any entry of the collection'}."
        )
    for edit in edits:
        for table_edit in edit.tables.values():
            if table_edit.table != "contigs":
                continue
            violations = [
                ("contigs", column, bad)
                for column, values in table_edit.values.items()
                if values.dtype.kind in "OSU" and (bad := non_ascii_values(values))
            ]
            if violations:
                raise MetaCoolerFormatError(
                    f"Refusing to write {edit.uri}: these /contigs values are not "
                    "ASCII:\n" + format_ascii_violations(violations)
                )


def _write_map(
    path: Path,
    group: str,
    edit: _MapEdit,
    metadata: dict | None,
    history: dict,
    *,
    always: bool = False,
) -> None:
    """Apply one map's edit, in place, and record it in the map's metadata."""
    if not edit.changed and not always:
        logger.info("%s already holds these values; nothing to write.", edit.uri)
        return
    with h5py.File(path, "r+") as f:
        grp = f[group.strip("/")] if group.strip("/") else f
        if _SCRATCH in grp:
            del grp[_SCRATCH]
        changed = [t for t in edit.tables.values() if t.changed]
        # Build every changed table before swapping any in: rebuilding them one
        # at a time left /genomes edited and /contigs not (and no history) when
        # the second build failed, plus the scratch group behind.
        try:
            for table_edit in changed:
                _build_table(grp, table_edit)
        except BaseException:
            if _SCRATCH in grp:
                del grp[_SCRATCH]
            raise
        for table_edit in changed:
            del grp[table_edit.table]
            grp.move(f"{_SCRATCH}/{table_edit.table}", table_edit.table)
        if _SCRATCH in grp:
            del grp[_SCRATCH]
        _record(grp, edit, metadata, history)
    _log_edit(edit)


def _build_table(grp: h5py.Group, edit: _TableEdit) -> None:
    """Build one table with its edit applied, in the scratch group.

    Untouched columns are copied as they are -- dtype, chunking and filters
    included -- and the column order is kept, a replaced column staying where
    it was and added ones going last.  The table is never visible half-built:
    the caller swaps it into place only once every changed table is built.
    """
    old = grp[edit.table]
    scratch = grp.require_group(_SCRATCH)
    new = scratch.create_group(edit.table, track_order=True)
    for key, value in old.attrs.items():
        new.attrs[key] = value
    for column in old:
        if column in edit.dropped:
            continue
        if column in edit.values:
            write_column(new, column, edit.values[column])
        else:
            old.copy(column, new, name=column)
    for column in edit.added:
        write_column(new, column, edit.values[column])


def _record(
    grp: h5py.Group, edit: _MapEdit, metadata: dict | None, history: dict
) -> None:
    """Append this run to the map's provenance (``commands`` and ``annotate``)."""
    raw = grp.attrs.get("metadata", "{}")
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    current = json.loads(raw) if isinstance(raw, str) else dict(raw)
    for key, value in (metadata or {}).items():
        if key == "commands":
            commands = list(current.get("commands", []))
            commands.extend(value if isinstance(value, list) else [value])
            current["commands"] = commands
        else:
            current[key] = value
    entry = dict(history)
    for kind in ("added", "replaced", "dropped"):
        entry[kind] = [
            f"{table_edit.table}.{column}"
            for table_edit in edit.tables.values()
            for column in getattr(table_edit, kind)
        ]
    current["annotate"] = [*current.get("annotate", []), entry]
    grp.attrs["metadata"] = json.dumps(current)


def _log_edit(edit: _MapEdit) -> None:
    """Say what was written to one map."""
    parts = []
    for kind in ("added", "replaced", "dropped"):
        names = [
            f"{t.table}.{c}" for t in edit.tables.values() for c in getattr(t, kind)
        ]
        if names:
            parts.append(f"{kind} {', '.join(names)}")
    logger.info("%s: %s", edit.uri, "; ".join(parts) if parts else "no change")


def _quote(names: Sequence[str]) -> str:
    """Quote the first few *names*, noting how many were cut."""
    shown = ", ".join(repr(n) for n in names[:_MAX_QUOTED])
    return shown + (
        f" (+{len(names) - _MAX_QUOTED} more)" if len(names) > _MAX_QUOTED else ""
    )
