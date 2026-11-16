"""Column specs: which columns of a user table hold its names.

Every option that tells metacooler how to read a user table --
``--contigs-format`` and ``--genomes-format`` -- takes the same grammar::

    SPEC   := ALIAS | COLUMN ["," GENOME] ["," "noheader"]
    COLUMN := a 1-based position (digits only) | a header name (anything else)
    GENOME := COLUMN | "NA"

A table's *key* columns are the ones holding names: the contig and the genome
of a contig-to-genome assignment table, the contig of a per-contig metadata
table, the genome of a genomes table.  Every other column is metadata.

Positions and names may be mixed.  A name needs a header line to be looked up
in, so ``Name,Final_bin,noheader`` is refused.  ``NA`` in the genome's place
says that the table has *no* genome column: no contig is assigned to a genome,
so ``create`` builds a map with every contig in ``__unbinned__`` (nothing is
binned yet) and ``reassign`` unbins the contigs the table lists.  It has to be
written out: a spec that simply leaves the genome column off is refused where
one is needed, so a forgotten column cannot pass for "no binning".  Names are only ever used to
*find* positions: the key columns are then relabelled by position (see
:func:`metacooler.io._utils.rename_key_columns`), never with ``df.rename``.

This module is pure Python with no heavy import, so the CLI can check a spec
while it parses its options, before any file is opened.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path

#: Trailing token of a spec saying that the table has no header line.
NOHEADER = "noheader"
#: Token, in the genome's place, saying that the table has no genome column
#: (any case).  The same ``NA`` that, as a genome *label*, means "unbinned".
NO_GENOME = "NA"

#: ``--contigs-format`` aliases: the layout of each contig-to-bin table a common
#: binner writes, taken from the tool's own writer rather than from memory --
#: an alias that names the wrong columns would misread a table without error.
#:
#: - ``metator``: MetaTOR's ``contig_data_final.txt``.
#: - ``dastool``: DAS Tool's ``*_DASTool_contigs2bin.tsv``, and the tables its
#:   ``Fasta_to_Contigs2Bin.sh`` builds from FASTA bins (MetaBAT2, MaxBin2, ...).
#: - ``binette``: Binette's ``final_contig_to_bin.tsv`` (same headerless layout).
#: - ``semibin``: SemiBin2's ``contig_bins.tsv``, whose unbinned contigs carry
#:   bin ``-1``.
#: - ``vamb``: VAMB's ``vae_clusters_unsplit.tsv`` / ``vae_clusters_split.tsv``,
#:   which put the cluster first.
#:
#: CONCOCT has no alias: its clustering tables are comma-separated.
CONTIGS_FORMAT_ALIASES: dict[str, str] = {
    "contigs": "contig_name,genome_name",
    "contig_genome": "1,2,noheader",
    "metator": "Name,Final_bin",
    "dastool": "1,2,noheader",
    "binette": "1,2,noheader",
    "semibin": "contig,bin",
    "vamb": "contigname,clustername",
}

#: Default ``--contigs-format`` of an assignment table (``create``, ``reassign``).
CONTIGS_FORMAT_DEFAULT = "contigs"
#: Default ``--contigs-format`` of a per-contig metadata table (``collect``,
#: ``annotate``).
CONTIG_METADATA_FORMAT_DEFAULT = "contig_name"
#: Default ``--genomes-format``.
GENOMES_FORMAT_DEFAULT = "genome_name"

# Keys of each kind of table, for messages: an assignment table needs both, a
# metadata table only the contig (its genome column is an optional check).
ASSIGNMENT_KEYS = ("contig", "genome")
GENOMES_KEYS = ("genome",)

# How many header names to quote when a named column is missing.
_MAX_QUOTED_COLUMNS = 20


@dataclasses.dataclass(frozen=True)
class ColumnSpec:
    """A parsed column spec.

    Attributes:
        columns: The key columns, in order: an ``int`` is a 1-based position,
            a ``str`` a header name, ``None`` (from ``NA``) a genome column the
            table does not have.
        header: Whether the table's first line names its columns.
        text: The spec as given -- an alias with its expansion -- for messages.
        option: The option the spec came from, for messages.
    """

    columns: tuple[int | str | None, ...]
    header: bool
    text: str
    option: str

    @property
    def no_genome(self) -> bool:
        """Whether the spec says, with ``NA``, that there is no genome column."""
        return len(self.columns) > 1 and self.columns[1] is None


def alias_listing(aliases: Mapping[str, str]) -> str:
    """Render *aliases* as ``name (expansion), ...`` for help texts and errors."""
    return ", ".join(f"{name} ({spec})" for name, spec in aliases.items())


def parse_column_spec(
    spec: str,
    *,
    keys: Sequence[str],
    required: int | None = None,
    aliases: Mapping[str, str] | None = None,
    option: str,
) -> ColumnSpec:
    """Parse a column spec, checking its grammar and how many columns it names.

    Args:
        spec: The spec as the user wrote it.
        keys: What the key columns hold, in order (e.g. ``("contig",
            "genome")``); the spec may name at most this many.
        required: How many key columns the spec must name at least.  Defaults
            to all of *keys*.
        aliases: Named specs accepted in place of a column list.  An alias wins
            over a header name it happens to equal; such a column can still be
            given by position.
        option: The option *spec* came from, for messages.

    Returns:
        The parsed spec.

    Raises:
        ValueError: If *spec* is malformed, names the wrong number of columns,
            gives a position below 1, gives a name with ``noheader``, puts
            ``NA`` anywhere but in the genome's place, or names one column
            twice.
    """
    aliases = aliases or {}
    required = len(keys) if required is None else required
    raw = spec.strip()
    expanded = aliases.get(raw, raw)
    text = f"{raw} (= {expanded})" if raw in aliases else raw

    tokens = [token.strip() for token in expanded.split(",")]
    header = True
    if len(tokens) > 1 and tokens[-1] == NOHEADER:
        header = False
        tokens = tokens[:-1]
    if any(not token for token in tokens):
        raise ValueError(f"{option} {text!r} has an empty column.")

    if not required <= len(tokens) <= len(keys):
        raise ValueError(_count_message(text, tokens, keys, required, aliases, option))

    columns: list[int | str | None] = []
    for index, token in enumerate(tokens):
        if token.upper() == NO_GENOME:
            if index == 0:
                raise ValueError(_no_genome_misplaced(text, keys, option))
            columns.append(None)
            continue
        if token.isdigit():
            position = int(token)
            if position < 1:
                raise ValueError(
                    f"{option} {text!r}: column positions are 1-based (>= 1), "
                    f"got {token}."
                )
            columns.append(position)
        elif not header:
            raise ValueError(
                f"{option} {text!r}: {token!r} is a column name, but "
                f"',{NOHEADER}' says the table has no header line to find it in. "
                "Give its 1-based position instead."
            )
        else:
            columns.append(token)
    named = [column for column in columns if column is not None]
    if len(set(named)) < len(named):
        raise ValueError(_same_column_message(text, named, keys, option))
    return ColumnSpec(columns=tuple(columns), header=header, text=text, option=option)


def _no_genome_misplaced(text: str, keys: Sequence[str], option: str) -> str:
    """Explain an ``NA`` given where a column of names is needed."""
    if len(keys) == 1:
        return (
            f"{option} {text!r}: {NO_GENOME} stands for a missing genome column, "
            f"but a {keys[0]}s table is nothing without its {keys[0]} names."
        )
    return (
        f"{option} {text!r}: {NO_GENOME} stands for a missing genome column, so "
        f"it only goes in the genome's place, after the {keys[0]} column -- e.g. "
        f"'1,{NO_GENOME}' or '{keys[0]}_name,{NO_GENOME}'."
    )


def _count_message(
    text: str,
    tokens: Sequence[str],
    keys: Sequence[str],
    required: int,
    aliases: Mapping[str, str],
    option: str,
) -> str:
    """Explain a spec that names too few or too many columns."""
    wanted = " and the ".join(keys) if required == len(keys) else keys[0]
    examples = (
        "'1,2' or 'contig_name,genome_name'"
        if tuple(keys) == ASSIGNMENT_KEYS and required == 2
        else f"'1' or '{'contig_name' if keys[0] == 'contig' else 'genome_name'}'"
    )
    if len(tokens) == 1 and len(keys) > 1 and required == len(keys):
        # A bare word given where two columns are needed is most likely a
        # misspelled alias, so say that first.
        listing = f"it is not an alias ({', '.join(aliases)}), and " if aliases else ""
        return (
            f"Unknown {option} {text!r}: {listing}a column spec here names two "
            f"columns, the {wanted}, each a 1-based position or a header name "
            f"-- e.g. {examples}."
        )
    count = f"{len(tokens)} column{'s' if len(tokens) != 1 else ''}"
    expected = f"{len(keys)}" if required == len(keys) else f"{required} to {len(keys)}"
    return (
        f"{option} {text!r} names {count}, but this table takes {expected}: the "
        f"{wanted} column{'s' if required > 1 else ''}, each a 1-based position "
        f"or a header name -- e.g. {examples}."
    )


def _same_column_message(
    text: str, columns: Sequence[int | str | None], keys: Sequence[str], option: str
) -> str:
    """Explain a spec that points two key columns at the same column."""
    repeated = next(c for c in columns if list(columns).count(c) > 1)
    return (
        f"{option} {text!r} points the {' and the '.join(keys[: len(columns)])} "
        f"columns at the same column ({repeated}). Give two different columns, "
        "e.g. '1,2'."
    )


def resolve_column_spec(
    spec: ColumnSpec,
    labels: Sequence[object],
    *,
    kind: str,
    path: str | Path,
    keys: Sequence[str],
    hint: str = "",
) -> list[int | None]:
    """Find the 0-based position of every key column of *spec* in a table.

    Args:
        spec: The parsed spec.
        labels: The table's column labels, as read (header names, or pandas'
            positional integers for a headerless table).
        kind: Human-readable table name, e.g. ``"Contig table"``.
        path: The table, for messages.
        keys: What the key columns hold, in the order of ``spec.columns``.
        hint: Extra advice appended when a named column is missing (e.g. the
            aliases on offer).

    Returns:
        One 0-based position per entry of ``spec.columns``, in the same order;
        ``None`` for an ``NA`` genome column.

    Raises:
        ValueError: If the table is too narrow for the spec, a named column is
            absent from its header, or two entries resolve to one column.
    """
    names = [str(label) for label in labels]
    wanted = sum(column is not None for column in spec.columns)
    if len(names) < wanted:
        what = " and the ".join(keys[:wanted])
        raise ValueError(
            f"{kind} needs at least {wanted} column{'s' if wanted > 1 else ''} "
            f"(the {what} name{'s' if wanted > 1 else ''}), but {path} has "
            f"{len(names)}."
        )
    positions: list[int | None] = []
    for column, key in zip(spec.columns, keys, strict=False):
        if column is None:
            positions.append(None)
            continue
        if isinstance(column, int):
            if column > len(names):
                raise ValueError(
                    f"{spec.option} {spec.text!r} requires at least {column} "
                    f"columns in the table, got {len(names)}: {path}"
                )
            positions.append(column - 1)
            continue
        if column not in names:
            quoted = ", ".join(names[:_MAX_QUOTED_COLUMNS])
            more = " ..." if len(names) > _MAX_QUOTED_COLUMNS else ""
            raise ValueError(
                f"{kind} has no column named {column!r}, which {spec.option} "
                f"{spec.text!r} asks for: {path}\n"
                f"  Its header is: {quoted}{more}\n"
                f"  Pass {spec.option} to say which column holds the {key} "
                f"names, by header name or 1-based position{hint}."
            )
        positions.append(names.index(column))
    present = [position for position in positions if position is not None]
    if len(set(present)) < len(present):
        raise ValueError(
            f"{spec.option} {spec.text!r} points the "
            f"{' and the '.join(keys[: len(present)])} columns at the same "
            f"column ({present[0] + 1}) of {path}. Give two different columns."
        )
    return positions
