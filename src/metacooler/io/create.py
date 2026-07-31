"""Create metacool files from pairs, bins, and MAG assignments."""

from __future__ import annotations

import contextlib
import dataclasses
import gzip
import logging
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import h5py
import numpy as np
import pandas as pd
from pyfaidx import Fasta

from metacooler.io._utils import (
    _finalize_metacool_file,
    check_ascii_names,
    check_output_path,
    resolve_count_dtype,
)
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.genomes import parse_genomes_table
from metacooler.schema import (
    CONTIGS_DTYPE,
    GENOMES_DTYPE,
    validate,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from typing import Literal


_UNBINNED_SENTINELS: frozenset[str] = frozenset(
    {"__unbinned__", "na", "nan", "n/a", "nd", "n/d", ""}
)
# Row budget for the in-memory pixel sort, shared by ``from_pairs`` and
# ``from_cool``. Beyond it, pixels go through cooler's external sort instead.
#
# This counts DEDUPLICATED pixels. ``_buffer_pixels`` folds its buffer as it
# fills (see that function), so the number compared here is the real pixel table
# being built, not the raw per-chunk row count. That distinction is what keeps
# the threshold honest: Hi-C pixel counts saturate — on a fixed bin grid, 20x
# more pairs yields well under 2x more pixels — so an unfolded row count
# overestimates the true memory need by an order of magnitude and spills to the
# external sort for datasets that would comfortably have sorted in RAM.
# At 3 int64 columns, 100 M pixels is ~2.4 GB resident plus the fold transient.
_SORT_IN_MEMORY_PIXEL_THRESHOLD = 100_000_000  # 100 M deduplicated pixels
# Fraction of the row budget that may accumulate unfolded before ``_buffer_pixels``
# re-aggregates. Folding on every chunk would re-scan the accumulator far too
# often; folding only at the budget would defeat the point. A quarter keeps the
# unfolded overhang bounded while amortising the fold cost across many chunks.
_PIXEL_FOLD_FRACTION = 4
# Ceiling on the raw text block pyarrow buffers per read batch. Binned pixels are
# coalesced back up to the requested chunksize afterwards, so raising chunksize
# costs chunk-granularity, not read memory.
_PAIRS_READ_BLOCK_MAX_BYTES = 128 * 1024**2  # 128 MiB (~1.2 M pairs)
# Column names metacooler writes itself into /genomes and /contigs.  A user table
# is free to carry extra columns, but one of these names would collide.
_GENOMES_RESERVED_COLUMNS = frozenset({"name", "length", "n_contigs", "n_bins"})
_CONTIGS_RESERVED_COLUMNS = frozenset({"name", "length", "genome_id"})
# How many distinct unknown contig names to remember for the dropped-pair report.
_MAX_REPORTED_UNKNOWN_CONTIGS = 10
# How many genome names to name when the contigs and genomes tables diverge.
_MAX_REPORTED_DIVERGENT_GENOMES = 10

logger = logging.getLogger("metacooler")


def _read_contig_lengths(
    bins: pd.DataFrame | None = None,
    fasta: str | Path | None = None,
    pairs: str | Path | None = None,
) -> dict[str, int]:
    """Extract contig name → length mapping from one of three sources.

    Exactly one of *bins*, *fasta*, or *pairs* should be provided.  When
    *fasta* is given, sequence lengths are read directly from the FASTA
    index.  When *bins* is given, the maximum ``end`` value per contig is
    used as the length.  When *pairs* is given, ``#chromsize:`` header
    lines are parsed.

    Args:
        bins: Pre-loaded bins DataFrame with columns ``chrom``, ``start``,
            ``end``.  Lengths are inferred as ``max(end)`` per contig.
        fasta: Path to a FASTA (or gzip-compressed ``.fa.gz``) file.  The
            file is indexed with pyfaidx to retrieve exact sequence lengths.
        pairs: Path to a ``.pairs`` (or ``.pairs.gz``) file whose header
            contains ``#chromsize: <name> <length>`` lines.

    Returns:
        Mapping of contig name to sequence length in base pairs.

    Raises:
        ValueError: If none of the three sources can provide contig lengths.
    """
    if fasta is not None:
        logging.info("Extracting contig lengths from FASTA index")
        fasta_path = Path(fasta)
        if fasta_path.suffix == ".gz":
            # pyfaidx doesn't handle plain gzip; stream-decompress to temp file
            logging.info("Unzipping FASTA to temporary file for pyfaidx indexing")
            tmp_fd = tempfile.NamedTemporaryFile(mode="wb", suffix=".fa", delete=False)  # noqa: SIM115
            fa_path = Path(tmp_fd.name)
            fai_path = fa_path.parent / (fa_path.name + ".fai")
            fa = None
            try:
                with gzip.open(fasta_path, "rb") as gz:
                    shutil.copyfileobj(gz, tmp_fd)
                tmp_fd.close()
                fa = Fasta(str(fa_path))
                result = {name: len(record) for name, record in fa.items()}
            finally:
                with contextlib.suppress(Exception):
                    if fa is not None:
                        fa.close()
                fa_path.unlink(missing_ok=True)
                fai_path.unlink(missing_ok=True)
            return result
        fa = Fasta(str(fasta))
        try:
            return {name: len(record) for name, record in fa.items()}
        finally:
            with contextlib.suppress(Exception):
                fa.close()
    if bins is not None:
        # Infer from max end per chrom
        return bins.groupby("chrom")["end"].max().to_dict()
    if pairs is not None:
        # Try to parse chromsizes from pairs header
        lengths = _parse_pairs_header(pairs)
        if lengths:
            return lengths
    raise ValueError(
        "Contig lengths required. Provide --fasta, a BED bins file, "
        "or a pairs file with #chromsize: header lines."
    )


def _apply_reference_order(
    contig_lengths: dict[str, int],
    fasta: str | Path | None = None,
    pairs: str | Path | None = None,
) -> dict[str, int]:
    """Reorder *contig_lengths* to follow the FASTA / pairs-header contig order.

    ``contig_ordering="fasta"`` simply preserves the insertion order of
    *contig_lengths*.  When those lengths were derived from a bins table, that
    order is whatever ``groupby("chrom")`` produced — i.e. alphabetical — yet the
    file would still record ``contig-ordering="fasta"``.  Two files could then
    advertise the same ordering with different physical layouts, which
    ``merge`` would happily scramble.

    Args:
        contig_lengths: Mapping of contig name to length, in arbitrary order.
        fasta: Assembly FASTA whose record order defines the reference order.
        pairs: Pairs file whose ``#chromsize:`` header defines the reference
            order when no *fasta* is available.

    Returns:
        The same mapping, reordered to follow the reference order.  Contigs
        absent from the reference are appended (with a warning).

    Raises:
        ValueError: If neither a FASTA nor a pairs header can supply an order.
    """
    reference: list[str] | None = None
    if fasta is not None:
        reference = list(_read_contig_lengths(fasta=fasta))
    elif pairs is not None:
        header = _parse_pairs_header(pairs)
        if header:
            reference = list(header)
    if reference is None:
        raise ValueError(
            "contig_ordering='fasta' needs a reference contig order, but no FASTA "
            "was given and the pairs file has no '#chromsize:' header lines. "
            "Contig lengths came from the bins table, whose order is alphabetical "
            "— writing it as 'fasta' would mislabel the file. Provide a FASTA, add "
            "'#chromsize:' lines to the pairs header, or choose another ordering "
            "(genome-then-length, genome-then-alpha, alpha)."
        )
    ordered = {c: contig_lengths[c] for c in reference if c in contig_lengths}
    leftover = [c for c in contig_lengths if c not in ordered]
    if leftover:
        logger.warning(
            "%d contig(s) are absent from the reference (FASTA / pairs header) "
            "order and were appended at the end: %s",
            len(leftover),
            leftover[:10],
        )
        ordered.update({c: contig_lengths[c] for c in leftover})
    return ordered


def _drop_nonpositive_contigs(
    contig_lengths: dict[str, int],
    bins_df: pd.DataFrame | None = None,
) -> tuple[dict[str, int], pd.DataFrame | None]:
    """Drop contigs with non-positive length (and their bins, if any).

    A contig whose length is ``<= 0`` (e.g. an empty FASTA record or a
    malformed ``#chromsize:`` header) tiles to zero bins.  Such a contig still
    occupies a slot in ``chrom_offset`` with an empty range, which makes the
    out-of-range bin clamp in :func:`_bin_pairs_chunk` fold that contig's
    contacts into a *neighbouring* contig's bin — silent cross-genome
    corruption.  Removing these contigs up front avoids the problem entirely.

    Args:
        contig_lengths: Mapping of contig name to length.
        bins_df: Optional bins DataFrame (string ``chrom`` names) to filter in
            lockstep.

    Returns:
        ``(contig_lengths, bins_df)`` with non-positive-length contigs removed.
    """
    bad = [c for c, ln in contig_lengths.items() if ln <= 0]
    if bad:
        logger.warning(
            "Dropping %d contig(s) with non-positive length: %s",
            len(bad),
            bad[:10],
        )
        contig_lengths = {c: ln for c, ln in contig_lengths.items() if ln > 0}
        if bins_df is not None:
            bins_df = bins_df[bins_df["chrom"].isin(contig_lengths)].copy()
    return contig_lengths, bins_df


def _parse_pairs_header(pairs_path: str | Path) -> dict[str, int] | None:
    """Parse ``#chromsize:`` lines from the header of a pairs file.

    Reads lines from the beginning of the file until the first non-comment
    line is encountered.  Handles plain-text and gzip-compressed files.

    Args:
        pairs_path: Path to the pairs (or ``.pairs.gz``) file.

    Returns:
        Mapping of contig name to length extracted from the header, or
        ``None`` if no ``#chromsize:`` lines were found.
    """
    path = Path(pairs_path)
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    lengths = {}
    with opener(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("#chromsize:"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        lengths[parts[1]] = int(parts[2])
                    except ValueError as exc:
                        raise ValueError(
                            f"Malformed #chromsize line in {path}: {line.strip()!r} "
                            f"(length {parts[2]!r} is not an integer)"
                        ) from exc
    return lengths if lengths else None


def _genome_id_order(genome_names: Iterable[str]) -> list[str]:
    """Return genome names in canonical ``genome_id`` order.

    The order is ``__unbinned__`` first (index 0, when present), then remaining
    names alphabetically.  This MUST match the row order produced by
    :func:`_build_genomes_table` so that, for genome-grouped contig orderings,
    the physical contig layout is grouped by genome in ascending ``genome_id``
    order — the precondition for ``genome_offset`` to be a valid offset array.

    Args:
        genome_names: Iterable of genome names.

    Returns:
        List of the unique names in canonical ``genome_id`` order.
    """
    names = set(genome_names)
    ordered = sorted(names - {"__unbinned__"})
    if "__unbinned__" in names:
        ordered = ["__unbinned__"] + ordered
    return ordered


def _build_contig_order(
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    ordering: Literal[
        "genome-then-alpha", "genome-then-length", "fasta", "alpha"
    ] = "genome-then-alpha",
) -> list[str]:
    """Return contig names sorted according to *ordering*.

    Args:
        contig_lengths: Mapping of contig name to length.  Keys determine
            the full set of contigs to order.
        contig_genome_map: Mapping of contig name to genome name.  Contigs
            not present default to ``"__unbinned__"``.
        ordering: Strategy used to sort contigs:

            * ``"genome-then-alpha"`` — group by genome (alphabetical),
              then sort contigs alphabetically within each group.
            * ``"genome-then-length"`` — group by genome (alphabetical),
              then sort contigs longest-first within each group.
            * ``"fasta"`` — preserve insertion order of *contig_lengths*
              (i.e. FASTA record order).
            * ``"alpha"`` — strict global alphabetical order.

    Returns:
        Ordered list of all contig names from *contig_lengths*.

    Raises:
        ValueError: If *ordering* is an unrecognised string.
    """
    all_contigs = set(contig_lengths.keys())

    if ordering == "fasta":
        # Preserve dict insertion order (Python 3.7+)
        return list(contig_lengths.keys())

    if ordering == "alpha":
        return sorted(all_contigs)

    if ordering == "genome-then-alpha":
        # Group by genome, then alphabetical order within each genome.
        # Genomes are visited in genome_id order (see _genome_id_order) so the
        # physical contig layout matches the genomes-table ordering.
        genome_groups: dict[str, list[str]] = {}
        for contig in all_contigs:
            genome = contig_genome_map.get(contig, "__unbinned__")
            genome_groups.setdefault(genome, []).append(contig)

        ordered = []
        for genome in _genome_id_order(genome_groups.keys()):
            contigs = genome_groups[genome]
            contigs.sort()
            ordered.extend(contigs)
        return ordered

    if ordering == "genome-then-length":
        # Group by genome (in genome_id order), then longest contig first
        # within each genome.
        genome_groups2: dict[str, list[str]] = {}
        for contig in all_contigs:
            genome = contig_genome_map.get(contig, "__unbinned__")
            genome_groups2.setdefault(genome, []).append(contig)

        ordered2 = []
        for genome in _genome_id_order(genome_groups2.keys()):
            contigs = genome_groups2[genome]
            contigs.sort(key=lambda c: (-contig_lengths[c], c))
            ordered2.extend(contigs)
        return ordered2

    raise ValueError(f"Unknown ordering: {ordering!r}")


def _normalize_genome_name(v: str) -> str:
    """Normalize sentinel values (empty, NA, nan, ND, n/d, …) to ``__unbinned__``."""
    return "__unbinned__" if v.strip().lower() in _UNBINNED_SENTINELS else v


def _protect_reserved_columns(
    df: pd.DataFrame,
    reserved: frozenset[str],
    keys: frozenset[str],
    context: str,
) -> pd.DataFrame:
    """Rename user columns that collide with metacooler's reserved column names.

    A user table is free to carry extra columns (``docs/format.md`` advertises
    this), but a column called ``length``/``name``/``n_bins`` would collide with
    the columns metacooler itself writes: the ``merge`` below would add ``_x``/
    ``_y`` suffixes (``KeyError: "['length'] not in index"``) or the later
    ``insert("n_bins")`` would fail outright.  Rename them out of the way
    instead, loudly.

    Args:
        df: User-supplied table.
        reserved: Column names metacooler writes itself.
        keys: Column names in *df* that are merge keys and must be left alone.
        context: Human-readable table name used in the warning.

    Returns:
        *df* unchanged, or a renamed copy.
    """
    renames: dict[str, str] = {}
    taken = set(df.columns)
    for col in df.columns:
        if col in keys or col not in reserved:
            continue
        new = f"user_{col}"
        while new in taken:
            new = f"user_{new}"
        renames[col] = new
        taken.add(new)
    if not renames:
        return df
    logger.warning(
        "%s: column(s) %s collide with metacooler's reserved column names; "
        "renamed to %s.",
        context,
        ", ".join(repr(k) for k in renames),
        ", ".join(repr(v) for v in renames.values()),
    )
    return df.rename(columns=renames)


def _report_genome_table_divergence(
    derived: list[str], extra_genomes_df: pd.DataFrame
) -> None:
    """Report genome names the contigs table and the genomes table disagree on.

    ``/genomes`` is *derived* from the contig→genome assignments; a ``--genomes``
    table only contributes extra metadata columns on top of it, via a left merge
    on the genome name.  Two consequences are easy to mistake for metacooler
    inventing or losing genomes, and both used to happen in complete silence:

    - A genome named in the contigs table but absent from the genomes table
      still gets a row, with ``length``/``n_contigs``/``n_bins`` computed from
      its contigs and empty metadata.  On a real MetaTOR dataset this added 66
      ``MetaTOR_MGE_*`` bins that the MAG-only genomes table did not list, and
      they looked fabricated because their metadata columns were blank.
    - A genome listed in the genomes table with no contigs assigned to it has
      nothing to derive a row from and never reaches the output at all.

    Both are legitimate — a genomes table may cover a superset or a subset of
    the bins — so this warns rather than raising.

    Args:
        derived: Genome names derived from the contig assignments.
        extra_genomes_df: Parsed ``--genomes`` table, keyed by ``genome_id``.
    """
    listed = set(extra_genomes_df["genome_id"])
    # __unbinned__ is metacooler's own synthetic catch-all, never a table row.
    unlisted = sorted(set(derived) - listed - {"__unbinned__"})
    unused = sorted(listed - set(derived))
    for names, message in (
        (
            unlisted,
            "%d genome(s) in the contigs table have no row in the genomes "
            "table: %s%s. They are kept, with length/n_contigs/n_bins computed "
            "from their contigs and empty metadata columns.",
        ),
        (
            unused,
            "%d genome(s) in the genomes table have no contig assigned to them "
            "in the contigs table: %s%s. They are absent from /genomes.",
        ),
    ):
        if names:
            logger.warning(
                message,
                len(names),
                ", ".join(names[:_MAX_REPORTED_DIVERGENT_GENOMES]),
                ", ..." if len(names) > _MAX_REPORTED_DIVERGENT_GENOMES else "",
            )


def _build_genomes_table(
    contig_order: list[str],
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    extra_genomes_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the ``/genomes`` DataFrame from contig data.

    Args:
        contig_order: Ordered list of contig names (typically the output of
            :func:`_build_contig_order`).
        contig_lengths: Mapping of contig name to length.
        contig_genome_map: Mapping of contig name to genome name.  Contigs
            absent from this map are assigned to ``"__unbinned__"``.
        extra_genomes_df: Optional DataFrame from a genomes-metadata TSV
            file (column ``genome_id`` plus any additional columns).  When
            provided, the extra columns are left-merged into the output by
            genome name.  Columns colliding with metacooler's own
            (``name``, ``length``, ``n_contigs``, ``n_bins``) are renamed with a
            ``user_`` prefix.

    Returns:
        DataFrame with columns ``name``, ``length``, ``n_contigs``, and
        any extra columns from *extra_genomes_df*.  Rows are ordered with
        ``"__unbinned__"`` first (when present), then alphabetically.
    """
    # Collect genome names
    genome_names = set()
    for contig in contig_order:
        genome_names.add(contig_genome_map.get(contig, "__unbinned__"))

    # Only inject __unbinned__ when at least one contig actually lands in it;
    # the loop above already adds "__unbinned__" to genome_names when any contig
    # lacks an assignment, so no extra injection is needed here.

    # Canonical genome_id order: __unbinned__ first (if present), then alphabetical.
    sorted_names = _genome_id_order(genome_names)

    # Compute lengths and n_contigs
    genome_lengths = {g: 0 for g in sorted_names}
    genome_ncontigs = {g: 0 for g in sorted_names}
    for contig in contig_order:
        genome = contig_genome_map.get(contig, "__unbinned__")
        genome_lengths[genome] += contig_lengths[contig]
        genome_ncontigs[genome] += 1

    df = pd.DataFrame(
        {
            "name": sorted_names,
            "length": [genome_lengths[g] for g in sorted_names],
            "n_contigs": [genome_ncontigs[g] for g in sorted_names],
        }
    )

    # The genomes table only decorates the derived rows; say so out loud
    # whenever the two tables cover different sets of genomes.
    if extra_genomes_df is not None:
        _report_genome_table_divergence(sorted_names, extra_genomes_df)

    # Merge extra columns from the genomes file if provided
    if extra_genomes_df is not None and len(extra_genomes_df.columns) > 1:
        extra_genomes_df = _protect_reserved_columns(
            extra_genomes_df,
            reserved=_GENOMES_RESERVED_COLUMNS,
            keys=frozenset({"genome_id"}),
            context="Genomes table",
        )
        extra_cols = [c for c in extra_genomes_df.columns if c != "genome_id"]
        extra = extra_genomes_df.rename(columns={"genome_id": "name"})[
            ["name"] + extra_cols
        ].copy()
        # De-duplicate on the merge key so a repeated genome name cannot explode
        # the genomes table into duplicate rows (last row wins).
        extra = extra.drop_duplicates(subset="name", keep="last")
        for col in extra_cols:
            if (
                pd.api.types.is_string_dtype(extra[col].dtype)
                or extra[col].dtype == object
            ):
                with contextlib.suppress(ValueError, TypeError):
                    extra[col] = pd.to_numeric(extra[col], errors="raise")
        df = df.merge(extra, on="name", how="left")

    # Ensure column order: name, length, n_contigs, then extras
    core = ["name", "length", "n_contigs"]
    extras = [c for c in df.columns if c not in core]
    df = df[core + extras]
    return df.astype({k: v for k, v in GENOMES_DTYPE.items() if k in df.columns})


def _build_contigs_table(
    contig_order: list[str],
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    genomes: pd.DataFrame,
    contig_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the ``/contigs`` DataFrame with genome IDs and optional extra columns.

    Args:
        contig_order: Ordered list of contig names.
        contig_lengths: Mapping of contig name to length in base pairs.
        contig_genome_map: Mapping of contig name to genome name.
        genomes: The already-built ``/genomes`` DataFrame; used to look up
            each genome's integer index (``genome_id``).
        contig_df: Optional DataFrame parsed from a contigs-metadata TSV
            file.  Extra columns beyond ``contig_name`` and ``genome_id``
            are left-merged into the output by contig name; columns colliding
            with metacooler's own (``name``, ``length``) are renamed with a
            ``user_`` prefix.

    Returns:
        DataFrame with columns ``name``, ``length``, ``genome_id``, and
        any extra columns from *contig_df*.  Row order matches
        *contig_order*.
    """
    genome_name_to_id = {name: i for i, name in enumerate(genomes["name"])}

    def _get_genome_id(c: str) -> int:
        gname = contig_genome_map.get(c, "__unbinned__")
        gid = genome_name_to_id.get(gname)
        if gid is None:
            gid = genome_name_to_id.get("__unbinned__")
        if gid is None:
            raise ValueError(
                f"Contig {c!r} has no genome assignment and __unbinned__ is absent"
            )
        return gid

    df = pd.DataFrame(
        {
            "name": contig_order,
            "length": [contig_lengths[c] for c in contig_order],
            "genome_id": [_get_genome_id(c) for c in contig_order],
        }
    )

    # Merge additional columns from the input contig table if provided
    if contig_df is not None and len(contig_df.columns) > 2:
        contig_df = _protect_reserved_columns(
            contig_df,
            reserved=_CONTIGS_RESERVED_COLUMNS,
            keys=frozenset({"contig_name", "genome_id"}),
            context="Contigs table",
        )
        # Identify extra columns beyond contig_name and genome_id
        extra_cols = [
            c for c in contig_df.columns if c not in ("contig_name", "genome_id")
        ]
        if extra_cols:
            rename_map = {"contig_name": "name"}
            extra = contig_df.rename(columns=rename_map)[["name"] + extra_cols].copy()
            # De-duplicate on the merge key so a repeated contig name cannot
            # explode the contigs table into phantom rows (which would also
            # desync attrs['nchroms'] from len(/contigs)).  Last row wins,
            # matching the last-wins semantics of the contig→genome mapping.
            extra = extra.drop_duplicates(subset="name", keep="last")
            # Convert types for extra columns
            for col in extra_cols:
                if (
                    pd.api.types.is_string_dtype(extra[col].dtype)
                    or extra[col].dtype == object
                ):
                    # Try numeric conversion
                    with contextlib.suppress(ValueError, TypeError):
                        extra[col] = pd.to_numeric(extra[col], errors="raise")
            df = df.merge(extra, on="name", how="left")

    # Ensure column order: name, length, genome_id, then extras
    core_cols = ["name", "length", "genome_id"]
    extra_cols = [c for c in df.columns if c not in core_cols]
    df = df[core_cols + extra_cols]
    return df.astype({k: v for k, v in CONTIGS_DTYPE.items() if k in df.columns})


def _validate_bin_coords(df: pd.DataFrame) -> None:
    """Sanity-check bin ``start``/``end`` coordinates before writing.

    Catches malformed hand-made BED files that would otherwise be written
    verbatim (and, for out-of-range values, silently wrap when cast to the
    int32 storage dtype).

    Args:
        df: Bins DataFrame with integer ``start`` and ``end`` columns.

    Raises:
        ValueError: If any coordinate is negative, exceeds the int32 range,
            or if any bin has ``end < start``.
    """
    start = np.asarray(df["start"], dtype=np.int64)
    end = np.asarray(df["end"], dtype=np.int64)
    i32 = np.iinfo(np.int32)
    if (start < 0).any() or (end < 0).any():
        raise ValueError("Bin coordinates must be non-negative.")
    if (start > i32.max).any() or (end > i32.max).any():
        raise ValueError(
            "Bin coordinates exceed the int32 storage range "
            f"({i32.max} bp); contigs larger than ~2.1 Gbp are not supported."
        )
    inverted = end < start
    if inverted.any():
        example = df.loc[inverted].head(3).to_dict("records")
        raise ValueError(
            f"{int(inverted.sum())} bin(s) have end < start, e.g. {example}."
        )
    zero_width = end == start
    if zero_width.any():
        logger.warning("%d zero-width bin(s) (end == start).", int(zero_width.sum()))


def _read_bins_bed(path: str | Path) -> pd.DataFrame:
    """Read a headerless BED file of bins into a ``chrom``/``start``/``end`` frame.

    Columns beyond the first three (BED4/BED6 name, score, strand, …) are read
    and discarded.  Reading with a three-name ``names=`` list instead — the same
    trap documented in :func:`metacooler.io.contigs.parse_contigs_table` — makes
    pandas promote the surplus leading column into the index, so ``chrom`` picks
    up the *start* coordinate, ``start`` the end and ``end`` the fourth column.
    Nothing downstream can tell: contig lengths are derived from this very frame,
    so the shifted table is self-consistent and creation proceeds with contigs
    named after coordinates.

    Args:
        path: Path to a tab-separated BED file without a header row.

    Returns:
        DataFrame with a string ``chrom`` column and int64 ``start``/``end``.

    Raises:
        ValueError: If the file has fewer than three tab-separated columns.
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    if df.shape[1] < 3:
        raise ValueError(
            "Bins BED file requires at least 3 tab-separated columns "
            f"(chrom, start, end); got {df.shape[1]} in {path}."
        )
    if df.shape[1] > 3:
        logger.debug(
            "Bins BED file has %d columns; using the first three as chrom/start/end.",
            df.shape[1],
        )
    df = df.iloc[:, :3].copy()
    df.columns = ["chrom", "start", "end"]
    return df.astype({"chrom": str, "start": np.int64, "end": np.int64})


def _build_bins(
    bins: pd.DataFrame,
    contig_name_to_id: dict[str, int],
) -> pd.DataFrame:
    """Validate a bins DataFrame and remap its contig names to integer IDs.

    Args:
        bins: Bins DataFrame with a string ``chrom`` column plus ``start`` and
            ``end`` (see :func:`_read_bins_bed` for the BED reader).
        contig_name_to_id: Mapping from contig name to integer ID.

    Returns:
        DataFrame with columns ``chrom`` (int32 contig ID), ``start``
        (int32), ``end`` (int32).

    Raises:
        ValueError: If any contig name in the bins file is not found in
            *contig_name_to_id*, or if bin coordinates are malformed.
    """
    df = bins.copy()
    _validate_bin_coords(df)
    # Map chrom name -> contig id
    mapped = df["chrom"].map(contig_name_to_id)
    if mapped.isna().any():
        missing = df["chrom"][mapped.isna()].unique()
        raise ValueError(f"Bins reference unknown contigs: {missing}")
    df["chrom"] = mapped
    return df[["chrom", "start", "end"]].astype(
        {"chrom": np.int32, "start": np.int32, "end": np.int32}
    )


def _tile_bins(
    contig_order: list[str],
    contig_lengths: dict[str, int],
    contig_name_to_id: dict[str, int],
    binsize: int,
) -> pd.DataFrame:
    """Generate fixed-size bins tiling every contig.

    Each contig is divided into non-overlapping windows of width *binsize*;
    the last bin per contig may be shorter than *binsize*.

    Args:
        contig_order: Ordered list of contig names; bin rows are emitted in
            this order.
        contig_lengths: Mapping of contig name to length.
        contig_name_to_id: Mapping of contig name to integer ID (used as
            the ``chrom`` column value).
        binsize: Bin width in base pairs.

    Returns:
        DataFrame with columns ``chrom`` (int32 contig ID), ``start``
        (int32, 0-based), and ``end`` (int32, clamped to contig length).

    Raises:
        ValueError: If *binsize* is not positive, or if a contig is longer than
            the int32 bin-coordinate storage range.
    """
    if binsize <= 0:
        raise ValueError(f"binsize must be a positive integer, got {binsize!r}.")
    i32_max = int(np.iinfo(np.int32).max)
    parts = []
    for contig in contig_order:
        cid = contig_name_to_id[contig]
        length = int(contig_lengths[contig])
        # Build in int64 and range-check *before* narrowing: an int32 ``arange``
        # silently wraps to negative starts on a contig longer than ~2.1 Gbp,
        # and schema.validate() does not catch negative bin coordinates.
        if length > i32_max:
            raise ValueError(
                f"Contig {contig!r} is {length} bp, which exceeds the int32 bin "
                f"coordinate range ({i32_max} bp); contigs larger than ~2.1 Gbp "
                "are not supported."
            )
        starts = np.arange(0, length, binsize, dtype=np.int64)
        ends = np.minimum(starts + binsize, length)
        parts.append(pd.DataFrame({"chrom": cid, "start": starts, "end": ends}))
    df = pd.concat(parts, ignore_index=True)
    return df.astype({"chrom": np.int32, "start": np.int32, "end": np.int32})


def _scan_pairs_header(path: Path) -> tuple[int, bool]:
    """Scan the top of a pairs file for its header size and whether data follows.

    Reads lines sequentially until the first non-comment line is found.
    Handles plain text, ``.gz``/``.bgz`` (gzip), and ``.bz2`` compressed
    files.

    Args:
        path: Path to the pairs file.

    Returns:
        ``(n_header, has_data_rows)`` — the count of consecutive leading ``#``
        lines, and whether at least one non-empty data row follows them.  The
        second element lets the caller recognise a header-only file *before*
        handing the stream to pyarrow, instead of pattern-matching pyarrow's
        error message (which would also swallow genuine parse failures).
    """
    import bz2 as _bz2
    import gzip as _gzip

    if path.suffix in (".gz", ".bgz"):
        open_fn = _gzip.open
    elif path.suffix == ".bz2":
        open_fn = _bz2.open
    else:
        open_fn = open
    n = 0
    has_data = False
    with open_fn(str(path), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                n += 1
                continue
            if not line.strip():
                # Tolerate stray blank lines between the header and the body.
                continue
            has_data = True
            break
    return n, has_data


def _stream_pairs(
    pairs_path: str | Path,
    c1: int = 2,
    p1: int = 3,
    c2: int = 4,
    p2: int = 5,
    chunksize: int = 1_000_000,
) -> Iterable[pd.DataFrame]:
    """Stream a pairs file in chunks via pyarrow for fast native gzip decoding.

    pyarrow uses multi-threaded C++ parsing and native gzip decompression,
    which is significantly faster than pandas/gzip for large compressed files.
    Comment lines (``#``) at the top of the file are counted once with a fast
    sequential scan and then skipped by pyarrow's ``skip_rows`` option.
    """
    import pyarrow as pa
    import pyarrow.csv as pa_csv

    path = Path(pairs_path)

    # Count header comment lines; pyarrow has no native comment-skip support.
    # The same scan tells us whether any data row follows at all.
    n_header, has_data = _scan_pairs_header(path)
    if not has_data:
        # A header-only pairs file (no data rows) is a legitimate input — a
        # sample with no surviving contacts. Detect it structurally here rather
        # than by pattern-matching pyarrow's error message, which would also
        # swallow a genuine parse failure and write a silently empty map.
        logger.warning("Pairs file %s has no data rows; writing empty map.", path)
        return

    # Column indices are 1-based; convert to 0-based.
    # pyarrow auto-generates column names "f0", "f1", ... when
    # autogenerate_column_names=True, so we address columns by those names.
    pa_col_names = [f"f{c1 - 1}", f"f{p1 - 1}", f"f{c2 - 1}", f"f{p2 - 1}"]
    final_col_names = ["chrom1", "pos1", "chrom2", "pos2"]

    # pyarrow buffers one whole block of raw text at a time, so this controls
    # peak memory rather than the caller's chunk granularity (``_stream_binned_pairs``
    # coalesces batches back up to *chunksize*). Aim for roughly chunksize rows
    # at ~120 bytes/row, but cap it so a large chunksize stays affordable.
    block_size = min(max(chunksize * 120, 8 * 1024 * 1024), _PAIRS_READ_BLOCK_MAX_BYTES)

    read_opts = pa_csv.ReadOptions(
        skip_rows=n_header,
        autogenerate_column_names=True,
        block_size=block_size,
    )
    parse_opts = pa_csv.ParseOptions(delimiter="\t")
    convert_opts = pa_csv.ConvertOptions(
        include_columns=pa_col_names,
        column_types={
            f"f{c1 - 1}": pa.string(),
            f"f{p1 - 1}": pa.int64(),
            f"f{c2 - 1}": pa.string(),
            f"f{p2 - 1}": pa.int64(),
        },
    )

    # pa.input_stream auto-detects a codec from the file extension, but its codec
    # table has no ``.bgz`` entry: a bgzipped file would be handed to the CSV
    # parser as raw binary, yielding zero recognisable contigs and a structurally
    # valid but empty map. bgzf is gzip-compatible, so name the codec explicitly.
    if path.suffix == ".bgz":
        logger.info("Reading %s as gzip (bgzf is gzip-compatible).", path)
        src = pa.input_stream(str(path), compression="gzip")
    else:
        src = pa.input_stream(str(path))
    with contextlib.closing(src):
        reader = pa_csv.open_csv(
            src,
            read_options=read_opts,
            parse_options=parse_opts,
            convert_options=convert_opts,
        )
        for batch in reader:
            df = batch.to_pandas()
            df.columns = final_col_names
            yield df


@dataclasses.dataclass
class _PairsDropStats:
    """Running tally of pairs discarded while binning, for an end-of-run report.

    Silently dropping pairs is the failure mode that hurts most in practice: a
    contig-name-format mismatch between the pairs file and the contig table
    discards 100 % of contacts, and without this tally the run still reports
    ``NNZ=0`` as if nothing were wrong.
    """

    total: int = 0
    unknown_contig: int = 0
    null_position: int = 0
    unknown_names: set[str] = dataclasses.field(default_factory=set)

    def note_unknown_names(self, names: Iterable[str]) -> None:
        """Remember up to a handful of example unrecognised contig names."""
        for name in names:
            if len(self.unknown_names) >= _MAX_REPORTED_UNKNOWN_CONTIGS:
                return
            self.unknown_names.add(str(name))

    def report(self) -> None:
        """Log a WARNING summary of everything that was dropped (if anything)."""
        dropped = self.unknown_contig + self.null_position
        if dropped == 0:
            return
        examples = sorted(self.unknown_names)
        if dropped == self.total:
            logger.warning(
                "ALL %d pair(s) were discarded (%d on unrecognised contigs, %d with "
                "a missing position) — the resulting map is empty. Do the contig "
                "names in the pairs file match the FASTA / contig table? "
                "Unrecognised names seen: %s",
                self.total,
                self.unknown_contig,
                self.null_position,
                examples,
            )
            return
        logger.warning(
            "Discarded %d of %d pair(s) (%.2f%%): %d on unrecognised contigs "
            "(e.g. %s), %d with a missing position.",
            dropped,
            self.total,
            100.0 * dropped / max(self.total, 1),
            self.unknown_contig,
            examples,
            self.null_position,
        )


def _positions_to_int64(values: np.ndarray, column: str) -> np.ndarray:
    """Return *values* as int64, refusing silently-lossy inputs.

    An empty ``pos`` field makes pyarrow produce a null, which pandas surfaces as
    a float64 ``NaN``; casting that straight to int64 yields ``INT64_MIN`` and a
    wildly out-of-range bin id.  Null rows are filtered out before this call, so
    anything non-integral reaching it is a genuinely malformed file.

    Args:
        values: Position array taken from the pairs chunk.
        column: Column name, used in the error message.

    Returns:
        int64 view/copy of *values*.

    Raises:
        ValueError: If *values* is neither integer nor an exactly-integral float.
    """
    if values.dtype.kind in ("i", "u"):
        return values.astype(np.int64, copy=False)
    if values.dtype.kind == "f":
        if values.size and not np.all(np.isfinite(values) & (np.mod(values, 1) == 0)):
            raise ValueError(f"Pairs column {column!r} contains non-integer positions.")
        return values.astype(np.int64)
    raise ValueError(
        f"Pairs column {column!r} has non-numeric dtype {values.dtype}; "
        "check the --pairs column indices."
    )


def _finite_position_mask(col: pd.Series) -> np.ndarray:
    """Boolean mask of rows whose position value is present (not null/NaN)."""
    arr = col.to_numpy()
    if arr.dtype.kind in ("i", "u"):
        return np.ones(len(arr), dtype=bool)
    return col.notna().to_numpy()


def _bin_pairs_chunk(
    chunk: pd.DataFrame,
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_starts: np.ndarray | None,
    stats: _PairsDropStats | None = None,
) -> pd.DataFrame | None:
    """Bin a single pairs chunk into aggregated ``(bin1_id, bin2_id, count)`` pixels.

    Supports both fixed-binsize mode (arithmetic) and fragment-level mode
    (binary search via ``np.searchsorted``).  Returns ``None`` when the chunk
    contains no valid pairs after contig filtering.

    Args:
        chunk: Raw pairs DataFrame with columns ``chrom1``, ``pos1``,
            ``chrom2``, ``pos2``.
        contig_name_to_id: Mapping from contig name to integer ID.
        chrom_offset: Cumulative bin counts per contig (length n_contigs+1).
        binsize: Fixed bin size, or ``None`` for fragment-level bins.
        bin_starts: Flattened array of bin start positions; required when
            ``binsize`` is ``None``.
        stats: Optional accumulator updated with the number of pairs seen and
            dropped, so the caller can report them once at end of stream.

    Returns:
        Aggregated pixel DataFrame or ``None`` if the chunk is empty after
        filtering.
    """
    # Fast contig-name → integer-ID encoding.  Unrecognised contigs map to
    # NaN and are coerced to -1 so they can be dropped in one vectorised pass.
    c1_codes = chunk["chrom1"].map(contig_name_to_id).fillna(-1).astype(np.int32).values
    c2_codes = chunk["chrom2"].map(contig_name_to_id).fillna(-1).astype(np.int32).values

    # Drop pairs referencing unknown contigs, and pairs whose position field was
    # empty (pyarrow null → float NaN → INT64_MIN bin id if left alone).
    contig_ok = (c1_codes >= 0) & (c2_codes >= 0)
    pos_ok = _finite_position_mask(chunk["pos1"]) & _finite_position_mask(chunk["pos2"])
    valid = contig_ok & pos_ok
    if stats is not None:
        stats.total += len(chunk)
        n_unknown = int((~contig_ok).sum())
        stats.unknown_contig += n_unknown
        stats.null_position += int((contig_ok & ~pos_ok).sum())
        if n_unknown and len(stats.unknown_names) < _MAX_REPORTED_UNKNOWN_CONTIGS:
            bad = ~contig_ok
            stats.note_unknown_names(
                pd.unique(
                    pd.concat(
                        [
                            chunk["chrom1"][bad & (c1_codes < 0)],
                            chunk["chrom2"][bad & (c2_codes < 0)],
                        ]
                    )
                )
            )
    if not valid.any():
        return None

    c1 = c1_codes[valid].astype(np.int32)
    c2 = c2_codes[valid].astype(np.int32)
    # Convert 1-based pairs positions to 0-based, clamp at 0
    raw_p1 = _positions_to_int64(chunk["pos1"].to_numpy()[valid], "pos1") - 1
    raw_p2 = _positions_to_int64(chunk["pos2"].to_numpy()[valid], "pos2") - 1
    neg_mask = (raw_p1 < 0) | (raw_p2 < 0)
    if neg_mask.any():
        logger.warning(
            "%d pairs clamped: negative 0-based positions encountered in chunk",
            int(neg_mask.sum()),
        )
    p1 = np.maximum(raw_p1, 0)
    p2 = np.maximum(raw_p2, 0)

    if binsize is not None:
        # Fixed-size bins: clamp so an out-of-range position never overflows
        # into the next contig's bins.
        bin1 = np.minimum(
            chrom_offset[c1] + p1 // binsize,
            chrom_offset[c1 + 1] - 1,
        )
        bin2 = np.minimum(
            chrom_offset[c2] + p2 // binsize,
            chrom_offset[c2 + 1] - 1,
        )
    else:
        # Fragment-level bins: binary search within each contig's bin starts
        assert bin_starts is not None, "bin_starts required for fragment-level binning"
        bin1 = np.empty(len(c1), dtype=np.int64)
        bin2 = np.empty(len(c2), dtype=np.int64)
        for cid in np.unique(c1):
            mask = c1 == cid
            lo, hi = int(chrom_offset[cid]), int(chrom_offset[cid + 1])
            bin1[mask] = np.maximum(
                lo,
                lo + np.searchsorted(bin_starts[lo:hi], p1[mask], side="right") - 1,
            )
        for cid in np.unique(c2):
            mask = c2 == cid
            lo, hi = int(chrom_offset[cid]), int(chrom_offset[cid + 1])
            bin2[mask] = np.maximum(
                lo,
                lo + np.searchsorted(bin_starts[lo:hi], p2[mask], side="right") - 1,
            )

    # Enforce upper triangle
    swap = bin1 > bin2
    bin1[swap], bin2[swap] = bin2[swap], bin1[swap]

    df = pd.DataFrame(
        {
            "bin1_id": bin1,
            "bin2_id": bin2,
            "count": np.ones(len(bin1), dtype=np.int32),
        }
    )
    return df.groupby(["bin1_id", "bin2_id"], sort=False, as_index=False)["count"].sum()


def _total_time_str(start_time: float) -> str:
    """Format elapsed time since *start_time* as a human-readable string.

    Args:
        start_time: Monotonic counter value captured at the beginning of the
            operation (e.g. ``time.perf_counter()``).

    Returns:
        Formatted string such as ``"3.2s"``, ``"1m:24s"``, or
        ``"1h:3m:47s"`` depending on the elapsed duration.
    """
    secs = time.perf_counter() - start_time
    if secs >= 3600:
        return f"{int(secs // 3600)}h:{int((secs % 3600) // 60)}m:{int(secs % 60)}s"
    elif secs >= 60:
        return f"{int(secs // 60)}m:{int(secs % 60)}s"
    else:
        return f"{secs:.1f}s"


def _parse_contig_map(
    contig_genome_map: dict[str, str] | str | Path,
    contig_table_format: str,
) -> tuple[dict[str, str], pd.DataFrame | None]:
    """Parse and normalize a contig-to-genome mapping from a dict or file path.

    Args:
        contig_genome_map: Either a pre-built ``{contig: genome}`` mapping or a
            path to a contig table file.
        contig_table_format: Format string passed to ``parse_contigs_table`` when
            *contig_genome_map* is a file path.

    Returns:
        Tuple of ``(cg_map, contig_df)`` where ``contig_df`` is ``None`` when
        *contig_genome_map* is already a dict.
    """
    contig_df: pd.DataFrame | None = None
    if isinstance(contig_genome_map, str | Path):
        cg_map, contig_df = parse_contigs_table(
            contig_genome_map, fmt=contig_table_format
        )
        logger.info(
            f"Parsed contig-genome map for {len(cg_map)} contigs from {contig_genome_map}"
        )
    else:
        cg_map = dict(contig_genome_map)
        logger.info(f"Parsed contig-genome map for {len(cg_map)} contigs from dict")
    cg_map = {k: _normalize_genome_name(v) for k, v in cg_map.items()}
    return cg_map, contig_df


def _load_extra_genomes_table(
    genomes_table: str | Path | None,
) -> pd.DataFrame | None:
    """Load an optional genomes metadata table from disk.

    Args:
        genomes_table: Path to the TSV file, or ``None``.

    Returns:
        Parsed DataFrame or ``None`` if *genomes_table* is ``None``.
    """
    if genomes_table is None:
        return None
    df = parse_genomes_table(genomes_table)
    logger.info(f"Parsed genomes table for {len(df)} genomes from {genomes_table}")
    return df


def _compute_all_offsets(
    bins_df: pd.DataFrame,
    contigs_df: pd.DataFrame,
    genomes_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute chrom_offset and genome_offset; insert ``n_bins`` into genomes_df.

    Args:
        bins_df: Bins DataFrame with an integer ``chrom`` column (contig IDs).
        contigs_df: Contigs DataFrame with an integer ``genome_id`` column.
        genomes_df: Genomes DataFrame; ``n_bins`` is inserted at column position 3.

    Returns:
        ``(chrom_offset, genome_offset)`` — int64 arrays of length ``n_contigs+1``
        and ``n_genomes+1`` respectively.
    """
    n_contigs = len(contigs_df)
    n_genomes = len(genomes_df)

    chrom_counts = np.bincount(
        bins_df["chrom"].values.astype(np.int64), minlength=n_contigs
    )
    chrom_offset = np.zeros(n_contigs + 1, dtype=np.int64)
    chrom_offset[1:] = np.cumsum(chrom_counts)

    genome_ids = contigs_df["genome_id"].values
    genome_counts = np.bincount(genome_ids.astype(np.int64), minlength=n_genomes)
    genome_offset = np.zeros(n_genomes + 1, dtype=np.int64)
    genome_offset[1:] = np.cumsum(genome_counts)

    genome_bin_counts = np.bincount(
        contigs_df["genome_id"].values.astype(np.int64),
        weights=chrom_counts.astype(np.float64),
        minlength=n_genomes,
    ).astype(np.int64)
    genomes_df.insert(3, "n_bins", genome_bin_counts)

    return chrom_offset, genome_offset


def _sort_and_dedup_pixels(
    pixel_iter: Iterable[pd.DataFrame],
    write_chunk: int = 5_000_000,
) -> Iterable[pd.DataFrame]:
    """Collect, globally sort, deduplicate, and yield pixel DataFrames.

    Accumulates all pixel DataFrames from *pixel_iter* into RAM, performs a
    single global groupby-sum + sort, then yields fixed-size sub-chunks ready
    for HDF5 appending.

    Args:
        pixel_iter: Iterable of unsorted pixel DataFrames with columns
            ``bin1_id``, ``bin2_id``, ``count``.
        write_chunk: Maximum rows per yielded chunk.

    Yields:
        Sorted, deduplicated pixel DataFrames of at most *write_chunk* rows.
    """
    frames = list(pixel_iter)
    if not frames:
        return
    logger.info("Sorting pixels in memory...")
    _t0 = time.time()
    all_pixels = pd.concat(frames, ignore_index=True)
    del frames
    all_pixels = all_pixels.groupby(["bin1_id", "bin2_id"], sort=True, as_index=False)[
        "count"
    ].sum()
    logger.info(f"Sorting pixels in memory... DONE ({time.time() - _t0:.1f}s elapsed)")
    logger.info(f"Total pixels after global dedup: {len(all_pixels):,}")
    for lo in range(0, len(all_pixels), write_chunk):
        yield all_pixels.iloc[lo : lo + write_chunk]


def _merge_pixel_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate binned pixel frames, re-summing bin pairs seen in more than one.

    Args:
        frames: Non-empty list of aggregated pixel DataFrames.

    Returns:
        A single aggregated pixel DataFrame.
    """
    if len(frames) == 1:
        return frames[0]
    return (
        pd.concat(frames, ignore_index=True)
        .groupby(["bin1_id", "bin2_id"], sort=False, as_index=False)["count"]
        .sum()
    )


def _buffer_pixels(
    pixel_iter: Iterable[pd.DataFrame],
    max_rows: int,
    fold_rows: int | None = None,
) -> tuple[list[pd.DataFrame], Iterable[pd.DataFrame] | None]:
    """Buffer pixel frames in RAM, folding as it goes, until *max_rows* is exceeded.

    Lets a caller pick a sort strategy from the pixel volume it actually sees
    rather than from a proxy such as the input file size.  This matters because
    the external-sort path calls ``cooler.create`` once per chunk, and each such
    call carries a fixed cost proportional to the number of contigs — several
    seconds for a metagenome — regardless of how many pixels the chunk holds.

    The buffer is periodically **folded**: the pending frames are aggregated
    into a single running frame, so the row count driving the decision is the
    deduplicated pixel table rather than the sum of per-chunk row counts.
    Without that, the two diverge badly.  Bin pairs recur across chunks, and
    Hi-C pixel counts saturate on a fixed bin grid, so a deep library can stream
    hundreds of millions of raw rows that collapse to single-digit millions of
    distinct pixels — enough to trip an unfolded budget and pay for an external
    sort that the data never needed.

    Args:
        pixel_iter: Iterable of unsorted pixel DataFrames.
        max_rows: Budget, in *deduplicated* pixels.  Buffering stops once the
            folded frame exceeds it.
        fold_rows: Rows that may accumulate unfolded before a fold is forced.
            Defaults to ``max_rows // _PIXEL_FOLD_FRACTION``.

    Returns:
        ``(frames, rest)``.  *frames* holds at most one (folded) DataFrame.
        *rest* is ``None`` when the whole stream fit within the budget;
        otherwise it is the unconsumed remainder of *pixel_iter*.
    """
    if fold_rows is None:
        fold_rows = max(max_rows // _PIXEL_FOLD_FRACTION, 1)
    folded: pd.DataFrame | None = None
    pending: list[pd.DataFrame] = []
    pending_rows = 0

    def _fold() -> None:
        nonlocal folded, pending, pending_rows
        if not pending:
            return
        group = pending if folded is None else [folded, *pending]
        # Drop our references before the merge so the concat inside it is not
        # competing with a second set of live frames.
        folded = None
        pending = []
        pending_rows = 0
        folded = _merge_pixel_frames(group)

    for df in pixel_iter:
        pending.append(df)
        pending_rows += len(df)
        if pending_rows >= fold_rows:
            _fold()
            assert folded is not None  # a fold with pending always sets it
            if len(folded) > max_rows:
                # Hand the accumulator back as ONE frame rather than slicing it:
                # in the external path every frame becomes its own
                # ``cooler.create`` call, and that call's cost is dominated by a
                # fixed per-call overhead that scales with contig count, not
                # with pixels. One large chunk is far cheaper than twenty small.
                return [folded], pixel_iter
    _fold()
    return ([] if folded is None else [folded]), None


def _chain_pixel_frames(
    frames: list[pd.DataFrame],
    rest: Iterable[pd.DataFrame],
) -> Iterable[pd.DataFrame]:
    """Yield the buffered *frames*, dropping each reference, then *rest*.

    Popping as it goes lets the consumed frames be freed while the rest of the
    stream is still being written, so the buffer does not stay resident for the
    whole external sort.
    """
    while frames:
        yield frames.pop(0)
    yield from rest


def _ensure_nonempty_pixel_stream(
    pixel_iter: Iterable[pd.DataFrame],
) -> Iterable[pd.DataFrame]:
    """Yield from *pixel_iter*, emitting one empty frame if it is empty.

    ``cooler.create_cooler`` writes a malformed pixel table (a dense block of
    phantom rows with ``nnz=0``) when handed a generator that yields *nothing*,
    but handles an empty DataFrame correctly. This guard turns a zero-pixel
    input (e.g. a header-only pairs file) into a valid empty ``.metacool``.
    """
    yielded = False
    for df in pixel_iter:
        yielded = True
        yield df
    if not yielded:
        yield pd.DataFrame(
            {
                "bin1_id": np.array([], dtype=np.int64),
                "bin2_id": np.array([], dtype=np.int64),
                "count": np.array([], dtype=np.int32),
            }
        )


def _stream_binned_pairs(
    pairs: str | Path,
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_starts: np.ndarray | None,
    chunksize: int,
    pairs_kwargs: dict,
) -> Iterable[pd.DataFrame]:
    """Stream a pairs file and yield binned pixel DataFrames with progress.

    Args:
        pairs: Path to the pairs file (plain text, gzipped, or bgzipped).
        contig_name_to_id: Mapping from contig name to integer ID.
        chrom_offset: Cumulative bin counts per contig (length n_contigs+1).
        binsize: Fixed bin size, or ``None`` for fragment-level bins.
        bin_starts: Flattened bin start positions; required when *binsize*
            is ``None``.
        chunksize: Rows per streaming batch.
        pairs_kwargs: Extra keyword arguments forwarded to ``_stream_pairs``.

    Yields:
        Aggregated pixel DataFrames (``bin1_id``, ``bin2_id``, ``count``).
    """
    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

    con = Console(stderr=True)
    with Progress(
        SpinnerColumn(),
        TextColumn("[cyan]Binning pairs"),
        TextColumn("{task.completed:,.0f} pairs processed"),
        TimeElapsedColumn(),
        console=con,
        transient=True,
        disable=not con.is_terminal,
    ) as progress:
        task = progress.add_task("", total=None)
        pending: list[pd.DataFrame] = []
        pending_pairs = 0
        stats = _PairsDropStats()
        for chunk in _stream_pairs(pairs, chunksize=chunksize, **pairs_kwargs):
            progress.update(task, advance=len(chunk))
            result = _bin_pairs_chunk(
                chunk, contig_name_to_id, chrom_offset, binsize, bin_starts, stats
            )
            if result is not None and not result.empty:
                pending.append(result)
            pending_pairs += len(chunk)
            # The reader may hand back smaller batches than requested (its block
            # size is capped); coalesce them so a consumer that pays a fixed cost
            # per chunk sees the chunk count the caller asked for.
            if pending_pairs >= chunksize:
                if pending:
                    yield _merge_pixel_frames(pending)
                pending = []
                pending_pairs = 0
        if pending:
            yield _merge_pixel_frames(pending)
    stats.report()


def _stream_remapped_pixels(
    cool_path_str: str,
    n_pixels: int,
    old_to_new_bin: np.ndarray,
    nproc: int,
    count_dtype: np.dtype,
    chunksize: int = 1_000_000,
) -> Iterable[pd.DataFrame]:
    """Stream pixels from a cool file with bin-ID remapping and rich progress.

    Args:
        cool_path_str: Filesystem path to the source cool file.
        n_pixels: Total pixel count (used for progress bar).
        old_to_new_bin: Array mapping old bin indices to new bin indices.
        nproc: Accepted for backwards compatibility but currently ignored — see
            the warning logged below.
        count_dtype: dtype for the output ``count`` column, as resolved by
            :func:`metacooler.io._utils.resolve_count_dtype`.
        chunksize: Pixels per read batch.

    Yields:
        Remapped pixel DataFrames with columns ``bin1_id``, ``bin2_id``,
        ``count``.
    """
    from rich.console import Console
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TaskProgressColumn,
        TextColumn,
        TimeRemainingColumn,
    )

    if nproc > 1:
        # The worker pool used to exist, but libhdf5 is not guaranteed to be
        # built thread-safe, so every read had to be taken under a single lock —
        # serialising the only expensive step. The threads bought nothing and
        # each leaked an open Cooler handle, so the read is now plainly serial.
        logger.warning(
            "nproc=%d has no effect on pixel remapping: the HDF5 reads are "
            "serialised, so the work is done in a single thread.",
            nproc,
        )

    clr = cooler.Cooler(cool_path_str)

    def _remap_chunk(lo: int, hi: int) -> pd.DataFrame:
        chunk = clr.pixels()[lo:hi].copy()
        chunk["bin1_id"] = old_to_new_bin[chunk["bin1_id"].values.astype(np.int64)]
        chunk["bin2_id"] = old_to_new_bin[chunk["bin2_id"].values.astype(np.int64)]
        # Drop pixels that touch a dropped contig (sentinel -1).
        chunk = chunk[(chunk["bin1_id"] >= 0) & (chunk["bin2_id"] >= 0)].copy()
        chunk["count"] = chunk["count"].astype(count_dtype)
        swap = chunk["bin1_id"] > chunk["bin2_id"]
        chunk.loc[swap, ["bin1_id", "bin2_id"]] = chunk.loc[
            swap, ["bin2_id", "bin1_id"]
        ].values
        return chunk[["bin1_id", "bin2_id", "count"]]

    ranges = [(i, min(i + chunksize, n_pixels)) for i in range(0, n_pixels, chunksize)]
    con = Console(stderr=True)
    with Progress(
        TextColumn("[cyan]Remapping pixels"),
        BarColumn(),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TextColumn("chunks"),
        TimeRemainingColumn(),
        console=con,
        transient=True,
        disable=not con.is_terminal,
    ) as progress:
        task = progress.add_task("", total=len(ranges))
        for lo, hi in ranges:
            yield _remap_chunk(lo, hi)
            progress.update(task, advance=1)


@contextlib.contextmanager
def _staged_output(final: Path) -> Iterator[Path]:
    """Build the output at a sibling temp path, moving it into place on success.

    Writing straight to the final path meant a partial or invalid file could sit
    at the name the user asked for, and the old ``except BaseException: unlink()``
    guard fired on ``KeyboardInterrupt`` — Ctrl-C during the final ``validate()``
    of a multi-hour run deleted the finished file. Staging solves both: nothing
    ever appears at the final name until validation has passed, and a failure
    leaves the artifact behind (under the temp name) so it can be diagnosed.

    Args:
        final: Path the caller wants the finished file to have.

    Yields:
        The temporary path to write to.
    """
    tmp = final.with_name(final.name + ".partial")
    check_output_path(tmp, force=True)
    try:
        yield tmp
    except Exception:
        logger.error(
            "Creation failed. The partial file was left at %s for inspection; "
            "%s was not written.",
            tmp,
            final,
        )
        raise
    os.replace(tmp, final)


def _remap_bins(
    bins: pd.DataFrame,
    contig_order: list[str],
) -> tuple[pd.DataFrame, np.ndarray]:
    """Remap a bins DataFrame to new contig IDs and return the old→new bin mapping.

    Bins are sorted by new chrom ID (stable sort) so that
    ``cooler.create_cooler`` receives them in the canonical order.

    Args:
        bins: Original bins DataFrame (``chrom`` column contains contig names).
        contig_order: List of contig names defining the new ordering.

    Returns:
        ``(bins_remapped, old_to_new_bin)`` where *bins_remapped* contains only
        bins whose contig appears in *contig_order*, sorted by new chrom ID, and
        *old_to_new_bin* is an array of length ``len(bins)`` mapping each original
        bin index to its new position.  Bins for contigs absent from
        *contig_order* receive a sentinel value of ``-1``; callers that stream
        pixels must filter those rows before writing.
    """
    n_original = len(bins)
    new_contig_name_to_id = {name: idx for idx, name in enumerate(contig_order)}
    bins_work = bins.copy()
    bins_work["chrom"] = bins_work["chrom"].map(new_contig_name_to_id)
    # Bins whose contig is absent from contig_order map to NaN; drop them and
    # record their original positions so the pixel remapper can use -1 as a
    # sentinel to discard pixels touching dropped contigs.
    kept_mask = bins_work["chrom"].notna().values
    original_kept_indices = np.where(kept_mask)[0]
    bins_remapped = bins_work.loc[kept_mask, ["chrom", "start", "end"]].copy()
    bins_remapped = bins_remapped.astype(
        {"chrom": np.int32, "start": np.int32, "end": np.int32}
    )
    sort_idx = np.argsort(bins_remapped["chrom"].values, kind="stable")
    bins_remapped = bins_remapped.iloc[sort_idx].reset_index(drop=True)
    n_kept = len(bins_remapped)
    # Build full-size mapping: dropped bins → -1, kept bins → new sorted index.
    old_to_new_bin = np.full(n_original, -1, dtype=np.int64)
    old_to_new_bin[original_kept_indices[sort_idx]] = np.arange(n_kept, dtype=np.int64)
    return bins_remapped, old_to_new_bin


def from_pairs(
    uri: str | Path,
    bins: str | Path | pd.DataFrame | None = None,
    pairs: str | Path | None = None,
    contig_genome_map: dict[str, str] | str | Path | None = None,
    contig_table_format: str = "contigs",
    fasta: str | Path | None = None,
    assembly_name: str | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "fasta", "alpha"
    ] = "genome-then-length",
    drop_unassigned: bool = False,
    binsize: int | None = None,
    genomes_table: str | Path | None = None,
    pairs_kwargs: dict | None = None,
    metadata: dict | None = None,
    chunksize: int = 50_000_000,
    sort_in_memory: bool | None = None,
) -> None:
    """Create a ``.metacool`` file from pairs, bins, and contig assignments.

    Args:
        uri: Output file path (``.metacool`` extension recommended).
        bins: Path to a BED file of bins (without header), or a pre-built DataFrame
            with columns ``chrom, start, end``.  If ``None``, *binsize* must be given
            and bins will be tiled from the fasta.
        pairs: Path to a pairs file (plain text, gzipped, or bgzipped).
        contig_genome_map: Either a ``dict[contig_name → genome_name]`` or a path
            to a contig table file (see *contig_table_format*).  Format-specific
            colspec strings are also accepted (see ``io/contigs.py`` for details).
        contig_table_format: Column spec when *contig_genome_map* is a file path.
            Built-in values are ``"contigs"``, ``"contig_genome"``, and ``"metator"``,
            but any valid colspec string is accepted.
        fasta: Assembly fasta (required if *bins* is ``None``).
        assembly_name: Name for this assembly / metagenome.
        contig_ordering: Canonical ordering strategy.
        drop_unassigned: If ``True``, drop contigs not present in the contig
            table instead of assigning them to ``__unbinned__``.
        binsize: Fixed bin size (used only when *bins* is ``None``).
        genomes_table: Optional path to a genomes metadata TSV **with a header
            row**. The first column must be named ``genome_name``; any additional
            columns are stored in the ``/genomes`` table.
        pairs_kwargs: Extra kwargs passed to the pairs reader (column
            indices, etc.).
        metadata: Additional metadata dict stored in file attrs.
        chunksize: Number of *pairs* to process per chunk.  This is chunk
            granularity, not read memory: pyarrow's read block is capped
            separately at ``_PAIRS_READ_BLOCK_MAX_BYTES``, and the binned
            batches are coalesced back up to *chunksize* afterwards.  Under
            ``sort_in_memory=False`` every chunk costs one ``cooler.create``
            call, whose runtime is dominated by a fixed overhead that scales
            with contig count rather than pixels (seconds per call on a
            metagenome), so small chunks are disproportionately expensive.
            The cost of a large value is the pixels held while coalescing one
            chunk — roughly ``chunksize`` rows at 24 bytes in the worst case.
        sort_in_memory: If ``True``, collect all pixels in RAM and sort before
            writing (fast for small files).  If ``False``, pass the unsorted
            stream to cooler's external sort (scalable for large files).  If
            ``None`` (default), buffer the binned pixels and use ``True`` unless
            they exceed 100 M rows, in which case fall back to ``False`` with a
            logged warning.
    """
    uri = Path(uri)
    if pairs_kwargs is None:
        pairs_kwargs = {}
    if pairs is None:
        raise ValueError("pairs file path required")
    if not isinstance(contig_genome_map, dict | str | Path):
        raise TypeError("contig_genome_map must be a dict or a path")
    if bins is not None and binsize is not None:
        # Passing both silently binned every contact arithmetically against a bin
        # table that was never tiled at *binsize*, scattering contacts across
        # unrelated fragments while still passing validate().
        raise ValueError(
            "binsize and bins are mutually exclusive: pass binsize to tile "
            "fixed-size bins, or bins to use a pre-defined bin table."
        )
    timer = time.perf_counter()

    # --- Parse contig map and optional genomes metadata table ---
    cg_map, contig_df = _parse_contig_map(contig_genome_map, contig_table_format)
    extra_genomes_df = _load_extra_genomes_table(genomes_table)

    # --- Resolve contig lengths and initial bins DataFrame ---
    if isinstance(bins, str | Path):
        bins_df = _read_bins_bed(bins)
        contig_lengths = _read_contig_lengths(bins=bins_df)
    elif isinstance(bins, pd.DataFrame):
        bins_df = bins.copy()
        contig_lengths = _read_contig_lengths(bins=bins_df)
    elif bins is None:
        # Integer binsize mode: need contig lengths from fasta or pairs header
        contig_lengths = _read_contig_lengths(fasta=fasta, pairs=pairs)
        bins_df = None
    else:
        raise TypeError("bins must be a path, DataFrame, or None")

    # --- Resolve the reference ("fasta") contig order, if one is requested ---
    # When lengths come from a bins table their dict order is whatever
    # groupby("chrom") produced (alphabetical), not the FASTA order.
    if contig_ordering == "fasta" and bins_df is not None:
        contig_lengths = _apply_reference_order(
            contig_lengths, fasta=fasta, pairs=pairs
        )

    contig_lengths, bins_df = _drop_nonpositive_contigs(contig_lengths, bins_df)
    logger.info(f"Determined lengths for {len(contig_lengths)} contigs")

    # --- Filter unassigned contigs ---
    if drop_unassigned:
        n_before = len(contig_lengths)
        contig_lengths = {
            ct: ls
            for ct, ls in contig_lengths.items()
            if ct in cg_map and cg_map[ct] != "__unbinned__"
        }
        if bins_df is not None:
            bins_df = bins_df[bins_df["chrom"].isin(contig_lengths)].copy()
        logger.info(
            f"Dropped {n_before - len(contig_lengths)} unassigned contigs; "
            f"{len(contig_lengths)} remain"
        )

    # --- Build canonical contig order and genome / contig tables ---
    contig_order = _build_contig_order(contig_lengths, cg_map, ordering=contig_ordering)
    logger.debug(f"Contig order: {contig_order[:10]} ... {contig_order[-10:]}")

    genomes_df = _build_genomes_table(
        contig_order,
        contig_lengths,
        cg_map,
        extra_genomes_df,
    )
    contig_name_to_id = {name: i for i, name in enumerate(contig_order)}
    contigs_df = _build_contigs_table(
        contig_order, contig_lengths, cg_map, genomes_df, contig_df
    )
    # Before the pairs are streamed: a name cooler cannot read costs a second
    # here and an entire binning pass if it is only caught at write time.
    check_ascii_names(genomes_df, contigs_df)
    logger.info(
        f"Built final genomes and contigs tables with {len(genomes_df)} genomes "
        f"and {len(contigs_df)} contigs"
    )

    # --- Build final bins table ---
    if bins_df is None:
        if binsize is None:
            raise ValueError("binsize required when bins is None")
        bins_df = _tile_bins(contig_order, contig_lengths, contig_name_to_id, binsize)
    else:
        bins_df = _build_bins(bins_df, contig_name_to_id)
    bins_df = bins_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    logger.info(
        f"Built final bins table with {len(bins_df)} bins tiling "
        f"{len(contig_lengths)} contigs"
    )

    # --- Compute chrom / genome offsets and insert n_bins into genomes_df ---
    chrom_offset, genome_offset = _compute_all_offsets(bins_df, contigs_df, genomes_df)

    # --- Stream pairs, bin, sort, and write the cooler file ---
    logger.info(f"Starting to create cooler file at {uri}...")
    bin_starts: np.ndarray | None = (
        None if binsize is not None else bins_df["start"].values.astype(np.int64)
    )
    raw_pixels: Iterable[pd.DataFrame] = _stream_binned_pairs(
        pairs,
        contig_name_to_id,
        chrom_offset,
        binsize,
        bin_starts,
        chunksize,
        pairs_kwargs,
    )
    pixel_source: Iterable[pd.DataFrame]
    if sort_in_memory is None:
        # Decide from the pixel volume actually produced, not from a proxy such
        # as the pairs file size: binning aggregates duplicates away, so even a
        # multi-GiB pairs file usually yields few enough pixels to sort in RAM —
        # which writes the cooler in a single pass instead of one temporary
        # cooler per chunk.
        frames, rest = _buffer_pixels(raw_pixels, _SORT_IN_MEMORY_PIXEL_THRESHOLD)
        if rest is None:
            pixel_source = _sort_and_dedup_pixels(pixel_iter=frames)
            ordered = True
        else:
            logger.warning(
                "Deduplicated pixels exceeded the %s-row in-memory budget "
                "(%s distinct pixels buffered so far) — switching to "
                "sort_in_memory=False (cooler external sort). Pass "
                "sort_in_memory=True to override.",
                f"{_SORT_IN_MEMORY_PIXEL_THRESHOLD:,}",
                f"{len(frames[0]):,}" if frames else "0",
            )
            pixel_source = _chain_pixel_frames(frames, rest)
            ordered = False
    elif sort_in_memory:
        pixel_source = _sort_and_dedup_pixels(pixel_iter=raw_pixels)
        ordered = True
    else:
        pixel_source = raw_pixels
        ordered = False
    pixel_source = _ensure_nonempty_pixel_stream(pixel_source)
    with _staged_output(uri) as staged:
        cooler.create_cooler(
            str(staged),
            bins=bins_df,
            pixels=pixel_source,
            ordered=ordered,
            symmetric_upper=True,
            assembly=assembly_name,
            metadata=metadata or {},
            dupcheck=False,
        )

        # --- Post-process: add metacooler-specific groups and attrs ---
        logger.info("Post-processing cooler file to add metacool structure")
        _finalize_metacool_file(
            staged,
            genomes_df,
            contigs_df,
            chrom_offset,
            genome_offset,
            contig_ordering=contig_ordering,
            assembly_name=assembly_name,
            binsize=binsize,
        )

        # --- Validate ---
        logger.info("Validating created metacool file")
        validate(staged)
    with h5py.File(uri, "r") as f:
        nnz = int(f.attrs.get("nnz"))
        total = int(f.attrs.get("sum", 0))
    n_bins = len(bins_df)
    n_contigs = len(contigs_df)
    n_genomes = len(genomes_df)
    logger.info(
        f"Created {n_bins} bins from {n_contigs} contigs in {n_genomes} genomes. "
        f"NNZ={nnz}, sum={total}.",
    )
    logger.info(f"Created metacool file: {uri} ({_total_time_str(timer)} elapsed)")


def from_cool(
    cool_uri: str | Path,
    output: str | Path,
    contig_genome_map: str | Path,
    contig_table_format: str = "contigs",
    genomes_table: str | Path | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "cool", "alpha"
    ] = "genome-then-length",
    assembly_name: str | None = None,
    drop_unassigned: bool = False,
    metadata: dict | None = None,
    nproc: int = 1,
    sort_in_memory: bool | None = None,
) -> None:
    """Create a metacool file from an existing cool file.

    Args:
        cool_uri: Path to the source cool file.
        output: Destination path for the output metacool file.
        contig_genome_map: Path to contig-to-genome assignment table
            (``str`` or ``Path`` only; dicts are not accepted here).
        contig_table_format: Column spec for the contig assignment table.
            Built-in values are ``"contigs"``, ``"contig_genome"``,
            ``"metator"``, and ``"I,J"``, but any valid colspec string is
            accepted.
        genomes_table: Optional path to genomes metadata table.
        contig_ordering: Strategy for ordering contigs within the metacool file.
            ``'cool'`` preserves the original contig order from the input cool file.
            Other choices: ``'genome-then-length'``, ``'genome-then-alpha'``,
            ``'alpha'``.
        assembly_name: Name for this assembly / metagenome, stored in file attrs.
        drop_unassigned: If ``True``, drop contigs absent from the contig table
            instead of collecting them under ``__unbinned__``.
        metadata: Optional metadata dictionary to store in the file attrs.
        nproc: Accepted for backwards compatibility, but currently has **no
            effect**: HDF5 reads must be serialised, so pixel remapping runs in
            a single thread.  Values greater than 1 log a warning.
        sort_in_memory: If ``True``, collect all pixels in RAM and sort before
            writing (fast for small files).  If ``False``, pass the unsorted
            stream to cooler's external sort (scalable for large files).  If
            ``None`` (default), auto-detect based on pixel count: files with
            more than 100 M pixels use ``False`` with a logged warning.

    Raises:
        ValueError: If ``cool_uri`` is not a valid cooler file.
    """
    cool_path = Path(cool_uri)
    output_path = Path(output)
    timer = time.perf_counter()

    # --- Validate and load cool file ---
    if not cooler.fileops.is_cooler(str(cool_path)):
        raise ValueError(f"Not a valid cooler file: {cool_path}")
    clr = cooler.Cooler(str(cool_path))
    bins = pd.DataFrame(clr.bins()[:])
    n_pixels = clr.info["nnz"]
    logger.info(f"Loaded cool file: {len(bins)} bins, {n_pixels} pixels")

    # --- Auto-detect sort_in_memory based on pixel count ---
    if sort_in_memory is None:
        if n_pixels > _SORT_IN_MEMORY_PIXEL_THRESHOLD:
            logger.warning(
                "Input has %d pixels — switching to sort_in_memory=False "
                "(cooler external sort). Pass sort_in_memory=True to override.",
                n_pixels,
            )
            sort_in_memory = False
        else:
            sort_in_memory = True

    # Build contig lengths from bins (use actual contigs in bins, not chroms table)
    chroms = clr.chroms()[:]
    contig_id_to_name = dict(enumerate(chroms["name"]))
    contig_lengths: dict[str, int] = {
        str(chrom): int(group["end"].max()) for chrom, group in bins.groupby("chrom")
    }
    contig_lengths, _ = _drop_nonpositive_contigs(contig_lengths)

    # --- Parse contig assignments ---
    contig_map, contig_df = _parse_contig_map(contig_genome_map, contig_table_format)

    # --- Filter unassigned contigs ---
    if drop_unassigned:
        contig_lengths = {
            ct: ln
            for ct, ln in contig_lengths.items()
            if ct in contig_map and contig_map[ct] != "__unbinned__"
        }

    # --- Build contig order ---
    if contig_ordering == "cool":
        # Preserve the original contig order from the cool file's chroms table,
        # filtered to only contigs present in contig_lengths (after drop_unassigned).
        contig_order = [
            name for name in contig_id_to_name.values() if name in contig_lengths
        ]
    else:
        contig_order = _build_contig_order(
            contig_lengths=contig_lengths,
            contig_genome_map=contig_map,
            ordering=contig_ordering,
        )

    # --- Build genome / contig tables ---
    extra_genomes_df = _load_extra_genomes_table(genomes_table)
    genomes = _build_genomes_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=contig_map,
        extra_genomes_df=extra_genomes_df,
    )
    contigs = _build_contigs_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=contig_map,
        genomes=genomes,
        contig_df=contig_df,
    )
    # Fail before the pixel stream, not after it (see from_pairs).
    check_ascii_names(genomes, contigs)

    # --- Remap bins to new contig IDs; compute offsets ---
    bins_remapped, old_to_new_bin = _remap_bins(bins, contig_order)
    chrom_offset, genome_offset = _compute_all_offsets(bins_remapped, contigs, genomes)

    # --- Stream, remap, sort, and write ---
    logger.info(
        f"Building metacool file: {len(genomes)} genomes, {len(contigs)} contigs"
    )
    # Preserve the source count width: hard-casting to int32 wraps the int64
    # counts ops/binnify.py deliberately writes, and floors float counts to 0.
    count_dtype = resolve_count_dtype(
        clr.pixels()[0:1]["count"].dtype, context="from_cool"
    )
    raw_pixels: Iterable[pd.DataFrame] = _stream_remapped_pixels(
        str(cool_path), clr.info["nnz"], old_to_new_bin, nproc, count_dtype
    )
    if sort_in_memory:
        pixel_source: Iterable[pd.DataFrame] = _sort_and_dedup_pixels(
            pixel_iter=raw_pixels
        )
        ordered = True
    else:
        pixel_source: Iterable[pd.DataFrame] = raw_pixels
        ordered = False
    pixel_source = _ensure_nonempty_pixel_stream(pixel_source)
    with _staged_output(output_path) as staged:
        cooler.create_cooler(
            str(staged),
            bins=bins_remapped,
            pixels=pixel_source,
            ordered=ordered,
            symmetric_upper=True,
            assembly=assembly_name or "imported",
            metadata=metadata or {},
            dtypes={"count": count_dtype},
        )

        # --- Add metacooler-specific groups and attributes ---
        logger.info("Post-processing cooler file to add metacool structure")
        _finalize_metacool_file(
            staged,
            genomes,
            contigs,
            chrom_offset,
            genome_offset,
            contig_ordering=contig_ordering,
            assembly_name=assembly_name,
        )

        # --- Validate ---
        logger.info("Validating created metacool file")
        validate(staged)
    with h5py.File(output_path, "r") as f:
        nnz = int(f.attrs.get("nnz"))
        total = int(f.attrs.get("sum"))
        n_bins = int(f.attrs.get("nbins"))
        n_contigs = int(f.attrs.get("nchroms"))
        n_genomes = int(f.attrs.get("ngenomes"))
    logger.info(
        f"Created {n_bins} bins from {n_contigs} contigs in {n_genomes} genomes. "
        f"NNZ={nnz}, sum={total}.",
    )
    logger.info(
        f"Created metacool file: {output_path} ({_total_time_str(timer)} elapsed)"
    )
