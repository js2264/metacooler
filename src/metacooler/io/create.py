"""Create metacool files from pairs, bins, and MAG assignments."""

from __future__ import annotations

import contextlib
import csv
import dataclasses
import gzip
import itertools
import logging
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar

import cooler
import h5py
import numpy as np
import pandas as pd
from pyfaidx import Fasta

from metacooler.io._utils import (
    _finalize_metacool_file,
    aggregate_pixels,
    check_ascii_names,
    check_output_path,
    infer_binsize,
    pixel_keys,
    pixels_from_keys,
    reduce_sorted_keys,
    resolve_count_dtype,
    store_qc,
    to_upper_triangle,
    write_cooler,
)
from metacooler.io.colspec import (
    ASSIGNMENT_KEYS,
    CONTIGS_FORMAT_ALIASES,
    CONTIGS_FORMAT_DEFAULT,
    GENOMES_FORMAT_DEFAULT,
    parse_column_spec,
)
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.genomes import parse_genomes_table
from metacooler.ps import (
    DEFAULT_BINS_PER_DECADE,
    PSCounts,
    assemble,
    pairs_distance_grid,
    position_pair_denominators,
    ps_table_path,
    tally_intra_contig,
    write_ps_table,
)
from metacooler.qc import (
    DEFAULT_DISTANCE_THRESHOLD,
    QC_COLUMN_PREFIX,
    QC_CONTIG_COLUMNS,
    QC_GENOME_COLUMNS,
    BinMaps,
    ContactQC,
    bin_maps,
    unbinned_genome_id,
)
from metacooler.schema import (
    CONTIGS_DTYPE,
    DEFAULT_CHUNKSIZE,
    GENOMES_DTYPE,
    validate,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Generator,
        Iterable,
        Iterator,
        Mapping,
        Sequence,
    )
    from typing import Literal


# Genome names a contig table uses to say "no genome": MetaTOR writes ``ND``,
# SemiBin2 ``-1``, and spreadsheets ``NA`` or an empty cell.
_UNBINNED_SENTINELS: frozenset[str] = frozenset(
    {"__unbinned__", "na", "nan", "n/a", "nd", "n/d", "-1", ""}
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
# is free to carry extra columns, but one of these names would collide.  The
# QC columns are reserved unconditionally, not only when QC is computed, so a
# user column called ``qc_inter_genome`` means the same thing on every file.
_GENOMES_RESERVED_COLUMNS = frozenset(
    {"name", "length", "n_contigs", "n_bins"} | set(QC_GENOME_COLUMNS)
)
_CONTIGS_RESERVED_COLUMNS = frozenset(
    {"name", "length", "genome_id"} | set(QC_CONTIG_COLUMNS)
)
# In /genomes the whole qc_* namespace is metacooler's (see qc.is_qc_column):
# every writer that derives a map from another strips it. A user column there
# would vanish at the first reassign/binnify/merge, so it is renamed on the way
# in instead. /contigs reserves only its two QC names; a user qc_* stays.
_GENOMES_RESERVED_PREFIXES = (QC_COLUMN_PREFIX,)
# How many distinct unknown contig names to remember for the dropped-pair report.
_MAX_REPORTED_UNKNOWN_CONTIGS = 10
# How many genome names to name when the contigs and genomes tables diverge.
_MAX_REPORTED_DIVERGENT_GENOMES = 10
# How many names to quote from each side when a contig table matches no contig
# of the assembly, and how many unmatched names to list when only some miss.
_MAX_SHOWN_TABLE_NAMES = 3
_MAX_REPORTED_STRAY_CONTIGS = 10

logger = logging.getLogger("metacooler")


def _read_contig_lengths(
    bins: pd.DataFrame | None = None,
    fasta: str | Path | None = None,
    pairs: str | Path | Sequence[str | Path] | None = None,
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
        pairs: One or several ``.pairs`` (or ``.pairs.gz``) files whose
            headers contain ``#chromsize: <name> <length>`` lines; several
            headers are united by :func:`_merge_pairs_headers`.

    Returns:
        Mapping of contig name to sequence length in base pairs.

    Raises:
        ValueError: If none of the three sources can provide contig lengths,
            or two pairs files give one contig different lengths.
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
        # Try to parse chromsizes from the pairs header(s)
        lengths = _merge_pairs_headers(_pairs_paths(pairs))
        if lengths:
            return lengths
    raise ValueError(
        "Contig lengths required. Provide --fasta, a BED bins file, "
        "or a pairs file with #chromsize: header lines."
    )


def _apply_reference_order(
    contig_lengths: dict[str, int],
    fasta: str | Path | None = None,
    pairs: str | Path | Sequence[str | Path] | None = None,
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
        pairs: Pairs file(s) whose ``#chromsize:`` headers define the
            reference order when no *fasta* is available: several are united
            by :func:`_merge_pairs_headers`, in file order.

    Returns:
        The same mapping, reordered to follow the reference order.  Contigs
        absent from the reference are appended (with a warning).

    Raises:
        ValueError: If neither a FASTA nor a pairs header can supply an order,
            or two pairs files give one contig different lengths.
    """
    reference: list[str] | None = None
    if fasta is not None:
        reference = list(_read_contig_lengths(fasta=fasta))
    elif pairs is not None:
        header = _merge_pairs_headers(_pairs_paths(pairs))
        if header:
            reference = list(header)
    if reference is None:
        raise ValueError(
            "contig_ordering='fasta' needs a reference contig order, but no FASTA "
            "was given and no pairs file has '#chromsize:' header lines. "
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


def _pairs_input_stream(path: Path):
    """Open a pairs file as a binary stream, decompressing by file extension.

    The single opener for every reader of a pairs file -- the header scans as
    well as the body -- so they cannot disagree about which codecs exist.  They
    used to: the body went through pyarrow (gzip, bz2, zstd, lz4, brotli) while
    the two header scans opened the file with their own narrower ``open``
    choices, so a ``.zst`` file -- advertised in the CLI reference -- died with
    ``UnicodeDecodeError`` before a single pair was read, and ``.bz2`` failed
    whenever ``#chromsize`` lines were needed.

    pyarrow's codec table has no ``.bgz`` entry, and would hand a bgzipped file
    to the CSV parser as raw binary (zero recognisable contigs, a structurally
    valid but empty map); bgzf is gzip-compatible, so it is named explicitly.

    Args:
        path: Pairs file path.

    Returns:
        A readable ``pyarrow.NativeFile``; close it when done.
    """
    import pyarrow as pa

    if path.suffix == ".bgz":
        return pa.input_stream(str(path), compression="gzip")
    return pa.input_stream(str(path))


def _iter_pairs_lines(path: Path, block_size: int = 1 << 20) -> Iterator[str]:
    """Yield the decoded text lines of a (possibly compressed) pairs file.

    pyarrow streams do not implement ``readline``, so lines are split out of
    fixed-size blocks.  Meant for header scans, which stop early.

    Args:
        path: Pairs file path.
        block_size: Bytes decompressed per read.

    Yields:
        Each line, with its trailing newline.
    """
    with contextlib.closing(_pairs_input_stream(path)) as src:
        carry = b""
        while True:
            block = src.read(block_size)
            if not block:
                break
            lines = (carry + block).split(b"\n")
            carry = lines.pop()
            for raw in lines:
                yield raw.decode("utf-8") + "\n"
        if carry:
            yield carry.decode("utf-8")


def _iter_pairs_raw_lines(path: Path, block_size: int = 1 << 20) -> Iterator[bytes]:
    """Yield the raw lines of a (possibly compressed) pairs file, without newlines.

    The header scans work on these bytes: a metagenome's pairs header carries
    one ``#chromsize:`` line per contig -- two million lines -- and decoding
    each one to ``str`` before looking at it doubled the cost of reading it.
    """
    with contextlib.closing(_pairs_input_stream(path)) as src:
        carry = b""
        while True:
            block = src.read(block_size)
            if not block:
                break
            lines = (carry + block).split(b"\n")
            carry = lines.pop()
            yield from lines
        if carry:
            yield carry


def _parse_pairs_header(pairs_path: str | Path) -> dict[str, int] | None:
    """Parse ``#chromsize:`` lines from the header of a pairs file.

    Reads lines from the beginning of the file until the first non-comment
    line is encountered.  Handles every codec :func:`_pairs_input_stream`
    does.

    Args:
        pairs_path: Path to the pairs (or ``.pairs.gz``) file.

    Returns:
        Mapping of contig name to length extracted from the header, or
        ``None`` if no ``#chromsize:`` lines were found.
    """
    path = Path(pairs_path)
    lengths = {}
    for raw in _iter_pairs_raw_lines(path):
        if not raw.startswith(b"#"):
            break
        if raw.startswith(b"#chromsize:"):
            parts = raw.split()
            if len(parts) >= 3:
                try:
                    lengths[parts[1].decode("utf-8")] = int(parts[2])
                except ValueError as exc:
                    line = raw.decode("utf-8", "replace").strip()
                    length = parts[2].decode("utf-8", "replace")
                    raise ValueError(
                        f"Malformed #chromsize line in {path}: {line!r} "
                        f"(length {length!r} is not an integer)"
                    ) from exc
    return lengths if lengths else None


def _pairs_paths(
    pairs: str | os.PathLike[str] | Iterable[str | os.PathLike[str]] | None,
) -> list[Path]:
    """Normalise *pairs* -- one path, or several -- to a non-empty list of paths.

    Several files are read back to back as one stream (see :func:`from_pairs`),
    so the same file given twice would count every one of its contacts twice;
    that is refused rather than summed.

    Args:
        pairs: A pairs file path, or an iterable of them, in reading order.

    Returns:
        The paths, in the given order.

    Raises:
        ValueError: If no path is given, or the same file is given twice.
    """
    if pairs is None:
        raise ValueError("pairs file path required")
    if isinstance(pairs, str | os.PathLike):
        paths = [Path(pairs)]
    else:
        paths = [Path(p) for p in pairs]
    if not paths:
        raise ValueError("pairs file path required: got an empty list")
    # Keyed on the file itself (device and inode), so an alias (``./x``,
    # ``sub/../x``, a symlink) or a hard link of a file already listed is caught
    # too; a path that does not exist falls back to its resolved form and is
    # reported missing by the reader.
    seen: dict[object, Path] = {}
    for path in paths:
        try:
            st = path.stat()
            key: object = (st.st_dev, st.st_ino)
        except OSError:
            key = path.resolve()
        first = seen.get(key)
        if first is not None:
            also = "" if first == path else f" (also as {str(first)!r})"
            raise ValueError(
                f"Pairs file {str(path)!r} is given twice{also}; every one of "
                "its contacts would be counted twice."
            )
        seen[key] = path
    return paths


def _merge_pairs_headers(paths: Sequence[Path]) -> dict[str, int] | None:
    """Unite the ``#chromsize:`` headers of several pairs files.

    Each file is one library aligned to the same assembly, so their headers
    normally agree line for line.  A file without a header, or with a subset of
    the contigs, adds nothing and is fine.  One contig given two different
    lengths, though, means the files were aligned to different assemblies: a
    pair would land in a different bin, or past the end of its contig,
    depending on which file it came from.  That is refused, naming both files.

    Args:
        paths: Pairs files, in reading order.

    Returns:
        Contig name → length in first-seen order (the first file's header
        order, then each contig new to a later file in that file's order), or
        ``None`` if no file has a ``#chromsize:`` line.

    Raises:
        ValueError: If two files give one contig different lengths.
    """
    merged: dict[str, int] = {}
    for i, path in enumerate(paths):
        header = _parse_pairs_header(path)
        if not header:
            continue
        if not merged:
            merged = header
            continue
        if header == merged:
            # The usual case -- every library carries the assembly's header --
            # decided by one C-level comparison instead of a loop per contig.
            continue
        conflicts = [c for c, n in header.items() if merged.get(c, n) != n]
        if conflicts:
            name = conflicts[0]
            source = next(
                p for p in paths[:i] if name in (_parse_pairs_header(p) or {})
            )
            raise ValueError(
                f"{len(conflicts)} contig(s) have different '#chromsize:' lengths "
                f"in different pairs files, e.g. {name!r}: {merged[name]} in "
                f"{str(source)!r} but {header[name]} in {str(path)!r}. Every "
                "pairs file must be aligned to the same assembly."
            )
        for name, length in header.items():
            merged.setdefault(name, length)
    return merged or None


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


def _genome_of(contigs: Iterable[str], contig_genome_map: dict[str, str]) -> np.ndarray:
    """Genome name of every contig, ``"__unbinned__"`` when unassigned (object array)."""
    get = contig_genome_map.get
    return np.array([get(c, "__unbinned__") for c in contigs], dtype=object)


def _genome_ranks(genome: np.ndarray) -> np.ndarray:
    """Position of each entry's genome in :func:`_genome_id_order` (int64 array)."""
    codes, uniques = pd.factorize(genome)
    position = {g: i for i, g in enumerate(_genome_id_order(uniques))}
    return np.array([position[g] for g in uniques], dtype=np.int64)[codes]


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

    if ordering in ("genome-then-alpha", "genome-then-length"):
        # Group by genome, visited in genome_id order (see _genome_id_order) so
        # the physical contig layout matches the genomes-table ordering; within
        # a genome, alphabetically or longest first (ties alphabetically).
        # One lexsort over all contigs: a Python sort per genome with a
        # per-contig key function took seconds on a metagenome.
        names = np.fromiter(contig_lengths.keys(), dtype=object, count=len(all_contigs))
        if names.size == 0:
            return []
        genome_rank = _genome_ranks(_genome_of(names, contig_genome_map))
        name_rank = np.empty(names.size, dtype=np.int64)
        name_rank[np.argsort(names.astype(str), kind="stable")] = np.arange(names.size)
        if ordering == "genome-then-alpha":
            order = np.lexsort((name_rank, genome_rank))
        else:
            lengths = np.fromiter(
                contig_lengths.values(), dtype=np.int64, count=names.size
            )
            order = np.lexsort((name_rank, -lengths, genome_rank))
        return names[order].tolist()

    raise ValueError(f"Unknown ordering: {ordering!r}")


def _normalize_genome_name(v: str) -> str:
    """Normalize sentinel values (empty, NA, nan, ND, n/d, …) to ``__unbinned__``."""
    return "__unbinned__" if v.strip().lower() in _UNBINNED_SENTINELS else v


def _protect_reserved_columns(
    df: pd.DataFrame,
    reserved: frozenset[str],
    keys: frozenset[str],
    context: str,
    reserved_prefixes: tuple[str, ...] = (),
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
        reserved_prefixes: Name prefixes metacooler owns outright, such as
            ``qc_`` in ``/genomes``: a user column with one is renamed too.

    Returns:
        *df* unchanged, or a renamed copy.
    """
    renames: dict[str, str] = {}
    taken = set(df.columns)
    for col in df.columns:
        if col in keys:
            continue
        if col not in reserved and not str(col).startswith(reserved_prefixes):
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


def _coerce_numeric_columns(df: pd.DataFrame, columns: Sequence[str]) -> None:
    """Store each text column of *df* in *columns* as numbers, when all parse.

    The one typing rule for user metadata: a column is numeric when every one
    of its values parses as a number (an empty cell is a missing value), and
    text otherwise -- so ``0.63`` sorts and filters as a number while a single
    ``N/A`` keeps its whole column text.  ``create``, ``collect`` and
    ``annotate`` all type their tables here, so one table gives one column type
    whichever command stores it.

    Args:
        df: Table to convert, in place.
        columns: Columns to consider; numeric ones are left alone.
    """
    for col in columns:
        if pd.api.types.is_string_dtype(df[col].dtype) or df[col].dtype == object:
            with contextlib.suppress(ValueError, TypeError):
                df[col] = pd.to_numeric(df[col], errors="raise")


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
    *,
    report_divergence: bool = True,
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
        report_divergence: Warn when *extra_genomes_df* and the contig
            assignments cover different genomes (see
            :func:`_report_genome_table_divergence`).  ``collect`` turns this
            off: its genomes come from the input files rather than from a
            contigs table, so it reports the mismatch itself, in those terms.

    Returns:
        DataFrame with columns ``name``, ``length``, ``n_contigs``, and
        any extra columns from *extra_genomes_df*.  Rows are ordered with
        ``"__unbinned__"`` first (when present), then alphabetically.
    """
    # One genome per contig ("__unbinned__" when unassigned, so it is only a
    # row when some contig lands in it), then per-genome sums in one pass.
    codes, uniques = pd.factorize(_genome_of(contig_order, contig_genome_map))
    lengths = np.fromiter(
        map(contig_lengths.__getitem__, contig_order),
        dtype=np.int64,
        count=len(contig_order),
    )
    genome_lengths = np.zeros(len(uniques), dtype=np.int64)
    np.add.at(genome_lengths, codes, lengths)
    genome_ncontigs = np.bincount(codes, minlength=len(uniques)).astype(np.int64)

    # Canonical genome_id order: __unbinned__ first (if present), then alphabetical.
    sorted_names = _genome_id_order(uniques)
    position = {g: i for i, g in enumerate(uniques)}
    rows = np.array([position[g] for g in sorted_names], dtype=np.int64)

    df = pd.DataFrame(
        {
            "name": sorted_names,
            "length": genome_lengths[rows],
            "n_contigs": genome_ncontigs[rows],
        }
    )

    # The genomes table only decorates the derived rows; say so out loud
    # whenever the two tables cover different sets of genomes.
    if extra_genomes_df is not None and report_divergence:
        _report_genome_table_divergence(sorted_names, extra_genomes_df)

    # Merge extra columns from the genomes file if provided
    if extra_genomes_df is not None and len(extra_genomes_df.columns) > 1:
        extra_genomes_df = _protect_reserved_columns(
            extra_genomes_df,
            reserved=_GENOMES_RESERVED_COLUMNS,
            keys=frozenset({"genome_id"}),
            context="Genomes table",
            reserved_prefixes=_GENOMES_RESERVED_PREFIXES,
        )
        extra_cols = [c for c in extra_genomes_df.columns if c != "genome_id"]
        extra = extra_genomes_df.rename(columns={"genome_id": "name"})[
            ["name"] + extra_cols
        ].copy()
        # De-duplicate on the merge key so a repeated genome name cannot explode
        # the genomes table into duplicate rows (last row wins).
        extra = extra.drop_duplicates(subset="name", keep="last")
        _coerce_numeric_columns(extra, extra_cols)
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
    unbinned_id = genome_name_to_id.get("__unbinned__")

    # Resolved once per distinct genome, then gathered for every contig.
    codes, uniques = pd.factorize(_genome_of(contig_order, contig_genome_map))
    unique_ids = np.empty(len(uniques), dtype=np.int64)
    for j, gname in enumerate(uniques):
        gid = genome_name_to_id.get(gname, unbinned_id)
        if gid is None:
            contig = contig_order[int(np.argmax(codes == j))]
            raise ValueError(
                f"Contig {contig!r} has no genome assignment and __unbinned__ is absent"
            )
        unique_ids[j] = gid

    df = pd.DataFrame(
        {
            "name": contig_order,
            "length": np.fromiter(
                map(contig_lengths.__getitem__, contig_order),
                dtype=np.int64,
                count=len(contig_order),
            ),
            "genome_id": unique_ids[codes],
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
            _coerce_numeric_columns(extra, extra_cols)
            df = df.merge(extra, on="name", how="left")

    # Ensure column order: name, length, genome_id, then extras
    core_cols = ["name", "length", "genome_id"]
    extra_cols = [c for c in df.columns if c not in core_cols]
    df = df[core_cols + extra_cols]
    return df.astype({k: v for k, v in CONTIGS_DTYPE.items() if k in df.columns})


def _apply_create_filters(
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    contig_df: pd.DataFrame | None,
    extra_genomes_df: pd.DataFrame | None,
    filters: Sequence[str] | None,
    bins_df: pd.DataFrame | None = None,
    *,
    genome_n_bins: Mapping[str, int] | None = None,
    report_divergence: bool = True,
    dropped_effect: str = (
        "are tiled into no bins, and their contacts are discarded while binning "
        "(reported in the drop summary below)"
    ),
) -> tuple[dict[str, int], pd.DataFrame | None]:
    """Drop the contigs that fail the ``--filter`` predicates, before binning.

    Filtering here rather than after the fact is the point of supporting it in
    ``create`` at all: a contig removed now is never tiled into bins and its
    contacts are never binned, so the sub-kilobase tail of a metagenomic
    assembly costs nothing instead of dominating the contig count (and with it
    the per-chunk ``cooler.create`` overhead).

    The ``/genomes`` and ``/contigs`` frames the expressions are evaluated
    against are built with the very same :func:`_build_genomes_table` /
    :func:`_build_contigs_table` the output file will use, so a filter written
    against a ``dump`` of an existing file means exactly the same thing here —
    including the ``user_`` renames applied to colliding column names and the
    numeric coercion of text columns.  The one absentee is ``genomes.n_bins``,
    which does not exist until the bins are tiled; it gets a specific error
    rather than a bare "unknown column".

    Args:
        contig_lengths: Mapping of contig name to length, already narrowed by
            ``drop_unassigned``.
        contig_genome_map: Mapping of contig name to genome name.
        contig_df: Optional parsed contig table supplying extra columns.
        extra_genomes_df: Optional parsed genomes table supplying extra columns.
        filters: Raw ``--filter`` expressions, or ``None``.
        bins_df: Optional bins DataFrame (string ``chrom`` names) to narrow in
            lockstep, so a pre-defined bin table cannot keep bins for a contig
            that was just dropped.
        genome_n_bins: Bin count of every genome, for a caller whose bins
            exist before anything is filtered (``collect``: they come with the
            input files).  When given, ``genomes.n_bins`` becomes filterable
            like any other column, instead of being refused with a hint.
        report_divergence: Forwarded to :func:`_build_genomes_table`.
        dropped_effect: What happens to a dropped contig, completing the
            sentence logged after filtering ("the N dropped contig(s) ...").

    Returns:
        ``(contig_lengths, bins_df)`` with the failing contigs removed.

    Raises:
        metacooler.filters.FilterError: If a filter is malformed, names a
            column the tables do not have, or leaves no contig at all.
    """
    from metacooler.filters import (
        FilterError,
        explain_empty,
        filter_mask,
        parse_filters,
    )

    parsed = parse_filters(filters)
    if not parsed:
        return contig_lengths, bins_df

    # Order is irrelevant for a row-wise predicate, so the tables are built
    # over the contigs as they stand rather than waiting for the canonical
    # order — which is computed from the *surviving* contigs anyway.
    order = list(contig_lengths)
    genomes_df = _build_genomes_table(
        order,
        contig_lengths,
        contig_genome_map,
        extra_genomes_df,
        report_divergence=report_divergence,
    )
    contigs_df = _build_contigs_table(
        order, contig_lengths, contig_genome_map, genomes_df, contig_df
    )
    contigs_df = contigs_df.assign(
        genome=contigs_df["genome_id"].map(genomes_df["name"])
    )
    hints: dict[str, str] = {}
    if genome_n_bins is not None:
        # Same position _compute_all_offsets gives the column in the output.
        genomes_df.insert(
            3, "n_bins", genomes_df["name"].map(genome_n_bins).astype(np.int64)
        )
    else:
        hints["n_bins"] = (
            "no bin has been tiled yet, so a genome's bin count is not known "
            "during 'create'. Filter on genomes.length or genomes.n_contigs "
            "instead, or filter the finished file with 'extract'"
        )

    keep_genomes = set(
        genomes_df.loc[
            filter_mask(genomes_df, parsed, table="genomes", hints=hints), "name"
        ]
    )
    keep_contigs = set(
        contigs_df.loc[
            filter_mask(contigs_df, parsed, table="contigs", hints=hints)
            & contigs_df["genome"].isin(keep_genomes),
            "name",
        ]
    )

    n_before, n_genomes_before = len(contig_lengths), len(genomes_df)
    filtered = {c: ln for c, ln in contig_lengths.items() if c in keep_contigs}
    if not filtered:
        raise FilterError(
            explain_empty(
                {"genomes": genomes_df, "contigs": contigs_df}, parsed, hints=hints
            )
        )
    n_genomes_after = len({contig_genome_map.get(c, "__unbinned__") for c in filtered})
    logger.info(
        "--filter kept %d of %d contig(s) in %d of %d genome(s); the %d dropped "
        "contig(s) %s.",
        len(filtered),
        n_before,
        n_genomes_after,
        n_genomes_before,
        n_before - len(filtered),
        dropped_effect,
    )
    if bins_df is not None:
        bins_df = bins_df[bins_df["chrom"].isin(filtered)].copy()
    return filtered, bins_df


def _check_no_binning(
    contig_genome_map: dict[str, str] | str | Path | None,
    contig_table_format: str,
    *,
    genomes_table: str | Path | None,
    drop_unassigned: bool,
) -> None:
    """Refuse what cannot go with a map that has no binning yet.

    Nothing is binned yet when there is no contig table at all (``None``: the
    map is a cooler plus the genome layer, with one ``__unbinned__`` genome),
    or when its spec says, with ``NA`` in the genome's place, that the table
    has no genome column.  Either way every contig lands in ``__unbinned__``,
    so two options cannot mean anything, and are refused before any work:
    ``drop_unassigned`` would drop every contig, and a genomes table would
    have no genome to describe.

    Args:
        contig_genome_map: The contig table's path, a ready-made mapping
            (which has no column spec), or ``None`` for no table.
        contig_table_format: The spec a table is read with.
        genomes_table: The genomes metadata table, if any.
        drop_unassigned: Whether unassigned contigs are to be dropped.

    Raises:
        ValueError: If *drop_unassigned* or *genomes_table* is given with no
            binning.
    """
    if contig_genome_map is None:
        cause = "no contig table was given"
    elif isinstance(contig_genome_map, dict):
        return
    else:
        spec = parse_column_spec(
            contig_table_format,
            keys=ASSIGNMENT_KEYS,
            aliases=CONTIGS_FORMAT_ALIASES,
            option="--contigs-format",
        )
        if not spec.no_genome:
            return
        cause = f"--contigs-format {spec.text!r} has no genome column"
    if drop_unassigned:
        raise ValueError(
            f"--drop-unassigned would drop every contig: {cause}, so no contig "
            "is assigned to a genome."
        )
    if genomes_table is not None:
        raise ValueError(
            f"A genomes table was given, but {cause}, so no contig is assigned "
            "to a genome and there is no genome for it to describe. Leave it "
            "out, and add it with 'metacooler annotate' once the contigs are "
            "binned ('metacooler reassign')."
        )
    logger.info(
        "%s: nothing is binned yet, and every contig is assigned to __unbinned__.",
        cause[0].upper() + cause[1:],
    )


def _quote_names(names: Iterable[str], limit: int) -> str:
    """Quote the first *limit* of *names*, comma-separated, ``...`` past them.

    ``repr`` makes whitespace visible -- a trailing space, a stray carriage
    return -- which is the usual reason two names that print alike differ.
    Only ``limit + 1`` names are drawn from *names*, so it may be a whole
    assembly's worth.
    """
    shown = list(itertools.islice(names, limit + 1))
    text = ", ".join(repr(str(name)) for name in shown[:limit])
    return text + (", ..." if len(shown) > limit else "")


def _check_contig_table_matches(
    contig_map: Mapping[str, str],
    assembly_contigs: Mapping[str, object],
    *,
    source: Mapping[str, str] | str | Path,
    contig_table_format: str,
    assembly_source: str,
) -> None:
    """Refuse a contig table that names no contig of the assembly.

    A contig the table does not list is collected under ``__unbinned__`` --
    a table may cover only part of an assembly -- so a table that matched
    *nothing* produced a map with a single ``__unbinned__`` genome holding every
    contig, exited 0, and logged nothing above INFO.  That is a wrong table, or
    the right one read with the wrong ``--contigs-format`` (``1,2`` over a
    table with a leading ID column takes the IDs for contig names and the
    contig names for genomes), never a metagenome anyone meant to build;
    ``reassign`` refuses the same table for the same reason.  So is a table
    with no row at all.

    A table that matches *some* contigs is accepted, and the names it holds
    beyond the assembly (another assembly's contigs, a header line read as
    data) are only warned about: one table covering several assemblies is
    legitimate.

    The callers run this against every contig the assembly names, before
    non-positive lengths, ``--drop-unassigned`` and ``--filter`` narrow them --
    the table is right or wrong about the assembly, whatever those then keep
    -- and before a single pair or pixel is read.

    Args:
        contig_map: Parsed ``{contig name: genome name}`` table.
        assembly_contigs: Every contig of the assembly, keyed by name, in the
            assembly's order (only the keys are read).
        source: What *contig_map* was parsed from: a path, which the messages
            name, or a mapping passed in directly.
        contig_table_format: The format a *source* path was read with.
        assembly_source: Where the assembly's contig names came from, for the
            messages, e.g. ``"the cool file x.cool"``.

    Raises:
        ValueError: If the table names no contig of the assembly, or no contig
            at all.
    """
    if isinstance(source, str | Path):
        table = f"the contig table {source}"
        hint = (
            f" and that --contigs-format (here {contig_table_format!r}) selects "
            "its contig-name and genome-name columns"
        )
    else:
        table, hint = "the contig-to-genome mapping", ""
    if not contig_map:
        raise ValueError(
            f"There is no row in {table}, so no contig would be assigned to a genome."
        )
    unmatched = [name for name in contig_map if name not in assembly_contigs]
    if len(unmatched) == len(contig_map):
        raise ValueError(
            f"None of the {len(contig_map)} contig name(s) in {table} matches a "
            f"contig in {assembly_source}, so no contig would be assigned to a "
            "genome.\n"
            f"  table names:      {_quote_names(contig_map, _MAX_SHOWN_TABLE_NAMES)}\n"
            "  assembly contigs: "
            f"{_quote_names(assembly_contigs, _MAX_SHOWN_TABLE_NAMES)}\n"
            f"Check that the table describes this assembly{hint}."
        )
    if unmatched:
        logger.warning(
            "%d of the %d name(s) in %s are not contigs in %s and were ignored: %s",
            len(unmatched),
            len(contig_map),
            table,
            assembly_source,
            _quote_names(unmatched, _MAX_REPORTED_STRAY_CONTIGS),
        )


def _require_assigned_contigs(contig_lengths: dict, n_before: int) -> None:
    """Refuse a ``--drop-unassigned`` that left no contig at all.

    Without this the run streamed every pair before cooler failed with "not
    enough values to unpack" -- hours later, and without naming the cause.  A
    table naming no contig of the assembly is refused before this, by
    :func:`_check_contig_table_matches`, so what reaches it is a table whose
    every match is unassigned.

    Args:
        contig_lengths: Contigs remaining after the drop.
        n_before: Contig count before it.

    Raises:
        ValueError: If no contig remains.
    """
    if not contig_lengths:
        raise ValueError(
            f"--drop-unassigned removed all {n_before} contig(s): the contig "
            "table assigns none of them to a genome (an empty, NA-like or "
            "__unbinned__ genome name counts as unassigned)."
        )


def _source_assembly(clr: cooler.Cooler) -> str | None:
    """Return a cooler's ``genome-assembly`` attribute, if it names one."""
    value = clr.info.get("genome-assembly")
    if value in (None, "", "unknown"):
        return None
    return str(value)


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
        # validate() rejects these too, but only after the whole pairs pass has
        # been paid for; a zero-width bin can never hold a contact anyway.
        example = df.loc[zero_width].head(3).to_dict("records")
        raise ValueError(
            f"{int(zero_width.sum())} bin(s) are zero-width (end == start), "
            f"e.g. {example}. Remove them from the bins file."
        )


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
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
        quoting=csv.QUOTE_NONE,
    )
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
    n = len(contig_order)
    # Tiled with three np.repeat calls rather than one DataFrame per contig plus
    # a final concat: the per-contig form costs ~200 us per contig (frame
    # construction, then a concat over as many blocks), which is invisible on a
    # genome assembly and becomes many minutes -- and gigabytes of transient
    # frame objects held at once -- on a metagenome with millions of contigs.
    ids = np.fromiter(
        (contig_name_to_id[c] for c in contig_order), dtype=np.int64, count=n
    )
    lengths = np.fromiter(
        (int(contig_lengths[c]) for c in contig_order), dtype=np.int64, count=n
    )
    # Build in int64 and range-check *before* narrowing: an int32 ``arange``
    # silently wraps to negative starts on a contig longer than ~2.1 Gbp, and
    # schema.validate() does not catch negative bin coordinates.
    over = np.flatnonzero(lengths > i32_max)
    if over.size:
        contig = contig_order[int(over[0])]
        raise ValueError(
            f"Contig {contig!r} is {lengths[over[0]]} bp, which exceeds the int32 "
            f"bin coordinate range ({i32_max} bp); contigs larger than ~2.1 Gbp "
            "are not supported."
        )
    # Ceiling division; clamped so a non-positive length tiles to no bins rather
    # than to a negative repeat count.
    n_bins = np.maximum(-(-lengths // binsize), 0)
    chrom = np.repeat(ids, n_bins)
    caps = np.repeat(lengths, n_bins)
    # Index of each bin within its own contig: a global arange minus the global
    # position at which that contig's run starts.
    run_starts = np.repeat(np.cumsum(n_bins) - n_bins, n_bins)
    starts = (np.arange(chrom.size, dtype=np.int64) - run_starts) * binsize
    ends = np.minimum(starts + binsize, caps)
    return pd.DataFrame({"chrom": chrom, "start": starts, "end": ends}).astype(
        {"chrom": np.int32, "start": np.int32, "end": np.int32}
    )


def _scan_pairs_header(path: Path) -> tuple[int, bool]:
    """Scan the top of a pairs file for its header size and whether data follows.

    Reads lines sequentially until the first non-comment line is found.
    Handles every codec :func:`_pairs_input_stream` does.

    Args:
        path: Path to the pairs file.

    Returns:
        ``(n_header, has_data_rows)`` — the count of consecutive leading ``#``
        lines, and whether at least one non-empty data row follows them.  The
        second element lets the caller recognise a header-only file *before*
        handing the stream to pyarrow, instead of pattern-matching pyarrow's
        error message (which would also swallow genuine parse failures).
    """
    n = 0
    has_data = False
    for raw in _iter_pairs_raw_lines(path):
        if raw.startswith(b"#"):
            n += 1
            continue
        if not raw.strip():
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
    chunksize: int = DEFAULT_CHUNKSIZE,
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
        logger.warning("Pairs file %s has no data rows.", path)
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
        # Contig names arrive dictionary-encoded: each batch holds every
        # distinct name once plus an int32 code per row, so the strings are
        # never materialised per pair and _decode_pairs_chunk maps only the
        # dictionary (a Python dict lookup per *distinct* name, not per row).
        column_types={
            f"f{c1 - 1}": pa.dictionary(pa.int32(), pa.string()),
            f"f{p1 - 1}": pa.int64(),
            f"f{c2 - 1}": pa.dictionary(pa.int32(), pa.string()),
            f"f{p2 - 1}": pa.int64(),
        },
    )

    with contextlib.closing(_pairs_input_stream(path)) as src:
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
    out_of_bounds: int = 0
    # An insertion-ordered set (a dict's keys): "the first names seen" must not
    # depend on string hashing, so that per-batch tallies merged in batch order
    # (create pairs --nproc) keep exactly the names a serial run keeps.
    unknown_names: dict[str, None] = dataclasses.field(default_factory=dict)

    def note_unknown_names(self, names: Iterable[str]) -> None:
        """Remember up to a handful of example unrecognised contig names."""
        for name in names:
            if len(self.unknown_names) >= _MAX_REPORTED_UNKNOWN_CONTIGS:
                return
            self.unknown_names[str(name)] = None

    def merge(self, other: _PairsDropStats) -> None:
        """Add a later batch's tally to this one, as if it had been seen here."""
        self.total += other.total
        self.unknown_contig += other.unknown_contig
        self.null_position += other.null_position
        self.out_of_bounds += other.out_of_bounds
        self.note_unknown_names(other.unknown_names)

    @property
    def dropped(self) -> int:
        """Total pairs discarded, across every reason."""
        return self.unknown_contig + self.null_position + self.out_of_bounds

    def report(self) -> None:
        """Log a WARNING summary of everything that was dropped (if anything)."""
        dropped = self.dropped
        if dropped == 0:
            return
        examples = sorted(self.unknown_names)
        if dropped == self.total:
            logger.warning(
                "ALL %d pair(s) were discarded (%d on unrecognised contigs, %d with "
                "a missing position, %d past the end of their contig) — the "
                "resulting map is empty. Do the contig names and lengths in the "
                "pairs file match the FASTA / contig table? "
                "Unrecognised names seen: %s",
                self.total,
                self.unknown_contig,
                self.null_position,
                self.out_of_bounds,
                examples,
            )
            return
        logger.warning(
            "Discarded %d of %d pair(s) (%.2f%%): %d on unrecognised contigs "
            "(e.g. %s), %d with a missing position, %d past the end of their "
            "contig (are the pairs aligned to this exact assembly?).",
            dropped,
            self.total,
            100.0 * dropped / max(self.total, 1),
            self.unknown_contig,
            examples,
            self.null_position,
            self.out_of_bounds,
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


def _contig_codes(names: pd.Series, contig_name_to_id: dict[str, int]) -> np.ndarray:
    """Map a column of contig names to ``int32`` ids, ``-1`` where unknown.

    A categorical column (what :func:`_stream_pairs` produces) is mapped through
    its categories, a few thousand names per batch, then gathered by code;
    anything else goes through ``Series.map``, one dict lookup per row.
    """
    if isinstance(names.dtype, pd.CategoricalDtype):
        # Over an object array: iterating the Arrow-backed categories pandas
        # builds costs as much again as the lookups, and on a metagenome each
        # batch's dictionary holds hundreds of thousands of names.
        cats = names.cat.categories.to_numpy(dtype=object)
        get = contig_name_to_id.get
        cat_ids = np.fromiter(
            (get(c, -1) for c in cats),
            dtype=np.int32,
            count=len(cats),
        )
        codes = names.cat.codes.to_numpy()
        # Code -1 is a missing name: unknown, as it was through Series.map.
        return np.where(codes >= 0, cat_ids[np.maximum(codes, 0)], -1).astype(np.int32)
    return names.map(contig_name_to_id).fillna(-1).astype(np.int32).to_numpy()


def _decode_pairs_chunk(
    chunk: pd.DataFrame,
    contig_name_to_id: dict[str, int],
    stats: _PairsDropStats | None = None,
    contig_length: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Turn a raw pairs chunk into usable ``(contig1, contig2, pos1, pos2)`` arrays.

    Everything that reads a pairs file goes through here — the binning pass of
    :func:`from_pairs` and the standalone ``ps --pairs`` — so "which pairs are
    usable" is decided once.  A second implementation would be free to disagree
    about an unrecognised contig name or an empty position field, and the two
    curves would then be normalised by different denominators without anything
    saying so.

    Args:
        chunk: Raw pairs DataFrame with columns ``chrom1``, ``pos1``,
            ``chrom2``, ``pos2``.
        contig_name_to_id: Mapping from contig name to integer ID.
        stats: Optional accumulator updated with the number of pairs seen and
            dropped, so the caller can report them once at end of stream.
        contig_length: Length of every contig, indexed by contig id.  When
            given, a pair with either end at or past its contig's end is
            dropped and counted, as ``cooler cload`` refuses it: binning used
            to fold it silently into the contig's last bin, and QC counted it
            as long-range signal.

    Returns:
        ``(c1, c2, p1, p2)`` for the usable rows — contig ids, and 0-based
        positions clamped at 0 — or ``None`` if the chunk has no usable row.
    """
    # Contig name → integer id; unrecognised (or missing) names become -1 so
    # they can be dropped in one vectorised pass.
    c1_codes = _contig_codes(chunk["chrom1"], contig_name_to_id)
    c2_codes = _contig_codes(chunk["chrom2"], contig_name_to_id)

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
    if contig_length is not None:
        inside = (raw_p1 < contig_length[c1]) & (raw_p2 < contig_length[c2])
        if not inside.all():
            if stats is not None:
                stats.out_of_bounds += int((~inside).sum())
            c1, c2, raw_p1, raw_p2 = (
                c1[inside],
                c2[inside],
                raw_p1[inside],
                raw_p2[inside],
            )
            if len(c1) == 0:
                return None
    neg_mask = (raw_p1 < 0) | (raw_p2 < 0)
    if neg_mask.any():
        logger.warning(
            "%d pairs clamped: negative 0-based positions encountered in chunk",
            int(neg_mask.sum()),
        )
    return c1, c2, np.maximum(raw_p1, 0), np.maximum(raw_p2, 0)


def _folded_bin_keys(chrom_offset: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Return ``(contig_id << 32) | start`` for every bin, as one sorted array.

    Bins are laid out contig by contig, each contig's sorted by start, so the
    folded key is non-decreasing and a single ``np.searchsorted`` over it finds
    the bin of any ``(contig, position)``.  Coordinates are guarded to fit
    ``int32`` at creation, so they never reach the contig bits.

    Args:
        chrom_offset: First bin of each contig, length ``n_contigs + 1``.
        starts: Start coordinate of every bin.

    Returns:
        ``int64`` array of folded keys, one per bin.
    """
    contig_of_bin = np.repeat(
        np.arange(len(chrom_offset) - 1, dtype=np.int64), np.diff(chrom_offset)
    )
    return (contig_of_bin << 32) | starts.astype(np.int64)


def _fragment_bin_ids(
    bin_keys: np.ndarray,
    chrom_offset: np.ndarray,
    contig: np.ndarray,
    pos: np.ndarray,
) -> np.ndarray:
    """Map ``(contig, 0-based position)`` pairs to variable-width bin ids.

    The bin holding *pos* is the last one of its contig starting at or before
    it.  A position before the contig's first bin start is clamped into that
    first bin; one past the last start lands in the last bin, because every key
    of the next contig is larger.

    Args:
        bin_keys: Output of :func:`_folded_bin_keys`.
        chrom_offset: First bin of each contig, length ``n_contigs + 1``.
        contig: Contig id of each position.
        pos: 0-based position, within ``[0, contig length)``.

    Returns:
        ``int64`` bin id of each position.
    """
    keys = (contig.astype(np.int64) << 32) | pos.astype(np.int64)
    idx = np.searchsorted(bin_keys, keys, side="right") - 1
    return np.maximum(idx, chrom_offset[contig]).astype(np.int64, copy=False)


def _bin_pairs_chunk(
    chunk: pd.DataFrame,
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_keys: np.ndarray | None,
    stats: _PairsDropStats | None = None,
    qc: ContactQC | None = None,
    ps: PSCounts | None = None,
    contig_genome_id: np.ndarray | None = None,
    contig_length: np.ndarray | None = None,
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
        bin_keys: Folded ``(contig, start)`` key of every bin, from
            :func:`_folded_bin_keys`; required when ``binsize`` is ``None``.
        stats: Optional accumulator updated with the number of pairs seen and
            dropped, so the caller can report them once at end of stream.
        qc: Optional contact-quality tally, updated from the *raw positions* of
            this chunk.  Doing it here is what makes the intra-contig
            short/long split exact: after binning, the only distances left are
            multiples of the bin size.
        ps: Optional P(s) tally, updated from the same raw positions and for
            the same reason — a curve built after binning cannot resolve
            anything below one bin.
        contig_genome_id: Genome id of every contig, required when *qc* or *ps*
            is given.
        contig_length: Length of every contig; pairs past a contig's end are
            dropped (see :func:`_decode_pairs_chunk`).

    Returns:
        Aggregated pixel DataFrame or ``None`` if the chunk is empty after
        filtering.
    """
    decoded = _decode_pairs_chunk(chunk, contig_name_to_id, stats, contig_length)
    if decoded is None:
        return None
    _tally_decoded(decoded, qc, ps, contig_genome_id)
    return _bin_decoded(decoded, chrom_offset, binsize, bin_keys)


_Decoded = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def _tally_decoded(
    decoded: _Decoded,
    qc: ContactQC | None,
    ps: PSCounts | None,
    contig_genome_id: np.ndarray | None,
) -> None:
    """Feed one decoded chunk's raw positions to the QC and P(s) tallies."""
    if qc is None and ps is None:
        return
    assert contig_genome_id is not None, "contig_genome_id required with qc/ps"
    c1, c2, p1, p2 = decoded
    # One subtraction serves both tallies; neither reads it where the two
    # ends sit on different contigs.
    distance = np.abs(p2 - p1)
    if qc is not None:
        qc.tally(c1, c2, contig_genome_id[c1], contig_genome_id[c2], distance)
    if ps is not None:
        tally_intra_contig(ps, contig_genome_id, c1, c2, distance)


def _bin_decoded(
    decoded: _Decoded,
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_keys: np.ndarray | None,
) -> pd.DataFrame:
    """Bin one decoded chunk into sorted, aggregated pixels (one contact per pair)."""
    c1, c2, p1, p2 = decoded
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
        # Fragment-level bins: one global binary search over keys that fold the
        # contig id in above the coordinate, so (contig, start) order is plain
        # integer order -- never a per-contig loop (see _folded_bin_keys).
        assert bin_keys is not None, "bin_keys required for fragment-level binning"
        bin1 = _fragment_bin_ids(bin_keys, chrom_offset, c1, p1)
        bin2 = _fragment_bin_ids(bin_keys, chrom_offset, c2, p2)

    # One contact per pair, in the upper triangle: sort the folded keys and
    # count the runs.  The result comes out sorted, which every later merge
    # exploits.
    bin1, bin2 = to_upper_triangle(bin1, bin2)
    keys = np.sort(pixel_keys(bin1, bin2))
    return pixels_from_keys(*reduce_sorted_keys(keys, None))


_T = TypeVar("_T")
_R = TypeVar("_R")


@dataclasses.dataclass
class _BinnedBatch:
    """What binning one pairs batch produced, before it is folded into the run."""

    pairs: int
    pixels: pd.DataFrame | None
    decoded: _Decoded | None
    stats: _PairsDropStats


def _decode_and_bin(
    chunk: pd.DataFrame,
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_keys: np.ndarray | None,
    contig_length: np.ndarray | None,
    keep_decoded: bool,
) -> _BinnedBatch:
    """The part of :func:`_bin_pairs_chunk` that touches no shared state.

    Safe to run on a worker thread: it reads the shared lookup tables and
    returns everything that must be accumulated -- the drop tally, and the
    decoded arrays the QC and P(s) tallies need -- for the caller to fold in,
    in batch order.
    """
    stats = _PairsDropStats()
    decoded = _decode_pairs_chunk(chunk, contig_name_to_id, stats, contig_length)
    if decoded is None:
        return _BinnedBatch(len(chunk), None, None, stats)
    pixels = _bin_decoded(decoded, chrom_offset, binsize, bin_keys)
    return _BinnedBatch(len(chunk), pixels, decoded if keep_decoded else None, stats)


def _map_in_order(
    function: Callable[[_T], _R],
    batches: Iterable[_T],
    nproc: int,
) -> Iterator[_R]:
    """``map(function, batches)``, on *nproc* worker threads when ``nproc > 1``.

    The thread that pulls *batches* (pyarrow's reader: gunzip and parse, which
    release the GIL) keeps reading while the workers bin (NumPy, which releases
    it too).  Results are yielded in batch order, with at most *nproc* batches
    in flight, so memory stays at *nproc* batches.
    """
    if nproc <= 1:
        for batch in batches:
            yield function(batch)
        return
    from collections import deque
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=nproc) as pool:
        in_flight: deque = deque()
        for batch in batches:
            in_flight.append(pool.submit(function, batch))
            del batch
            if len(in_flight) >= nproc:
                yield in_flight.popleft().result()
        while in_flight:
            yield in_flight.popleft().result()


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
    contig_genome_map: dict[str, str] | str | Path | None,
    contig_table_format: str,
) -> tuple[dict[str, str], pd.DataFrame | None]:
    """Parse and normalize a contig-to-genome mapping from a dict or file path.

    Args:
        contig_genome_map: Either a pre-built ``{contig: genome}`` mapping, a
            path to a contig table file, or ``None`` for no table (nothing is
            binned: every contig goes to ``__unbinned__``).
        contig_table_format: Format string passed to ``parse_contigs_table`` when
            *contig_genome_map* is a file path.

    Returns:
        Tuple of ``(cg_map, contig_df)`` where ``contig_df`` is ``None`` when
        *contig_genome_map* is not a file path.
    """
    contig_df: pd.DataFrame | None = None
    if contig_genome_map is None:
        return {}, None
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
    # Normalise each distinct genome name once (there are thousands, against
    # millions of contigs), and rebuild the map only if a sentinel occurs.
    renamed = {
        v: n for v in set(cg_map.values()) if (n := _normalize_genome_name(v)) != v
    }
    if renamed:
        cg_map = {k: renamed.get(v, v) for k, v in cg_map.items()}
    return cg_map, contig_df


def _load_extra_genomes_table(
    genomes_table: str | Path | None,
    fmt: str = GENOMES_FORMAT_DEFAULT,
) -> pd.DataFrame | None:
    """Load an optional genomes metadata table from disk.

    Args:
        genomes_table: Path to the TSV file, or ``None``.
        fmt: Which column holds the genome names (see
            :func:`metacooler.io.genomes.parse_genomes_table`).

    Returns:
        Parsed DataFrame or ``None`` if *genomes_table* is ``None``.
    """
    if genomes_table is None:
        return None
    df = parse_genomes_table(genomes_table, fmt)
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
    all_pixels = aggregate_pixels(frames, consume=True)
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
    return aggregate_pixels(frames)


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
) -> tuple[Iterable[pd.DataFrame], bool]:
    """Return *pixel_iter* made safe to write, and whether it was empty.

    ``cooler.create_cooler`` writes a malformed pixel table (a dense block of
    phantom rows with ``nnz=0``) when handed a generator that yields *nothing*,
    but handles an empty DataFrame correctly -- on the *ordered* path only:
    with ``ordered=False`` its unordered writer fails with "No objects to
    concatenate".  So an empty stream becomes one empty frame, and the caller
    must write it ordered (an empty table is trivially sorted).  This turns a
    zero-pixel input (a header-only pairs file, every pair on an unknown or
    filtered contig, an ``nnz=0`` cooler) into a valid empty ``.metacool``
    under every sort strategy.

    Only the first frame is pulled to decide, so a lazy stream stays lazy.

    Args:
        pixel_iter: The pixel frames to write.

    Returns:
        ``(frames, was_empty)``.
    """
    it = iter(pixel_iter)
    first = next(it, None)
    if first is None:
        empty = pd.DataFrame(
            {
                "bin1_id": np.array([], dtype=np.int64),
                "bin2_id": np.array([], dtype=np.int64),
                "count": np.array([], dtype=np.int32),
            }
        )
        return iter([empty]), True
    return itertools.chain([first], it), False


def _stream_pairs_files(
    paths: Sequence[Path],
    chunksize: int,
    pairs_kwargs: dict,
) -> Iterator[pd.DataFrame]:
    """Stream several pairs files back to back, as one sequence of raw batches.

    No batch spans two files, but nothing downstream cares where one file ends:
    :func:`_stream_binned_pairs` coalesces batches across the boundary, and the
    QC, P(s) and drop tallies are sums.  So several libraries of one assembly
    give exactly the map, tallies and P(s) table of their concatenation.

    Args:
        paths: Pairs files, in reading order.
        chunksize: Rows per batch, forwarded to :func:`_stream_pairs`.
        pairs_kwargs: Column indices, forwarded to :func:`_stream_pairs`.

    Yields:
        Raw pairs DataFrames (``chrom1``, ``pos1``, ``chrom2``, ``pos2``).
    """
    for i, path in enumerate(paths, start=1):
        if len(paths) > 1:
            logger.info("Reading pairs file %d of %d: %s", i, len(paths), path)
        yield from _stream_pairs(path, chunksize=chunksize, **pairs_kwargs)


def _stream_binned_pairs(
    pairs: str | Path | Sequence[str | Path],
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_keys: np.ndarray | None,
    chunksize: int,
    pairs_kwargs: dict,
    stats: _PairsDropStats | None = None,
    qc: ContactQC | None = None,
    ps: PSCounts | None = None,
    contig_genome_id: np.ndarray | None = None,
    contig_length: np.ndarray | None = None,
    nproc: int = 1,
) -> Iterable[pd.DataFrame]:
    """Stream pairs files and yield binned pixel DataFrames with progress.

    Args:
        pairs: Path to a pairs file (plain text, gzipped, or bgzipped), or a
            sequence of them, read back to back in the given order (see
            :func:`_stream_pairs_files`).
        contig_name_to_id: Mapping from contig name to integer ID.
        chrom_offset: Cumulative bin counts per contig (length n_contigs+1).
        binsize: Fixed bin size, or ``None`` for fragment-level bins.
        bin_keys: Folded ``(contig, start)`` bin keys from
            :func:`_folded_bin_keys`; required when *binsize* is ``None``.
        chunksize: Rows per streaming batch.
        pairs_kwargs: Extra keyword arguments forwarded to ``_stream_pairs``.
        stats: Drop tally to fill in.  Passed in (rather than created here)
            when the caller needs the pair counts afterwards — the QC block
            reports them as the denominator it did *not* use.
        qc: Optional contact-quality tally, forwarded to
            :func:`_bin_pairs_chunk`.
        ps: Optional P(s) tally, forwarded likewise.
        contig_genome_id: Genome id of every contig, required when *qc* or *ps*
            is given.
        nproc: Worker threads binning batches while the pairs file is read
            (see :func:`_map_in_order`); 1 bins on the reading thread.

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
        stats = stats if stats is not None else _PairsDropStats()
        if qc is not None or ps is not None:
            assert contig_genome_id is not None, "contig_genome_id required with qc/ps"

        def bin_batch(chunk: pd.DataFrame) -> _BinnedBatch:
            return _decode_and_bin(
                chunk,
                contig_name_to_id,
                chrom_offset,
                binsize,
                bin_keys,
                contig_length,
                keep_decoded=qc is not None or ps is not None,
            )

        batches = _stream_pairs_files(_pairs_paths(pairs), chunksize, pairs_kwargs)
        for batch in _map_in_order(bin_batch, batches, nproc):
            progress.update(task, advance=batch.pairs)
            # Folded in batch order, so every tally is the serial one.
            stats.merge(batch.stats)
            if batch.decoded is not None:
                _tally_decoded(batch.decoded, qc, ps, contig_genome_id)
            if batch.pixels is not None and not batch.pixels.empty:
                pending.append(batch.pixels)
            pending_pairs += batch.pairs
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


def _tally_pixel_qc(
    stream: Iterable[pd.DataFrame],
    tally: ContactQC,
    maps: BinMaps,
) -> Iterable[pd.DataFrame]:
    """Pass a pixel stream through unchanged, tallying QC on the way.

    ``from_cool`` has no read positions to classify, but it is already
    streaming every pixel to the writer; teeing that stream costs one pass over
    each chunk in NumPy and saves the caller a second full read of the file
    afterwards.

    Args:
        stream: Pixel chunks headed for the writer.
        tally: Accumulator to update.
        maps: Per-bin lookup arrays for the *output* bin table.

    Yields:
        Each chunk of *stream*, untouched.
    """
    for chunk in stream:
        tally.tally_pixels(
            chunk["bin1_id"].to_numpy(),
            chunk["bin2_id"].to_numpy(),
            chunk["count"].to_numpy(),
            maps,
        )
        yield chunk


def _stream_remapped_pixels(
    cool_path_str: str,
    n_pixels: int,
    old_to_new_bin: np.ndarray,
    nproc: int,
    count_dtype: np.dtype,
    chunksize: int = DEFAULT_CHUNKSIZE,
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
        chunk = clr.pixels()[["bin1_id", "bin2_id", "count"]][lo:hi]
        b1 = old_to_new_bin[chunk["bin1_id"].to_numpy(dtype=np.int64)]
        b2 = old_to_new_bin[chunk["bin2_id"].to_numpy(dtype=np.int64)]
        counts = chunk["count"].to_numpy().astype(count_dtype)
        # Drop pixels that touch a dropped contig (sentinel -1).
        keep = (b1 >= 0) & (b2 >= 0)
        if not keep.all():
            b1, b2, counts = b1[keep], b2[keep], counts[keep]
        b1, b2 = to_upper_triangle(b1, b2)
        return pd.DataFrame({"bin1_id": b1, "bin2_id": b2, "count": counts}, copy=False)

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


def _partial_path(path: Path) -> Path:
    """Where a writer builds *path* before moving it into place."""
    return path.with_name(path.name + ".partial")


def _refuse_staging_an_input(
    final: str | Path,
    inputs: Iterable[str | Path | None],
    companions: Sequence[Path] = (),
) -> None:
    """Refuse, up front, a run whose staged output would be one of its inputs.

    :func:`_staged_output` refuses it too, but only once the run gets there --
    after every pair is binned, or the whole QC pass of ``reassign``.  Each
    writer calls this where it checks its output path, before it reads
    anything.

    Args:
        final: The output path.
        inputs: What the run reads, as paths or URIs.
        companions: Files staged beside *final* (see :func:`_staged_output`).

    Raises:
        ValueError: If ``<final>.partial`` (or a companion's) is an input.
    """
    for path in (Path(final), *companions):
        check_output_path(
            _partial_path(path), force=True, inputs=inputs, defer_removal=True
        )


def _move_into_place(
    companions: Sequence[tuple[Path, Path]], partial: Path, final: Path
) -> None:
    """Move staged *companions*, then *partial*, into place -- all or none.

    A companion's previous file is moved aside, to ``<name>.previous``, before
    its new one takes its place, and if a later move fails -- *final*'s
    included -- every new companion goes back to its partial name and every
    previous file comes back, so a failed run leaves the previous files as
    they were.  Once *final* is in place the previous files are removed; one
    that cannot be is left where it is, with a WARNING (and an interrupt is
    re-raised), never put back: the new files are complete and match.  Only a
    process killed during the moves, where nothing can undo them, leaves them
    half done: a previous companion still at ``<name>.previous``, or a new
    companion beside the previous *final*.  Either way, a ``<name>.previous``
    left behind makes the next run refuse to start
    (:func:`_refuse_a_set_aside_file`).

    Args:
        companions: ``(partial, final)`` pairs, moved first.
        partial: The main output, as staged.
        final: Where it goes, last.
    """
    aside: list[tuple[Path, Path]] = []
    placed: list[tuple[Path, Path]] = []
    try:
        for staged, path in companions:
            if path.exists():
                previous = _set_aside_path(path)
                os.replace(path, previous)
                aside.append((previous, path))
            os.replace(staged, path)
            placed.append((staged, path))
        os.replace(partial, final)
    except BaseException:
        for staged, path in reversed(placed):
            os.replace(path, staged)
        for previous, path in reversed(aside):
            os.replace(previous, path)
        raise
    for previous, path in aside:
        try:
            previous.unlink(missing_ok=True)
        except BaseException as error:
            logger.warning(
                "%s is in place with %s, but the previous %s, set aside as %s, "
                "could not be removed (%r). Remove it: the next run refuses to "
                "start while it is there.",
                final,
                path,
                path.name,
                previous,
                error,
            )
            if not isinstance(error, Exception):
                raise


def _set_aside_path(path: Path) -> Path:
    """Where :func:`_move_into_place` keeps *path*'s previous file meanwhile."""
    return path.with_name(path.name + ".previous")


def _refuse_a_set_aside_file(final: Path, companions: Iterable[Path]) -> None:
    """Refuse to start while a companion's previous file is still set aside.

    :func:`_move_into_place` leaves ``<name>.previous`` behind only when it was
    killed during its moves, or could not remove that file afterwards.  Which
    of the two files then matches *final* cannot be told from here, and
    setting the companion aside again would overwrite one of them, so the run
    stops before it reads anything.

    Args:
        final: The output path.
        companions: Files moved into place with *final*.

    Raises:
        FileExistsError: If ``<companion>.previous`` exists.
    """
    for path in companions:
        previous = _set_aside_path(path)
        if previous.exists():
            raise FileExistsError(
                f"{str(previous)!r} exists. A run leaves such a file when it "
                "was killed while it moved its map and table into place, or "
                "could not remove the previous table afterwards. If so, check "
                f"which table matches {str(final)!r}, {str(path)!r} or "
                f"{str(previous)!r}: keep that one as {str(path)!r} and remove "
                "the other. If the file is one of yours, rename it."
            )


@contextlib.contextmanager
def _staged_output(
    final: Path,
    *,
    keep_partial: bool = True,
    inputs: Iterable[str | Path | None] = (),
    companions: Sequence[Path] = (),
) -> Generator[Path, None, None]:
    """Build the output at a sibling temp path, moving it into place on success.

    Writing straight to the final path meant a partial or invalid file could sit
    at the name the user asked for, and the old ``except BaseException: unlink()``
    guard fired on ``KeyboardInterrupt`` — Ctrl-C during the final ``validate()``
    of a multi-hour run deleted the finished file. Staging solves both: nothing
    ever appears at the final name until validation has passed, and a failure
    leaves the artifact behind (under the temp name) so it can be diagnosed.

    Args:
        final: Path the caller wants the finished file to have.
        keep_partial: What a failure does with the temp file.  ``True`` leaves
            it for inspection and says so at ERROR (an exception) or silently
            (``KeyboardInterrupt``) -- right for a person at a terminal.
            ``False`` removes it whatever stopped the run, ``KeyboardInterrupt``
            included, and says so at INFO -- for an API caller that stages its
            own output and only wants the finished file or nothing.  The final
            path is never written either way, and a failure before the temp
            file exists has nothing to keep or remove, so it says nothing.
        inputs: What the caller reads, as paths or URIs.  A temp file left by
            an earlier run is removed before writing, unless it is one of
            these -- a partial map given as input, to be rescued into the very
            path it was staged for -- which is refused instead.
        companions: Files the caller writes beside *final*, each at its
            :func:`_partial_path`, and moved into place with it, before it
            (:func:`_move_into_place`): ``from_pairs``' ``.ps.tsv`` table.
            On failure their partial files are kept or removed like
            *final*'s.

    Yields:
        The temporary path to write to.

    Raises:
        ValueError: If a temp path is one of *inputs*.
    """
    tmp = _partial_path(final)
    staged = [(tmp, final), *((_partial_path(c), c) for c in companions)]
    for partial, _ in staged:
        check_output_path(partial, force=True, inputs=inputs)
    try:
        yield tmp
        _move_into_place(staged[1:], tmp, final)
    except BaseException as error:
        left = [(partial, path) for partial, path in staged if partial.exists()]
        for partial, path in left:
            if not keep_partial:
                partial.unlink()
                logger.info(
                    "Creation failed; the partial file %s was removed and %s "
                    "was not written.",
                    partial,
                    path,
                )
            elif isinstance(error, Exception):
                logger.error(
                    "Creation failed. The partial file was left at %s for "
                    "inspection; %s was not written.",
                    partial,
                    path,
                )
        raise


def _remap_bins(
    bins: pd.DataFrame,
    contig_order: list[str],
    *,
    keep_extra_columns: bool = False,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Remap a bins DataFrame to new contig IDs and return the old→new bin mapping.

    Bins are sorted by new chrom ID (stable sort) so that
    ``cooler.create_cooler`` receives them in the canonical order.

    Args:
        bins: Original bins DataFrame (``chrom`` column contains contig names).
        contig_order: List of contig names defining the new ordering.
        keep_extra_columns: Carry every column beyond ``chrom``/``start``/``end``
            through the permutation instead of dropping it.  ``from_cool`` leaves
            this off — it is *importing* a foreign cooler, whose per-bin columns
            (a ``weight`` from a balancing run on the source's own bin ids) it
            makes no claim about.  ``reassign`` turns it on, because there the
            permutation is the only thing that happens: the bins themselves are
            unchanged, so a balancing weight still describes its bin and dropping
            it would silently unbalance the file.

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
    core = ["chrom", "start", "end"]
    extras = (
        [c for c in bins_work.columns if c not in core] if keep_extra_columns else []
    )
    bins_remapped = bins_work.loc[kept_mask, core + extras].copy()
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
    pairs: str | Path | Sequence[str | Path] | None = None,
    contig_genome_map: dict[str, str] | str | Path | None = None,
    contig_table_format: str = CONTIGS_FORMAT_DEFAULT,
    fasta: str | Path | None = None,
    assembly_name: str | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "fasta", "alpha"
    ] = "genome-then-length",
    drop_unassigned: bool = False,
    filters: Sequence[str] | None = None,
    binsize: int | None = None,
    genomes_table: str | Path | None = None,
    genomes_table_format: str = GENOMES_FORMAT_DEFAULT,
    pairs_kwargs: dict | None = None,
    metadata: dict | None = None,
    chunksize: int = DEFAULT_CHUNKSIZE,
    sort_in_memory: bool | None = None,
    qc: bool = True,
    qc_distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
    ps: bool = True,
    ps_bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
    force: bool = False,
    nproc: int = 1,
    keep_partial: bool = True,
) -> None:
    """Create a ``.metacool`` file from pairs, bins, and contig assignments.

    Args:
        uri: Output file path (``.metacool`` extension recommended).
        bins: Path to a BED file of bins (without header), or a pre-built DataFrame
            with columns ``chrom, start, end``.  If ``None``, *binsize* must be given
            and bins will be tiled from the fasta.
        pairs: Path to a pairs file (plain text, gzipped, or bgzipped), or a
            sequence of them -- typically one per Hi-C library of the same
            assembly.  Several files are read back to back in the given order,
            in one pass, as if they were one concatenated file: their contacts
            are summed, and the QC block, the P(s) table and the drop report
            cover them all, exactly -- where per-library maps joined by
            ``merge`` lose the QC block and have no P(s) table.  Their
            ``#chromsize:`` headers are united, and one contig given two
            different lengths is refused, as is the same file given twice.
        contig_genome_map: Either a ``dict[contig_name → genome_name]`` or a path
            to a contig table file (see *contig_table_format*).  Contigs it
            does not list are collected under ``__unbinned__``, and names that
            are not contigs of the assembly are ignored with a warning, but it
            must name at least one contig of the assembly.  ``None`` builds a
            map with no binning yet: a cooler plus the genome layer, every
            contig in ``__unbinned__`` (bin them later with ``reassign``).
        contig_table_format: Which columns of a table *contig_genome_map* hold
            the contig and the genome names: ``CONTIG,GENOME[,noheader]``, each
            a 1-based position or a header name, or an alias (see
            :func:`metacooler.io.contigs.parse_contigs_table`).
        fasta: Assembly fasta (required if *bins* is ``None``).
        assembly_name: Name for this assembly / metagenome.
        contig_ordering: Canonical ordering strategy.
        drop_unassigned: If ``True``, drop contigs not present in the contig
            table instead of assigning them to ``__unbinned__``.
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`)
            applied to the derived ``/genomes`` and ``/contigs`` tables before
            anything is binned, e.g. ``["contigs.length >= 5kb"]``.  A dropped
            contig gets no bins, and its pairs are discarded while binning.
            ``genomes.n_bins`` is not available here.
        binsize: Fixed bin size (used only when *bins* is ``None``).
        genomes_table: Optional path to a genomes metadata TSV.  One column
            holds the genome names (see *genomes_table_format*); every other
            column is stored in the ``/genomes`` table.
        genomes_table_format: Which column of *genomes_table* holds the genome
            names: a 1-based position or a header name, plus ``,noheader`` for
            a table without a header line.  Defaults to ``genome_name``.
        pairs_kwargs: Extra kwargs passed to the pairs reader (column
            indices, etc.), the same for every pairs file.
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
        qc: Tally contact-quality statistics while the pairs are being binned
            and store them in ``metadata["qc"]`` and the ``qc_*`` columns of
            ``/genomes`` and ``/contigs`` (see :mod:`metacooler.qc`).  Free
            apart from a few vectorised operations per chunk, since the pass
            over the pairs is happening anyway — and only possible *here*,
            because binning discards the sub-bin distances the short/long
            split needs.  The per-contig short/long counts are what lets a
            later :func:`metacooler.io.reassign.reassign` keep that split.
        qc_distance_threshold: Intra-contig short/long split for *qc*, in bp.
        ps: Accumulate a per-genome P(s) curve from the same pass and write it
            beside the output as ``<name>.ps.tsv`` (see :mod:`metacooler.ps`).
            Measured on read positions, so it resolves distances the binned map
            cannot — which is the whole reason to do it here rather than
            afterwards.  Costs one ``searchsorted`` and one ``bincount`` per
            chunk over the intra-contig rows.
        ps_bins_per_decade: Log-grid resolution of the P(s) curve.
        force: Overwrite an existing ``<name>.ps.tsv``.  The ``.metacool``
            itself is guarded by the caller (the CLI checks it before any work
            starts).
        nproc: Threads binning pairs batches while the file is being read and
            decompressed (which stays on one thread, being one gzip stream).
            Results are folded in file order, so the output and every tally
            are those of ``nproc=1``; memory grows by one batch per thread.
        keep_partial: On a failure while the output is being written, leave
            ``<uri>.partial`` behind for inspection, as the default does
            (logged at ERROR).  ``False`` removes it instead, whatever stopped
            the run (``KeyboardInterrupt`` included), and logs that at INFO.
            The output path itself is never written by a failed run, nor is
            the ``.ps.tsv`` table, staged as ``<name>.ps.tsv.partial`` and
            moved into place just before the map: the previous table is set
            aside (as ``<name>.ps.tsv.previous``) until the map is in place,
            and put back if the map cannot be.  Once the map is in place the
            set-aside table is removed; if that fails it is left, with a
            WARNING, and not put back: the new map and table match, and the
            run succeeds unless it was interrupted.  Only a process killed
            during those moves can leave them half done: the previous table
            still set aside, or the new table beside the previous map.

    Raises:
        ValueError: If *nproc* is below 1, no pairs file is given, the pairs
            files are not compatible (see *pairs*), or *contig_genome_map*
            names no contig of the assembly (checked before any pair is read).
        FileNotFoundError: If a pairs file does not exist -- checked before
            any is read, so a typo in the last of several libraries does not
            cost a binning pass over the others.
        FileExistsError: If ``<name>.ps.tsv.previous`` exists, with or
            without *ps*.  A run leaves it when it was killed while it moved
            its map and table into place, or could not remove that table
            afterwards, and which table matches the map is for the user to
            check; a file of the user's own of that name is refused alike.
            Checked before anything is read.
    """
    if nproc < 1:
        raise ValueError(f"nproc must be at least 1, got {nproc}")
    uri = Path(uri)
    if pairs_kwargs is None:
        pairs_kwargs = {}
    pairs_paths = _pairs_paths(pairs)
    missing = [str(p) for p in pairs_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Pairs file(s) not found: {', '.join(missing)}")
    if len(pairs_paths) > 1:
        logger.info(
            "Reading %d pairs files in one pass: %s",
            len(pairs_paths),
            ", ".join(map(str, pairs_paths)),
        )
    if contig_genome_map is not None and not isinstance(
        contig_genome_map, dict | str | Path
    ):
        raise TypeError("contig_genome_map must be a dict, a path, or None")
    if bins is not None and binsize is not None:
        # Passing both silently binned every contact arithmetically against a bin
        # table that was never tiled at *binsize*, scattering contacts across
        # unrelated fragments while still passing validate().
        raise ValueError(
            "binsize and bins are mutually exclusive: pass binsize to tile "
            "fixed-size bins, or bins to use a pre-defined bin table."
        )
    _check_no_binning(
        contig_genome_map,
        contig_table_format,
        genomes_table=genomes_table,
        drop_unassigned=drop_unassigned,
    )
    # A P(s) table still set aside: an earlier run was killed while it moved its
    # map and table into place, and which table matches the map is unknown.
    _refuse_a_set_aside_file(uri, (ps_table_path(uri),))
    # Guard the P(s) sidecar before the streaming pass, not after it: finding
    # out that the table cannot be written is worth seconds, not hours.
    ps_path = ps_table_path(uri) if ps else None
    if ps_path is not None:
        check_output_path(ps_path, force, defer_removal=True)
    inputs = [
        *pairs_paths,
        bins if isinstance(bins, str | Path) else None,
        fasta,
        contig_genome_map if isinstance(contig_genome_map, str | Path) else None,
        genomes_table,
    ]
    companions = () if ps_path is None else (ps_path,)
    _refuse_staging_an_input(uri, inputs, companions)
    timer = time.perf_counter()

    # --- Parse contig map and optional genomes metadata table ---
    cg_map, contig_df = _parse_contig_map(contig_genome_map, contig_table_format)
    extra_genomes_df = _load_extra_genomes_table(genomes_table, genomes_table_format)

    # --- Resolve contig lengths and initial bins DataFrame ---
    if isinstance(bins, str | Path):
        bins_df = _read_bins_bed(bins)
        contig_lengths = _read_contig_lengths(bins=bins_df)
        assembly_source = f"the bins file {bins}"
    elif isinstance(bins, pd.DataFrame):
        bins_df = bins.copy()
        contig_lengths = _read_contig_lengths(bins=bins_df)
        assembly_source = "the bins table"
    elif bins is None:
        # Integer binsize mode: need contig lengths from fasta or pairs header
        contig_lengths = _read_contig_lengths(fasta=fasta, pairs=pairs_paths)
        assembly_source = (
            f"the FASTA {fasta}" if fasta is not None else "the pairs header(s)"
        )
        bins_df = None
    else:
        raise TypeError("bins must be a path, DataFrame, or None")

    # --- Check the contig table against the assembly, before a pair is read ---
    # Against every contig the assembly names: non-positive lengths,
    # --drop-unassigned and --filter only narrow the set afterwards.  No table
    # at all is no binning yet, and has nothing to check.
    if contig_genome_map is not None:
        _check_contig_table_matches(
            cg_map,
            contig_lengths,
            source=contig_genome_map,
            contig_table_format=contig_table_format,
            assembly_source=assembly_source,
        )

    # --- Resolve the reference ("fasta") contig order, if one is requested ---
    # When lengths come from a bins table their dict order is whatever
    # groupby("chrom") produced (alphabetical), not the FASTA order.
    if contig_ordering == "fasta" and bins_df is not None:
        contig_lengths = _apply_reference_order(
            contig_lengths, fasta=fasta, pairs=pairs_paths
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
        _require_assigned_contigs(contig_lengths, n_before)

    # --- Apply --filter predicates before anything is tiled or binned ---
    contig_lengths, bins_df = _apply_create_filters(
        contig_lengths,
        cg_map,
        contig_df,
        extra_genomes_df,
        filters,
        bins_df,
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

    # Everything the tables were built from is dead now, and the binning pass
    # below is the memory-hungry one: `_buffer_pixels` accumulates up to
    # _SORT_IN_MEMORY_PIXEL_THRESHOLD rows and folds them, which transiently
    # doubles. Holding these on for its whole duration costs gigabytes on a
    # metagenome -- `contig_df` alone is the user's text table with one Python
    # str per cell (~1.5 GiB at 1.9 M contigs x 10 columns, measured).
    del contig_df, cg_map, extra_genomes_df, contig_lengths

    # --- Stream pairs, bin, sort, and write the cooler file ---
    logger.info(f"Starting to create cooler file at {uri}...")
    bin_keys: np.ndarray | None = (
        None
        if binsize is not None
        else _folded_bin_keys(chrom_offset, bins_df["start"].to_numpy())
    )
    drop_stats = _PairsDropStats()
    qc_tally = (
        ContactQC(
            len(genomes_df),
            qc_distance_threshold,
            source="pairs",
            unbinned_genome=unbinned_genome_id(genomes_df["name"]),
            n_contigs=len(contigs_df),
        )
        if qc
        else None
    )
    contig_genome_id = contigs_df["genome_id"].to_numpy(dtype=np.int64)
    contig_length = contigs_df["length"].to_numpy(dtype=np.int64)
    ps_edges = pairs_distance_grid(contig_length, ps_bins_per_decade) if ps else None
    ps_tally = PSCounts(len(genomes_df), ps_edges) if ps_edges is not None else None
    raw_pixels: Iterable[pd.DataFrame] = _stream_binned_pairs(
        pairs_paths,
        contig_name_to_id,
        chrom_offset,
        binsize,
        bin_keys,
        chunksize,
        pairs_kwargs,
        drop_stats,
        qc_tally,
        ps_tally,
        contig_genome_id if (qc or ps) else None,
        contig_length,
        nproc=nproc,
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
    pixel_source, was_empty = _ensure_nonempty_pixel_stream(pixel_source)
    ordered = ordered or was_empty
    # The map and its P(s) sidecar are staged together, and both moved into
    # place only once both are written and the map has validated.
    with _staged_output(
        uri, keep_partial=keep_partial, inputs=inputs, companions=companions
    ) as staged:
        write_cooler(
            staged,
            bins_df,
            pixel_source,
            binsize=binsize,
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

        # --- Contact-quality statistics ---
        # Written here rather than through write_cooler's `metadata`: the tally
        # is only complete once the pixel stream above has been consumed, which
        # happens inside that very call.
        if qc_tally is not None:
            store_qc(
                staged,
                qc_tally.summary(
                    binsize=binsize,
                    n_pairs_read=drop_stats.total,
                    n_pairs_discarded=drop_stats.dropped,
                ),
                qc_tally.genome_table(genomes_df["name"].to_numpy()),
                qc_tally.contig_table(contigs_df["name"].to_numpy()),
            )

        # --- Validate ---
        logger.info("Validating created metacool file")
        validate(staged)

        # --- P(s) sidecar: a companion file, not part of the .metacool ---
        if ps_tally is not None and ps_edges is not None and ps_path is not None:
            write_ps_table(
                _partial_path(ps_path),
                assemble(
                    genomes_df["name"].to_numpy(),
                    ps_edges,
                    ps_tally.matrix,
                    position_pair_denominators(
                        contig_length, contig_genome_id, len(genomes_df), ps_edges
                    ),
                ),
                force=True,
            )

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
    contig_genome_map: str | Path | None = None,
    contig_table_format: str = CONTIGS_FORMAT_DEFAULT,
    genomes_table: str | Path | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "cool", "alpha"
    ] = "genome-then-length",
    assembly_name: str | None = None,
    drop_unassigned: bool = False,
    filters: Sequence[str] | None = None,
    metadata: dict | None = None,
    nproc: int = 1,
    sort_in_memory: bool | None = None,
    qc: bool = True,
    qc_distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
    genomes_table_format: str = GENOMES_FORMAT_DEFAULT,
    keep_partial: bool = True,
) -> None:
    """Create a metacool file from an existing cool file.

    Args:
        cool_uri: Path to the source cool file.
        output: Destination path for the output metacool file.
        contig_genome_map: Path to contig-to-genome assignment table
            (``str`` or ``Path`` only; dicts are not accepted here).  Contigs
            it does not list are collected under ``__unbinned__``, and names
            that are not contigs of the cool file are ignored with a warning,
            but it must name at least one contig of the cool file.  ``None``
            (the default) builds a map with no binning yet: the cool file plus
            the genome layer, every contig in ``__unbinned__`` (bin them later
            with ``reassign``; ``contig_ordering="cool"`` keeps the cool
            file's own contig order, and makes that reassignment cheap).
        contig_table_format: Which columns of the table hold the contig and
            the genome names: ``CONTIG,GENOME[,noheader]``, each a 1-based
            position or a header name, or an alias (see
            :func:`metacooler.io.contigs.parse_contigs_table`).
        genomes_table: Optional path to genomes metadata table.
        contig_ordering: Strategy for ordering contigs within the metacool file.
            ``'cool'`` preserves the original contig order from the input cool file.
            Other choices: ``'genome-then-length'``, ``'genome-then-alpha'``,
            ``'alpha'``.
        assembly_name: Name for this assembly / metagenome, stored in file attrs.
        drop_unassigned: If ``True``, drop contigs absent from the contig table
            instead of collecting them under ``__unbinned__``.
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`)
            applied to the derived ``/genomes`` and ``/contigs`` tables before
            the bins are remapped.  A dropped contig keeps none of its bins,
            and every pixel touching one is discarded.  ``genomes.n_bins`` is
            not available here.
        metadata: Optional metadata dictionary to store in the file attrs.
        nproc: Accepted for backwards compatibility, but currently has **no
            effect**: HDF5 reads must be serialised, so pixel remapping runs in
            a single thread.  Values greater than 1 log a warning.
        sort_in_memory: If ``True``, collect all pixels in RAM and sort before
            writing (fast for small files).  If ``False``, pass the unsorted
            stream to cooler's external sort (scalable for large files).  If
            ``None`` (default), auto-detect based on pixel count: files with
            more than 100 M pixels use ``False`` with a logged warning.
        qc: Tally contact-quality statistics off the pixel stream that is
            being copied anyway, and store them in ``metadata["qc"]`` and the
            ``qc_*`` columns of ``/genomes`` and ``/contigs`` (see
            :mod:`metacooler.qc`).  A cool file has no read positions left, so
            intra-contig distances are measured between bin midpoints and the
            block records ``source="pixels"``.
        qc_distance_threshold: Intra-contig short/long split for *qc*, in bp.
        genomes_table_format: Which column of *genomes_table* holds the genome
            names: a 1-based position or a header name, plus ``,noheader`` for
            a table without a header line.  Defaults to ``genome_name``.
        keep_partial: On a failure while the output is being written, leave
            ``<output>.partial`` behind for inspection, as the default does
            (logged at ERROR).  ``False`` removes it instead, whatever stopped
            the run (``KeyboardInterrupt`` included), and logs that at INFO.
            The output path itself is never written by a failed run.

    Raises:
        ValueError: If ``cool_uri`` is not a valid cooler file, or
            *contig_genome_map* names no contig of it (checked before any
            pixel is read).
    """
    cool_path = Path(cool_uri)
    output_path = Path(output)
    inputs = [cool_uri, contig_genome_map, genomes_table]
    _refuse_staging_an_input(output_path, inputs)
    timer = time.perf_counter()
    _check_no_binning(
        contig_genome_map,
        contig_table_format,
        genomes_table=genomes_table,
        drop_unassigned=drop_unassigned,
    )

    # --- Validate and load cool file ---
    if not cooler.fileops.is_cooler(str(cool_path)):
        raise ValueError(f"Not a valid cooler file: {cool_path}")
    clr = cooler.Cooler(str(cool_path))
    # A metacool is symmetric-upper, and the remap below folds every pixel into
    # the upper triangle. On a "square" cooler -- which stores both triangles --
    # that sums each off-diagonal contact with its mirror, doubling it, and the
    # result still validates. Folding is only correct for a symmetric matrix,
    # and square storage is chosen precisely when one is not guaranteed.
    storage_mode = clr.info.get("storage-mode", "symmetric-upper")
    if storage_mode != "symmetric-upper":
        raise ValueError(
            f"{cool_path} uses storage-mode {storage_mode!r}; metacooler only "
            "imports symmetric-upper coolers (the cooler default). Re-create "
            "it with symmetric_upper=True, keeping the upper triangle."
        )
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

    # Build contig lengths from bins (use actual contigs in bins, not chroms table).
    # One vectorised max per contig, not one materialised sub-frame per contig:
    # ``clr.bins()[:]`` hands back ``chrom`` as a Categorical, so iterating the
    # groupby builds a DataFrame per contig (~14 us each, measured) and costs
    # half a minute on a metagenome before a single pixel has been read.
    # ``observed=True`` also drops the empty groups a zero-length chrom would
    # otherwise contribute as a NaN max, which ``int()`` cannot convert.
    chroms = clr.chroms()[:]
    # Every contig of the cool file, keyed by name in the file's order: the
    # "cool" ordering preserves it, and the contig table is checked against it.
    cool_contigs = dict.fromkeys(chroms["name"])
    _max_end = bins.groupby("chrom", observed=True)["end"].max()
    contig_lengths: dict[str, int] = {
        str(name): int(end)
        for name, end in zip(_max_end.index, _max_end.to_numpy(), strict=True)
    }
    del _max_end
    contig_lengths, _ = _drop_nonpositive_contigs(contig_lengths)

    # --- Parse contig assignments ---
    contig_map, contig_df = _parse_contig_map(contig_genome_map, contig_table_format)
    # Against every contig of the cool file, before --drop-unassigned and
    # --filter narrow them, and before a single pixel is read.  No table at all
    # is no binning yet, and has nothing to check.
    if contig_genome_map is not None:
        _check_contig_table_matches(
            contig_map,
            cool_contigs,
            source=contig_genome_map,
            contig_table_format=contig_table_format,
            assembly_source=f"the cool file {cool_path}",
        )
    extra_genomes_df = _load_extra_genomes_table(genomes_table, genomes_table_format)

    # --- Filter unassigned contigs ---
    if drop_unassigned:
        n_before = len(contig_lengths)
        contig_lengths = {
            ct: ln
            for ct, ln in contig_lengths.items()
            if ct in contig_map and contig_map[ct] != "__unbinned__"
        }
        _require_assigned_contigs(contig_lengths, n_before)

    # --- Apply --filter predicates before the bins are remapped ---
    contig_lengths, _ = _apply_create_filters(
        contig_lengths, contig_map, contig_df, extra_genomes_df, filters
    )

    # --- Build contig order ---
    if contig_ordering == "cool":
        # Preserve the original contig order from the cool file's chroms table,
        # filtered to only contigs present in contig_lengths (after drop_unassigned).
        contig_order = [name for name in cool_contigs if name in contig_lengths]
    else:
        contig_order = _build_contig_order(
            contig_lengths=contig_lengths,
            contig_genome_map=contig_map,
            ordering=contig_ordering,
        )

    # --- Build genome / contig tables ---
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

    # Free the sources the tables were derived from before the pixel stream, for
    # the reason given in from_pairs: the streaming pass is where the memory is
    # wanted, and none of these is read again.
    del bins, chroms, cool_contigs, contig_df, contig_map, extra_genomes_df
    del contig_lengths, contig_order

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
    qc_tally = (
        ContactQC(
            len(genomes),
            qc_distance_threshold,
            source="pixels",
            unbinned_genome=unbinned_genome_id(genomes["name"]),
            n_contigs=len(contigs),
        )
        if qc
        else None
    )
    if qc_tally is not None:
        raw_pixels = _tally_pixel_qc(
            raw_pixels, qc_tally, bin_maps(bins_remapped, contigs)
        )
    if sort_in_memory:
        # Contigs are only reordered -- bins map one-to-one -- so the remapped
        # pixels stay unique and this is one argsort on the pixel key: the
        # reduction after it finds no repeated key and returns the counts as
        # they are.  It is kept because a cooler written without cooler's
        # duplicate check can repeat a pixel, which is then summed.
        pixel_source: Iterable[pd.DataFrame] = _sort_and_dedup_pixels(
            pixel_iter=raw_pixels
        )
        ordered = True
    else:
        pixel_source: Iterable[pd.DataFrame] = raw_pixels
        ordered = False
    pixel_source, was_empty = _ensure_nonempty_pixel_stream(pixel_source)
    ordered = ordered or was_empty
    with _staged_output(
        output_path,
        keep_partial=keep_partial,
        inputs=inputs,
    ) as staged:
        write_cooler(
            staged,
            bins_remapped,
            pixel_source,
            ordered=ordered,
            symmetric_upper=True,
            # Keep the source's assembly name unless the caller overrides it.
            assembly=assembly_name or _source_assembly(clr) or "imported",
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

        # --- Contact-quality statistics ---
        # Complete only now: the tally rides on the pixel stream that
        # write_cooler above consumed.
        if qc_tally is not None:
            store_qc(
                staged,
                qc_tally.summary(binsize=infer_binsize(bins_remapped)),
                qc_tally.genome_table(genomes["name"].to_numpy()),
                qc_tally.contig_table(contigs["name"].to_numpy()),
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
