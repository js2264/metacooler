"""Create metacool files from pairs, bins, and MAG assignments."""

from __future__ import annotations

import contextlib
import gzip
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import h5py
import numpy as np
import pandas as pd
from pyfaidx import Fasta

from metacooler.io._utils import _finalize_metacool_file
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.genomes import parse_genomes_table
from metacooler.schema import (
    CONTIGS_DTYPE,
    GENOMES_DTYPE,
    validate,
)

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal


_UNBINNED_SENTINELS: frozenset[str] = frozenset(
    {"__unbinned__", "na", "nan", "n/a", "nd", "n/d", ""}
)

logger = logging.getLogger("metacooler")


def _read_contig_lengths(
    bins: pd.DataFrame | None = None,
    fasta: str | Path | None = None,
    pairs: str | Path | None = None,
) -> dict[str, int]:
    """Extract contig lengths from bins, fasta, or pairs header."""
    if fasta is not None:
        fasta_path = Path(fasta)
        if fasta_path.suffix == ".gz":
            # pyfaidx doesn't handle plain gzip; decompress to temp
            import tempfile

            with (
                tempfile.NamedTemporaryFile(
                    mode="wb", suffix=".fa", delete=False
                ) as tmp,
                gzip.open(fasta_path, "rb") as gz,
            ):
                tmp.write(gz.read())
            fa_path = Path(tmp.name)
            fai_path = fa_path.parent / (fa_path.name + ".fai")
            try:
                fa = Fasta(tmp.name)
                result = {name: len(record) for name, record in fa.items()}
            finally:
                with contextlib.suppress(Exception):
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


def _parse_pairs_header(pairs_path: str | Path) -> dict[str, int] | None:
    """Parse #chromsize: lines from a pairs file header.

    Returns a dict {contig_name: length} or None if no chromsizes found.
    """
    path = Path(pairs_path)
    opener = gzip.open if path.suffix == ".gz" else open
    lengths = {}
    with opener(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("#chromsize:"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    lengths[parts[1]] = int(parts[2])
    return lengths if lengths else None


def _build_contig_order(
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    ordering: Literal[
        "genome-then-alpha", "genome-then-length", "fasta", "alpha", "custom"
    ] = "genome-then-alpha",
    custom_order: list[str] | None = None,
) -> list[str]:
    """Return contig names in canonical order."""
    all_contigs = set(contig_lengths.keys())

    if ordering == "fasta":
        # Preserve dict insertion order (Python 3.7+)
        return list(contig_lengths.keys())

    if ordering == "alpha":
        return sorted(all_contigs)

    if ordering == "custom":
        if custom_order is None:
            raise ValueError("custom_order required when ordering='custom'")
        missing = all_contigs - set(custom_order)
        if missing:
            raise ValueError(f"custom_order missing contigs: {missing}")
        return custom_order

    if ordering == "genome-then-alpha":
        # Group by genome, then alphabetical order within each genome
        genome_groups: dict[str, list[str]] = {}
        for contig in all_contigs:
            genome = contig_genome_map.get(contig, "__unbinned__")
            genome_groups.setdefault(genome, []).append(contig)

        ordered = []
        for genome in sorted(genome_groups.keys()):
            contigs = genome_groups[genome]
            contigs.sort()
            ordered.extend(contigs)
        return ordered

    if ordering == "genome-then-length":
        # Group by genome, then longest contig first within each genome
        genome_groups2: dict[str, list[str]] = {}
        for contig in all_contigs:
            genome = contig_genome_map.get(contig, "__unbinned__")
            genome_groups2.setdefault(genome, []).append(contig)

        ordered2 = []
        for genome in sorted(genome_groups2.keys()):
            contigs = genome_groups2[genome]
            contigs.sort(key=lambda c: (-contig_lengths[c], c))
            ordered2.extend(contigs)
        return ordered2

    raise ValueError(f"Unknown ordering: {ordering!r}")


def _normalize_genome_name(v: str) -> str:
    """Normalize sentinel values (empty, NA, nan, ND, n/d, …) to ``__unbinned__``."""
    return "__unbinned__" if v.strip().lower() in _UNBINNED_SENTINELS else v


def _build_genomes_table(
    contig_order: list[str],
    contig_lengths: dict[str, int],
    contig_genome_map: dict[str, str],
    extra_genomes_df: pd.DataFrame | None = None,
    keep_unbinned: bool = True,
) -> pd.DataFrame:
    """Build the genomes DataFrame from contig data.

    Args:
        contig_order: Ordered list of contig names.
        contig_lengths: Mapping of contig name → length.
        contig_genome_map: Mapping of contig name → genome name.
        extra_genomes_df: Optional DataFrame from a genomes table file
            (column ``genome_id`` plus any extra columns).  Extra columns
            are merged into the output by genome name.
    """
    # Collect genome names
    genome_names = set()
    for contig in contig_order:
        genome_names.add(contig_genome_map.get(contig, "__unbinned__"))

    # Only inject __unbinned__ when at least one contig actually lands in it
    has_unbinned = any(
        contig_genome_map.get(c, "__unbinned__") == "__unbinned__" for c in contig_order
    )
    if keep_unbinned and has_unbinned and "__unbinned__" not in genome_names:
        genome_names.add("__unbinned__")

    # Sort: __unbinned__ first (if present), then alphabetical
    sorted_names = sorted(genome_names - {"__unbinned__"})
    if "__unbinned__" in genome_names:
        sorted_names = ["__unbinned__"] + sorted_names

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

    # Merge extra columns from the genomes file if provided
    if extra_genomes_df is not None and len(extra_genomes_df.columns) > 1:
        extra_cols = [c for c in extra_genomes_df.columns if c != "genome_id"]
        extra = extra_genomes_df.rename(columns={"genome_id": "name"})[
            ["name"] + extra_cols
        ].copy()
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
    """Build the contigs DataFrame with genome_id and optional extra columns."""
    genome_name_to_id = {name: i for i, name in enumerate(genomes["name"])}

    df = pd.DataFrame(
        {
            "name": contig_order,
            "length": [contig_lengths[c] for c in contig_order],
            "genome_id": [
                genome_name_to_id.get(
                    contig_genome_map.get(c, "__unbinned__"),
                    genome_name_to_id.get("__unbinned__", 0),
                )
                for c in contig_order
            ],
        }
    )

    # Merge additional columns from the input contig table if provided
    if contig_df is not None and len(contig_df.columns) > 2:
        # Identify extra columns beyond contig_name and genome_id
        extra_cols = [
            c for c in contig_df.columns if c not in ("contig_name", "genome_id")
        ]
        if extra_cols:
            rename_map = {"contig_name": "name"}
            extra = contig_df.rename(columns=rename_map)[["name"] + extra_cols].copy()
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


def _build_bins(
    bins: str | Path | pd.DataFrame,
    contig_order: list[str],
    contig_name_to_id: dict[str, int],
) -> pd.DataFrame:
    """Read a BED file (or DataFrame) of bins and remap chrom names to contig ids."""
    if isinstance(bins, str | Path):
        df = pd.read_csv(
            bins,
            sep="\t",
            header=None,
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": np.int32, "end": np.int32},
        )
    else:
        df = bins.copy()
    # Map chrom name -> contig id
    df["chrom"] = df["chrom"].map(contig_name_to_id)
    if df["chrom"].isna().any():
        missing = df[df["chrom"].isna()]["chrom"].unique()
        raise ValueError(f"Bins reference unknown contigs: {missing}")
    return df[["chrom", "start", "end"]]


def _tile_bins(
    contig_order: list[str],
    contig_lengths: dict[str, int],
    contig_name_to_id: dict[str, int],
    binsize: int,
) -> pd.DataFrame:
    """Generate fixed-size bins tiling each contig."""
    parts = []
    for contig in contig_order:
        cid = contig_name_to_id[contig]
        length = contig_lengths[contig]
        starts = np.arange(0, length, binsize, dtype=np.int32)
        ends = np.minimum(starts + binsize, length).astype(np.int32)
        parts.append(pd.DataFrame({"chrom": cid, "start": starts, "end": ends}))
    df = pd.concat(parts, ignore_index=True)
    return df.astype({"chrom": np.int32, "start": np.int32, "end": np.int32})


def _count_pairs_header(path: Path) -> int:
    """Count the number of ``#``-prefixed header lines at the top of a pairs file."""
    import bz2 as _bz2
    import gzip as _gzip

    if path.suffix in (".gz", ".bgz"):
        open_fn = _gzip.open
    elif path.suffix == ".bz2":
        open_fn = _bz2.open
    else:
        open_fn = open
    n = 0
    with open_fn(str(path), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                n += 1
            else:
                break
    return n


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

    # Count header comment lines; pyarrow has no native comment-skip support
    n_header = _count_pairs_header(path)

    # Column indices are 1-based; convert to 0-based.
    # pyarrow auto-generates column names "f0", "f1", ... when
    # autogenerate_column_names=True, so we address columns by those names.
    pa_col_names = [f"f{c1 - 1}", f"f{p1 - 1}", f"f{c2 - 1}", f"f{p2 - 1}"]
    final_col_names = ["chrom1", "pos1", "chrom2", "pos2"]

    # Aim for approximately chunksize rows per batch; assume ~120 bytes/row.
    block_size = max(chunksize * 120, 8 * 1024 * 1024)

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

    # pa.input_stream auto-decompresses .gz/.bgz via native C++ gzip
    src = pa.input_stream(str(path))
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


def _bin_pairs_chunk(
    chunk: pd.DataFrame,
    contig_name_to_id: dict[str, int],
    chrom_offset: np.ndarray,
    binsize: int | None,
    bin_starts: np.ndarray | None,
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

    Returns:
        Aggregated pixel DataFrame or ``None`` if the chunk is empty after
        filtering.
    """
    # Use pd.Categorical for fast contig-name → integer-ID encoding.
    # This mirrors cooler's _sanitize_records approach: unrecognised contigs
    # receive code == -1 and are dropped in one vectorised pass.
    contig_list = list(contig_name_to_id.keys())
    c1_cat = pd.Categorical(chunk["chrom1"], categories=contig_list, ordered=True)
    c2_cat = pd.Categorical(chunk["chrom2"], categories=contig_list, ordered=True)
    c1_codes = c1_cat.codes  # int8/int16 ndarray, -1 for unknown
    c2_codes = c2_cat.codes

    # Drop pairs referencing unknown contigs in one boolean mask
    valid = (c1_codes >= 0) & (c2_codes >= 0)
    if not valid.any():
        return None

    c1 = c1_codes[valid].astype(np.int32)
    c2 = c2_codes[valid].astype(np.int32)
    # Convert 1-based pairs positions to 0-based, clamp at 0
    p1 = np.maximum(chunk["pos1"].values[valid] - 1, 0)
    p2 = np.maximum(chunk["pos2"].values[valid] - 1, 0)

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


def _total_time_str(start_time: pd.Timestamp) -> str:
    tottime = pd.Timestamp.now() - start_time
    if tottime.total_seconds() >= 3600:
        return f"{int(tottime.total_seconds() // 3600)}h:{int((tottime.total_seconds() % 3600) // 60)}m:{int(tottime.total_seconds() % 60)}s"
    elif tottime.total_seconds() >= 60:
        return f"{int(tottime.total_seconds() // 60)}m:{int(tottime.total_seconds() % 60)}s"
    else:
        return f"{tottime.total_seconds():.1f}s"


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
    all_pixels = pd.concat(frames, ignore_index=True)
    del frames
    all_pixels = all_pixels.groupby(["bin1_id", "bin2_id"], sort=True, as_index=False)[
        "count"
    ].sum()
    logger.info(f"Total pixels after global dedup: {len(all_pixels):,}")
    for lo in range(0, len(all_pixels), write_chunk):
        yield all_pixels.iloc[lo : lo + write_chunk]


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
        for chunk in _stream_pairs(pairs, chunksize=chunksize, **pairs_kwargs):
            progress.update(task, advance=len(chunk))
            result = _bin_pairs_chunk(
                chunk, contig_name_to_id, chrom_offset, binsize, bin_starts
            )
            if result is not None and not result.empty:
                yield result


def _stream_remapped_pixels(
    cool_path_str: str,
    n_pixels: int,
    old_to_new_bin: np.ndarray,
    nproc: int,
    chunksize: int = 1_000_000,
) -> Iterable[pd.DataFrame]:
    """Stream pixels from a cool file with bin-ID remapping and rich progress.

    Args:
        cool_path_str: Filesystem path to the source cool file.
        n_pixels: Total pixel count (used for progress bar).
        old_to_new_bin: Array mapping old bin indices to new bin indices.
        nproc: Number of parallel reader threads (1 = serial).
        chunksize: Pixels per read batch.

    Yields:
        Remapped pixel DataFrames with columns ``bin1_id``, ``bin2_id``,
        ``count``.
    """
    import threading
    from concurrent.futures import ThreadPoolExecutor

    from rich.console import Console
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TaskProgressColumn,
        TextColumn,
        TimeRemainingColumn,
    )

    thread_local = threading.local()

    def _remap_chunk(lo: int, hi: int) -> pd.DataFrame:
        if not hasattr(thread_local, "clr"):
            thread_local.clr = cooler.Cooler(cool_path_str)
        chunk = thread_local.clr.pixels()[lo:hi].copy()
        chunk["bin1_id"] = old_to_new_bin[chunk["bin1_id"].values.astype(np.int64)]
        chunk["bin2_id"] = old_to_new_bin[chunk["bin2_id"].values.astype(np.int64)]
        # Drop pixels that touch a dropped contig (sentinel -1).
        chunk = chunk[(chunk["bin1_id"] >= 0) & (chunk["bin2_id"] >= 0)].copy()
        chunk["count"] = chunk["count"].astype(np.int32)
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
        if nproc == 1:
            for lo, hi in ranges:
                yield _remap_chunk(lo, hi)
                progress.update(task, advance=1)
        else:
            with ThreadPoolExecutor(max_workers=nproc) as pool:
                futures = [pool.submit(_remap_chunk, lo, hi) for lo, hi in ranges]
                for f in futures:
                    yield f.result()
                    progress.update(task, advance=1)


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
    contig_table_format: Literal["contig_genome", "metator", "contigs"] = "contigs",
    fasta: str | Path | None = None,
    assembly_name: str | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "fasta", "alpha", "custom"
    ] = "genome-then-length",
    custom_order: list[str] | None = None,
    drop_unassigned: bool = False,
    binsize: int | None = None,
    genomes_table: str | Path | None = None,
    pairs_kwargs: dict | None = None,
    metadata: dict | None = None,
    chunksize: int = 1_000_000,
) -> None:
    """Create a ``.metacool`` file from pairs, bins, and contig assignments.

    Args:
        uri: Output file path (``.metacool`` extension recommended).
        bins: Path to a BED file of bins, or a pre-built DataFrame
            (``chrom, start, end``).  If ``None``, *binsize* must be given
            and bins will be tiled from the fasta.
        pairs: Path to a pairs file (plain text, gzipped, or bgzipped).
        contig_genome_map: Either a ``dict[contig_name → genome_name]`` or
            a path to a contig table file.
        contig_table_format: Format when *contig_genome_map* is a file path.
        fasta: Assembly fasta (required if *bins* is ``None``).
        assembly_name: Name for this assembly / metagenome.
        contig_ordering: Canonical ordering strategy.
        custom_order: Explicit contig order when *contig_ordering* is
            ``"custom"``.
        drop_unassigned: If ``True``, drop contigs not present in the contig
            table instead of assigning them to ``__unbinned__``.
        binsize: Fixed bin size (used only when *bins* is ``None``).
        genomes_table: Optional path to a genomes metadata TSV **with a header
            row**. The first column must be named ``genome_name``; any additional
            columns are stored in the ``/genomes`` table.
        pairs_kwargs: Extra kwargs passed to the pairs reader (column
            indices, etc.).
        metadata: Additional metadata dict stored in file attrs.
        chunksize: Number of pairs to process per chunk.
    """
    uri = Path(uri)
    if pairs_kwargs is None:
        pairs_kwargs = {}
    if pairs is None:
        raise ValueError("pairs file path required")
    if not isinstance(contig_genome_map, dict | str | Path):
        raise TypeError("contig_genome_map must be a dict or a path")
    timer = pd.Timestamp.now()

    # --- Parse contig map and optional genomes metadata table ---
    cg_map, contig_df = _parse_contig_map(contig_genome_map, contig_table_format)
    extra_genomes_df = _load_extra_genomes_table(genomes_table)

    # --- Resolve contig lengths and initial bins DataFrame ---
    if isinstance(bins, str | Path):
        bins_df = pd.read_csv(
            bins,
            sep="\t",
            header=None,
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": np.int32, "end": np.int32},
        )
        contig_lengths = bins_df.groupby("chrom")["end"].max().to_dict()
    elif isinstance(bins, pd.DataFrame):
        bins_df = bins.copy()
        contig_lengths = bins_df.groupby("chrom")["end"].max().to_dict()
    elif bins is None:
        # Integer binsize mode: need contig lengths from fasta or pairs header
        contig_lengths = _read_contig_lengths(fasta=fasta, pairs=pairs)
        bins_df = None
    else:
        raise TypeError("bins must be a path, DataFrame, or None")
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
    contig_order = _build_contig_order(
        contig_lengths, cg_map, ordering=contig_ordering, custom_order=custom_order
    )
    logger.debug(f"Contig order: {contig_order[:10]} ... {contig_order[-10:]}")

    genomes_df = _build_genomes_table(
        contig_order,
        contig_lengths,
        cg_map,
        extra_genomes_df,
        keep_unbinned=not drop_unassigned,
    )
    contig_name_to_id = {name: i for i, name in enumerate(contig_order)}
    contigs_df = _build_contigs_table(
        contig_order, contig_lengths, cg_map, genomes_df, contig_df
    )
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
        bins_df = _build_bins(bins_df, contig_order, contig_name_to_id)
    bins_df = bins_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    logger.info(
        f"Built final bins table with {len(bins_df)} bins tiling "
        f"{len(contig_lengths)} contigs"
    )

    # --- Compute chrom / genome offsets and insert n_bins into genomes_df ---
    chrom_offset, genome_offset = _compute_all_offsets(bins_df, contigs_df, genomes_df)

    # --- Stream pairs, bin, sort, and write the cooler file ---
    bin_starts: np.ndarray | None = (
        None if binsize is not None else bins_df["start"].values.astype(np.int64)
    )
    raw_pixels = _stream_binned_pairs(
        pairs,
        contig_name_to_id,
        chrom_offset,
        binsize,
        bin_starts,
        chunksize,
        pairs_kwargs,
    )

    logger.info(f"Starting to create cooler file at {uri}...")
    cooler.create_cooler(
        str(uri),
        bins=bins_df,
        pixels=_sort_and_dedup_pixels(raw_pixels),
        ordered=True,
        symmetric_upper=True,
        assembly=assembly_name,
        metadata=metadata or {},
        dupcheck=False,
    )

    # --- Post-process: add metacooler-specific groups and attrs ---
    logger.info("Post-processing cooler file to add metacool structure")
    _finalize_metacool_file(
        uri,
        genomes_df,
        contigs_df,
        chrom_offset,
        genome_offset,
        contig_ordering=contig_ordering,
        assembly_name=assembly_name,
    )

    # --- Validate ---
    logger.info("Validating created metacool file")
    validate(uri)
    with h5py.File(uri, "r") as f:
        nnz = len(f["pixels"]["count"])
        total = int(f["pixels"]["count"][:].sum())
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
    contig_table_format: Literal[
        "contigs", "contig_genome", "metator", "I,J"
    ] = "contigs",
    genomes_table: str | Path | None = None,
    contig_ordering: Literal[
        "genome-then-length", "genome-then-alpha", "cool", "alpha", "custom"
    ] = "genome-then-length",
    assembly_name: str | None = None,
    drop_unassigned: bool = False,
    metadata: dict | None = None,
    nproc: int = 1,
) -> None:
    """Create a metacool file from an existing cool file.

    Args:
        cool_uri: Path to the source cool file.
        output: Destination path for the output metacool file.
        contig_genome_map: Path to contig-to-genome assignment table
            (``str`` or ``Path`` only; dicts are not accepted here).
        contig_table_format: Format of the contig assignment table.
        genomes_table: Optional path to genomes metadata table.
        contig_ordering: Strategy for ordering contigs within the metacool file.
            ``'cool'`` preserves the original contig order from the input cool file.
            Other choices: ``'genome-then-length'``, ``'genome-then-alpha'``,
            ``'alpha'``, ``'custom'``.
        assembly_name: Name for this assembly / metagenome, stored in file attrs.
        drop_unassigned: If ``True``, drop contigs absent from the contig table
            instead of collecting them under ``__unbinned__``.
        metadata: Optional metadata dictionary to store in the file attrs.
        nproc: Number of parallel reader threads for pixel remapping.
            Values greater than 1 use a ``ThreadPoolExecutor`` to overlap
            HDF5 reads across chunks.

    Raises:
        ValueError: If ``cool_uri`` is not a valid cooler file.
    """
    cool_path = Path(cool_uri)
    output_path = Path(output)
    timer = pd.Timestamp.now()

    # --- Validate and load cool file ---
    if not cooler.fileops.is_cooler(str(cool_path)):
        raise ValueError(f"Not a valid cooler file: {cool_path}")
    clr = cooler.Cooler(str(cool_path))
    bins = pd.DataFrame(clr.bins()[:])
    logger.info(f"Loaded cool file: {len(bins)} bins, {clr.info['nnz']} pixels")

    # Build contig lengths from bins (use actual contigs in bins, not chroms table)
    chroms = clr.chroms()[:]
    contig_id_to_name = dict(enumerate(chroms["name"]))
    contig_lengths: dict[str, int] = {
        str(chrom): int(group["end"].max()) for chrom, group in bins.groupby("chrom")
    }

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
        cool_order = [
            name for name in contig_id_to_name.values() if name in contig_lengths
        ]
        contig_order = _build_contig_order(
            contig_lengths=contig_lengths,
            contig_genome_map=contig_map,
            ordering="custom",
            custom_order=cool_order,
        )
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

    # --- Remap bins to new contig IDs; compute offsets ---
    bins_remapped, old_to_new_bin = _remap_bins(bins, contig_order)
    chrom_offset, genome_offset = _compute_all_offsets(bins_remapped, contigs, genomes)

    # --- Stream, remap, sort, and write ---
    logger.info(
        f"Building metacool file: {len(genomes)} genomes, {len(contigs)} contigs"
    )
    raw_pixels = _stream_remapped_pixels(
        str(cool_path), clr.info["nnz"], old_to_new_bin, nproc
    )
    cooler.create_cooler(
        str(output_path),
        bins=bins_remapped,
        pixels=_sort_and_dedup_pixels(raw_pixels),
        ordered=True,
        symmetric_upper=True,
        assembly=assembly_name or "imported",
        metadata=metadata or {},
    )

    # --- Add metacooler-specific groups and attributes ---
    _finalize_metacool_file(
        output_path,
        genomes,
        contigs,
        chrom_offset,
        genome_offset,
        contig_ordering=contig_ordering,
        assembly_name=assembly_name,
    )

    logger.info(
        f"Created metacool file: {output_path} ({_total_time_str(timer)} elapsed)"
    )
