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

from metacooler._version import __version__
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.genomes import parse_genomes_table
from metacooler.schema import (
    CONTIGS_DTYPE,
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    GENOMES_DTYPE,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
    validate,
)

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal


_UNBINNED_SENTINELS: frozenset[str] = frozenset(
    {"__unbinned__", "na", "nan", "n/a", "nd", "n/d", ""}
)

logger = logging.getLogger("metacooler")


def _str_to_hdf5(arr: np.ndarray) -> np.ndarray:
    """Convert a string/object array to fixed-length bytes for HDF5 storage."""
    if arr.dtype.kind in ("O", "U", "S") or str(arr.dtype) == "string":
        str_arr = [
            "" if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v)
            for v in arr
        ]
        max_len = max((len(s) for s in str_arr), default=1)
        return np.array([s.encode("utf-8") for s in str_arr], dtype=f"S{max_len}")
    return np.asarray(arr)


def _finalize_metacool_file(
    path: Path,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    chrom_offset: np.ndarray | None,
    genome_offset: np.ndarray,
    *,
    contig_ordering: str = "genome-then-alpha",
    assembly_name: str | None = None,
) -> None:
    """Add metacooler-specific groups and attributes to a cooler file.

    Handles: chroms soft-link, genomes/contigs HDF5 groups (with byte-string
    encoding), stripping the enum dtype from ``bins/chrom``, writing index
    offsets, and setting all format attributes.

    Args:
        path: Path to the cooler file to finalize (opened in append mode).
        genomes: ``/genomes`` DataFrame.
        contigs: ``/contigs`` DataFrame.
        chrom_offset: Pre-computed ``chrom_offset`` array (length n_contigs+1).
            If ``None``, computed from the ``bins/chrom`` dataset.
        genome_offset: Pre-computed ``genome_offset`` array (length n_genomes+1).
        contig_ordering: Value for the ``contig-ordering`` HDF5 attribute.
        assembly_name: Optional assembly name stored in ``assembly-name`` attr.
    """
    with h5py.File(path, "a") as f:
        # Recreate chroms as a soft link to contigs (cooler compat)
        if "chroms" in f:
            del f["chroms"]
        f["chroms"] = h5py.SoftLink("contigs")

        # Write /genomes
        if "genomes" in f:
            del f["genomes"]
        g = f.create_group("genomes", track_order=True)
        for col in genomes.columns:
            g.create_dataset(col, data=_str_to_hdf5(genomes[col].values))

        # Write /contigs
        if "contigs" in f:
            del f["contigs"]
        c = f.create_group("contigs", track_order=True)
        for col in contigs.columns:
            c.create_dataset(col, data=_str_to_hdf5(contigs[col].values))

        # Strip cooler's enum dtype from bins/chrom so it stores plain int32
        chrom_data = f["bins"]["chrom"][:]
        del f["bins"]["chrom"]
        f["bins"].create_dataset("chrom", data=chrom_data, dtype=np.int32)

        # Compute chrom_offset from bins if not provided
        if chrom_offset is None:
            n_contigs = len(contigs)
            contig_counts = np.bincount(
                chrom_data.astype(np.int64), minlength=n_contigs
            )
            chrom_offset = np.zeros(n_contigs + 1, dtype=np.int64)
            chrom_offset[1:] = np.cumsum(contig_counts)

        # Write index offsets
        idx = f.require_group("indexes")
        for key, data in (
            ("chrom_offset", chrom_offset),
            ("genome_offset", genome_offset),
        ):
            if key in idx:
                del idx[key]
            idx.create_dataset(key, data=data.astype(np.int64))

        # Format attributes
        f.attrs["format"] = COOLER_FORMAT
        f.attrs["format-url"] = COOLER_FORMAT_URL
        f.attrs["metacooler-format"] = METACOOLER_FORMAT
        f.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        f.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        f.attrs["generated-by"] = f"metacooler-{__version__}"
        f.attrs["ngenomes"] = len(genomes)
        f.attrs["contig-ordering"] = contig_ordering
        if assembly_name:
            f.attrs["assembly-name"] = assembly_name


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


def _stream_pairs(
    pairs_path: str | Path,
    c1: int = 2,
    p1: int = 3,
    c2: int = 4,
    p2: int = 5,
    chunksize: int = 1_000_000,
) -> Iterable[pd.DataFrame]:
    """Stream a pairs file in chunks.

    Uses ``usecols`` to avoid parsing unneeded columns and reads position
    columns directly as ``int64`` to skip the ``str → int`` conversion.
    """
    path = Path(pairs_path)
    compression = "infer"
    if path.suffix == ".gz":
        compression = "gzip"
    elif path.suffix == ".bz2":
        compression = "bz2"

    # Column indices are 1-based in CLI, convert to 0-based
    usecols_0based = [c1 - 1, p1 - 1, c2 - 1, p2 - 1]
    col_names = ["chrom1", "pos1", "chrom2", "pos2"]

    # Per-column dtypes: chrom columns as str, position columns as int64.
    # This avoids reading everything as str and converting afterwards.
    col_dtypes = {
        c1 - 1: str,
        p1 - 1: np.int64,
        c2 - 1: str,
        p2 - 1: np.int64,
    }

    reader = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=usecols_0based,
        dtype=col_dtypes,
        comment="#",
        compression=compression,
        chunksize=chunksize,
        on_bad_lines="warn",
    )
    for chunk in reader:
        chunk.columns = col_names
        yield chunk


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
        # Fixed-size bins: direct arithmetic
        bin1 = chrom_offset[c1] + p1 // binsize
        bin2 = chrom_offset[c2] + p2 // binsize
    else:
        # Fragment-level bins: binary search within each contig's bin starts
        assert bin_starts is not None, "bin_starts required for fragment-level binning"
        bin1 = np.empty(len(c1), dtype=np.int64)
        bin2 = np.empty(len(c2), dtype=np.int64)
        for cid in np.unique(c1):
            mask = c1 == cid
            lo, hi = int(chrom_offset[cid]), int(chrom_offset[cid + 1])
            bin1[mask] = (
                lo + np.searchsorted(bin_starts[lo:hi], p1[mask], side="right") - 1
            )
        for cid in np.unique(c2):
            mask = c2 == cid
            lo, hi = int(chrom_offset[cid]), int(chrom_offset[cid + 1])
            bin2[mask] = (
                lo + np.searchsorted(bin_starts[lo:hi], p2[mask], side="right") - 1
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


def from_pairs(
    uri: str | Path,
    bins: str | Path | pd.DataFrame | None = None,
    pairs: str | Path | None = None,
    contig_genome_map: dict[str, str] | str | Path | None = None,
    contig_table_format: Literal["contig_genome", "metator", "contigs"] = "contigs",
    fasta: str | Path | None = None,
    assembly_name: str | None = None,
    contig_ordering: Literal[
        "genome-then-alpha", "genome-then-length", "fasta", "alpha", "custom"
    ] = "genome-then-alpha",
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

    # --- Parse contig table ---
    contig_df: pd.DataFrame | None = None
    if isinstance(contig_genome_map, str | Path):
        cg_map, contig_df = parse_contigs_table(
            contig_genome_map, fmt=contig_table_format
        )
        logger.info(
            f"Parsed contig-genome map for {len(cg_map)} contigs from {contig_genome_map}"
        )
    elif isinstance(contig_genome_map, dict):
        cg_map: dict[str, str] = contig_genome_map
        logger.info(f"Parsed contig-genome map for {len(cg_map)} contigs from dict")

    # Normalize empty/NA genome names to __unbinned__
    n_contigs0: int = len(cg_map)
    cg_map = {k: _normalize_genome_name(v) for k, v in cg_map.items()}
    n_contigs1: int = len(cg_map)
    if n_contigs0 != n_contigs1:
        logger.info(
            f"Normalized {n_contigs0 - n_contigs1} contig genome names to __unbinned__"
        )

    # --- Parse optional genomes metadata table ---
    extra_genomes_df: pd.DataFrame | None = None
    if genomes_table is not None:
        extra_genomes_df = parse_genomes_table(genomes_table)
        logger.info(
            f"Parsed genomes table for {len(extra_genomes_df)} genomes from {genomes_table}"
        )

    # --- Get contig lengths and build bins ---
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

    # --- Filter unassigned ---
    if drop_unassigned:
        n_contigs0: int = len(contig_lengths)
        contig_lengths = {
            ct: ls
            for ct, ls in contig_lengths.items()
            if ct in cg_map and cg_map[ct] != "__unbinned__"
        }
        if bins_df is not None:
            bins_df = bins_df[bins_df["chrom"].isin(contig_lengths)].copy()
        n_contigs1: int = len(contig_lengths)
        logger.info(
            f"Dropped {n_contigs0 - n_contigs1} unassigned contigs; {n_contigs1} remain"
        )

    # --- Build canonical order ---
    contig_order = _build_contig_order(
        contig_lengths, cg_map, ordering=contig_ordering, custom_order=custom_order
    )
    logger.debug(f"Contig order: {contig_order[:10]} ... {contig_order[-10:]}")

    # --- Build genomes / contigs tables ---
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
        f"Built final genomes and contigs tables with {len(genomes_df)} genomes and {len(contigs_df)} contigs"
    )

    # --- Build bins ---
    if bins_df is None:
        if binsize is None:
            raise ValueError("binsize required when bins is None")
        bins_df = _tile_bins(contig_order, contig_lengths, contig_name_to_id, binsize)
    else:
        bins_df = _build_bins(bins_df, contig_order, contig_name_to_id)
    bins_df = bins_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    logger.info(
        f"Built final bins table with {len(bins_df)} bins tiling {len(contig_lengths)} contigs"
    )

    # --- Compute offsets ---
    n_contigs = len(contigs_df)
    n_bins = len(bins_df)

    # chrom_offset: first bin index for each contig (O(n_bins) via bincount)
    chrom_counts = np.bincount(
        bins_df["chrom"].values.astype(np.int64), minlength=n_contigs
    )
    chrom_offset = np.zeros(n_contigs + 1, dtype=np.int64)
    chrom_offset[1:] = np.cumsum(chrom_counts)
    logger.debug(f"Chrom offsets: {chrom_offset[:10]} ... {chrom_offset[-10:]}")

    # genome_offset: first contig index for each genome (O(n_contigs) via bincount)
    n_genomes = len(genomes_df)
    genome_ids = contigs_df["genome_id"].values
    genome_counts = np.bincount(genome_ids.astype(np.int64), minlength=n_genomes)
    genome_offset = np.zeros(n_genomes + 1, dtype=np.int64)
    genome_offset[1:] = np.cumsum(genome_counts)
    logger.debug(f"Genome offsets: {genome_offset[:10]} ... {genome_offset[-10:]}")

    # n_bins per genome: sum chrom_counts by genome_id
    genome_bin_counts = np.bincount(
        contigs_df["genome_id"].values.astype(np.int64),
        weights=chrom_counts.astype(np.float64),
        minlength=n_genomes,
    ).astype(np.int64)
    genomes_df.insert(3, "n_bins", genome_bin_counts)
    logger.debug(
        f"Genome bin counts: {genome_bin_counts[:10]} ... {genome_bin_counts[-10:]}"
    )

    # --- Stream pairs and bin them ---
    # Flattened bin start positions — needed for fragment-level binning.
    _bin_starts: np.ndarray | None = (
        None if binsize is not None else bins_df["start"].values.astype(np.int64)
    )

    def pixel_stream():
        """Yield aggregated pixel DataFrames one chunk at a time."""
        from rich.console import Console
        from rich.progress import (
            Progress,
            SpinnerColumn,
            TextColumn,
            TimeElapsedColumn,
        )

        _con = Console(stderr=True)
        with Progress(
            SpinnerColumn(),
            TextColumn("[cyan]Binning pairs"),
            TextColumn("{task.completed:,.0f} pairs processed"),
            TimeElapsedColumn(),
            console=_con,
            transient=True,
            disable=not _con.is_terminal,
        ) as _progress:
            _task = _progress.add_task("", total=None)
            for _chunk in _stream_pairs(pairs, chunksize=chunksize, **pairs_kwargs):
                _progress.update(_task, advance=len(_chunk))
                _result = _bin_pairs_chunk(
                    _chunk, contig_name_to_id, chrom_offset, binsize, _bin_starts
                )
                if _result is not None and not _result.empty:
                    yield _result

    # --- Create the cooler file ---
    logger.info(f"Starting to create cooler file at {uri}...")
    cooler.create_cooler(
        str(uri),
        bins=bins_df,
        pixels=pixel_stream(),
        ordered=False,
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

    # Validate
    logger.info("Validating created metacool file")
    validate(uri)
    with h5py.File(uri, "r") as f:
        nnz = len(f["pixels"]["count"])
        total = int(f["pixels"]["count"][:].sum())
    logger.info(
        f"Created {n_bins} bins from {n_contigs} contigs in {n_genomes} genomes. NNZ={nnz}, sum={total}.",
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

    # Validate input cool file
    if not cooler.fileops.is_cooler(str(cool_path)):
        raise ValueError(f"Not a valid cooler file: {cool_path}")

    # Load cool file data (bins are small; pixels are streamed below)
    clr = cooler.Cooler(str(cool_path))
    bins = clr.bins()[:]

    logger.info(f"Loaded cool file: {len(bins)} bins, {clr.info['nnz']} pixels")

    # Get contig names from cool file's chroms table
    chroms = clr.chroms()[:]
    contig_id_to_name = dict(enumerate(chroms["name"]))

    # Build contig lengths from bins (use actual contigs in bins, not chroms table)
    contig_lengths: dict[str, int] = {
        str(chrom): int(group["end"].max()) for chrom, group in bins.groupby("chrom")
    }

    # Parse contig assignments
    contig_map, contig_df = parse_contigs_table(contig_genome_map, contig_table_format)

    # Normalize empty/NA genome names to __unbinned__
    contig_map = {k: _normalize_genome_name(v) for k, v in contig_map.items()}

    # Drop contigs not present in the contig table (or assigned to __unbinned__)
    if drop_unassigned:
        contig_lengths = {
            ct: ln
            for ct, ln in contig_lengths.items()
            if ct in contig_map and contig_map[ct] != "__unbinned__"
        }

    # Build contig order
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

    # Build genomes table
    if genomes_table is not None:
        genomes_df_extra = parse_genomes_table(genomes_table)
    else:
        genomes_df_extra = None

    genomes = _build_genomes_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=contig_map,
        extra_genomes_df=genomes_df_extra,
    )

    # Build contigs table
    contigs = _build_contigs_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=contig_map,
        genomes=genomes,
        contig_df=contig_df,
    )

    # Remap bins to new contig IDs (based on reordered contigs)
    new_contig_name_to_id = {name: idx for idx, name in enumerate(contig_order)}

    # Remap bins chrom column (bins have contig names, not IDs)
    bins_remapped = bins.copy()
    bins_remapped["chrom"] = bins_remapped["chrom"].map(new_contig_name_to_id)
    if bins_remapped["chrom"].isna().any():
        missing = bins_remapped[bins_remapped["chrom"].isna()]["chrom"].unique()
        raise ValueError(f"Bins reference contigs not in contig_order: {missing}")

    bins_remapped = bins_remapped[["chrom", "start", "end"]].copy()
    bins_remapped["chrom"] = bins_remapped["chrom"].astype(np.int32)
    bins_remapped["start"] = bins_remapped["start"].astype(np.int32)
    bins_remapped["end"] = bins_remapped["end"].astype(np.int32)

    # Sort bins by new chrom ID so cooler.create_cooler stores them in the correct
    # order.  Without this, cooler renumbers enum codes by first-appearance order,
    # decoupling the code→name mapping from the physical bin data (start/end).
    sort_idx = np.argsort(bins_remapped["chrom"].values, kind="stable")
    bins_remapped = bins_remapped.iloc[sort_idx].reset_index(drop=True)

    # Build old-bin-ID → new-bin-ID mapping for pixel remapping
    n_bins = len(bins_remapped)
    old_to_new_bin = np.empty(n_bins, dtype=np.int64)
    old_to_new_bin[sort_idx] = np.arange(n_bins, dtype=np.int64)

    # n_bins per genome: sum bin counts by genome_id via contigs
    _chrom_counts = np.bincount(
        bins_remapped["chrom"].values.astype(np.int64), minlength=len(contigs)
    )
    _genome_bin_counts = np.bincount(
        contigs["genome_id"].values.astype(np.int64),
        weights=_chrom_counts.astype(np.float64),
        minlength=len(genomes),
    ).astype(np.int64)
    genomes.insert(3, "n_bins", _genome_bin_counts)

    # Stream and remap pixel bin IDs in chunks using a generator to avoid
    # loading the full pixel table into memory.
    _chunksize = 1_000_000
    _n_pixels = clr.info["nnz"]
    _cool_path_str = str(cool_path)  # pickling-safe string for thread-local opens

    def _remap_chunk(lo: int, hi: int) -> pd.DataFrame:
        """Read and remap a single pixel range; opens a thread-local Cooler handle."""
        import threading

        tl = _remap_chunk.__dict__.setdefault("_tl", threading.local())
        if not hasattr(tl, "clr"):
            tl.clr = cooler.Cooler(_cool_path_str)
        _chunk = tl.clr.pixels()[lo:hi].copy()
        _chunk["bin1_id"] = old_to_new_bin[_chunk["bin1_id"].values.astype(np.int64)]
        _chunk["bin2_id"] = old_to_new_bin[_chunk["bin2_id"].values.astype(np.int64)]
        _chunk["count"] = _chunk["count"].astype(np.int32)
        _swap = _chunk["bin1_id"] > _chunk["bin2_id"]
        _chunk.loc[_swap, ["bin1_id", "bin2_id"]] = _chunk.loc[
            _swap, ["bin2_id", "bin1_id"]
        ].values
        return _chunk[["bin1_id", "bin2_id", "count"]]

    def _remap_pixel_stream():
        from rich.console import Console
        from rich.progress import (
            BarColumn,
            MofNCompleteColumn,
            Progress,
            TaskProgressColumn,
            TextColumn,
            TimeRemainingColumn,
        )

        _ranges = [
            (i, min(i + _chunksize, _n_pixels)) for i in range(0, _n_pixels, _chunksize)
        ]
        _con = Console(stderr=True)
        with Progress(
            TextColumn("[cyan]Remapping pixels"),
            BarColumn(),
            TaskProgressColumn(),
            MofNCompleteColumn(),
            TextColumn("chunks"),
            TimeRemainingColumn(),
            console=_con,
            transient=True,
            disable=not _con.is_terminal,
        ) as _progress:
            _task = _progress.add_task("", total=len(_ranges))

            if nproc == 1:
                for _lo, _hi in _ranges:
                    yield _remap_chunk(_lo, _hi)
                    _progress.update(_task, advance=1)
            else:
                # Parallel path: submit all ranges; h5py releases the GIL on reads
                # so ThreadPoolExecutor gives genuine I/O overlap.
                from concurrent.futures import ThreadPoolExecutor

                with ThreadPoolExecutor(max_workers=nproc) as _pool:
                    _futures = [
                        _pool.submit(_remap_chunk, _lo, _hi) for _lo, _hi in _ranges
                    ]
                    for _f in _futures:
                        yield _f.result()
                        _progress.update(_task, advance=1)

    logger.info(
        f"Building metacool file: {len(genomes)} genomes, {len(contigs)} contigs"
    )

    # Create cooler file from remapped bins and pixel stream.
    # Use ordered=False so cooler performs an external sort, avoiding the need
    # to materialise the full pixel table in memory before writing.
    cooler.create_cooler(
        str(output_path),
        bins=bins_remapped,
        pixels=_remap_pixel_stream(),
        ordered=False,
        symmetric_upper=True,
        assembly=assembly_name or "imported",
        metadata=metadata or {},
    )

    # Compute genome_offset from contigs
    genome_ids = contigs["genome_id"].values
    n_genomes = len(genomes)
    genome_counts = np.bincount(genome_ids.astype(np.int64), minlength=n_genomes)
    genome_offset = np.zeros(n_genomes + 1, dtype=np.int64)
    genome_offset[1:] = np.cumsum(genome_counts)

    # Add metacooler-specific groups and attributes
    _finalize_metacool_file(
        output_path,
        genomes,
        contigs,
        None,
        genome_offset,
        contig_ordering=contig_ordering,
        assembly_name=assembly_name,
    )

    logger.info(f"Successfully created metacool file: {output_path}")
