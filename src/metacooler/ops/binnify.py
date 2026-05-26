"""Re-binning operation for metacool files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")


def binnify(
    uri: str | Path,
    output: str | Path,
    binsize: int,
    level: Literal["contig", "genome"] = "contig",
) -> None:
    """Re-bin a metacool file to a new fixed bin size.

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.metacool`` file.
        binsize: New bin size in base pairs.
        level: Tiling level — ``contig`` (tile per contig) or
            ``genome`` (tile per genome, contigs concatenated).

    Raises:
        ValueError: If the input has a fixed binsize that is >= *binsize*
            (new binsize must be *strictly* larger).
    """
    from metacooler.api import MetaCooler

    mc = MetaCooler(uri)

    # Validate binsize
    current_binsize = mc.binsize
    if current_binsize is not None and binsize <= current_binsize:
        raise ValueError(
            f"New binsize ({binsize}) must be strictly larger than current binsize ({current_binsize})."
        )

    contigs = mc.contigs
    genomes = mc.genomes
    old_bins = mc.bins()
    old_pixels = mc.pixels()

    # Get contig codes for each bin
    bin_chrom_codes = old_bins["chrom"].cat.codes.values

    if level == "contig":
        # Tile per contig
        new_bins_rows = []
        for cid in range(len(contigs)):
            contig_length = int(contigs.iloc[cid]["length"])
            for start in range(0, contig_length, binsize):
                end = min(start + binsize, contig_length)
                new_bins_rows.append({"chrom": cid, "start": start, "end": end})
        new_bins = pd.DataFrame(new_bins_rows)
    elif level == "genome":
        # Tile per genome (contigs concatenated within each genome).
        # Bins that span contig boundaries are split into one sub-bin per contig.
        new_bins_rows = []
        for gid in range(len(genomes)):
            genome_contigs = contigs[contigs["genome_id"] == gid]
            offset = 0
            contig_offsets = {}
            for cid in genome_contigs.index:
                contig_offsets[cid] = offset
                offset += int(contigs.iloc[cid]["length"])
            genome_length = offset
            for start in range(0, genome_length, binsize):
                end = min(start + binsize, genome_length)
                # Walk through all contigs this [start, end) range touches
                # and emit one sub-bin per contig segment.
                for cid in genome_contigs.index:
                    clen = int(contigs.iloc[cid]["length"])
                    contig_start = contig_offsets[cid]
                    contig_end = contig_start + clen
                    if end <= contig_start or start >= contig_end:
                        continue  # this contig not touched by [start, end)
                    local_start = max(start, contig_start) - contig_start
                    local_end = min(end, contig_end) - contig_start
                    new_bins_rows.append(
                        {"chrom": cid, "start": local_start, "end": local_end}
                    )
        new_bins = pd.DataFrame(new_bins_rows)
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    if new_bins.empty:
        raise ValueError("No bins generated. Check input data.")

    # Build mapping: old bin → new bin (using midpoint assignment)
    old_midpoints = (old_bins["start"].values + old_bins["end"].values) / 2
    old_chroms = bin_chrom_codes

    # Build new bin lookup: for each contig, sorted start positions and indices
    new_bin_lookup: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for cid, grp in new_bins.groupby("chrom", sort=False):
        new_bin_lookup[int(cid)] = (  # ty:ignore[invalid-argument-type]
            grp["start"].values.astype(np.int64),
            grp.index.values.astype(np.int64),
        )

    # Map old bins to new bins (vectorised per-contig searchsorted)
    old_to_new = np.full(len(old_bins), -1, dtype=np.int64)
    for cid in range(len(contigs)):
        old_mask = old_chroms == cid
        if not old_mask.any():
            continue
        new_starts, new_indices = new_bin_lookup.get(
            cid, (np.array([], dtype=np.int64), np.array([], dtype=np.int64))
        )
        if len(new_starts) == 0:
            continue
        mids = old_midpoints[old_mask]
        idx = np.searchsorted(new_starts, mids, side="right") - 1
        valid = idx >= 0
        positions = np.where(old_mask)[0]
        old_to_new[positions[valid]] = new_indices[idx[valid]]

    # Filter out unmapped bins
    valid_mask = old_to_new >= 0

    # Remap pixels
    b1 = old_pixels["bin1_id"].values
    b2 = old_pixels["bin2_id"].values
    counts = old_pixels["count"].values

    # Filter pixels where both bins are mapped
    pixel_mask = valid_mask[b1] & valid_mask[b2]
    new_b1 = old_to_new[b1[pixel_mask]]
    new_b2 = old_to_new[b2[pixel_mask]]
    new_counts = counts[pixel_mask]

    # Enforce upper triangle
    swap = new_b1 > new_b2
    new_b1[swap], new_b2[swap] = new_b2[swap].copy(), new_b1[swap].copy()

    # Aggregate
    pixel_df = pd.DataFrame({"bin1_id": new_b1, "bin2_id": new_b2, "count": new_counts})
    pixel_df = pixel_df.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()

    # Build output bins with contig names
    contig_names = contigs["name"].values
    new_bins_out = pd.DataFrame(
        {
            "chrom": contig_names[new_bins["chrom"].values.astype(int)],
            "start": new_bins["start"].values.astype(np.int32),
            "end": new_bins["end"].values.astype(np.int32),
        }
    )

    logger.info(
        "Re-binned %s → %s: binsize=%d, %d bins, %d pixels",
        uri,
        output,
        binsize,
        len(new_bins_out),
        len(pixel_df),
    )

    # Write as a new metacool file using a temporary pairs approach
    # We'll write directly using cooler + h5py to preserve metacooler structure
    _write_rebinned_metacool(
        output=Path(output),
        bins=new_bins_out,
        pixels=pixel_df,
        genomes=genomes,
        contigs=contigs,
        binsize=binsize,
        source_uri=uri,
    )


def _write_rebinned_metacool(
    output: Path,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    binsize: int,
    source_uri: str | Path,
) -> None:
    """Write a re-binned metacool file preserving the full metacooler schema."""
    import cooler
    import h5py

    from metacooler._version import __version__
    from metacooler.schema import (
        COOLER_FORMAT,
        COOLER_FORMAT_URL,
        METACOOLER_FORMAT,
        METACOOLER_FORMAT_URL,
        METACOOLER_FORMAT_VERSION,
    )

    # Create the cooler first
    cooler.create_cooler(
        str(output),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
    )

    # Now add metacooler-specific tables and attributes
    with h5py.File(output, "a") as f:
        # Set metacooler format attributes
        f.attrs["format"] = COOLER_FORMAT
        f.attrs["format-url"] = COOLER_FORMAT_URL
        f.attrs["metacooler-format"] = METACOOLER_FORMAT
        f.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        f.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        f.attrs["generated-by"] = f"metacooler-{__version__}"
        f.attrs["ngenomes"] = len(genomes)
        f.attrs["contig-ordering"] = "rebinned"

        def _to_hdf5_data(arr):
            """Convert array to HDF5-compatible dtype."""
            if (
                arr.dtype == object
                or arr.dtype.kind in ("U", "S")
                or str(arr.dtype) == "string"
            ):
                # Convert strings to fixed-length bytes for HDF5
                str_arr = [str(x) for x in arr]
                max_len = max(len(s) for s in str_arr) if len(str_arr) > 0 else 1
                return np.array(
                    [s.encode("utf-8") for s in str_arr], dtype=f"S{max_len}"
                )
            return np.asarray(arr)

        # Write genomes table
        g = f.require_group("genomes")
        for col in genomes.columns:
            data = _to_hdf5_data(genomes[col].values)
            if col in g:
                del g[col]
            g.create_dataset(col, data=data)

        # Write contigs table (overwrite chroms-based one)
        if "contigs" in f:
            del f["contigs"]
        c = f.create_group("contigs")
        for col in contigs.columns:
            data = _to_hdf5_data(contigs[col].values)
            c.create_dataset(col, data=data)

        # Create soft link chroms -> contigs (cooler compat)
        if "chroms" in f:
            del f["chroms"]
        f["chroms"] = h5py.SoftLink("contigs")

        # Build indexes
        idx = f.require_group("indexes")

        # chrom_offset: cumulative bin count per contig (O(n_bins) via bincount)
        bin_chroms = bins["chrom"].values
        contig_names = contigs["name"].values
        bin_chrom_codes = pd.Categorical(
            bin_chroms, categories=contig_names
        ).codes.astype(np.int64)
        chrom_counts = np.bincount(bin_chrom_codes, minlength=len(contig_names))
        chrom_offset = np.concatenate([[0], np.cumsum(chrom_counts)])
        if "chrom_offset" in idx:
            del idx["chrom_offset"]
        idx.create_dataset("chrom_offset", data=chrom_offset.astype(np.int64))

        # genome_offset: cumulative bin count per genome (O(n_bins) via bincount)
        # Map each bin's contig code to a genome id, then bincount
        contig_genome_ids = contigs["genome_id"].values  # shape: (n_contigs,)
        bin_genome_ids = contig_genome_ids[bin_chrom_codes]
        genome_counts = np.bincount(bin_genome_ids, minlength=len(genomes))
        genome_offset = np.concatenate([[0], np.cumsum(genome_counts)])
        if "genome_offset" in idx:
            del idx["genome_offset"]
        idx.create_dataset("genome_offset", data=genome_offset.astype(np.int64))

    logger.info("Written rebinned metacool to %s", output)
