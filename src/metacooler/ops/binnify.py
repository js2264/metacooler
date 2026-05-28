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

    # Remap pixels in chunks and accumulate into a dict to handle the many-to-one
    # bin mapping: multiple source pixels can land on the same coarser output pair.
    def _remap_pixel_chunks(
        chunksize: int = 1_000_000,
    ):
        """Remap old bin IDs to new bin IDs and accumulate counts.

        Streams pixels from the source file in chunks, remaps each
        ``bin1_id`` / ``bin2_id`` via *old_to_new*, enforces the upper
        triangle, and accumulates pixel counts in a dict.  Yields a
        single sorted DataFrame at the end.

        Args:
            chunksize: Number of pixels to load per chunk.

        Yields:
            One ``pd.DataFrame`` with columns ``bin1_id``, ``bin2_id``,
            ``count`` containing all remapped pixels.
        """
        accumulator: dict[tuple[int, int], int] = {}
        n_pix = len(mc._clr.pixels())
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            b1 = chunk["bin1_id"].values
            b2 = chunk["bin2_id"].values
            counts = chunk["count"].values

            pixel_mask = valid_mask[b1] & valid_mask[b2]
            if not pixel_mask.any():
                continue
            new_b1 = old_to_new[b1[pixel_mask]]
            new_b2 = old_to_new[b2[pixel_mask]]
            new_counts = counts[pixel_mask]

            # Enforce upper triangle
            swap = new_b1 > new_b2
            new_b1[swap], new_b2[swap] = new_b2[swap].copy(), new_b1[swap].copy()

            # Accumulate counts — many old pairs can map to the same new pair
            for i1, i2, c in zip(
                new_b1.tolist(), new_b2.tolist(), new_counts.tolist(), strict=False
            ):
                key = (int(i1), int(i2))
                accumulator[key] = accumulator.get(key, 0) + int(c)

        if accumulator:
            keys_sorted = sorted(accumulator)
            b1_out = np.array([k[0] for k in keys_sorted], dtype=np.int64)
            b2_out = np.array([k[1] for k in keys_sorted], dtype=np.int64)
            counts_out = np.array([accumulator[k] for k in keys_sorted], dtype=np.int32)
            yield pd.DataFrame(
                {"bin1_id": b1_out, "bin2_id": b2_out, "count": counts_out}
            )

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
        "Re-binned %s → %s: binsize=%d, %d bins",
        uri,
        output,
        binsize,
        len(new_bins_out),
    )

    # Recompute n_bins per genome from the new bins
    new_bin_genome_ids = contigs["genome_id"].values[
        pd.Categorical(new_bins_out["chrom"], categories=contigs["name"].values).codes
    ]
    new_genome_bin_counts = np.bincount(
        new_bin_genome_ids.astype(np.int64), minlength=len(genomes)
    ).astype(np.int64)
    genomes = genomes.copy()
    if "n_bins" in genomes.columns:
        genomes["n_bins"] = new_genome_bin_counts
    else:
        genomes.insert(3, "n_bins", new_genome_bin_counts)

    # Write as a new metacool file using a temporary pairs approach
    # We'll write directly using cooler + h5py to preserve metacooler structure
    _write_rebinned_metacool(
        output=Path(output),
        bins=new_bins_out,
        pixels=_remap_pixel_chunks(),
        genomes=genomes,
        contigs=contigs,
        binsize=binsize,
        source_uri=uri,
    )


def _write_rebinned_metacool(
    output: Path,
    bins: pd.DataFrame,
    pixels,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    binsize: int,
    source_uri: str | Path,
) -> None:
    """Write a re-binned metacool file preserving the full metacooler schema.

    Args:
        output: Destination path for the output ``.metacool`` file.
        bins: Bins DataFrame with contig names in the ``chrom`` column.
        pixels: Iterable of pixel DataFrames (streamed; never fully in RAM).
        genomes: Genomes DataFrame.
        contigs: Contigs DataFrame.
        binsize: New bin size in base pairs.
        source_uri: Source file URI (used for logging only).
    """
    import cooler

    from metacooler.io._utils import _finalize_metacool_file
    from metacooler.schema import validate

    # Pre-compute index offsets from the bins DataFrame
    contig_names = contigs["name"].values
    bin_chrom_codes = pd.Categorical(
        bins["chrom"].values, categories=contig_names
    ).codes.astype(np.int64)
    chrom_counts = np.bincount(bin_chrom_codes, minlength=len(contig_names))
    chrom_offset = np.concatenate([[0], np.cumsum(chrom_counts)]).astype(np.int64)

    contig_genome_ids = contigs["genome_id"].values
    bin_genome_ids = contig_genome_ids[bin_chrom_codes]
    genome_counts = np.bincount(bin_genome_ids.astype(np.int64), minlength=len(genomes))
    genome_offset = np.concatenate([[0], np.cumsum(genome_counts)]).astype(np.int64)

    # Create the cooler first; the generator yields one pre-sorted, aggregated DataFrame
    cooler.create_cooler(
        str(output),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
    )

    _finalize_metacool_file(
        output,
        genomes,
        contigs,
        chrom_offset,
        genome_offset,
        contig_ordering="rebinned",
    )

    validate(output)
    logger.info("Written rebinned metacool to %s", output)
