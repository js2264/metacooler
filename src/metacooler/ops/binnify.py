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

    Supports both fixed-width and variable-width (custom) source binnings;
    the new *binsize* is applied uniformly across contigs or genomes
    regardless of how the source was originally binned.

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.metacool`` file.
        binsize: New bin size in base pairs.
        level: Tiling level. ``contig`` (default) tiles each contig
            independently, preserving the contig hierarchy. ``genome``
            artificially ligates each genome's contigs (ordered by decreasing
            length) into one sequence and tiles that uniformly; the output's
            contigs become the ligated genomes (one per genome), so bins may
            span the original contig boundaries.

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
        # Tile each contig independently.
        new_bins_rows = []
        for cid in range(len(contigs)):
            contig_length = int(contigs.iloc[cid]["length"])
            for start in range(0, contig_length, binsize):
                end = min(start + binsize, contig_length)
                new_bins_rows.append({"chrom": cid, "start": start, "end": end})
        new_bins = pd.DataFrame(new_bins_rows)
        if new_bins.empty:
            raise ValueError("No bins generated. Check input data.")
        # Canonical (contig-id, start) order so bin rows match the /contigs table.
        new_bins = new_bins.sort_values(["chrom", "start"], kind="stable").reset_index(
            drop=True
        )
        old_to_new = _map_old_bins_per_contig(
            old_bins, bin_chrom_codes, len(contigs), new_bins
        )
    elif level == "genome":
        # Ligate each genome's contigs (ordered by decreasing length) into a
        # single artificial sequence, then tile that sequence uniformly. The
        # output's "contigs" ARE the ligated genomes (one contig per genome), so
        # bins genuinely span contig boundaries within a genome.
        new_bins, old_to_new, contigs, genomes = _tile_per_genome_ligated(
            contigs, genomes, old_bins, bin_chrom_codes, binsize
        )
        if new_bins.empty:
            raise ValueError("No bins generated. Check input data.")
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

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
        partial_frames: list[pd.DataFrame] = []
        n_pix = len(mc._clr.pixels())
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            b1 = chunk["bin1_id"].values
            b2 = chunk["bin2_id"].values
            # int64 so the per-chunk and cross-chunk count sums cannot overflow int32
            counts = chunk["count"].values.astype(np.int64)

            pixel_mask = valid_mask[b1] & valid_mask[b2]
            if not pixel_mask.any():
                continue
            new_b1 = old_to_new[b1[pixel_mask]]
            new_b2 = old_to_new[b2[pixel_mask]]
            new_counts = counts[pixel_mask]

            # Enforce upper triangle
            swap = new_b1 > new_b2
            new_b1[swap], new_b2[swap] = new_b2[swap].copy(), new_b1[swap].copy()

            chunk_df = pd.DataFrame(
                {"bin1_id": new_b1, "bin2_id": new_b2, "count": new_counts}
            )
            partial_frames.append(
                chunk_df.groupby(["bin1_id", "bin2_id"], sort=False, as_index=False)[
                    "count"
                ].sum()
            )

        if partial_frames:
            combined = pd.concat(partial_frames, ignore_index=True)
            result = combined.groupby(
                ["bin1_id", "bin2_id"], sort=True, as_index=False
            )["count"].sum()
            yield result

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
        assembly_name=mc.info.get("assembly"),
        contig_ordering=mc.info.get("ordering", "genome-then-alpha"),
    )


def _map_old_bins_per_contig(
    old_bins: pd.DataFrame,
    bin_chrom_codes: np.ndarray,
    n_contigs: int,
    new_bins: pd.DataFrame,
) -> np.ndarray:
    """Map each old bin to a new bin within the SAME contig (midpoint rule)."""
    old_midpoints = (old_bins["start"].values + old_bins["end"].values) / 2
    new_bin_lookup: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for cid, grp in new_bins.groupby("chrom", sort=False):
        new_bin_lookup[int(cid)] = (  # ty:ignore[invalid-argument-type]
            grp["start"].values.astype(np.int64),
            grp.index.values.astype(np.int64),
        )

    old_to_new = np.full(len(old_bins), -1, dtype=np.int64)
    for cid in range(n_contigs):
        old_mask = bin_chrom_codes == cid
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
    return old_to_new


def _tile_per_genome_ligated(
    contigs: pd.DataFrame,
    genomes: pd.DataFrame,
    old_bins: pd.DataFrame,
    bin_chrom_codes: np.ndarray,
    binsize: int,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Ligate contigs per genome (decreasing length) and tile at genome level.

    Returns ``(new_bins, old_to_new, out_contigs, out_genomes)`` where the
    output contigs are the ligated genomes (one per genome, named after the
    genome) so that a re-binned bin can span the original contig boundaries.
    """
    n_genomes = len(genomes)
    contig_len = contigs["length"].values.astype(np.int64)
    contig_gid = contigs["genome_id"].values.astype(np.int64)

    # Ligated offset of each original contig within its genome. Contigs are
    # laid out by DECREASING length (ties broken by contig id for determinism).
    contig_offset = np.zeros(len(contigs), dtype=np.int64)
    genome_length = np.zeros(n_genomes, dtype=np.int64)
    for gid in range(n_genomes):
        cids = np.where(contig_gid == gid)[0]
        order = sorted(cids.tolist(), key=lambda c: (-int(contig_len[c]), int(c)))
        off = 0
        for c in order:
            contig_offset[c] = off
            off += int(contig_len[c])
        genome_length[gid] = off

    # Uniform tiling of each ligated genome (already in (gid, start) order).
    rows = []
    genome_bin_start = np.zeros(n_genomes + 1, dtype=np.int64)
    for gid in range(n_genomes):
        glen = int(genome_length[gid])
        n_new = 0
        for start in range(0, glen, binsize):
            end = min(start + binsize, glen)
            rows.append({"chrom": gid, "start": start, "end": end})
            n_new += 1
        genome_bin_start[gid + 1] = genome_bin_start[gid] + n_new
    new_bins = pd.DataFrame(rows)
    if new_bins.empty:
        return new_bins, np.array([], dtype=np.int64), contigs, genomes

    # Map each old bin to the ligated-genome bin containing its midpoint.
    old_chrom = bin_chrom_codes.astype(np.int64)
    old_mid = (
        old_bins["start"].values.astype(np.int64)
        + old_bins["end"].values.astype(np.int64)
    ) // 2
    ligated_pos = contig_offset[old_chrom] + old_mid
    gid_of_old = contig_gid[old_chrom]
    old_to_new = genome_bin_start[gid_of_old] + (ligated_pos // binsize)

    # Output "contigs" are the ligated genomes (one contig per genome).
    out_contigs = pd.DataFrame(
        {
            "name": genomes["name"].values,
            "length": genome_length,
            "genome_id": np.arange(n_genomes, dtype=np.int32),
        }
    )
    out_genomes = genomes.copy()
    out_genomes["length"] = genome_length
    out_genomes["n_contigs"] = 1
    return new_bins, old_to_new, out_contigs, out_genomes


def _write_rebinned_metacool(
    output: Path,
    bins: pd.DataFrame,
    pixels,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    binsize: int,
    source_uri: str | Path,
    assembly_name: str | None = None,
    contig_ordering: str = "genome-then-alpha",
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
        assembly_name: Optional assembly name to propagate from source file.
        contig_ordering: Contig ordering string to propagate from source file.
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

    # genome_offset indexes the /contigs table (CSR over contigs), so it must
    # count CONTIGS per genome — not bins. Counting bins would produce offsets
    # running up to n_bins, i.e. out-of-range indices into the contigs table.
    contig_genome_ids = contigs["genome_id"].values
    genome_counts = np.bincount(
        contig_genome_ids.astype(np.int64), minlength=len(genomes)
    )
    genome_offset = np.concatenate([[0], np.cumsum(genome_counts)]).astype(np.int64)

    # Create the cooler first; the generator yields one pre-sorted, aggregated DataFrame.
    # int64 count so a coarse bin-pair summing >2**31 source contacts is not
    # saturated by cooler's default int32.
    cooler.create_cooler(
        str(output),
        bins=bins,
        pixels=pixels,
        dtypes={"count": np.int64},
        ordered=True,
        symmetric_upper=True,
    )

    _finalize_metacool_file(
        output,
        genomes,
        contigs,
        chrom_offset,
        genome_offset,
        contig_ordering=contig_ordering,
        assembly_name=assembly_name,
    )

    validate(output)
    logger.info("Written rebinned metacool to %s", output)
