"""Extract a genome-level cool file from a metacool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import _filter_by_genomes_and_contigs

logger = logging.getLogger("metacooler")


def extract(
    uri: str | Path,
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
    output: str | Path = "extracted.cool",
    split: bool = False,
) -> None:
    """Extract selected genomes and/or contigs from a metacool file.

    ``--genomes`` and ``--contigs`` can be combined: the union of all
    contigs from specified genomes plus directly specified contigs is
    extracted.

    Args:
        uri: Path to the source ``.metacool`` file.
        genomes: List of genome names to extract.
        contigs: List of contig names to extract.
        output: Destination path. If *split* is True, treated as a
            directory prefix (one ``.cool`` per genome).
        split: If True, write one ``.cool`` per genome instead of a
            single merged file.

    Raises:
        KeyError: If any genome or contig is not found.
        ValueError: If neither genomes nor contigs are specified.
    """
    from metacooler.api import MetaCooler

    if not genomes and not contigs:
        raise ValueError("At least one of --genomes or --contigs must be specified.")

    mc = MetaCooler(uri)
    selected_contig_ids = _filter_by_genomes_and_contigs(
        mc, genomes=genomes, contigs=contigs
    )

    # Build genome_contig_map only when split mode needs it
    genome_contig_map: dict[str, list[int]] = {}
    if split and genomes:
        all_genomes_df = mc.genomes
        all_contigs = mc.contigs
        if contigs:
            logger.warning(
                "--contigs specified together with --split: extra contigs will be "
                "appended to every genome's output file. "
                "Consider running without --split or omitting --contigs."
            )
        for g in genomes:
            gid = int(all_genomes_df.loc[all_genomes_df["name"] == g].index[0])
            genome_contig_map[g] = all_contigs.index[
                all_contigs["genome_id"] == gid
            ].tolist()

        # Write one cool per genome
        output_path = Path(output)
        for g in genomes:
            cids = genome_contig_map[g]
            if contigs:
                # Also include directly specified contigs
                extra = set()
                for c in contigs:
                    cid = int(all_contigs.loc[all_contigs["name"] == c].index[0])
                    extra.add(cid)
                cids = list(set(cids) | extra)
            out_file = output_path.parent / f"{output_path.stem}_{g}.cool"
            _extract_contig_ids(mc, sorted(cids), out_file, assembly=g)
            logger.info("Extracted genome %r → %s", g, out_file)
    else:
        if split and not genomes and contigs:
            logger.warning(
                "--split requires at least one --genomes argument to split by genome. "
                "Since only --contigs were specified, writing a single merged output file.",
            )
        # Single output file with all selected contigs
        _extract_contig_ids(
            mc, sorted(selected_contig_ids), Path(output), assembly="multi"
        )


def _extract_contig_ids(
    mc,
    contig_ids: list[int],
    output: Path,
    assembly: str = "unknown",
) -> None:
    """Extract bins/pixels for a set of contig IDs and write a cool file."""
    all_bins = mc.bins()
    bin_chrom_codes = all_bins["chrom"].cat.codes.values

    # Boolean mask: bins belonging to selected contigs
    genome_bin_mask = np.isin(bin_chrom_codes, contig_ids)
    genome_bin_indices = np.where(genome_bin_mask)[0]

    global_to_local_bin = {
        int(gbin): local for local, gbin in enumerate(genome_bin_indices)
    }

    genome_bins = all_bins.iloc[genome_bin_indices].copy()
    genome_bins_chrom_names = genome_bins["chrom"].values

    bins_df = pd.DataFrame(
        {
            "chrom": genome_bins_chrom_names,
            "start": genome_bins["start"].values.astype(np.int32),
            "end": genome_bins["end"].values.astype(np.int32),
        }
    )

    genome_bin_set = set(genome_bin_indices.tolist())
    genome_bin_array = np.array(sorted(genome_bin_set), dtype=np.int64)
    valid_min, valid_max = genome_bin_array[0], genome_bin_array[-1]
    chunksize = 1_000_000

    # Use bin1_offset index to find pixel range for this genome's bins
    with h5py.File(str(mc.store), "r") as _f:
        _root = _f[mc._group_path.lstrip("/")] if mc._group_path != "/" else _f
        bin1_offset = _root["indexes"]["bin1_offset"][:]
    pix_lo = int(bin1_offset[valid_min])
    n_pix = len(mc._clr.pixels())
    pix_hi = (
        int(bin1_offset[valid_max + 1]) if valid_max + 1 < len(bin1_offset) else n_pix
    )

    def _pixel_stream():
        for _i in range(pix_lo, pix_hi, chunksize):
            chunk = mc._clr.pixels()[_i : _i + chunksize].copy()
            mask = chunk["bin1_id"].isin(genome_bin_set) & chunk["bin2_id"].isin(
                genome_bin_set
            )
            if not mask.any():
                continue
            chunk = chunk[mask].copy()
            chunk["bin1_id"] = chunk["bin1_id"].map(global_to_local_bin)
            chunk["bin2_id"] = chunk["bin2_id"].map(global_to_local_bin)
            chunk["bin1_id"] = chunk["bin1_id"].astype(np.int64)
            chunk["bin2_id"] = chunk["bin2_id"].astype(np.int64)
            chunk["count"] = chunk["count"].astype(np.int32)
            swap = chunk["bin1_id"] > chunk["bin2_id"]
            chunk.loc[swap, ["bin1_id", "bin2_id"]] = chunk.loc[
                swap, ["bin2_id", "bin1_id"]
            ].values
            yield chunk[["bin1_id", "bin2_id", "count"]]

    cooler.create_cooler(
        str(output),
        bins=bins_df,
        pixels=_pixel_stream(),
        ordered=False,
        symmetric_upper=True,
        assembly=assembly,
        metadata={"source-file": str(mc.store)},
    )
    logger.info("Written %d bins to %s", len(bins_df), output)
