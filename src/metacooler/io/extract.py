"""Extract a genome-level cool file from a metacool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

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
    all_genomes = mc.genomes
    all_contigs = mc.contigs

    # Resolve contig IDs from genomes
    selected_contig_ids: set[int] = set()
    genome_contig_map: dict[str, list[int]] = {}

    if genomes:
        for g in genomes:
            if g not in all_genomes["name"].values:
                raise KeyError(
                    f"Genome {g!r} not found. Available: {all_genomes['name'].tolist()}"
                )
            gid = int(all_genomes.loc[all_genomes["name"] == g].index[0])
            cids = all_contigs.index[all_contigs["genome_id"] == gid].tolist()
            selected_contig_ids.update(cids)
            genome_contig_map[g] = cids

    if contigs:
        for c in contigs:
            if c not in all_contigs["name"].values:
                raise KeyError(
                    f"Contig {c!r} not found. Available: {all_contigs['name'].tolist()}"
                )
            cid = int(all_contigs.loc[all_contigs["name"] == c].index[0])
            selected_contig_ids.add(cid)

    if split and genomes:
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
    pixels = mc.pixels()

    mask = pixels["bin1_id"].isin(genome_bin_set) & pixels["bin2_id"].isin(
        genome_bin_set
    )
    genome_pixels = pixels[mask].copy()

    genome_pixels["bin1_id"] = genome_pixels["bin1_id"].map(global_to_local_bin)
    genome_pixels["bin2_id"] = genome_pixels["bin2_id"].map(global_to_local_bin)
    genome_pixels = genome_pixels.dropna(subset=["bin1_id", "bin2_id"])
    genome_pixels["bin1_id"] = genome_pixels["bin1_id"].astype(np.int64)
    genome_pixels["bin2_id"] = genome_pixels["bin2_id"].astype(np.int64)
    genome_pixels["count"] = genome_pixels["count"].astype(np.int32)

    swap = genome_pixels["bin1_id"] > genome_pixels["bin2_id"]
    genome_pixels.loc[swap, ["bin1_id", "bin2_id"]] = genome_pixels.loc[
        swap, ["bin2_id", "bin1_id"]
    ].values
    genome_pixels = genome_pixels.groupby(["bin1_id", "bin2_id"], as_index=False)[
        "count"
    ].sum()

    cooler.create_cooler(
        str(output),
        bins=bins_df,
        pixels=genome_pixels,
        ordered=True,
        symmetric_upper=True,
        assembly=assembly,
        metadata={"source-file": str(mc.store)},
    )
    logger.info(
        "Written %d bins, %d pixels to %s", len(bins_df), len(genome_pixels), output
    )
