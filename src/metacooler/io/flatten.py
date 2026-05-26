"""Flatten a metacool file into a standard cool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np

logger = logging.getLogger("metacooler")


def flatten(
    uri: str | Path,
    output: str | Path,
) -> None:
    """Flatten a metacool file into a standard cool file.

    Removes the genome hierarchy: all contigs become "chromosomes" in
    the output cool file.  Bins and pixels are preserved as-is (bin IDs
    do not change since contigs are already the ``chrom`` column).

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.cool`` file.
    """
    from metacooler.api import MetaCooler

    mc = MetaCooler(uri)

    # Get contigs as the chromosomes for the output cool
    contigs = mc.contigs

    # Build bins DataFrame — chrom column already uses contig names
    bins = mc.bins().copy()
    # Ensure chrom is string, not categorical codes
    if hasattr(bins["chrom"], "cat"):
        bins["chrom"] = bins["chrom"].astype(str)

    bins = bins[["chrom", "start", "end"]].copy()
    bins["start"] = bins["start"].astype(np.int32)
    bins["end"] = bins["end"].astype(np.int32)

    # Get pixels — bin IDs are unchanged
    pixels = mc.pixels().copy()
    pixels["bin1_id"] = pixels["bin1_id"].astype(np.int64)
    pixels["bin2_id"] = pixels["bin2_id"].astype(np.int64)
    pixels["count"] = pixels["count"].astype(np.int32)

    # Enforce upper triangle
    swap = pixels["bin1_id"] > pixels["bin2_id"]
    pixels.loc[swap, ["bin1_id", "bin2_id"]] = pixels.loc[
        swap, ["bin2_id", "bin1_id"]
    ].values

    logger.info(
        "Flattening %s → %s: %d contigs, %d bins, %d pixels",
        uri,
        output,
        len(contigs),
        len(bins),
        len(pixels),
    )

    cooler.create_cooler(
        str(output),
        bins=bins,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
        metadata={"source-file": str(uri), "flattened": True},
    )
    logger.info("Written flattened cool to %s", output)
