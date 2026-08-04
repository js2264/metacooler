"""Flatten a metacool file into a standard cool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np

from metacooler.io._utils import enforce_upper_triangle, resolve_count_dtype
from metacooler.schema import DEFAULT_CHUNKSIZE

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

    # Stream pixels in chunks to avoid loading the full table into memory.
    chunksize = DEFAULT_CHUNKSIZE
    n_pix = mc._clr.info["nnz"]

    # Preserve the width of the source ``count`` column: a binnified source
    # writes int64 and would silently wrap if narrowed back to int32 here.
    count_dtype = resolve_count_dtype(
        mc._clr.pixels()[:0]["count"].dtype, context="flatten"
    )

    def _pixel_stream():
        for _i in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[_i : _i + chunksize].copy()
            chunk["bin1_id"] = chunk["bin1_id"].astype(np.int64)
            chunk["bin2_id"] = chunk["bin2_id"].astype(np.int64)
            chunk["count"] = chunk["count"].astype(count_dtype)
            enforce_upper_triangle(chunk)
            yield chunk[["bin1_id", "bin2_id", "count"]]

    logger.info(
        "Flattening %s → %s: %d contigs, %d bins",
        uri,
        output,
        len(contigs),
        len(bins),
    )

    cooler.create_cooler(
        str(output),
        bins=bins,
        pixels=_pixel_stream(),
        ordered=True,
        symmetric_upper=True,
        dtypes={"count": count_dtype},
        metadata={"source-file": str(uri), "flattened": True},
    )
    logger.info("Written flattened cool to %s", output)
