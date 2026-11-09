"""Contact-quality statistics for a metacool file that already exists.

The classification itself lives in :mod:`metacooler.qc`; this module is the
file-facing half — it opens a map, walks its pixel table in chunks, and
optionally writes the result back.

What it measures is necessarily the *binned* map, so an intra-contig distance
here is the gap between two bin midpoints.  ``create pairs`` computes the same
block from raw read positions while it is already reading the pairs file, and
records ``source="pairs"``; prefer that one when both exist.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

from metacooler.api import MetaCooler
from metacooler.qc import (
    DEFAULT_DISTANCE_THRESHOLD,
    QC_METADATA_KEY,
    ContactQC,
    bin_maps,
    unbinned_genome_id,
)
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger("metacooler")


class QCResult(NamedTuple):
    """What :func:`qc` computed.

    Attributes:
        summary: The file-level block, as stored under ``metadata["qc"]``.
        per_genome: One row per genome, with the ``qc_*`` columns.
    """

    summary: dict
    per_genome: pd.DataFrame


def qc(
    uri: str | Path,
    *,
    distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
    chunksize: int = DEFAULT_CHUNKSIZE,
    store: bool = False,
) -> QCResult:
    """Compute contact-quality statistics from a metacool file's pixel table.

    Args:
        uri: ``file.metacool`` path or ``file.mmetacool::entries/name`` URI.
        distance_threshold: Intra-contig short/long split, in bp.
        chunksize: Pixels read per iteration.
        store: If ``True``, write the result back into the file —
            ``metadata["qc"]`` plus the per-genome ``qc_*`` columns of
            ``/genomes``.

    Returns:
        A :class:`QCResult`.

    Raises:
        ValueError: If *distance_threshold* is negative.
        OSError: If *store* is set and the file cannot be opened for writing.
    """
    mc = MetaCooler(uri)
    genomes = mc.genomes
    tally = ContactQC(
        n_genomes=len(genomes),
        distance_threshold=distance_threshold,
        source="pixels",
        unbinned_genome=unbinned_genome_id(genomes["name"]),
    )
    maps = bin_maps(mc.bins(), mc.contigs)

    # cooler hands `binsize` back as whatever the attribute holds — a NumPy
    # integer normally, but a string on a file another tool wrote — and this is
    # both compared with an int and JSON-serialised below.
    binsize = mc.binsize
    binsize = None if binsize is None else int(binsize)
    if binsize is not None and binsize >= tally.distance_threshold:
        logger.warning(
            "Bin size (%d bp) is at least the short/long distance threshold "
            "(%d bp), so every intra-contig contact except the diagonal is "
            "reported as long-range. Re-run 'create pairs' (which measures "
            "distances between read positions) or lower --threshold to split "
            "the intra-contig signal meaningfully.",
            binsize,
            tally.distance_threshold,
        )

    pixels = mc._clr.pixels()
    n_pixels = len(pixels)
    for lo in range(0, n_pixels, chunksize):
        chunk = pixels[lo : lo + chunksize]
        tally.tally_pixels(
            chunk["bin1_id"].to_numpy(),
            chunk["bin2_id"].to_numpy(),
            chunk["count"].to_numpy(),
            maps,
        )

    summary = tally.summary(binsize=binsize)
    per_genome = tally.genome_table(genomes["name"].to_numpy())
    logger.info(
        "QC over %d pixel(s): %s contacts, trans-contig ratio %s",
        n_pixels,
        f"{summary['n_contacts']:,}",
        summary["trans_contig_ratio"],
    )

    if store:
        from metacooler.io._utils import store_qc

        existing = mc.metadata.get(QC_METADATA_KEY)
        if isinstance(existing, dict) and existing.get("source") == "pairs":
            logger.warning(
                "Replacing a block measured on read positions with one measured "
                "between bin midpoints: the intra-contig short/long split will "
                "be coarser than the one 'create pairs' recorded.",
            )
        store_qc(mc.store, summary, per_genome, group_path=mc._group_path)
        logger.info("Stored QC block and per-genome columns in %s", mc.uri)
    return QCResult(summary, per_genome)
