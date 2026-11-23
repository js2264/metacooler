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
    BinMaps,
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
        distance_threshold: Intra-contig short/long split, in bp: long-range
            from this distance on, measured between bin midpoints.  Warns
            when it is no greater than the bin size, where little but the
            diagonal is short-range.
        chunksize: Pixels read per iteration.
        store: If ``True``, write the result back into the file —
            ``metadata["qc"]`` plus the per-genome ``qc_*`` columns of
            ``/genomes`` and the per-contig ones of ``/contigs``.

    Returns:
        A :class:`QCResult`.

    Raises:
        ValueError: If *distance_threshold* is negative.
        OSError: If *store* is set and the file cannot be written: another
            process holds its lock (``MapLockedError``), it changed while it
            was read (``MapChangedError``), or it is read-only.
    """
    if not store:
        return _qc(uri, distance_threshold=distance_threshold, chunksize=chunksize)
    from metacooler.io.staging import lock, map_file

    # Stored in place: under the map's lock from before it is read, through a
    # staged copy (see metacooler.io.staging).
    with lock(map_file(uri)):
        return _qc(
            uri, distance_threshold=distance_threshold, chunksize=chunksize, store=True
        )


def _qc(
    uri: str | Path,
    *,
    distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
    chunksize: int = DEFAULT_CHUNKSIZE,
    store: bool = False,
) -> QCResult:
    """:func:`qc`'s work, under the map's lock when it stores the result."""
    mc = MetaCooler(uri)
    genomes = mc.genomes
    contigs = mc.contigs
    tally = ContactQC(
        n_genomes=len(genomes),
        distance_threshold=distance_threshold,
        source="pixels",
        unbinned_genome=unbinned_genome_id(genomes["name"]),
        n_contigs=len(contigs) if store else None,
    )

    # cooler hands `binsize` back as whatever the attribute holds — a NumPy
    # integer normally, but a string on a file another tool wrote — and this is
    # both compared with an int and JSON-serialised below.
    binsize = mc.binsize
    binsize = None if binsize is None else int(binsize)
    tally.warn_if_split_within_bin(
        binsize,
        "A useful split needs a --threshold above the bin size, or a map made "
        "from pairs: 'create pairs' measures the split on read positions, and "
        "'metacooler metadata get -k qc' prints the block it stored.",
    )

    n_pixels = tally_map(mc, tally, chunksize=chunksize)

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
        from metacooler.io.staging import staged_inplace
        from metacooler.qc import QC_CONTIG_COLUMNS

        with staged_inplace(
            mc.uri,
            planned=mc,
            touches=[f"contigs/{column}" for column in QC_CONTIG_COLUMNS],
        ) as staged:
            store_qc(
                staged,
                summary,
                per_genome,
                tally.contig_table(contigs["name"].to_numpy()),
                group_path=mc._group_path,
            )
        logger.info(
            "Stored QC block and per-genome and per-contig columns in %s", mc.uri
        )
    return QCResult(summary, per_genome)


def tally_map(
    mc: MetaCooler,
    tally: ContactQC,
    *,
    chunksize: int = DEFAULT_CHUNKSIZE,
    maps: BinMaps | None = None,
) -> int:
    """Feed every pixel of an open map to a tally, chunk by chunk.

    Args:
        mc: The open map, whose pixels are read through cooler's selector.
        tally: Accumulator to update.
        chunksize: Pixels read per iteration.
        maps: The contig, genome and midpoint of every bin, which classify the
            pixels.  Defaults to the map's own (:func:`bin_maps` of its bins
            and ``/contigs``); ``reassign`` passes the genomes of a new
            assignment.

    Returns:
        The number of pixels read.
    """
    if maps is None:
        maps = bin_maps(mc.bins(), mc.contigs)
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
    return n_pixels
