"""Network edge list export from metacool files."""

from __future__ import annotations

import contextlib
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.io._utils import (
    aggregate_pixels,
    check_output_path,
    chrom_codes,
    to_upper_triangle,
)
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")

#: Number of per-chunk groupby frames buffered before they are folded into the
#: running aggregate.  Keeps peak memory at O(FOLD_EVERY x chunk edges) instead
#: of growing one frame per chunk for the whole pixel table.
_FOLD_EVERY = 32


def _fold(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge per-chunk edge frames, re-summing node pairs seen more than once.

    The frames are pixel-shaped (``bin1_id``/``bin2_id``/``count`` standing for
    ``node1``/``node2``/``value``), so this is :func:`aggregate_pixels`: a
    key-based merge whose output is sorted by node pair.

    Args:
        frames: Non-empty list of sorted, duplicate-free edge frames.

    Returns:
        One frame with a single row per node pair, sorted.
    """
    return aggregate_pixels(frames)


def _aggregate_pairs(
    mc,
    bin_to_node: np.ndarray,
    weight: str,
    n_pix: int,
    chunksize: int,
) -> pd.DataFrame | None:
    """Sum pixel weights per node pair, folding partial results as we go.

    Args:
        mc: Open :class:`~metacooler.api.MetaCooler`.
        bin_to_node: Array mapping each bin ID to its node ID (contig or genome).
        weight: Pixels column used as the edge weight.
        n_pix: Total number of pixels in the source table.
        chunksize: Number of pixels read per chunk.

    Returns:
        A frame with columns ``node1``/``node2``/``value`` sorted by node pair,
        where ``node1 <= node2``, or ``None`` if the pixel table is empty.
    """
    pending: list[pd.DataFrame] = []
    running: pd.DataFrame | None = None

    for lo in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[lo : lo + chunksize]
        n1, n2 = to_upper_triangle(
            bin_to_node[chunk["bin1_id"].to_numpy()],
            bin_to_node[chunk["bin2_id"].to_numpy()],
        )
        values = chunk[weight].to_numpy()
        # A node pair sums every pixel between two contigs or genomes, and
        # aggregate_pixels keeps the dtype it is given: an int32 count column
        # (what create writes) wrapped past 2**31 to a negative edge weight.
        # Widen integers as ops/aggregate does; float weights pass unchanged.
        if np.issubdtype(values.dtype, np.integer):
            values = values.astype(np.int64, copy=False)
        df = pd.DataFrame(
            {"bin1_id": n1, "bin2_id": n2, "count": values},
            copy=False,
        )
        pending.append(aggregate_pixels([df]))
        if len(pending) >= _FOLD_EVERY:
            running = _fold(pending if running is None else [running, *pending])
            pending = []

    if running is None and not pending:
        return None
    edges = _fold(pending if running is None else [running, *pending])
    return edges.rename(
        columns={"bin1_id": "node1", "bin2_id": "node2", "count": "value"}
    )


def network(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "contig",
    weight: str = "count",
    out: str | Path | None = None,
    header: bool = False,
    chunksize: int = DEFAULT_CHUNKSIZE,
    force: bool = False,
) -> None:
    """Export a network edge list from a metacool file.

    Args:
        uri: Path to the source ``.metacool`` file.
        level: Aggregation level — ``bin``, ``contig``, or ``genome``.
        weight: Column name in the pixels table to use as edge weight.
        out: Output file path. If None, writes to stdout.
        header: If ``True``, write a ``node1\\tnode2\\tvalue`` header line.
        chunksize: Number of pixels to process per chunk.
        force: Overwrite *out* if it already exists.

    Raises:
        ValueError: If *weight* column does not exist in the pixels table, or
            if *level* is not one of ``bin``/``contig``/``genome``.
        FileExistsError: If *out* exists and *force* is False.
    """
    from metacooler.api import MetaCooler

    mc = MetaCooler(uri)

    # Validate weight column
    with mc.open() as f:
        pixel_cols = list(f["pixels"].keys())
    if weight not in pixel_cols:
        raise ValueError(
            f"Weight column {weight!r} not found in pixels table. "
            f"Available columns: {pixel_cols}"
        )
    if level not in ("bin", "contig", "genome"):
        raise ValueError(f"level must be 'bin', 'contig', or 'genome', got {level!r}")

    if out is not None:
        check_output_path(out, force, inputs=[uri])

    n_pix = len(mc._clr.pixels())

    if level == "bin":
        # Raw pixel edge list — stream chunks directly to output file.
        # Wrap stdout in nullcontext so the ``with`` block does not close it
        # (closing sys.stdout would break any later output in the process).
        cm = open(out, "w") if out else contextlib.nullcontext(sys.stdout)  # noqa: SIM115
        with cm as fout:
            if header:
                fout.write("node1\tnode2\tvalue\n")
            for lo in range(0, n_pix, chunksize):
                chunk = mc._clr.pixels()[lo : lo + chunksize]
                df = pd.DataFrame(
                    {
                        "node1": chunk["bin1_id"].values,
                        "node2": chunk["bin2_id"].values,
                        "value": chunk[weight].values,
                    }
                )
                fout.write(df.to_csv(sep="\t", index=False, header=False))
        return

    # Contig and genome levels: map every bin to its node, then reduce.
    bin_to_contig = chrom_codes(mc.bins())
    if level == "contig":
        bin_to_node = bin_to_contig
        node_names = mc.contigs["name"].values
    else:
        contig_to_genome = mc.contigs["genome_id"].values
        bin_to_node = contig_to_genome[bin_to_contig]
        node_names = mc.genomes["name"].values

    agg = _aggregate_pairs(mc, bin_to_node, weight, n_pix, chunksize)
    if agg is None:
        edges = pd.DataFrame(columns=["node1", "node2", "value"])
    else:
        edges = pd.DataFrame(
            {
                "node1": node_names[agg["node1"].values],
                "node2": node_names[agg["node2"].values],
                "value": agg["value"].values,
            }
        )

    logger.info("Generated %d edges at %s level", len(edges), level)

    # Write TSV (header optional for compatibility with igraph, networkx, etc.)
    if out:
        edges.to_csv(out, sep="\t", index=False, header=header)
    else:
        edges.to_csv(sys.stdout, sep="\t", index=False, header=header)
