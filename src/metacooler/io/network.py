"""Network edge list export from metacool files."""

from __future__ import annotations

import contextlib
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.io._utils import check_output_path, chrom_codes

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")

#: Number of per-chunk groupby frames buffered before they are folded into the
#: running aggregate.  Keeps peak memory at O(FOLD_EVERY x chunk edges) instead
#: of growing one frame per chunk for the whole pixel table.
_FOLD_EVERY = 32


def _fold(frames: list[pd.DataFrame], *, sort: bool = False) -> pd.DataFrame:
    """Concatenate per-chunk edge frames and re-sum duplicated node pairs.

    Args:
        frames: Non-empty list of frames with ``node1``/``node2``/``value``.
        sort: Sort the output by node pair (only worth paying once, at the end).

    Returns:
        One frame with a single row per ``(node1, node2)`` pair.
    """
    combined = frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True)
    return combined.groupby(["node1", "node2"], sort=sort, as_index=False)[
        "value"
    ].sum()


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
        n1 = bin_to_node[chunk["bin1_id"].values].copy()
        n2 = bin_to_node[chunk["bin2_id"].values].copy()

        swap = n1 > n2
        n1[swap], n2[swap] = n2[swap].copy(), n1[swap].copy()

        df = pd.DataFrame({"node1": n1, "node2": n2, "value": chunk[weight].values})
        pending.append(
            df.groupby(["node1", "node2"], sort=False, as_index=False)["value"].sum()
        )
        if len(pending) >= _FOLD_EVERY:
            running = _fold(pending if running is None else [running, *pending])
            pending = []

    if running is None and not pending:
        return None
    return _fold(pending if running is None else [running, *pending], sort=True)


def network(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "contig",
    weight: str = "count",
    out: str | Path | None = None,
    header: bool = False,
    chunksize: int = 1_000_000,
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
        check_output_path(out, force)

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
