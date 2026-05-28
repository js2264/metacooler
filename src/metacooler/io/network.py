"""Network edge list export from metacool files."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")


def network(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "contig",
    weight: str = "count",
    out: str | Path | None = None,
    header: bool = False,
    chunksize: int = 1_000_000,
) -> None:
    """Export a network edge list from a metacool file.

    Args:
        uri: Path to the source ``.metacool`` file.
        level: Aggregation level — ``bin``, ``contig``, or ``genome``.
        weight: Column name in the pixels table to use as edge weight.
        out: Output file path. If None, writes to stdout.
        header: If ``True``, write a ``node1\\tnode2\\tvalue`` header line.
        chunksize: Number of pixels to process per chunk.

    Raises:
        ValueError: If *weight* column does not exist in the pixels table.
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

    n_pix = len(mc._clr.pixels())

    if level == "bin":
        # Raw pixel edge list — stream chunks directly
        rows = []
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            rows.append(
                pd.DataFrame(
                    {
                        "node1": chunk["bin1_id"].values,
                        "node2": chunk["bin2_id"].values,
                        "value": chunk[weight].values,
                    }
                )
            )
        edges = (
            pd.concat(rows, ignore_index=True)
            if rows
            else pd.DataFrame(columns=["node1", "node2", "value"])
        )

    elif level == "contig":
        bins = mc.bins()
        bin_to_contig = bins["chrom"].cat.codes.values

        # Incremental aggregation using a dict accumulator to avoid OOM
        accumulator: dict[tuple[int, int], int] = {}
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            c1 = bin_to_contig[chunk["bin1_id"].values].copy()
            c2 = bin_to_contig[chunk["bin2_id"].values].copy()
            vals = chunk[weight].values

            swap = c1 > c2
            c1[swap], c2[swap] = c2[swap].copy(), c1[swap].copy()

            for i in range(len(c1)):
                key = (int(c1[i]), int(c2[i]))
                accumulator[key] = accumulator.get(key, 0) + int(vals[i])

        if accumulator:
            keys = np.array(list(accumulator.keys()), dtype=np.int64)
            vals_arr = np.array(list(accumulator.values()), dtype=np.int64)
            contigs = mc.contigs
            contig_names = contigs["name"].values
            edges = pd.DataFrame(
                {
                    "node1": contig_names[keys[:, 0]],
                    "node2": contig_names[keys[:, 1]],
                    "value": vals_arr,
                }
            )
        else:
            edges = pd.DataFrame(columns=["node1", "node2", "value"])

    elif level == "genome":
        bins = mc.bins()
        contigs = mc.contigs
        bin_to_contig = bins["chrom"].cat.codes.values
        contig_to_genome = contigs["genome_id"].values

        accumulator2: dict[tuple[int, int], int] = {}
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            g1 = contig_to_genome[bin_to_contig[chunk["bin1_id"].values]].copy()
            g2 = contig_to_genome[bin_to_contig[chunk["bin2_id"].values]].copy()
            vals = chunk[weight].values

            swap = g1 > g2
            g1[swap], g2[swap] = g2[swap].copy(), g1[swap].copy()

            for i in range(len(g1)):
                key = (int(g1[i]), int(g2[i]))
                accumulator2[key] = accumulator2.get(key, 0) + int(vals[i])

        if accumulator2:
            keys = np.array(list(accumulator2.keys()), dtype=np.int64)
            vals_arr = np.array(list(accumulator2.values()), dtype=np.int64)
            genomes = mc.genomes
            genome_names = genomes["name"].values
            edges = pd.DataFrame(
                {
                    "node1": genome_names[keys[:, 0]],
                    "node2": genome_names[keys[:, 1]],
                    "value": vals_arr,
                }
            )
        else:
            edges = pd.DataFrame(columns=["node1", "node2", "value"])

    else:
        raise ValueError(f"level must be 'bin', 'contig', or 'genome', got {level!r}")

    logger.info("Generated %d edges at %s level", len(edges), level)

    # Write TSV (header optional for compatibility with igraph, networkx, etc.)
    if out:
        edges.to_csv(out, sep="\t", index=False, header=header)
    else:
        edges.to_csv(sys.stdout, sep="\t", index=False, header=header)
