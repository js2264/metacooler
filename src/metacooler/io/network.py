"""Network edge list export from metacool files."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

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
) -> None:
    """Export a network edge list from a metacool file.

    Args:
        uri: Path to the source ``.metacool`` file.
        level: Aggregation level — ``bin``, ``contig``, or ``genome``.
        weight: Column name in the pixels table to use as edge weight.
        out: Output file path. If None, writes to stdout.
        header: If ``True``, write a ``node1\\tnode2\\tvalue`` header line.

    Raises:
        ValueError: If *weight* column does not exist in the pixels table.
    """
    from metacooler.api import MetaCooler

    mc = MetaCooler(uri)
    pixels = mc.pixels()

    # Validate weight column
    if weight not in pixels.columns:
        raise ValueError(
            f"Weight column {weight!r} not found in pixels table. "
            f"Available columns: {list(pixels.columns)}"
        )

    if level == "bin":
        # Raw pixel edge list
        edges = pd.DataFrame(
            {
                "node1": pixels["bin1_id"].values,
                "node2": pixels["bin2_id"].values,
                "value": pixels[weight].values,
            }
        )
    elif level == "contig":
        bins = mc.bins()
        bin_to_contig = bins["chrom"].cat.codes.values

        c1 = bin_to_contig[pixels["bin1_id"].values]
        c2 = bin_to_contig[pixels["bin2_id"].values]
        vals = pixels[weight].values

        # Enforce upper-triangular at the aggregated level
        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap].copy(), c1[swap].copy()

        edge_df = pd.DataFrame({"node1": c1, "node2": c2, "value": vals})
        edges = edge_df.groupby(["node1", "node2"], as_index=False)["value"].sum()

        # Map contig IDs to names
        contigs = mc.contigs
        contig_names = contigs["name"].values
        edges["node1"] = contig_names[edges["node1"].values]
        edges["node2"] = contig_names[edges["node2"].values]

    elif level == "genome":
        bins = mc.bins()
        contigs = mc.contigs
        bin_to_contig = bins["chrom"].cat.codes.values
        contig_to_genome = contigs["genome_id"].values

        g1 = contig_to_genome[bin_to_contig[pixels["bin1_id"].values]]
        g2 = contig_to_genome[bin_to_contig[pixels["bin2_id"].values]]
        vals = pixels[weight].values

        # Enforce upper-triangular
        swap = g1 > g2
        g1[swap], g2[swap] = g2[swap].copy(), g1[swap].copy()

        edge_df = pd.DataFrame({"node1": g1, "node2": g2, "value": vals})
        edges = edge_df.groupby(["node1", "node2"], as_index=False)["value"].sum()

        # Map genome IDs to names
        genomes = mc.genomes
        genome_names = genomes["name"].values
        edges["node1"] = genome_names[edges["node1"].values]
        edges["node2"] = genome_names[edges["node2"].values]
    else:
        raise ValueError(f"level must be 'bin', 'contig', or 'genome', got {level!r}")

    logger.info("Generated %d edges at %s level", len(edges), level)

    # Write TSV (header optional for compatibility with igraph, networkx, etc.)
    if out:
        edges.to_csv(out, sep="\t", index=False, header=header)
    else:
        edges.to_csv(sys.stdout, sep="\t", index=False, header=header)
