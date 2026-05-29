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
        # Raw pixel edge list — stream chunks directly to output file
        with open(out, "w") if out else sys.stdout as fout:
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

    elif level == "contig":
        bins = mc.bins()
        bin_to_contig = bins["chrom"].cat.codes.values

        partial_frames: list[pd.DataFrame] = []
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            c1 = bin_to_contig[chunk["bin1_id"].values].copy()
            c2 = bin_to_contig[chunk["bin2_id"].values].copy()
            vals = chunk[weight].values

            swap = c1 > c2
            c1[swap], c2[swap] = c2[swap].copy(), c1[swap].copy()

            df = pd.DataFrame({"node1": c1, "node2": c2, "value": vals})
            partial_frames.append(
                df.groupby(["node1", "node2"], sort=False, as_index=False)[
                    "value"
                ].sum()
            )

        if partial_frames:
            combined = pd.concat(partial_frames, ignore_index=True)
            agg = combined.groupby(["node1", "node2"], sort=True, as_index=False)[
                "value"
            ].sum()
            contigs = mc.contigs
            contig_names = contigs["name"].values
            edges = pd.DataFrame(
                {
                    "node1": contig_names[agg["node1"].values],
                    "node2": contig_names[agg["node2"].values],
                    "value": agg["value"].values,
                }
            )
        else:
            edges = pd.DataFrame(columns=["node1", "node2", "value"])

    elif level == "genome":
        bins = mc.bins()
        contigs = mc.contigs
        bin_to_contig = bins["chrom"].cat.codes.values
        contig_to_genome = contigs["genome_id"].values

        partial_frames2: list[pd.DataFrame] = []
        for lo in range(0, n_pix, chunksize):
            chunk = mc._clr.pixels()[lo : lo + chunksize]
            g1 = contig_to_genome[bin_to_contig[chunk["bin1_id"].values]].copy()
            g2 = contig_to_genome[bin_to_contig[chunk["bin2_id"].values]].copy()
            vals = chunk[weight].values

            swap = g1 > g2
            g1[swap], g2[swap] = g2[swap].copy(), g1[swap].copy()

            df = pd.DataFrame({"node1": g1, "node2": g2, "value": vals})
            partial_frames2.append(
                df.groupby(["node1", "node2"], sort=False, as_index=False)[
                    "value"
                ].sum()
            )

        if partial_frames2:
            combined2 = pd.concat(partial_frames2, ignore_index=True)
            agg2 = combined2.groupby(["node1", "node2"], sort=True, as_index=False)[
                "value"
            ].sum()
            genomes = mc.genomes
            genome_names = genomes["name"].values
            edges = pd.DataFrame(
                {
                    "node1": genome_names[agg2["node1"].values],
                    "node2": genome_names[agg2["node2"].values],
                    "value": agg2["value"].values,
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
