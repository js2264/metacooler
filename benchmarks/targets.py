#!/usr/bin/env python
"""Pick what the query benchmarks ask for: one contig, one MAG and one MAG pair.

    python targets.py --contigs contigs.tsv --sizes contigs.sizes --genomes genomes.tsv

prints ``CONTIG GENOME GENOME_A GENOME_B`` on one line, for a bash ``read``.

The MAG is drawn uniformly from the binned genomes of at least 500 kb, and the
contig from that MAG. The pair is drawn in proportion to abundance x length, as
the simulated library joins two genomes, so that it shares contacts: two MAGs
drawn uniformly from a thousand share one contact in the median, and a query
that returns nothing checks nothing. The seed is fixed, so every render asks
the same questions of the same dataset.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def pick(table, sizes, abundance, seed=0):
    """One contig, one MAG (>= 500 kb) and one pair of MAGs."""
    rng = np.random.default_rng(seed)
    t = table.join(sizes, on="contig_name")
    t = t[t["genome_name"] != "__unbinned__"]
    length = t.groupby("genome_name")["length"].sum()
    mags = length.index[length >= 500_000].to_numpy()
    if len(mags) < 2:  # a tiny community: its two largest genomes
        mags = length.sort_values().index[-2:].to_numpy()
    genome = rng.choice(mags)
    contig = rng.choice(t.loc[t["genome_name"] == genome, "contig_name"].to_numpy())
    weight = abundance.reindex(mags).to_numpy() * length[mags].to_numpy()
    a, b = rng.choice(mags, 2, replace=False, p=weight / weight.sum())
    return contig, genome, a, b


def main(argv=None):
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--contigs", required=True, type=Path, help="binning table")
    parser.add_argument("--sizes", required=True, type=Path, help="contig lengths")
    parser.add_argument(
        "--genomes", required=True, type=Path, help="genome table, with abundance"
    )
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    table = pd.read_csv(args.contigs, sep="\t", usecols=[0, 1], dtype=str)
    sizes = pd.read_csv(
        args.sizes, sep="\t", header=None, names=["contig_name", "length"]
    ).set_index("contig_name")
    abundance = pd.read_csv(args.genomes, sep="\t", index_col="genome_name")[
        "abundance"
    ]
    print(*pick(table, sizes, abundance, args.seed))


if __name__ == "__main__":
    main()
