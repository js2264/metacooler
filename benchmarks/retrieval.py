#!/usr/bin/env python
"""Time contig, genome and genome-pair queries in-process, three ways.

    python retrieval.py --cool M.cool --metacool M.metacool \
        --contigs contigs.tsv --sizes contigs.sizes --genomes genomes.tsv

One query is one contact block, returned as the stored (upper-triangle) pixels:
a contig against itself, a genome against itself, or two genomes against each
other. Three ways of answering it:

``cool``
    The genome-agnostic ``.cool`` and the binning table beside it. A genome's
    contigs are found by scanning the table; each contig's bins by cooler's
    ``chrom_offset`` index; its pixels by ``bin1_offset``, one row block per
    contig.
``metacooler``
    The ``metacooler.MetaCooler`` API, as released: ``pixels_for_genome_pair``
    for genomes and genome pairs, ``contig_extent`` and a matrix slice for
    contigs.
``index``
    The ``.metacool``'s own indexes, read with h5py: ``genome_offset`` gives a
    genome's contig range, ``chrom_offset`` its bin range, ``bin1_offset`` the
    one row block of the pixel table that holds it.

``cool`` and ``index`` share the code that reads row blocks, so they differ in
one thing only: how a genome's bins are found. Every query is answered by all
three, which must agree on its pixel count and contact sum. The first query of
each kind is run once untimed (lazy tables, caches), then every query is timed
once. Writes one TSV row per query and way, plus an ``open`` row per way: what
it takes to be ready for the first query (the median of three).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from metacooler.api import MetaCooler


def blocks(grp, bin_ranges, keep):
    """Stored pixels of the rows in ``bin_ranges`` whose bin2 passes ``keep``."""
    bin1_offset = grp["indexes/bin1_offset"]
    nnz, total = 0, 0
    for lo, hi in bin_ranges:
        p0, p1 = int(bin1_offset[lo]), int(bin1_offset[hi])
        if p0 == p1:
            continue
        bin2 = grp["pixels/bin2_id"][p0:p1]
        sel = keep(bin2)
        nnz += int(sel.sum())
        total += int(grp["pixels/count"][p0:p1][sel].sum())
    return nnz, total


def in_ranges(bin_ranges):
    """Membership test for a union of bin ranges, as a vectorised function."""
    if not bin_ranges:
        return lambda b: np.zeros(len(b), bool)
    starts, ends = np.array(bin_ranges).T
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]

    def keep(b):
        i = np.searchsorted(starts, b, side="right") - 1
        return (i >= 0) & (b < ends[np.maximum(i, 0)])

    return keep


class Cool:
    """A .cool and its binning table: the genome is looked up in the table."""

    name = "cool"

    def __init__(self, path, table):
        self.f = h5py.File(path, "r")
        names = self.f["chroms/name"][:].astype(str)
        self.chrom_id = dict(zip(names, range(len(names)), strict=True))
        t = pd.read_csv(table, sep="\t", usecols=[0, 1], dtype=str)
        self.table_contig = t.iloc[:, 0].to_numpy()
        self.table_genome = t.iloc[:, 1].to_numpy()

    def _contig_bins(self, contig_ids):
        offset = self.f["indexes/chrom_offset"]
        return [(int(offset[c]), int(offset[c + 1])) for c in sorted(contig_ids)]

    def _genome_bins(self, genome):
        contigs = self.table_contig[self.table_genome == genome]  # linear scan
        return self._contig_bins(self.chrom_id[c] for c in contigs)

    def contig(self, name):
        ranges = self._contig_bins([self.chrom_id[name]])
        return blocks(self.f, ranges, in_ranges(ranges))

    def genome(self, name):
        ranges = self._genome_bins(name)
        return blocks(self.f, ranges, in_ranges(ranges))

    def genome_pair(self, a, b):
        ra, rb = self._genome_bins(a), self._genome_bins(b)
        na, ta = blocks(self.f, ra, in_ranges(rb))
        nb, tb = blocks(self.f, rb, in_ranges(ra))
        return na + nb, ta + tb


class Index:
    """A .metacool read through its own indexes."""

    name = "index"

    def __init__(self, path, table):
        self.f = h5py.File(path, "r")
        names = self.f["genomes/name"][:].astype(str)
        self.genome_id = dict(zip(names, range(len(names)), strict=True))
        self._contig_id = None

    def _genome_bins(self, name):
        g = self.genome_id[name]
        c0, c1 = (int(x) for x in self.f["indexes/genome_offset"][g : g + 2])
        offset = self.f["indexes/chrom_offset"]
        return int(offset[c0]), int(offset[c1])

    def contig(self, name):
        if self._contig_id is None:  # contig names, read once (as cooler does)
            names = self.f["contigs/name"][:].astype(str)
            self._contig_id = dict(zip(names, range(len(names)), strict=True))
        c = self._contig_id[name]
        lo, hi = (int(x) for x in self.f["indexes/chrom_offset"][c : c + 2])
        return blocks(self.f, [(lo, hi)], lambda b: b < hi)

    def genome(self, name):
        lo, hi = self._genome_bins(name)
        return blocks(self.f, [(lo, hi)], lambda b: b < hi)

    def genome_pair(self, a, b):
        (a0, a1), (b0, b1) = sorted([self._genome_bins(a), self._genome_bins(b)])
        return blocks(self.f, [(a0, a1)], lambda x: (x >= b0) & (x < b1))


class API:
    """metacooler.MetaCooler, as released."""

    name = "metacooler"

    def __init__(self, path, table):
        self.mc = MetaCooler(path)

    @staticmethod
    def _sum(px):
        return len(px), int(px["count"].sum())

    def contig(self, name):
        lo, hi = self.mc.contig_extent(name)
        px = self.mc.matrix(balance=False, as_pixels=True)[lo:hi, lo:hi]
        return self._sum(px[px["bin1_id"] <= px["bin2_id"]])

    def genome(self, name):
        return self._sum(self.mc.pixels_for_genome_pair(name, name))

    def genome_pair(self, a, b):
        return self._sum(self.mc.pixels_for_genome_pair(a, b))


def pick_queries(table, sizes, abundance, n, seed):
    """``n`` genomes of at least 500 kb, one contig of each, and ``n`` pairs.

    The genomes are drawn uniformly. The pairs are drawn in proportion to
    abundance x length, as the simulated library joins two genomes, so that a
    pair shares contacts: two MAGs drawn uniformly from a thousand have one in
    the median, and a query that returns nothing checks nothing.
    """
    rng = np.random.default_rng(seed)
    t = table.join(sizes, on="contig_name")
    t = t[t["genome_name"] != "__unbinned__"]
    length = t.groupby("genome_name")["length"].sum()
    mags = length.index[length >= 500_000].to_numpy()
    if len(mags) == 0:
        mags = length.sort_values().index[-n:].to_numpy()
    genomes = rng.choice(mags, size=min(n, len(mags)), replace=False)
    contigs = [rng.choice(t.loc[t["genome_name"] == g, "contig_name"]) for g in genomes]
    weight = abundance.reindex(mags).to_numpy() * length[mags].to_numpy()
    weight = weight / weight.sum()
    pairs = [tuple(rng.choice(mags, 2, replace=False, p=weight)) for _ in range(n)]
    return (
        [("contig", c) for c in contigs]
        + [("genome", g) for g in genomes]
        + [("genome_pair", p) for p in pairs]
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--cool", type=Path)
    parser.add_argument("--metacool", type=Path)
    parser.add_argument("--contigs", required=True, type=Path, help="binning table")
    parser.add_argument("--sizes", required=True, type=Path, help="contig lengths")
    parser.add_argument(
        "--genomes", required=True, type=Path, help="genome table, with abundance"
    )
    parser.add_argument("--dataset", default="")
    parser.add_argument("--queries", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--targets",
        action="store_true",
        help="only print one contig, one genome and a genome pair, and exit",
    )
    args = parser.parse_args(argv)
    if not args.targets and not (args.cool and args.metacool):
        parser.error("--cool and --metacool are required, unless --targets")

    table = pd.read_csv(args.contigs, sep="\t", usecols=[0, 1], dtype=str)
    sizes = pd.read_csv(
        args.sizes, sep="\t", header=None, names=["contig_name", "length"]
    ).set_index("contig_name")
    abundance = pd.read_csv(args.genomes, sep="\t", index_col="genome_name")[
        "abundance"
    ]
    queries = pick_queries(table, sizes, abundance, args.queries, args.seed)
    if args.targets:
        first = {}
        for query, target in queries:
            first.setdefault(query, target)
        print(first["contig"], first["genome"], *first["genome_pair"])
        return
    n_contigs = len(table)
    n_genomes = table.loc[
        table["genome_name"] != "__unbinned__", "genome_name"
    ].nunique()
    del table, sizes, abundance

    rows, answers = [], {}
    paths = {"cool": args.cool, "index": args.metacool, "metacooler": args.metacool}
    for kind in (Cool, Index, API):
        opened = []
        for _ in range(3):  # the median of three opens
            t0 = time.perf_counter()
            way = kind(paths[kind.name], args.contigs)
            opened.append(time.perf_counter() - t0)
        rows.append((kind.name, "open", "", np.nan, np.nan, float(np.median(opened))))
        warm = set()
        for query, target in queries:
            args_ = target if isinstance(target, tuple) else (target,)
            run = getattr(way, query)
            if query not in warm:  # once per kind of query: lazy tables, caches
                warm.add(query)
                run(*args_)
            t0 = time.perf_counter()
            nnz, total = run(*args_)
            seconds = time.perf_counter() - t0
            label = "|".join(args_)
            rows.append((kind.name, query, label, nnz, total, seconds))
            expected = answers.setdefault((query, label), (nnz, total))
            if (nnz, total) != expected:
                sys.exit(
                    f"{kind.name} disagrees on {query} {label}: "
                    f"{(nnz, total)} != {expected}"
                )

    frame = pd.DataFrame(
        rows, columns=["way", "query", "target", "nnz", "sum", "seconds"]
    )
    frame.insert(0, "dataset", args.dataset)
    frame.insert(1, "contigs", n_contigs)
    frame.insert(2, "genomes", n_genomes)
    frame.to_csv(args.output or sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
