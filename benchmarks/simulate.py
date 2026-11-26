#!/usr/bin/env python
"""Simulate a metagenomic Hi-C dataset for the metacooler benchmarks.

One community of genomes, one assembly of it cut into contigs, one binning of
those contigs into genomes, and one Hi-C library sequenced from it — all drawn
from a preset and a seed, so the same command always writes the same data.

    python simulate.py --preset M --outdir data/M [--keep]

writes, in OUTDIR (``--keep`` leaves a directory alone when it already holds
the preset, made by this version of the generator):

    contigs.sizes      contig, length; in assembly order (cooler's chromsizes)
    contigs.tsv        contig_name, genome_name, gc, coverage  (``--contigs``)
    genomes.tsv        genome_name and per-genome metadata     (``--genomes``)
    contacts.pairs.gz  4DN pairs: upper triangle, sorted, ``#chromsize:`` lines
    summary.json       what the dataset is: sizes, N50s, contact categories

The model, in three parts:

* **Genomes.** ``mags`` bacterial genomes (2.5 Mb median, log-normal) and
  ``genomes - mags`` small bins standing for plasmids, phages and other mobile
  elements (20 kb median) are *binned*. The rest of the assembly belongs to
  low-abundance genomes no binner recovered: their contigs are left out of the
  binning and land in ``__unbinned__``.
* **Contigs.** Every genome is cut into contigs with log-normal lengths around a
  per-genome median, so some MAGs come out in a few long contigs and others in
  hundreds of short ones; unbinned contigs are short (~0.9 kb median), as they
  are in real assemblies. Nothing is shorter than 500 bp, the usual assembly
  cut-off. The assembly lists the contigs in random order, so a genome's
  contigs are scattered through it, as in any assembler's output.
* **Contacts.** A read pair picks a genome in proportion to abundance x length.
  96 % of pairs stay in that genome: 60 % at a separation drawn from
  P(s) ~ 1/s (200 bp up to half the chromosome), 40 % uniformly along it, on
  a circular chromosome whose contigs lie in random order. The other 4 % join
  two independently drawn genomes (random ligation). The share of intra-contig,
  inter-contig and inter-genome pairs then follows from the contig lengths, as
  it does in a real library.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pa_csv

#: Bumped whenever the model or the files change, so that a cached dataset made
#: by an older generator is never mistaken for the current one.
GENERATOR_VERSION = 1

#: name: binned genomes, of which MAGs; contigs in the assembly; read pairs;
#: range of the per-MAG median contig length (log-uniform); seed.
PRESETS: dict[str, dict] = {
    "S": dict(genomes=10, mags=10, contigs=1_000, pairs=1_000_000,
              mag_contig_median=(20_000, 100_000), seed=1),
    "M": dict(genomes=100, mags=50, contigs=20_000, pairs=10_000_000,
              mag_contig_median=(4_000, 40_000), seed=2),
    "L": dict(genomes=2_000, mags=200, contigs=200_000, pairs=30_000_000,
              mag_contig_median=(2_000, 30_000), seed=3),
    "XL": dict(genomes=20_000, mags=1_000, contigs=2_000_000, pairs=100_000_000,
               mag_contig_median=(2_000, 30_000), seed=4),
}  # fmt: skip

MIN_CONTIG = 500  # bp, the assembly's length cut-off
CONTIG_SIGMA = 1.0  # log-normal spread of contig lengths within a genome
MAG_LENGTH = (2.5e6, 0.5, 4e5, 8e6)  # median, sigma, min, max (bp)
MGE_LENGTH = (2e4, 1.0, 1e3, 5e5)
MGE_CONTIG_MEDIAN = (1_000, 10_000)
UNBINNED_LENGTH = (400, 0.9, 2e5)  # median above MIN_CONTIG, sigma, max (bp)
HIDDEN_GENOME_CONTIGS = 40  # mean contigs per unbinned (hidden) genome
P_TRANS = 0.04  # pairs joining two independently drawn genomes
P_LOCAL = 0.60  # cis pairs at a P(s) ~ 1/s separation; the others uniform
S_MIN = 200  # bp, shortest separation of the 1/s component
CHUNK = 5_000_000  # pairs drawn per batch (fixed: the draws depend on it)
WRITE_CHUNK = 4_000_000  # pairs formatted per batch


def lognormal(rng, median, sigma, n, lo=None, hi=None):
    x = rng.lognormal(np.log(median), sigma, n)
    return np.clip(x, lo, hi) if lo is not None or hi is not None else x


def contig_lengths(rng, target, median, n_contigs=None):
    """Cut each genome of ``target`` bp into contigs around its ``median``."""
    excess = np.maximum(median - MIN_CONTIG, 100.0)
    mean = MIN_CONTIG + excess * np.exp(CONTIG_SIGMA**2 / 2)
    if n_contigs is None:
        n_contigs = np.maximum(1, np.rint(target / mean)).astype(np.int64)
    owner = np.repeat(np.arange(len(target)), n_contigs)
    length = MIN_CONTIG + np.rint(
        rng.lognormal(np.log(excess[owner]), CONTIG_SIGMA)
    ).astype(np.int64)
    return owner, length


def draw_community(rng, preset):
    """The genomes of the community and the contigs of its assembly."""
    n_mag, n_mge = preset["mags"], preset["genomes"] - preset["mags"]

    mag_target = lognormal(rng, *MAG_LENGTH[:2], n_mag, *MAG_LENGTH[2:])
    mag_median = np.exp(rng.uniform(*np.log(preset["mag_contig_median"]), n_mag))
    mge_target = lognormal(rng, *MGE_LENGTH[:2], n_mge, *MGE_LENGTH[2:])
    mge_median = np.exp(rng.uniform(*np.log(MGE_CONTIG_MEDIAN), n_mge))
    mag_owner, mag_len = contig_lengths(rng, mag_target, mag_median)
    mge_owner, mge_len = contig_lengths(rng, mge_target, mge_median)

    n_binned = len(mag_len) + len(mge_len)
    n_unbinned = preset["contigs"] - n_binned
    if n_unbinned < 0:
        raise ValueError(
            f"the binned genomes alone make {n_binned:,} contigs, more than the "
            f"{preset['contigs']:,} asked for: raise 'contigs' or the MAG medians"
        )
    # Unbinned contigs belong to low-abundance genomes: group them into hidden
    # genomes of ~HIDDEN_GENOME_CONTIGS contigs, for the contact model only.
    sizes = 1 + rng.poisson(HIDDEN_GENOME_CONTIGS - 1, n_unbinned // 8 + 1)
    sizes = sizes[: np.searchsorted(np.cumsum(sizes), n_unbinned) + 1]
    sizes[-1] -= sizes.sum() - n_unbinned
    sizes = sizes[sizes > 0]
    n_hidden = len(sizes)
    hid_owner = np.repeat(np.arange(n_hidden), sizes)
    median, sigma, longest = UNBINNED_LENGTH
    hid_len = MIN_CONTIG + np.rint(
        lognormal(rng, median, sigma, n_unbinned, None, longest - MIN_CONTIG)
    ).astype(np.int64)

    width = len(str(max(n_mag, n_mge, 1)))
    genomes = pd.DataFrame(
        {
            "name": [f"MAG_{i + 1:0{width}d}" for i in range(n_mag)]
            + [f"MGE_{i + 1:0{width}d}" for i in range(n_mge)]
            + [f"hidden_{i + 1}" for i in range(n_hidden)],
            "kind": ["MAG"] * n_mag + ["MGE"] * n_mge + ["unbinned"] * n_hidden,
            "abundance": np.concatenate(
                [
                    lognormal(rng, 1.0, 1.2, n_mag),
                    lognormal(rng, 1.0, 1.5, n_mge),
                    lognormal(rng, 0.2, 1.0, n_hidden),
                ]
            ),
            "gc": rng.uniform(0.28, 0.72, n_mag + n_mge + n_hidden),
        }
    )
    contigs = pd.DataFrame(
        {
            "genome": np.concatenate(
                [mag_owner, n_mag + mge_owner, n_mag + n_mge + hid_owner]
            ),
            "length": np.concatenate([mag_len, mge_len, hid_len]),
        }
    )
    genomes["length"] = np.bincount(
        contigs["genome"], weights=contigs["length"], minlength=len(genomes)
    ).astype(np.int64)
    genomes["n_contigs"] = np.bincount(contigs["genome"], minlength=len(genomes))
    return genomes, contigs


def lay_out(rng, genomes, contigs):
    """Place contigs on their chromosomes (random order) and in the assembly."""
    n = len(contigs)
    # The chromosome order of a genome's contigs is unknown: draw it.
    chrom_order = np.lexsort((rng.random(n), contigs["genome"].to_numpy()))
    contigs = contigs.iloc[chrom_order].reset_index(drop=True)
    contigs["start"] = np.r_[0, np.cumsum(contigs["length"].to_numpy())[:-1]]
    genomes["start"] = np.r_[0, np.cumsum(genomes["length"].to_numpy())[:-1]]
    # The assembly lists contigs in an order unrelated to the genomes.
    contigs["asm"] = rng.permutation(n)
    contigs["name"] = [f"k141_{i + 1}" for i in contigs["asm"]]
    return contigs


def annotate(rng, genomes, contigs):
    """Per-genome and per-contig metadata, as a binner and an assembler report."""
    g = genomes
    binned = g["kind"] != "unbinned"
    completeness = np.where(
        g["kind"] == "MAG",
        100 * rng.beta(8, 1.2, len(g)),
        100 * rng.beta(1, 6, len(g)),
    )
    contamination = np.minimum(rng.exponential(2.0, len(g)), 40.0)
    g["completeness"] = completeness.round(1)
    g["contamination"] = contamination.round(1)
    g["quality"] = np.select(
        [(completeness > 90) & (contamination < 5),
         (completeness >= 50) & (contamination < 10)],
        ["high", "medium"],
        "low",
    )  # fmt: skip
    gid = contigs["genome"].to_numpy()
    contigs["gc"] = np.clip(
        rng.normal(g["gc"].to_numpy()[gid], 0.015), 0.15, 0.85
    ).round(3)
    contigs["coverage"] = (
        20 * g["abundance"].to_numpy()[gid] * rng.lognormal(0, 0.25, len(contigs))
    ).round(2)
    contigs["genome_name"] = np.where(
        binned.to_numpy()[gid], g["name"].to_numpy()[gid], "__unbinned__"
    )
    return g, contigs


def draw_pairs(rng, genomes, contigs, n_pairs):
    """Read pairs as (contig1, pos1, contig2, pos2, strands), contigs by rank."""
    weight = genomes["abundance"].to_numpy() * genomes["length"].to_numpy()
    cumw = np.cumsum(weight)
    glen = genomes["length"].to_numpy().astype(np.float64)
    gstart = genomes["start"].to_numpy()
    cstart = contigs["start"].to_numpy()
    asm = contigs["asm"].to_numpy().astype(np.int32)

    out = {k: np.empty(n_pairs, np.int32) for k in ("c1", "p1", "c2", "p2")}
    strands = np.empty((2, n_pairs), np.int8)
    for lo in range(0, n_pairs, CHUNK):
        n = min(CHUNK, n_pairs - lo)
        g1 = np.searchsorted(cumw, rng.random(n) * cumw[-1], side="right")
        g2 = np.searchsorted(cumw, rng.random(n) * cumw[-1], side="right")
        trans = rng.random(n) < P_TRANS
        local = rng.random(n) < P_LOCAL
        half = glen[g1] / 2
        log_lo = np.log(S_MIN)
        log_hi = np.log(np.maximum(half, S_MIN))
        sep = np.where(
            local,
            np.exp(log_lo + rng.random(n) * (log_hi - log_lo)),
            rng.random(n) * half,
        )
        sign = np.where(rng.random(n) < 0.5, -1.0, 1.0)
        x1 = np.floor(rng.random(n) * glen[g1])
        x2 = np.where(
            trans,
            np.floor(rng.random(n) * glen[g2]),
            np.mod(np.floor(x1 + sign * sep), glen[g1]),
        )
        g2 = np.where(trans, g2, g1)
        for end, g, x in (("1", g1, x1), ("2", g2, x2)):
            pos = gstart[g] + x.astype(np.int64)
            idx = np.searchsorted(cstart, pos, side="right") - 1
            out["c" + end][lo : lo + n] = asm[idx]
            out["p" + end][lo : lo + n] = pos - cstart[idx] + 1  # 1-based
        strands[:, lo : lo + n] = rng.integers(0, 2, (2, n), dtype=np.int8)

    # Upper triangle in assembly order, then sorted as `pairtools sort` would.
    c1, p1, c2, p2 = out.pop("c1"), out.pop("p1"), out.pop("c2"), out.pop("p2")
    flip = (c1 > c2) | ((c1 == c2) & (p1 > p2))
    c1[flip], c2[flip] = c2[flip], c1[flip]
    p1[flip], p2[flip] = p2[flip], p1[flip]
    strands[0, flip], strands[1, flip] = strands[1, flip], strands[0, flip]
    del flip
    order = np.lexsort((p2, p1, c2, c1))
    c1, p1, c2, p2 = c1[order], p1[order], c2[order], p2[order]
    return c1, p1, c2, p2, strands[:, order]


def write_pairs(path, name, sizes, pairs):
    """4DN pairs, gzipped with pigz when there is one; md5 of the plain text."""
    c1, p1, c2, p2, strands = pairs
    names = pa.array(sizes.index.to_numpy(dtype=object), pa.string())
    header = [
        "## pairs format v1.0",
        "#sorted: chr1-chr2-pos1-pos2",
        "#shape: upper triangle",
        f"#genome_assembly: {name}",
        *(f"#chromsize: {c} {n}" for c, n in sizes.items()),
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2",
    ]
    md5 = hashlib.md5()
    options = pa_csv.WriteOptions(
        include_header=False, delimiter="\t", quoting_style="none"
    )
    gzip = ["pigz", "-p", "4", "-c"] if shutil.which("pigz") else ["gzip", "-c"]
    with open(path, "wb") as handle:
        proc = subprocess.Popen(gzip, stdin=subprocess.PIPE, stdout=handle)
        assert proc.stdin is not None
        head = ("\n".join(header) + "\n").encode()
        md5.update(head)
        proc.stdin.write(head)
        dot = pa.array(["."])
        sign = pa.array(["+", "-"])
        for lo in range(0, len(c1), WRITE_CHUNK):
            hi = min(lo + WRITE_CHUNK, len(c1))
            n = hi - lo
            table = pa.table(
                {
                    "readID": pa.DictionaryArray.from_arrays(
                        pa.array(np.zeros(n, np.int8)), dot
                    ),
                    "chrom1": names.take(pa.array(c1[lo:hi])),
                    "pos1": pa.array(p1[lo:hi]),
                    "chrom2": names.take(pa.array(c2[lo:hi])),
                    "pos2": pa.array(p2[lo:hi]),
                    "strand1": pa.DictionaryArray.from_arrays(
                        pa.array(strands[0, lo:hi]), sign
                    ),
                    "strand2": pa.DictionaryArray.from_arrays(
                        pa.array(strands[1, lo:hi]), sign
                    ),
                }
            )
            sink = pa.BufferOutputStream()
            pa_csv.write_csv(table, sink, write_options=options)
            block = sink.getvalue()
            md5.update(memoryview(block))
            proc.stdin.write(memoryview(block))
        proc.stdin.close()
        if proc.wait():
            raise RuntimeError(f"{gzip[0]} failed with status {proc.returncode}")
    return md5.hexdigest()


def n50(lengths):
    x = np.sort(np.asarray(lengths))[::-1]
    c = np.cumsum(x)
    i = int(np.searchsorted(c, c[-1] / 2))
    return int(x[i]), i + 1


def per_genome_n50(genome, length):
    """N50 of each genome's contigs, as a Series indexed by genome."""
    frame = pd.DataFrame({"genome": genome, "length": length})
    frame = frame.sort_values(["genome", "length"], ascending=[True, False])
    by = frame.groupby("genome")["length"]
    reached = by.cumsum() >= by.transform("sum") / 2
    return frame[reached.to_numpy()].groupby("genome")["length"].first()


def summarize(preset_name, preset, genomes, contigs, pairs):
    c1, _, c2, _, _ = pairs
    kind = genomes["kind"].to_numpy()
    genome_of = contigs.sort_values("asm")["genome"].to_numpy()  # by rank
    binned_of = kind[genome_of] != "unbinned"
    g1, g2, b1, b2 = genome_of[c1], genome_of[c2], binned_of[c1], binned_of[c2]
    n = len(c1)
    binned = genomes[kind != "unbinned"]
    mags = genomes[kind == "MAG"]
    lengths = contigs["length"].to_numpy()
    in_mag = kind[contigs["genome"].to_numpy()] == "MAG"
    mag_n50 = per_genome_n50(contigs["genome"].to_numpy()[in_mag], lengths[in_mag])
    contig_n50, contig_l50 = n50(lengths)
    return {
        "dataset": preset_name,
        "seed": preset["seed"],
        "generator_version": GENERATOR_VERSION,
        "genomes": len(binned),
        "mags": len(mags),
        "mge_bins": int((kind == "MGE").sum()),
        "contigs": len(contigs),
        "binned_contigs": int(binned_of.sum()),
        "assembly_bp": int(lengths.sum()),
        "binned_bp": int(
            lengths[kind[contigs["genome"].to_numpy()] != "unbinned"].sum()
        ),
        "contig_n50": contig_n50,
        "contig_l50": contig_l50,
        "contig_median": int(np.median(lengths)),
        "contig_max": int(lengths.max()),
        "contigs_per_genome_median": float(binned["n_contigs"].median()),
        "contigs_per_genome_mean": float(binned["n_contigs"].mean()),
        "contigs_per_genome_max": int(binned["n_contigs"].max()),
        "contigs_per_mag_median": float(mags["n_contigs"].median()),
        "mag_bp_median": int(mags["length"].median()),
        "mag_bp_max": int(mags["length"].max()),
        "mag_contig_n50_median": int(mag_n50.median()),
        "pairs": int(n),
        "pairs_intra_contig": float((c1 == c2).mean()),
        "pairs_inter_contig_intra_genome": float(
            ((c1 != c2) & (g1 == g2) & b1 & b2).sum() / n
        ),
        "pairs_inter_genome": float(((g1 != g2) & b1 & b2).sum() / n),
        "pairs_unbinned": float(((c1 != c2) & ~(b1 & b2)).sum() / n),
    }


FILES = ("contigs.sizes", "contigs.tsv", "genomes.tsv", "contacts.pairs.gz")


def up_to_date(out, name, preset):
    """Whether ``out`` holds this preset, made by this generator, files intact."""
    try:
        summary = json.loads((out / "summary.json").read_text())
    except (OSError, ValueError):
        return False
    made = {k: summary.get(k) for k in ("dataset", "generator_version", "seed")}
    if made != {"dataset": name, "generator_version": GENERATOR_VERSION,
                "seed": preset["seed"]}:  # fmt: skip
        return False
    if (summary.get("contigs"), summary.get("pairs")) != (
        preset["contigs"],
        preset["pairs"],
    ):
        return False
    sizes = summary.get("files_bytes", {})
    return all(
        (out / f).is_file() and (out / f).stat().st_size == sizes.get(f) for f in FILES
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--preset", choices=PRESETS, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pairs", type=int, help="override the preset's pair count")
    parser.add_argument(
        "--keep", action="store_true", help="leave an up-to-date OUTDIR alone"
    )
    args = parser.parse_args(argv)

    preset = dict(PRESETS[args.preset])
    if args.pairs:
        preset["pairs"] = args.pairs
    out = args.outdir
    if args.keep and up_to_date(out, args.preset, preset):
        print(f"{args.preset}: {out} is up to date (generator v{GENERATOR_VERSION})")
        return
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(preset["seed"])
    t0 = time.perf_counter()

    genomes, contigs = draw_community(rng, preset)
    contigs = lay_out(rng, genomes, contigs)
    genomes, contigs = annotate(rng, genomes, contigs)
    pairs = draw_pairs(rng, genomes, contigs, preset["pairs"])

    asm = contigs.sort_values("asm")
    sizes = asm.set_index("name")["length"]
    sizes.to_csv(out / "contigs.sizes", sep="\t", header=False)
    asm[["name", "genome_name", "gc", "coverage"]].rename(
        columns={"name": "contig_name"}
    ).to_csv(out / "contigs.tsv", sep="\t", index=False)
    binned = genomes[genomes["kind"] != "unbinned"]
    binned[
        ["name", "completeness", "contamination", "quality", "gc", "abundance"]
    ].rename(columns={"name": "genome_name"}).round({"gc": 3, "abundance": 4}).to_csv(
        out / "genomes.tsv", sep="\t", index=False
    )
    md5 = write_pairs(out / "contacts.pairs.gz", args.preset, sizes, pairs)

    summary = summarize(args.preset, preset, genomes, contigs, pairs)
    summary["pairs_md5_uncompressed"] = md5
    summary["files_bytes"] = {f: (out / f).stat().st_size for f in FILES}
    summary["seconds"] = round(time.perf_counter() - t0, 1)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"{args.preset}: wrote {out} in {summary['seconds']} s")


if __name__ == "__main__":
    main()
