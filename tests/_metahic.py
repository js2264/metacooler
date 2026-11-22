"""Synthetic metagenomic Hi-C maps for the normalization tests.

:func:`make_metahic` simulates a small community -- a few genomes of several
contigs each, plus unbinned contigs -- and writes its assembly, its contig
table (with a ``cov`` coverage column and an ``RS`` site-count column) and a
pairs file in which most contacts stay inside a contig or a genome, as real
metaHi-C contacts do, and a few join random contigs (the "spurious" ones).
It then builds the map with :func:`metacooler.from_pairs`, measuring the QC
block on read positions at a 1 kb threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from metacooler.io.create import from_pairs


@dataclass(frozen=True)
class MetaHiC:
    """Paths of one simulated community and its map."""

    fasta: Path
    contigs: Path
    pairs: Path
    metacool: Path
    contig_genome: dict[str, str]


def make_metahic(
    directory: Path,
    *,
    genomes: int = 3,
    contigs_per_genome: int = 3,
    unbinned: int = 1,
    n_pairs: int = 8_000,
    min_length: int = 2_000,
    max_length: int = 8_000,
    binsize: int = 1_000,
    seed: int = 0,
) -> MetaHiC:
    """Simulate a community and build its ``.metacool`` map in *directory*.

    Args:
        directory: Where the files go.
        genomes: Number of genomes (``G0``, ``G1``, ...).
        contigs_per_genome: Contigs of each genome.
        unbinned: Contigs assigned to no genome.
        n_pairs: Read pairs simulated.
        min_length: Shortest contig, in bp.
        max_length: Longest contig, in bp (exclusive).
        binsize: Bin size of the map.
        seed: Seed of the simulation.

    Returns:
        A :class:`MetaHiC`.
    """
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    n_binned = genomes * contigs_per_genome
    n_contigs = n_binned + unbinned
    names = [f"c{i}" for i in range(n_contigs)]
    lengths = rng.integers(min_length, max_length, size=n_contigs)
    genome_of = np.concatenate(
        [np.repeat(np.arange(genomes), contigs_per_genome), np.full(unbinned, -1)]
    )
    abundance = rng.uniform(0.5, 4.0, size=genomes + 1)
    contig_abundance = abundance[genome_of]  # unbinned contigs share the last

    fasta = directory / "assembly.fa"
    with open(fasta, "w") as f:
        for name, length in zip(names, lengths, strict=True):
            seq = "".join(rng.choice(list("ACGT"), size=int(length)))
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i : i + 70] + "\n")

    contig_genome = {
        name: f"G{genome_of[i]}" for i, name in enumerate(names) if genome_of[i] >= 0
    }
    table = directory / "contigs.tsv"
    with open(table, "w") as f:
        f.write("contig_name\tgenome_name\tcov\tRS\n")
        for i, name in enumerate(names):
            genome = contig_genome.get(name, "NA")
            cov = contig_abundance[i] * rng.uniform(0.8, 1.25) * 10
            f.write(f"{name}\t{genome}\t{cov:.4f}\t{int(lengths[i] // 250) + 1}\n")

    p_first = lengths * contig_abundance
    first = rng.choice(n_contigs, size=n_pairs, p=p_first / p_first.sum())
    kind = rng.random(n_pairs)
    second = first.copy()
    same_genome = (kind >= 0.55) & (kind < 0.93)
    for i in np.flatnonzero(same_genome):
        g = genome_of[first[i]]
        pool = np.flatnonzero(genome_of == g) if g >= 0 else np.array([first[i]])
        second[i] = rng.choice(pool)
    noise = kind >= 0.93
    second[noise] = rng.choice(
        n_contigs, size=int(noise.sum()), p=p_first / p_first.sum()
    )
    # Positions in [1, length - 1]: well inside the contig, never dropped.
    pos1 = (rng.random(n_pairs) * (lengths[first] - 1)).astype(np.int64) + 1
    pos2 = (rng.random(n_pairs) * (lengths[second] - 1)).astype(np.int64) + 1
    intra = second == first
    # Contacts inside a contig decay with distance.
    offset = np.minimum(rng.geometric(1 / 1500, size=n_pairs), lengths[first] - 2)
    sign = np.where(rng.random(n_pairs) < 0.5, -1, 1)
    near = np.clip(pos1 + sign * offset, 1, lengths[first] - 1)
    pos2 = np.where(intra, near, pos2)

    pairs = directory / "contacts.pairs"
    with open(pairs, "w") as f:
        f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
        for i in range(n_pairs):
            f.write(
                f"r{i}\t{names[first[i]]}\t{pos1[i]}\t{names[second[i]]}\t{pos2[i]}\t+\t+\n"
            )

    metacool = directory / "map.metacool"
    from_pairs(
        uri=metacool,
        bins=None,
        binsize=binsize,
        pairs=pairs,
        contig_genome_map=table,
        contig_table_format="contigs",
        fasta=fasta,
        chunksize=max(n_pairs // 3, 1),
        qc_distance_threshold=1_000,
        ps=False,
    )
    return MetaHiC(fasta, table, pairs, metacool, contig_genome)
