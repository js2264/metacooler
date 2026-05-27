#! /usr/bin/env python

"""
These synthetic fixtures are generated programmatically and checked into the repo for reproducible testing.  They can be regenerated with the Python script below.

Files

| File | Description |
|---|---|
| `assembly_1.fa` | 44 contigs (~0.66 Mbp), 10 MAGs + 5 unbinned, lengths 2–50 kb |
| `contigs_table_1.tsv` | contigs-level metadata table (with headers)
| `mag_table_1.tsv` | MAG-level metadata table (with headers)
| `contacts_1_rep1.pairs.gz` | 50,000 synthetic contacts over assembly_1 (80% intra-MAG, 20% inter-MAG), witg proper header |
| `contacts_1_rep2.pairs.gz` | 40,000 synthetic contacts over assembly_1 (70% intra-MAG, 30% inter-MAG), witg proper header |
| `assembly_2.fa` | 30 contigs (~0.45 Mbp), 8 MAGs + 2 unbinned, lengths 2–50 kb |
| `contigs_table_2.tsv` | contigs-level metadata table (with headers) |
| `mag_table_2.tsv` | MAG-level metadata table (with headers) |
| `contacts_2_rep1.pairs.gz` | 30,000 synthetic contacts over assembly_2 (80% intra-MAG, 20% inter-MAG), witg proper header |
| `contacts_2_rep2.pairs.gz` | 70,000 synthetic contacts over assembly_2 (70% intra-MAG, 30% inter-MAG), witg proper header |
"""

import gzip

import cooler
import numpy as np
import pyfaidx

np.random.seed(seed=42)

# For first assembly
n_mags = 10
contigs_per_mag: list[int] = [3, 4, 5, 3, 4, 5, 3, 4, 5, 3]
n_unbinned = 5
min_len = 2000
max_len = 50000
binsize = 1000
n_pairs = 50000

# Generate contigs
# Collect (length, mag) without baking the MAG index into the name,
# then shuffle before assigning opaque sequential names so alphabetical
# sort of contig names reveals no MAG grouping.
mag_names: list[str] = [f"MAG_{i:02d}" for i in range(1, n_mags + 1)]
_raw: list = []

for mag_idx, n_contigs in enumerate(iterable=contigs_per_mag):
    for _ in range(n_contigs):
        length = int(np.random.lognormal(np.log(15000), 0.8))
        length = max(min_len, min(max_len, length))
        length = (length // binsize) * binsize
        if length < min_len:
            length = min_len
        _raw.append((length, mag_names[mag_idx]))

for _ in range(n_unbinned):
    length = int(np.random.lognormal(np.log(8000), 0.6))
    length = max(min_len, min(max_len, length))
    length = (length // binsize) * binsize
    if length < min_len:
        length = min_len
    _raw.append((length, "__unbinned__"))

# Shuffle before naming so sequential names don't encode MAG membership
np.random.shuffle(_raw)
contigs: list = [
    (f"contig_1_{i:03d}", length, mag) for i, (length, mag) in enumerate(_raw)
]

# Write FASTA
with open("assembly_1.fa", "w") as f:
    for name, length, _ in contigs:
        f.write(f">{name}\n")
        seq = "".join(np.random.choice(list("ACGT"), size=length))
        for i in range(0, length, 80):
            f.write(seq[i : i + 80] + "\n")

# Write contigs table
# Unbinned contigs use a mix of "__unbinned__", "" (empty), and "NA" labels
# to exercise all sentinel-value normalisation paths.
_UNBINNED_LABEL_VARIANTS = ["__unbinned__", "", "NA"]
_unbinned_counter_1 = 0
with open("contigs_table_1.tsv", "w") as f:
    # columns should be with header: contig_name, genome_name, gc_content, coverage
    f.write("contig_name\tgenome_name\tgc_content\tcoverage\n")
    for name, _, mag in contigs:
        if mag == "__unbinned__":
            label = _UNBINNED_LABEL_VARIANTS[
                _unbinned_counter_1 % len(_UNBINNED_LABEL_VARIANTS)
            ]
            _unbinned_counter_1 += 1
        else:
            label = mag
        f.write(
            f"{name}\t{label}\t{np.random.uniform(0.4, 0.6):.2f}\t{np.random.uniform(1, 100):.2f}\n"
        )

# Write MAG table
# Shuffle MAG order so they are not listed sequentially
shuffled_mag_names_1 = list(mag_names)
np.random.shuffle(shuffled_mag_names_1)
with open("mag_table_1.tsv", "w") as f:
    # columns should be with header: genome_name, phylum, completeness, contamination, quality
    # Random values for phylum, completeness, contamination, quality
    f.write("genome_name\tphylum\tcompleteness\tcontamination\tquality\n")
    for name in shuffled_mag_names_1:
        phylum = np.random.choice(["Proteobacteria", "Bacteroidetes", "Firmicutes"])
        completeness = np.random.uniform(0.5, 1.0)
        contamination: int | float = np.random.uniform(low=0, high=0.5)
        quality = np.random.choice(["high", "medium", "low"])
        f.write(
            f"{name}\t{phylum}\t{completeness:.2f}\t{contamination:.2f}\t{quality}\n"
        )

# Write pairs
contig_names: list = [c[0] for c in contigs]
contig_lengths: dict = {c[0]: c[1] for c in contigs}

with gzip.open(filename="contacts_1_rep1.pairs.gz", mode="wt") as f:
    f.write("#pairs format v1.0\n")
    f.write("#sorted: chr1-chr2-pos1-pos2\n")
    f.write("#shape: upper triangle\n")
    f.write("#genome_assembly: synthetic_assembly_1\n")
    for name, length, _ in contigs:
        f.write(f"#chromsize: {name} {length}\n")
    f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    for i in range(n_pairs):
        if np.random.rand() < 0.8:
            mag = np.random.choice(mag_names + ["__unbinned__"])
            mag_contigs = [c for c in contigs if c[2] == mag]
            if len(mag_contigs) < 2:
                mag_contigs = contigs
            c1, c2 = np.random.choice([c[0] for c in mag_contigs], size=2, replace=True)
        else:
            c1, c2 = np.random.choice(contig_names, size=2, replace=False)
        p1 = np.random.randint(1, contig_lengths[c1])
        p2 = np.random.randint(1, contig_lengths[c2])
        f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

with gzip.open("contacts_1_rep2.pairs.gz", "wt") as f:
    f.write("#pairs format v1.0\n")
    f.write("#sorted: chr1-chr2-pos1-pos2\n")
    f.write("#shape: upper triangle\n")
    f.write("#genome_assembly: synthetic_assembly_1\n")
    for name, length, _ in contigs:
        f.write(f"#chromsize: {name} {length}\n")
    f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    for i in range(n_pairs):
        if np.random.rand() < 0.7:
            mag = np.random.choice(mag_names + ["__unbinned__"])
            mag_contigs = [c for c in contigs if c[2] == mag]
            if len(mag_contigs) < 2:
                mag_contigs = contigs
            c1, c2 = np.random.choice([c[0] for c in mag_contigs], size=2, replace=True)
        else:
            c1, c2 = np.random.choice(contig_names, size=2, replace=False)
        p1 = np.random.randint(1, contig_lengths[c1])
        p2 = np.random.randint(1, contig_lengths[c2])
        f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

## Now for second assembly
n_mags = 8
contigs_per_mag = [3, 4, 5, 3, 4, 5, 1, 3]
n_unbinned = 2
min_len = 2000
max_len = 50000
binsize = 1000
n_pairs = 50000

# Generate contigs
mag_names = [f"MAG_{i:02d}" for i in range(1, n_mags + 1)]
_raw = []

for mag_idx, n_contigs in enumerate(contigs_per_mag):
    for _ in range(n_contigs):
        length = int(np.random.lognormal(np.log(15000), 0.8))
        length = max(min_len, min(max_len, length))
        length = (length // binsize) * binsize
        if length < min_len:
            length = min_len
        _raw.append((length, mag_names[mag_idx]))

for _ in range(n_unbinned):
    length = int(np.random.lognormal(np.log(8000), 0.6))
    length = max(min_len, min(max_len, length))
    length = (length // binsize) * binsize
    if length < min_len:
        length = min_len
    _raw.append((length, "__unbinned__"))

# Shuffle before naming so sequential names don't encode MAG membership
np.random.shuffle(_raw)
contigs = [(f"contig_2_{i:03d}", length, mag) for i, (length, mag) in enumerate(_raw)]

# Write FASTA
with open("assembly_2.fa", "w") as f:
    for name, length, _ in contigs:
        f.write(f">{name}\n")
        seq = "".join(np.random.choice(list("ACGT"), size=length))
        for i in range(0, length, 80):
            f.write(seq[i : i + 80] + "\n")

# Write contigs table
_unbinned_counter_2 = 0
with open("contigs_table_2.tsv", "w") as f:
    # columns should be with header: contig_name, genome_name, gc_content, coverage
    f.write("contig_name\tgenome_name\tgc_content\tcoverage\n")
    for name, _, mag in contigs:
        if mag == "__unbinned__":
            label = _UNBINNED_LABEL_VARIANTS[
                _unbinned_counter_2 % len(_UNBINNED_LABEL_VARIANTS)
            ]
            _unbinned_counter_2 += 1
        else:
            label = mag
        f.write(
            f"{name}\t{label}\t{np.random.uniform(0.4, 0.6):.2f}\t{np.random.uniform(1, 100):.2f}\n"
        )

# Write MAG table
# Shuffle MAG order so they are not listed sequentially
shuffled_mag_names_2 = list(mag_names)
np.random.shuffle(shuffled_mag_names_2)
with open("mag_table_2.tsv", "w") as f:
    f.write("genome_name\tphylum\tcompleteness\tcontamination\tquality\n")
    for name in shuffled_mag_names_2:
        phylum = np.random.choice(["Proteobacteria", "Bacteroidetes", "Firmicutes"])
        completeness = np.random.uniform(0.5, 1.0)
        contamination = np.random.uniform(0, 0.5)
        quality = np.random.choice(["high", "medium", "low"])
        f.write(
            f"{name}\t{phylum}\t{completeness:.2f}\t{contamination:.2f}\t{quality}\n"
        )

# Write pairs
contig_names = [c[0] for c in contigs]
contig_lengths = {c[0]: c[1] for c in contigs}

with gzip.open("contacts_2_rep1.pairs.gz", "wt") as f:
    f.write("#pairs format v1.0\n")
    f.write("#sorted: chr1-chr2-pos1-pos2\n")
    f.write("#shape: upper triangle\n")
    f.write("#genome_assembly: synthetic_assembly_1\n")
    for name, length, _ in contigs:
        f.write(f"#chromsize: {name} {length}\n")
    f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    for i in range(n_pairs):
        if np.random.rand() < 0.8:
            mag = np.random.choice(mag_names + ["__unbinned__"])
            mag_contigs = [c for c in contigs if c[2] == mag]
            if len(mag_contigs) < 2:
                mag_contigs = contigs
            c1, c2 = np.random.choice([c[0] for c in mag_contigs], size=2, replace=True)
        else:
            c1, c2 = np.random.choice(contig_names, size=2, replace=False)
        p1 = np.random.randint(1, contig_lengths[c1])
        p2 = np.random.randint(1, contig_lengths[c2])
        f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

with gzip.open("contacts_2_rep2.pairs.gz", "wt") as f:
    f.write("#pairs format v1.0\n")
    f.write("#sorted: chr1-chr2-pos1-pos2\n")
    f.write("#shape: upper triangle\n")
    f.write("#genome_assembly: synthetic_assembly_1\n")
    for name, length, _ in contigs:
        f.write(f"#chromsize: {name} {length}\n")
    f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    for i in range(n_pairs):
        if np.random.rand() < 0.7:
            mag = np.random.choice(mag_names + ["__unbinned__"])
            mag_contigs = [c for c in contigs if c[2] == mag]
            if len(mag_contigs) < 2:
                mag_contigs = contigs
            c1, c2 = np.random.choice([c[0] for c in mag_contigs], size=2, replace=True)
        else:
            c1, c2 = np.random.choice(contig_names, size=2, replace=False)
        p1 = np.random.randint(1, contig_lengths[c1])
        p2 = np.random.randint(1, contig_lengths[c2])
        f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

fa_records = pyfaidx.Fasta("assembly_1.fa")
mat = cooler.util.digest(fa_records, enzyme="DpnII")  # type: ignore
mat.to_csv("assembly_1_DpnII-digested.bed", sep="\t", header=False, index=False)
with open("assembly_1.chrom.sizes", "w") as f:
    for name, seq in fa_records.items():
        f.write(f"{name}\t{len(seq)}\n")
fa_records = pyfaidx.Fasta("assembly_2.fa")
mat = cooler.util.digest(fa_records, enzyme="DpnII")  # type: ignore
mat.to_csv("assembly_2_DpnII-digested.bed", sep="\t", header=False, index=False)
with open("assembly_2.chrom.sizes", "w") as f:
    for name, seq in fa_records.items():
        f.write(f"{name}\t{len(seq)}\n")
