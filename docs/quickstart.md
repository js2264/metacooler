# Quickstart

This guide takes you from raw MetaTOR output to a genome-level contact heatmap in five commands.

## Input files

You need three files (a fourth is optional):

| File | Description |
|------|-------------|
| `contacts.pairs.gz` | 4DN-compliant pairs file — contacts between contig positions |
| `fragments.bed` | BED file of bins (restriction fragments or fixed-size windows) |
| `contig_data_final.txt` | MetaTOR contig-to-MAG assignment table |
| `mag_table.tsv` | *(optional)* MAG metadata (completeness, taxonomy…) |

!!! tip "Generating input files"
    Pairs files can be produced from raw FASTQ reads with
    [`metator pipeline`](https://github.com/koszullab/metator) or
    [`pairtools`](https://github.com/open2c/pairtools).
    Bin BED files can be generated with `cooler makebins` (fixed-size) or
    `cooler digest` (restriction fragments).

## Step 1 — Create a metacool file

```bash
metacooler create pairs \
    --pairs   contacts.pairs.gz \
    --bins    fragments.bed \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_table.tsv \
    --output  gut_T0.metacool
```

This streams the pairs file in chunks, bins contacts, and writes an HDF5 file
with the full metacooler schema. The `--contigs-format metator` flag tells
metacooler to look for `Name` and `Final_bin` columns in the contig table.

!!! note "No MetaTOR?"
    If your contig table is a simple two-column TSV with no header
    (`contig_name<TAB>genome_name`), use `--contigs-format contig_genome` instead.
    See the [Creating from pairs](guides/from-pairs.md) guide for all supported formats.

## Step 2 — Inspect the file

```bash
metacooler info gut_T0.metacool
```

```json
{
  "bin-size": 2000,
  "bin-type": "fixed",
  "creation-date": "2026-05-12T09:41:00+00:00",
  "format": "HDF5::MetaCooler",
  "format-url": "https://github.com/js2264/metacooler",
  "format-version": 1,
  "generated-by": "metacooler-0.1.0",
  "assembly": null,
  "nbins": 48291,
  "ncontigs": 1203,
  "ngenomes": 47,
  "nnz": 2847193,
  "sum": 18043712,
  "storage-mode": "symmetric-upper",
  "ordering": "genome-then-length"
}
```

```bash
# Rich table format
metacooler info gut_T0.metacool --table

# Inspect the HDF5 hierarchy
metacooler tree gut_T0.metacool
```

## Step 3 — Browse tables

```bash
# List all genomes / MAGs
metacooler dump gut_T0.metacool -t genomes -H

# List contigs belonging to a specific MAG
metacooler dump gut_T0.metacool -t contigs --genomes MAG_01 -H
```

## Step 4 — Visualise the contact map

```bash
# Genome-level heatmap (default level). Note: -f for input, -o for output.
metacooler show -f gut_T0.metacool --level genome -o genome_heatmap.png

# Intra-genome heatmap at bin resolution
metacooler show -f gut_T0.metacool --level bin --genome MAG_01 -o MAG_01_heatmap.png

# A specific genomic range (implies --level bin)
metacooler show -f gut_T0.metacool -r contig_0_0:0-50000 -o range.png
```

## Step 5 — Export for downstream tools

```bash
# Aggregate to genome level → standard .cool (works with cooltools, HiGlass)
metacooler aggregate gut_T0.metacool --level genome -o genome_level.cool

# Extract one MAG as a standalone cool file
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool

# Export the genome-level network edge list
metacooler network gut_T0.metacool --level genome -o genome_network.tsv
```

## Python API equivalent

The same workflow from Python:

```python
import metacooler

# Create
metacooler.from_pairs(
    uri="gut_T0.metacool",
    pairs="contacts.pairs.gz",
    bins="fragments.bed",
    contig_genome_map="contig_data_final.txt",
    contig_table_format="metator",
    genomes_table="mag_table.tsv",
)

# Open and explore
mc = metacooler.MetaCooler("gut_T0.metacool")

print(mc.info)
print(mc.genomes)
print(mc.contigs.head())

# Genome-level contact map: write a standard cooler, then read it back
mc.to_cooler("genome_level.cool", level="genome")

# Bins for a specific MAG
bins = mc.bins_for_genome("MAG_01")
```
