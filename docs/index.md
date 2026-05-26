# metacooler

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue)](https://www.python.org/downloads/)
[![PyPI version](https://img.shields.io/pypi/v/metacooler.svg)](https://pypi.org/project/metacooler/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/js2264/metacooler/branch/devel/graph/badge.svg?token=XMYgQxb7KH)](https://codecov.io/gh/js2264/metacooler)

**Hi-C contact maps for metagenomes.**

metacooler extends the [cooler](https://github.com/open2c/cooler) HDF5 format to handle multi-genome (metagenomic) Hi-C data. It adds a `genomes/` hierarchy on top of cooler's `chroms/bins/pixels` model while preserving **full backward compatibility**: any `.metacool` file opens with `cooler.Cooler()` out of the box.

---

## Why metacooler?

Metagenomic Hi-C binning tools such as [MetaTOR](https://github.com/koszullab/metator) and [bin3C](https://github.com/cerebis/bin3C) produce per-contig contact maps across dozens or hundreds of MAGs simultaneously. Standard cooler files have no concept of genome identity — contigs from different MAGs appear as independent "chromosomes" with no shared metadata.

metacooler solves this by:

- Storing a **`/genomes` table** linking contigs to MAGs, with optional quality metadata (completeness, contamination, taxonomy…)
- Computing **genome-level and contig-level aggregated contact matrices** without re-creating the file
- Providing a **`/contigs` table** as the canonical chromosome table, with a soft link at `/chroms` so cooler compatibility is transparent
- Offering a **collection format** (`.mmetacool`) to group multiple assemblies or samples in a single file

## Quick look

```bash
# Create from MetaTOR output
metacooler create pairs \
    --pairs  contacts.pairs.gz \
    --bins   fragments.bed \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_table.tsv \
    --output  gut_T0.metacool

# Inspect
metacooler info gut_T0.metacool

# Genome-level heatmap
metacooler show gut_T0.metacool --level genome --out heatmap.png

# Export one MAG as a standard cool file
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool
```

## Compatibility

| Tool | Compatibility |
|------|--------------|
| `cooler.Cooler()` | ✅ Full — reads `.metacool` as a standard cooler |
| cooltools | ✅ Works on files produced by `metacooler extract` or `aggregate` |
| HiGlass | ✅ Works on `.cool` files exported from metacooler |
| pairtools | ✅ Pairs files are the primary input format |
| MetaTOR | ✅ Native support via `--contigs-format metator` |

## Overview

```
.metacool file (HDF5)
├── /genomes     ← genome/MAG table  [name, length, n_contigs, …metadata]
├── /contigs     ← contig table      [name, length, genome_id, …metadata]
├── /bins        ← genomic bins      [chrom (=contig_id), start, end]
├── /pixels      ← contacts          [bin1_id, bin2_id, count]
├── /indexes     ← fast lookup       [chrom_offset, bin1_offset, genome_offset]
└── /chroms      → soft link to /contigs   (cooler compatibility)
```
