# metacooler

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue)](https://www.python.org/downloads/)
[![PyPI version](https://img.shields.io/pypi/v/metacooler.svg)](https://pypi.org/project/metacooler/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/js2264/metacooler/branch/devel/graph/badge.svg?token=XMYgQxb7KH)](https://codecov.io/gh/js2264/metacooler)

**Hi-C contact maps for metagenomes.**

metacooler extends the [cooler](https://github.com/open2c/cooler) HDF5 format to
handle multi-genome (metagenomic) Hi-C data. It adds a `genomes/` hierarchy on
top of cooler's `chroms/bins/pixels` model while preserving **full backward
compatibility**: any `.metacool` file opens with `cooler.Cooler()` out of the
box.

```bash
pip install metacooler
```

<div class="grid cards" markdown>

- :material-rocket-launch: **[Quickstart](quickstart.md)** — from raw pairs to a genome heatmap in five commands
- :material-lightbulb-on: **[Concepts](concepts.md)** — the mental model: genomes, contigs, bins, pixels
- :material-console: **[CLI reference](cli.md)** — every command and option
- :material-language-python: **[Python API](api/metacooler.md)** — the `MetaCooler` class and functions

</div>

---

## Why metacooler?

Metagenomic Hi-C binning tools such as
[MetaTOR](https://github.com/koszullab/metator) and
[bin3C](https://github.com/cerebis/bin3C) produce per-contig contact maps across
dozens or hundreds of MAGs simultaneously. Standard cooler files have no concept
of genome identity — contigs from different MAGs appear as independent
"chromosomes" with no shared metadata.

metacooler solves this by:

- Storing a **`/genomes` table** linking contigs to MAGs, with optional quality
  metadata (completeness, contamination, taxonomy…).
- Computing **genome-level and contig-level aggregated contact matrices** without
  re-creating the file.
- Recording **contact-quality statistics** — the 3D ratio, and how much signal is
  intra-contig, intra-MAG and inter-MAG — at creation time, per file *and* per
  MAG, plus **P(s)** curves for judging the library itself.
- Providing a **`/contigs` table** as the canonical chromosome table, with a soft
  link at `/chroms` so cooler compatibility is transparent.
- Offering a **collection format** (`.mmetacool`) to group multiple assemblies,
  samples, or resolutions in a single file.

See [Concepts](concepts.md) for the full data model.

## Quick look

```bash
# Create from MetaTOR output
metacooler create pairs \
    --pairs   contacts.pairs.gz \
    --bins    fragments.bed \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_table.tsv \
    -o        gut_T0.metacool

# Inspect
metacooler info gut_T0.metacool

# Contact quality: the 3D ratio and its breakdown
metacooler qc gut_T0.metacool --table

# P(s), one curve per MAG — create pairs already wrote gut_T0.ps.tsv at
# base-pair resolution; this re-reads it from the binned map, or from the
# original pairs with --pairs
metacooler ps gut_T0.metacool --by genome --plot ps.png

# Genome-level heatmap
metacooler show gut_T0.metacool --level genome -o heatmap.png

# Export one MAG as a standard cool file
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool
```

## Compatibility

| Tool | Compatibility |
|------|--------------|
| `cooler.Cooler()` | ✅ Full — reads `.metacool` as a standard cooler |
| cooltools | ✅ Works on files produced by `metacooler extract` or `aggregate` |
| HiGlass | ✅ Via `cooler zoomify` on a `.metacool` (or an exported `.cool`), which writes the multi-resolution `.mcool` HiGlass ingests — see [Multi-resolution files](guides/multireso.md#higlass) |
| pairtools | ✅ Pairs files are the primary input format |
| MetaTOR | ✅ Native support via `--contigs-format metator` |

## File layout at a glance

```
.metacool file (HDF5)
├── /genomes     ← genome/MAG table  [name, length, n_contigs, n_bins, …metadata]
├── /contigs     ← contig table      [name, length, genome_id, …metadata]
├── /bins        ← genomic bins      [chrom (=contig_id), start, end]
├── /pixels      ← contacts          [bin1_id, bin2_id, count]
├── /indexes     ← fast lookup       [chrom_offset, bin1_offset, genome_offset]
└── /chroms      → soft link to /contigs   (cooler compatibility)
```

The full byte-level specification lives in the
[format reference](format.md).

## Documentation map

| If you want to… | Go to |
|------------------|-------|
| Install the package | [Installation](installation.md) |
| Follow an end-to-end example | [Quickstart](quickstart.md) |
| Understand the data model | [Concepts](concepts.md) |
| Build a file from pairs or a cooler | [Creating from pairs](guides/from-pairs.md), [from a cooler](guides/from-cool.md) |
| Combine one cooler per genome into one file | [Collecting single-genome cools](guides/collect.md) |
| Get data back out | [Exporting data](guides/export.md) |
| Work at multiple resolutions | [Multi-resolution files](guides/multireso.md) |
| Look up a command | [CLI reference](cli.md) |
| Call metacooler from Python | [Python API](api/metacooler.md) |
| Read the on-disk format | [Format specification](format.md) |
