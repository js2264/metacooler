# Why metacooler?

Metagenomic Hi-C binning tools such as
[MetaTOR](https://github.com/koszullab/metator) and
[bin3C](https://github.com/cerebis/bin3C) produce per-contig contact maps across
dozens or hundreds of MAGs simultaneously. Standard cooler files have no concept
of genome identity — contigs from different MAGs appear as independent
"chromosomes" with no shared metadata.

## What metacooler adds

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

## A quick look

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

The [Quickstart](quickstart.md) walks through the same session step by step.

## Compatibility

| Tool | Compatibility |
|------|--------------|
| `cooler.Cooler()` | ✅ Full — reads `.metacool` as a standard cooler |
| cooltools | ✅ Works on files produced by `metacooler extract` or `aggregate` |
| HiGlass | ✅ Via `cooler zoomify` on a `.metacool` (or an exported `.cool`), which writes the multi-resolution `.mcool` HiGlass ingests — see [Multi-resolution files](guides/multireso.md#higlass) |
| pairtools | ✅ Pairs files are the primary input format |
| MetaTOR | ✅ Native support via `--contigs-format metator` |

How a single file can be both a cooler and a metacooler is explained in
[Cooler compatibility](concepts.md#cooler-compatibility).

## Where next

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
