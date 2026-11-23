# metacooler: Hi-C contact maps for metagenomes


[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue)](https://www.python.org/downloads/) [![PyPI version](https://img.shields.io/pypi/v/metacooler.svg)](https://pypi.org/project/metacooler/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![codecov](https://codecov.io/gh/js2264/metacooler/branch/devel/graph/badge.svg?token=XMYgQxb7KH)](https://codecov.io/gh/js2264/metacooler)

An extension of the [cooler](https://github.com/open2c/cooler) format for multi-genome (metagenomic) Hi-C data. `metacooler` adds a `genomes/` hierarchy on top of cooler’s `chroms/bins/pixels` model, while preserving **full backward compatibility**: any `.metacool` file can be opened with `cooler.Cooler()` out of the box.

`metacooler` builds on `cooler` to provide a framework for storing and analysing metagenomic contact maps: every `.metacool` file is also a valid cooler, so cooler, cooltools and the rest of the ecosystem keep working on it.

## Installation

``` console
pip install metacooler
```

For development, clone the repository and run `uv sync`, which installs the runtime dependencies plus the `dev` dependency group (pytest, ruff, ty, pre-commit).

## Quick start

``` console
metacooler create pairs --help
```

``` console
metacooler create pairs \
    --bins fragments.bed \
    --contigs contigs2mag.tsv \
    --pairs contacts.pairs \
    --output output.metacool \
    --genomes mag_table.tsv       ## <-- Optional

metacooler info output.metacool
metacooler dump -t genomes output.metacool
metacooler aggregate --level genome output.metacool -o genome.cool
metacooler show output.metacool --level genome -o heatmap.png
metacooler qc output.metacool --table
metacooler ps output.metacool --pairs contacts.pairs --by genome --plot ps.png
```

Already have one `.cool` file per genome? `collect` zips them into a single metacool file, one genome per input:

``` console
metacooler collect E_coli.cool B_subtilis.cool S_aureus.cool -o species.metacool
```

Every command takes its input as a positional path (or an `file.mmetacool::entries/<name>` URI), and shares `-o/--output`, `-f/--force`, `-v/--verbose`, and `-q/--quiet`.

                                                                                    
     Usage: metacooler create pairs [OPTIONS]                                       
                                                                                    
     Create a metacool file from pairs and contig assignments.                      
                                                                                    
    ╭─ Options ────────────────────────────────────────────────────────────────────╮
    │ *  --pairs                   Pairs file: plain text, or compressed with      │
    │                              gzip, bgzip, bzip2, zstd or lz4, read from its  │
    │                              content. Repeat the option to read several      │
    │                              files — one per Hi-C library of the same        │
    │                              assembly — in one pass, as if concatenated:     │
    │                              their contacts are summed, and the QC block,    │
    │                              the P(s) table and the dropped-pairs report     │
    │                              cover them all. Not comma-split: a path may     │
    │                              contain a comma. [required]                     │
    │ *  --bins                    BED file of bins (without column names), or an  │
    │                              integer bin size in bp (e.g. 1000). If an       │
    │                              integer, contig lengths must be available from  │
    │                              the pairs header or a --fasta should be         │
    │                              provided. [required]                            │
    │ *  --output              -o  Output .metacool file path. [required]          │
    │    --contigs                 Contig-to-genome assignment table; which        │
    │                              columns hold the names is set by                │
    │                              --contigs-format. Leave it out to build a map   │
    │                              with no binning yet: every contig in            │
    │                              __unbinned__, to be binned later with reassign. │
    │    --contigs-format          Which columns of the --contigs table hold the   │
    │                              contig and the genome names: CONTIG,GENOME,     │
    │                              each a 1-based position or a header name — e.g. │
    │                              2,5 or Name,Final_cluster — plus ,noheader when │
    │                              the table has no header line. NA as the genome  │
    │                              column (1,NA,noheader, contig_name,NA) says the │
    │                              table assigns no genome: create bins nothing    │
    │                              yet (every contig in __unbinned__), reassign    │
    │                              unbins the contigs it lists. [default: contigs] │
    │                              Aliases for the tables of common binners:       │
    │                              contigs (contig_name,genome_name),              │
    │                              contig_genome (1,2,noheader), metator           │
    │                              (Name,Final_cluster; Final_bin if the table has │
    │                              no Final_cluster), dastool (1,2,noheader),      │
    │                              binette (1,2,noheader), semibin (contig,bin),   │
    │                              vamb (contigname,clustername).                  │
    │    --genomes                 Optional genomes metadata table. Its            │
    │                              genome-name column (see --genomes-format) is    │
    │                              matched against the genomes; every other column │
    │                              is stored in the /genomes table.                │
    │    --genomes-format          Which column of the --genomes table holds the   │
    │                              genome names: a 1-based position or a header    │
    │                              name — e.g. Name (CheckM2) or user_genome       │
    │                              (GTDB-Tk) — plus ,noheader when the table has   │
    │                              no header line. [default: genome_name]          │
    │    --filter                  Keep only rows matching <table>.<column> <op>   │
    │                              <value>, where <table> is genomes or contigs —  │
    │                              e.g. --filter 'genomes.length > 1Mb'.           │
    │                              Operators: > >= < <= == != =~ !~ (=~ is a regex │
    │                              match). Numbers accept bp/kb/Mb/Gb suffixes.    │
    │                              Repeat the option to require several conditions │
    │                              at once.                                        │
    │    --metadata                Path to a JSON or YAML file containing metadata │
    │                              to ingest into the file.                        │
    │    --qc/--no-qc              Tally contact-quality statistics during         │
    │                              creation and store them in metadata.qc and the  │
    │                              qc_* columns of /genomes and /contigs. Read     │
    │                              them back with metacooler qc. [default: qc]     │
    │    --qc-threshold            Distance splitting short- from long-range       │
    │                              intra-contig contacts, measured between read    │
    │                              positions: long-range from this distance on.    │
    │                              The default counts as long-range the pairs more │
    │                              than 1 kb apart, as MetaTOR's mge does. Accepts │
    │                              bp/kb/Mb/Gb suffixes. [default: 1001]           │
    │    --force               -f  Overwrite output file if it already exists.     │
    │    --verbose             -v  Increase verbosity (debug logging).             │
    │    --quiet               -q  Suppress informational logging (warnings and    │
    │                              errors still shown).                            │
    │    --contig-order            Strategy for ordering contigs in the output     │
    │                              metacool file.                                  │
    │                              - genome-then-length: group by genome,          │
    │                              longest-first within each genome.               │
    │                              - genome-then-alpha: group by genome,           │
    │                              alphabetical within each genome.                │
    │                              - fasta: preserve order from the input FASTA    │
    │                              (or pairs header).                              │
    │                              - alpha: sort all contigs alphabetically.       │
    │                              [default: genome-then-length]                   │
    │    --fasta                   Assembly FASTA. Used to extract contig lengths  │
    │                              when --bins is an integer and contig sizes are  │
    │                              absent from the pairs header.                   │
    │    --assembly-name           Name for this assembly / metagenome.            │
    │    --drop-unassigned         Drop contigs not assigned to any genome.        │
    │    --chunksize               Number of pairs to process per chunk. This is   │
    │                              chunk granularity, not read memory (the pyarrow │
    │                              read block is capped separately). Under         │
    │                              'external' sorting each chunk is written as a   │
    │                              temporary cooler at a fixed per-call cost that  │
    │                              scales with contig count, so smaller chunks     │
    │                              cost disproportionately more. [default:         │
    │                              50000000]                                       │
    │    --nproc                   Threads binning pairs while the file is read    │
    │                              and decompressed (which stays on one thread).   │
    │                              The output is identical whatever the value;     │
    │                              each thread holds one more batch of pairs in    │
    │                              memory. [default: 1]                            │
    │    --ps/--no-ps              Accumulate a per-MAG P(s) curve from the same   │
    │                              pass over the pairs and write it beside the     │
    │                              output as <output>.ps.tsv. Measured on read     │
    │                              positions, so it resolves distances the binned  │
    │                              map cannot. No plot is drawn — pass the table   │
    │                              to metacooler ps or your own plotting.          │
    │                              [default: ps]                                   │
    │    --ps-bins-per-decade      Number of log-spaced distance bins per decade   │
    │                              in the P(s) table. [default: 10]                │
    │    --sort-strategy           Pixel sorting strategy. 'auto': use 'memory'    │
    │                              unless the binned pixels exceed 100 M rows.     │
    │                              'memory': collect all pixels in RAM and sort    │
    │                              before writing (fast for small files).          │
    │                              'external': stream unsorted pixels to cooler's  │
    │                              external sort (scalable for large files).       │
    │                              [default: auto]                                 │
    │    --help                -h  Show this message and exit.                     │
    ╰──────────────────────────────────────────────────────────────────────────────╯

## Input files

**Bins file**: BED-like file with `chrom`, `start`, `end` columns defining the genomic bins.

- Can be generated from a FASTA assembly with `cooler makebins` to get fixed-size bins
- Can be generated from a FASTA assembly with `cooler digest` to get restriction fragments.
- Otherwise an arbitrary resolution can be used by specifying `--bins 1000` for 1kb bins, etc.

``` txt
contig_1_000    0   721
contig_1_000    721 747
contig_1_000    747 767
contig_1_000    767 1013
contig_1_000    1013    1099
...
contig_1_043    8869    8959
contig_1_043    8959    9000
```

**Contig-genome map**: TSV file mapping contigs to genomes (MAGs). This is used to populate the `contigs` table and map contigs to genomes. Must be either:

- A two-column file with NO HEADER: first column is treated as contig name, second column as genome/MAG name.
- A TSV file with header, whose first column is `contig_name` and second column `genome_name`. Additional metadata columns will be stored in the `contigs/` table.

``` txt
contig_name genome_name gc_content  coverage
contig_1_000    MAG_04  0.51    99.40
contig_1_001    MAG_02  0.59    75.43
contig_1_002    MAG_04  0.45    8.81
contig_1_003    __unbinned__    0.52    8.28
...
contig_1_042    MAG_05  0.52    18.45
contig_1_043    MAG_05  0.42    97.78
```

**Pairs file**: [4DN-compliant](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) tab-delimited file with `readID`, `chrom1`, `pos1`, `chrom2`, `pos2`, `strand1`, `strand2` columns defining the Hi-C contacts. Can be generated from raw FASTQ reads with [`metator pipeline`](https://github.com/koszullab/metator).

``` txt
#pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#genome_assembly: synthetic_assembly_1
#chromsize: contig_1_000 18000
#chromsize: contig_1_001 12000
#chromsize: contig_1_002 3000
#chromsize: contig_1_003 6000
...
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
read_0  contig_1_007    41718   contig_1_006    6634    +   +
read_1  contig_1_037    2208    contig_1_000    12785   +   +
read_2  contig_1_003    5884    contig_1_028    2231    +   +
read_3  contig_1_014    3164    contig_1_021    47093   +   +
read_4  contig_1_008    138 contig_1_038    6042    +   +
...
```

*\[Optional\]* **Genomes table**: TSV file with at least a column named `genome_name`, and optional metadata columns describing the genomes/MAGs.

``` txt
genome_name phylum  completeness    contamination   quality
MAG_02  Firmicutes  0.92    0.14    low
MAG_03  Proteobacteria  0.73    0.03    low
MAG_01  Bacteroidetes   0.63    0.10    high
MAG_04  Firmicutes  0.99    0.24    medium
...
MAG_09  Proteobacteria  0.95    0.04    low
MAG_05  Proteobacteria  0.89    0.27    high
```

*\[Optional\]* **Metadata file**: a YAML/JSON file with additional metadata.

``` yaml
sample:
  name: gut_microbiome_T0
  description: This is a gut microbiome sample at time point 0
  organism: human gut metagenome
  enzyme: HpaII
  enzyme_site: "C^CGG"
sequencing:
  platform: Illumina NovaSeq 6000
  read_length: 150
  n_read_pairs: 245000000
  n_valid_pairs: 180000000
assembly:
  name: gut_microbiome_T0_assembly.fa
  description: metagenome assembly for the gut microbiome sample at time point 0
  assembler: MEGAHIT v1.2.9
  assembly_date: 2024-06-01
  assembly_stats:
    total_length: 3500000
    n_contigs: 120
    n50: 50000
```

## Design

- **Metacool identity**: `metacooler-format="HDF5::MetaCooler"` attribute + `genomes/` group
- **Implement genomes table**: `name`, `length`, `n_contigs`, `n_bins`, plus optional metadata columns
- **Support multi-contigs and multi-genomes visualization**: `show` command at `--level bin/contig/genome`
- **Normalize metagenome-aware**: `normalize ice` (balancing genome by genome), `normcc` (MetaCC) and `hiczin` (HiCzin) store cooler-style weights in `/bins`, beside those any other tool stores with `metacooler.store_weights`; `show`/`dump`/`network --balance` read them, and so does cooler itself
- **Report contact quality**: the 3D ratio and its intra-contig / intra-MAG / inter-MAG breakdown, tallied at creation time into `metadata.qc` and per-MAG `/genomes` columns, recomputable with `qc`
- **Export P(s)**: contact probability vs genomic distance, pooled or per MAG, as a table and a log-log plot (`ps`). `create pairs` writes the per-MAG curve at base-pair resolution beside the output as `<output>.ps.tsv`
- **Cooler compatibility**: soft link `chroms → contigs`, `format="HDF5::Cooler"`, standard `indexes/chrom_offset`
- **Canonical contig ordering**: configurable (`genome-then-length`, `genome-then-alpha`, `fasta`, `alpha`; `cool` for `create cool` and `collect`)
- **Export sequences**: `fasta` writes the contigs of each genome, or of each cluster of a clustering round, from the assembly: each record as it is, or stitched into one sequence per genome
- **Collect single-genome maps**: `collect` stacks one `.cool` per genome into a metacool file, renaming contigs `<genome>__<contig>` so shared names cannot collide
- **Unassigned contigs**: default to `__unbinned__` genome
- **Re-attributable bins**: `reassign` swaps in a new contig-to-genome table on an existing file — the `/genomes` table and the contig layout are rebuilt, the contacts are not recomputed
- **Multi-resolution**: `binnify` (exact multiples of the source bin size only) and `zoomify` → a `.mmetacool` collection

## Releasing

Changelog entries for in-progress work go under `## Unreleased` in `docs/changelog.md` as they land on `devel`. To cut a release:

``` console
uv run bump-my-version bump minor   # or patch / major
git push --follow-tags
```

[`bump-my-version`](https://github.com/callowayproject/bump-my-version) (configured in `pyproject.toml`) promotes `## Unreleased` to a dated `## X.Y.Z - YYYY-MM-DD` section, bumps `__version__` in `src/metacooler/_version.py` — the single source of truth for the package version — and commits and tags the result, all in one step. Pushing the tag triggers [`release.yml`](.github/workflows/release.yml): it builds the sdist/wheel, fails if the tag doesn’t match `__version__`, publishes to PyPI via [trusted publishing](https://docs.pypi.org/trusted-publishers/) (no stored token), and opens a GitHub release with that version’s changelog section and the built artifacts attached. The documentation site is redeployed by [`docs.yml`](.github/workflows/docs.yml) from `devel`’s push of the release commit, not from the tag, which the `github-pages` environment refuses.
