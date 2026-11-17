# Benchmarks: metacooler vs cooler

`benchmarks.qmd` times metacooler against cooler **from the command line**, on
four reproducible **simulated** metagenomes, from 10³ to 2 × 10⁶ contigs. No real
reads, assembly or binning is used (there are no accessions): `simulate.py`
draws each one from a preset and a fixed seed.

| | genomes (MAGs) | contigs | assembly | read pairs |
|---|---|---|---|---|
| S | 10 (10) | 1,000 | 27 Mb | 1 M |
| M | 100 (50) | 20,000 | 157 Mb | 10 M |
| L | 2,000 (200) | 200,000 | 764 Mb | 30 M |
| XL | 20,000 (1,000) | 2,000,000 | 5.1 Gb | 100 M |

It covers building a map from pairs (`cooler cload pairs` vs `metacooler create
pairs`), its size on disk (against the `.cool` with its binning tables, and one
`.cool` per MAG), converting between the formats, multi-resolution pyramids
(`zoomify`), querying a contig, a genome and a genome pair, and genome-level
aggregation. Each command runs as a fresh process under GNU `time` and
`timeout`, 3 to 5 times (once if a run takes over 10 min), and records wall
clock, CPU and peak memory. A run stopped at the time limit is reported where
it was stopped, as a lower bound. `cooler cload pairs` is not run on XL, where
it has never finished within the limit. The results are one row per run in
`_results/timings.tsv`, and every table and figure is drawn from that file.

Every timing is against **cooler 0.12.1**, the version `uv.lock` pins.
`bench.sh` names it (`COOLER_VERSION`) and the notebook refuses to time any
other, so a lock that moves cooler on has to move that line too.

| file | what |
|---|---|
| `benchmarks.qmd` | the notebook: bash cells that run and time the commands; Python for tables and figures |
| `bench.sh` | settings, and `bench`, the timing function every bash cell sources |
| `simulate.py` | the dataset generator: presets, seeds, model |
| `targets.py` | what the queries ask for: one contig, one MAG and one MAG pair per dataset, from a fixed seed |

## Render

```bash
uv sync && source .venv/bin/activate   # metacooler, cooler and jupyter on PATH
cd benchmarks
quarto render benchmarks.qmd                         # the whole ladder: ~2.5 h on 2 cores, 8 GB RAM
BENCH_DATASETS="S M" quarto render benchmarks.qmd    # a quick look: minutes
```

Needs GNU `time` (`apt install time`); `pigz` speeds up dataset generation.

| variable | default | |
|---|---|---|
| `BENCH_DATASETS` | `S M L XL` | which datasets, smallest first |
| `BENCH_TIMEOUT` | `1800` | seconds a run may take before it is stopped (and reported as a lower bound) |
| `BENCH_MIN_REPS`, `BENCH_MAX_REPS` | `3`, `5` | runs per command, at least and at most |
| `BENCH_BUDGET` | `120` | seconds: once a command's runs add up to this, no more than `BENCH_MIN_REPS` |
| `BENCH_LONG` | `600` | seconds: a command whose run takes longer runs once |
| `BENCH_REPS` | | every command exactly this many times instead (`1` for a quick look) |
| `BENCH_SKIP_COOLER_CREATE` | `XL` | datasets on which `cooler cload pairs` is not run; empty runs it on all |
| `BENCH_BINSIZE` | `1000` | bin size of every map |
| `BENCH_WORK`, `BENCH_RESULTS` | `_work/`, `_results/` | where data and maps, and tables and figures, go |
| `BENCH_REDO` | `0` | `1` times again commands already in `timings.tsv` |

A render re-uses what is already in `timings.tsv`, so the page can be edited
and re-rendered without re-running the benchmarks.

## In CI

[`.github/workflows/benchmarks.yml`](../.github/workflows/benchmarks.yml)
renders the notebook for every `v*` tag and on demand (Actions → Benchmarks →
Run workflow, to pick the datasets and the time limit). The
page, the tables, the figures and the command logs are kept as a workflow
artifact. On a release, the page and a zip of `_results/` are attached to it.

A run on a `v*` tag, or on `devel` with *publish-page* ticked (the default),
also publishes the page on the documentation site, at
<https://js2264.github.io/metacooler/benchmarks/>, with `_results/` beside it.
Each Pages deployment replaces the whole site, so the run does not deploy the
page itself. It replaces the two assets of the `benchmark-page` release and
starts [`docs.yml`](../.github/workflows/docs.yml) on `devel` (the
`github-pages` environment refuses a deployment from a tag), which fetches
them on every build. A run in which a command failed publishes nothing. Untick
*publish-page* for a run that should not replace the published page.

The datasets are deterministic, so they are only generated once per generator
version. They are kept in the release `benchmark-data-v<GENERATOR_VERSION>`,
one tarball per dataset, and fetched from there. Bump `GENERATOR_VERSION` in
`simulate.py` whenever the model or a preset changes.
