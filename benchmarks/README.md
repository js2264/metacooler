# Benchmarks: metacooler vs cooler

`benchmarks.qmd` times metacooler against cooler **from the command line**, on
four reproducible synthetic metagenomes, from 10³ to 2 × 10⁶ contigs:

| | genomes (MAGs) | contigs | assembly | read pairs |
|---|---|---|---|---|
| S | 10 (10) | 1,000 | 27 Mb | 1 M |
| M | 100 (50) | 20,000 | 157 Mb | 10 M |
| L | 2,000 (200) | 200,000 | 764 Mb | 30 M |
| XL | 20,000 (1,000) | 2,000,000 | 5.1 Gb | 100 M |

It covers building a map from pairs (`cooler cload pairs` vs `metacooler create
pairs`), converting between the formats, multi-resolution pyramids (`zoomify`),
querying a contig, a genome and a genome pair, and genome-level aggregation.
Each command runs as a fresh process under GNU `time` and `timeout`. The results
are one row per run in `_results/timings.tsv`, and every table and figure is
drawn from that file.

Every timing is against **cooler 0.12.1**, pinned as `COOLER_VERSION` in
`bench.sh`: the notebook refuses to run with any other cooler, and the workflow
installs that one. (`uv.lock` pins an older cooler for the test suite.)

| file | what |
|---|---|
| `benchmarks.qmd` | the notebook: bash cells that run and time the commands; Python for tables and figures |
| `bench.sh` | settings, and `bench`, the timing function every bash cell sources |
| `simulate.py` | the dataset generator: presets, seeds, model |
| `targets.py` | what the queries ask for: one contig, one MAG and one MAG pair per dataset, from a fixed seed |

## Render

```bash
uv sync && uv pip install cooler==0.12.1   # the cooler the benchmarks are against
source .venv/bin/activate                    # metacooler, cooler and jupyter on PATH
cd benchmarks
quarto render benchmarks.qmd                         # the whole ladder: ~4 h on 4 cores, 16 GB RAM
BENCH_DATASETS="S M" quarto render benchmarks.qmd    # a quick look: minutes
```

Needs GNU `time` (`apt install time`); `pigz` speeds up dataset generation.
Call the tools from the activated environment, not through `uv run`: `uv run`
and `uv sync` put back the cooler that `uv.lock` pins.

| variable | default | |
|---|---|---|
| `BENCH_DATASETS` | `S M L XL` | which datasets, smallest first |
| `BENCH_TIMEOUT` | `1800` | seconds a run may take before it is stopped (and reported as a lower bound) |
| `BENCH_REPS` | 5 / 3 / 2 / 1 | runs per command (S / M / L / XL) |
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

The datasets are deterministic, so they are only generated once per generator
version. They are kept in the release `benchmark-data-v<GENERATOR_VERSION>`,
one tarball per dataset, and fetched from there. Bump `GENERATOR_VERSION` in
`simulate.py` whenever the model or a preset changes.
