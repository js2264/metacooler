# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**metacooler** extends the [cooler](https://github.com/open2c/cooler) HDF5 format for **multi-genome (metagenomic) Hi-C contact maps**. It layers a `genomes/` → `contigs/` hierarchy on top of cooler's `chroms/bins/pixels` model. The central invariant: **every `.metacool` file must remain openable by `cooler.Cooler()`**. This is achieved by writing `format="HDF5::Cooler"`, a soft-link `chroms → contigs`, and standard `indexes/chrom_offset`, alongside the metacooler-specific `metacooler-format="HDF5::MetaCooler"` attribute and `genomes/` group. Do not break this dual identity when editing the format.

Two file types exist:
- `.metacool` — a single metagenomic contact map (a cooler + `genomes/` hierarchy).
- `.mmetacool` — a collection: `@format="HDF5::MMetaCooler"` with an `entries/<name>/` group per metacool. Addressed via cooler-style URIs: `file.mmetacool::entries/<name>`.

## Commands

Uses `uv` (not pip directly) and `ty` as the type checker (**not mypy**).

```bash
uv sync --dev                                  # install with dev deps

uv run pytest tests/ -x -v                     # run tests (fail fast); --cov is on by default via pyproject
uv run pytest tests/test_create_correctness.py -x       # single test file
uv run pytest tests/test_create.py::test_name -x        # single test

uv run ruff check src/ tests/                  # lint
uv run ruff format src/ tests/                 # format (line-length 88, double quotes; no black)
uv run ty check src/ tests/                    # type check

make fix          # ty check + ruff format + ruff check --fix + regenerate README.md from README.qmd
make regenerate   # rebuild tests/data/ fixtures (generate-raw-files.py + raw-to-cool.sh)
make cli-tests    # run tests/data/cli_tests.sh end-to-end CLI smoke tests

uv run metacooler --help                       # run CLI
```

**README.md is generated** from `README.qmd` via `quarto render` (see the `make fix` target). Edit `README.qmd`, never `README.md` directly — it is also excluded from ruff.

## Architecture

```
src/metacooler/
├── api.py       # MetaCooler class — the primary Python API
├── schema.py    # format constants, dtypes, is_metacooler/is_mmetacooler, validate(), read_table()
├── uri.py       # parse_uri("file::group") → (Path, group_path)
├── cli/         # one file per command: cmd_<name>.py
├── io/          # file read/write: create, extract, dump, flatten, network, merge, mmetacool, contigs, genomes
└── ops/         # computation: aggregate, binnify, zoomify, show, range
```

**Layer separation is intentional**: `io/` reads and writes files; `ops/` does computation on already-open files. Keep new code on the correct side.

**`MetaCooler` (api.py)** wraps a `cooler.Cooler` instance and delegates standard reads (`bins`, `pixels`, `matrix`, `chroms`) to it, adding genome/contig-aware accessors (`genomes`, `contigs`, `bins_for_genome`, `pixels_for_genome_pair`, `genome_extent`, …). The constructor accepts a URI: a bare path opens a `.metacool`; a `file.mmetacool::entries/<name>` URI opens one entry. Opening the root of an `.mmetacool` raises `KeyError` listing the available entries.

**On-disk layout** (`.metacool`):
```
/genomes  [name, length, n_contigs, n_bins]
/contigs  [name, length, genome_id]
/bins     [chrom(=contig id), start, end]
/pixels   [bin1_id, bin2_id, count, ...]      # sorted by (bin1_id, bin2_id)
/indexes  [chrom_offset, bin1_offset, genome_offset]
/chroms   → soft link to /contigs             # cooler compatibility
```
`schema.validate()` enforces this (required groups/columns, offset lengths = count+1, pixel sort order, genome_id range). Run it mentally when changing anything that writes the format.

## Conventions (non-obvious)

- **Never add barrel imports to `__init__.py`.** Public symbols are exposed via a lazy `__getattr__` dispatch (one `if name == ...` branch per symbol) to keep import time low and avoid circular imports. Add new public API by extending `__all__` and adding a branch.
- **CLI commands self-register**: each `cli/cmd_<name>.py` imports `cli` from `cli/main.py` and decorates with `@cli.command()`; `main.py` imports every `cmd_*` module at the bottom to trigger registration, and lists it in `COMMAND_GROUPS`. Adding a command = new file + one import line + one entry in the right group.
- **URI arguments** use the `MetaCoolURI` click param type (in `cli/main.py`), which validates the file/group exists and passes the raw string through for downstream `parse_uri()`.
- **HDF5 string writes**: h5py cannot write pandas `object`/`StringDtype` columns — convert to fixed-length bytes (`S`) before `create_dataset()`. On read, `schema.read_table()` decodes `S`/bytes back to `str` (matching cooler).
- **Logging**: use `logging.getLogger("metacooler")`, never `print`. Every CLI command exposes `-v/--verbose` for debug logging.
- **Streaming creation**: `io/create.from_pairs` processes pairs in chunks (`--chunksize`) with an `auto`/`memory`/`external` sort strategy (2 GiB threshold) so large files don't blow memory. Preserve this when touching the create path.
- **Python ≥3.10**: use `match`, `X | Y` unions, `from __future__ import annotations`.
- **Comma-separated multi-value options** (`--genomes`, `--contigs`) use the `_parse_comma_list` click callback pattern (see `cli/cmd_extract.py`, `cli/cmd_dump.py`) so a value can be repeated or given as a comma-joined list.

## Testing

- Fixtures live in `tests/conftest.py`: `tiny_fasta`, `tiny_contigs_table`, `tiny_bins_bed`, `tiny_pairs`, `tiny_genomes_table`, `tiny_metacooler`, `tiny_cool`, `tiny_mmetacool`, plus variants for each contig-table format (`contigs`, `contig_genome`, `metator`, column-index) and edge cases (variable bins, reversed cool order).
- Static test data is in `tests/data/`, regenerated with `make regenerate`.
- `@pytest.mark.xfail` marks known issues — do not remove a marker without fixing the underlying bug.

## Reference

- `COMMANDS.md` — full CLI reference (every command, options, and the API entry point each maps to).
- `docs/format.md` — the on-disk format specification.
- `reviewed/` — dated code-review notes from prior passes.
