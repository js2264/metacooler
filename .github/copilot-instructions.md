# Copilot Instructions — metacooler

## Project Overview

**metacooler** is a Python library and CLI extending the [cooler](https://github.com/open2c/cooler) HDF5 format for **multi-genome (metagenomic) Hi-C contact maps**. It adds a `genomes/` hierarchy on top of cooler's `chroms/bins/pixels` model while preserving full backward compatibility — any `.metacool` file opens with `cooler.Cooler()`.

- **Status**: Alpha (v0.1.0)
- **Python**: ≥3.10 (use `match`, `X | Y` union types, `from __future__ import annotations`)
- **Build**: hatchling
- **Package manager**: `uv` (not pip directly)
- **CLI framework**: rich-click (wraps click)
- **Linting**: ruff (check + format). **No black.**
- **Testing**: pytest with pytest-cov. Run: `uv run pytest tests/ -x -v`
- **Format**: HDF5 via h5py + cooler

## Architecture

```
src/metacooler/
├── api.py          # MetaCooler class — primary Python API
├── schema.py       # HDF5 format constants, dtypes, validators
├── uri.py          # URI parsing (file::group_path syntax)
├── _version.py     # Version string
├── cli/            # rich-click CLI commands (one file per command)
│   ├── main.py     # Root CLI group + command registration
│   ├── cmd_*.py    # Individual commands
│   └── logging_utils.py
├── io/             # I/O operations (create, extract, dump, flatten, network, merge, mmetacool)
└── ops/            # Computational operations (aggregate, binnify, zoomify, show)
```

### Key patterns

- **One CLI command per file**: `cli/cmd_<name>.py` imports `cli` from `main.py` and decorates with `@cli.command()`
- **I/O vs Ops separation**: `io/` handles file reading/writing, `ops/` handles computation
- **MetaCooler class**: wraps `cooler.Cooler` with genome/contig-aware accessors
- **Lazy imports in `__init__.py`**: use `__getattr__` pattern, not barrel imports
- **Streaming creation**: pairs processed in chunks for large files
- **URI syntax**: `file.mmetacool::groups/entry1` for addressing entries in multi-metacool files

## Data Format

### `.metacool` (single file)
```
/genomes:   [name, length, n_contigs]
/contigs:   [name, length, genome_id]
/bins:      [chrom(=contig_id), start, end]
/pixels:    [bin1_id, bin2_id, count, ...]
/indexes:   [chrom_offset, bin1_offset, genome_offset]
/chroms:    → soft link to /contigs (cooler compat)
```

### `.mmetacool` (collection)
```
@attrs: format=HDF5::MMetaCooler, format-version=1
/groups/<name>: each is a full metacooler
```

## Code Style

- **Ruff config**: line-length=88, select=E,F,I,UP,B,SIM, ignore=E501
- **Quotes**: double
- **Type hints**: use `X | Y` union syntax, `from __future__ import annotations`
- **Docstrings**: Google-style, include Args/Returns/Raises sections
- **Logging**: use `logging.getLogger("metacooler")`, not print
- **HDF5 string handling**: convert object/str dtypes to fixed-length bytes (`S`) for h5py writes

## Build & Test

```bash
# Install
uv sync

# Test
uv run pytest tests/ -x -v
uv run pytest tests/ --cov=metacooler --cov-report=term-missing

# Lint
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/

# Run CLI
uv run metacooler --help
uv run metacooler info out.metacool
```

## Conventions

- **Never use `__init__.py` barrel imports** — import from submodules directly
- **CLI commands always have `-v/--verbose`** for debug logging via `setup_logging()`
- **Backward compat for CLI options**: when renaming (e.g. `--genome` → `--genomes`), keep the old flag as a deprecated alias with `warnings.warn()`
- **Comma-separated multi-value options**: use `callback=_parse_comma_list` pattern for `--genomes`, `--contigs`
- **HDF5 writes**: always convert pandas string columns to bytes before `create_dataset()` — h5py cannot write `object` or `StringDtype` arrays
- **Test fixtures**: defined in `tests/conftest.py` — `tiny_fasta`, `tiny_mag_table`, `tiny_bins_bed`, `tiny_pairs`, `tiny_metacooler`, `tiny_mmetacool`
- **Test data**: static files in `tests/data/` (assembly.fa, contacts.pairs, bins.bed, etc.)
- **xfail tests**: known issues marked with `@pytest.mark.xfail` — do not remove without fixing

## Key Dependencies

| Package | Role |
|---------|------|
| `cooler` | Core HDF5 contact map format |
| `h5py` | Low-level HDF5 access |
| `pandas` | DataFrames for tables |
| `numpy` | Numeric arrays |
| `scipy` | Sparse matrices for aggregation |
| `pyfaidx` | FASTA parsing for contig lengths |
| `rich-click` | CLI with rich formatting |
| `bioframe` | Genomic interval operations |
