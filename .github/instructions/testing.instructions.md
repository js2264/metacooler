---
description: "Use when writing or modifying tests for metacooler. Covers fixtures, test patterns, CLI testing."
applyTo: "tests/**"
---
# Testing Conventions

## Running tests
```bash
uv run pytest tests/ -x -v                    # stop on first failure
uv run pytest tests/test_<name>.py -v          # single file
uv run pytest tests/ --cov=metacooler          # with coverage
```

## Fixtures (in `tests/conftest.py`)
- `tiny_fasta` — 4-contig FASTA (A=1000bp, B=600bp, C=400bp, D=200bp)
- `tiny_contigs_table` — 2-MAG assignment TSV (contig_D → `__unbinned__`)
- `tiny_bins_bed` — fixed 100bp bins for 4 contigs
- `tiny_pairs` — 200 synthetic contacts (seed=42)
- `tiny_cool` — `.cool` file built from `tiny_pairs` + `tiny_bins_bed` via `cooler.create_cooler()`
- `tiny_metacooler` — complete metacool created via `create_metacooler()` API
- `tiny_mmetacool` — 2-entry mmetacool created via `create_mmetacool()` API
- `tiny_genomes_table` — genomes metadata with extra columns

## Test patterns
- **I/O tests**: create temp files, verify with `MetaCooler()` or `cooler.Cooler()`
- **CLI tests**: ALWAYS use `click.testing.CliRunner` — never use `subprocess` or `sys.executable` to invoke the metacooler CLI:
  ```python
  from click.testing import CliRunner
  from metacooler.cli.main import cli

  runner = CliRunner()
  result = runner.invoke(cli, ["command", "--option", "value"])
  assert result.exit_code == 0, result.output
  ```
  Check `result.exit_code` (not `.returncode`) and `result.output` (not `.stdout`/`.stderr`).
  Unhandled exceptions appear in `result.exception`, not `result.output`.
- **External tool tests** (e.g. the `cooler` CLI): subprocess is acceptable when explicitly testing third-party tool compatibility.
- **Error tests**: `pytest.raises(ExceptionType, match="pattern")`
- **Known failures**: `@pytest.mark.xfail(reason="...")` — do not remove without fixing
- **Expensive fixtures**: use `scope="module"` with `tmp_path_factory` when a fixture is shared across many tests to avoid repeated setup cost.

## Static test data
Located in `tests/data/`: assembly.fa, contacts.pairs, bins.bed, contigs_info.tsv, MAG_info.tsv, toy.metacool, toy.mmetacool, metadata.json/yaml.

## Test file naming
- `test_<module>.py` for unit tests of a specific module
- `test_cli.py` / `test_cli_edge_cases.py` for CLI tests
- `test_cli_mmetacool.py` for mmetacool CLI commands
- `test_info_mmetacool.py` for info/tree/attrs mmetacool support
