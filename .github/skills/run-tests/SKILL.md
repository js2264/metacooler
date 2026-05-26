---
name: run-tests
description: "Run metacooler tests and interpret results. Use when: running tests, checking coverage, debugging test failures, validating changes."
---
# Run metacooler Tests

## When to Use
- After making code changes
- Debugging test failures
- Checking code coverage
- Validating a new feature

## Quick Commands

### Run all tests
```bash
uv run pytest tests/ -x -v
```

### Run specific test file
```bash
uv run pytest tests/test_<name>.py -x -v
```

### Run with coverage
```bash
uv run pytest tests/ --cov=metacooler --cov-report=term-missing
```

### Run specific test
```bash
uv run pytest tests/test_<name>.py::TestClass::test_method -v
```

### Lint check
```bash
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
```

### Auto-fix lint
```bash
uv run ruff check --fix src/ tests/
uv run ruff format src/ tests/
```

## Test Structure

| Test file | Covers |
|-----------|--------|
| `test_api.py` | MetaCooler class |
| `test_schema.py`, `test_schema_validation.py` | Format validation |
| `test_schema_mmetacool.py` | MMetaCooler schema |
| `test_cli.py`, `test_cli_edge_cases.py` | CLI commands |
| `test_cli_mmetacool.py` | ls/cp/mv/rm commands |
| `test_info_mmetacool.py` | info/tree/attrs for mmetacool |
| `test_extract.py` | Genome extraction |
| `test_dump.py` | Table dumping |
| `test_flatten.py` | Flatten to cool |
| `test_network.py` | Network edge list |
| `test_binnify.py` | Re-binning |
| `test_zoomify.py` | Multi-resolution |
| `test_mmetacool_io.py` | MMetaCooler I/O |
| `test_uri.py` | URI parsing |
| `test_aggregate.py` | Aggregation |
| `test_merge.py` | File merging |
| `test_show.py` | Visualization |
| `test_mag.py` | MAG table parsing |
| `test_cload.py` | Basic loading |
| `test_create_edge_cases.py` | Creation edge cases |

## Interpreting Results
- **xfail**: Known failures marked with `@pytest.mark.xfail` — expected to fail
- **PendingDeprecationWarning**: rich-click parameter deprecations — safe to ignore
- **Coverage target**: aim for ≥80% on new code
