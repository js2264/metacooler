# Installation

## Requirements

- Python ≥ 3.10, on Linux or macOS
- Nothing to compile: every dependency publishes wheels for both

## From PyPI

```bash
pip install metacooler
```

## With uv (recommended for development)

```bash
uv add metacooler
```

## Development install

Clone the repository and install with all development dependencies:

```bash
git clone https://github.com/js2264/metacooler
cd metacooler
uv sync
```

The development toolchain (`ruff`, `pytest`, `pre-commit`, `ty`) is declared as
a [dependency group](https://packaging.python.org/en/latest/specifications/dependency-groups/),
which `uv sync` installs by default. With pip ≥ 25.1, the equivalent is
`pip install -e . --group dev`.

## Building the documentation

The documentation site is built with [Zensical](https://zensical.org) and
[mkdocstrings](https://mkdocstrings.github.io), which the `docs` dependency
group installs on its own: mkdocstrings reads metacooler's sources without
importing them, so neither metacooler nor its dependencies are needed.

```bash
uv run --only-group docs zensical serve   # live preview at http://localhost:8000
uv run --only-group docs zensical build   # render the static site into ./site
```

With pip ≥ 25.1, `pip install --group docs` in the clone installs the same.

## Verifying the install

```bash
metacooler --version
```
