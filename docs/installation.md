# Installation

## Requirements

- Python ≥ 3.10
- A C compiler (required by `h5py`) — present by default on Linux and macOS

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
uv sync --dev
```

This installs `ruff`, `pytest`, `pre-commit`, and `ty` in addition to the runtime dependencies.

## Optional dependencies

### Plotly backend

An optional Plotly backend for interactive heatmaps (not yet exposed in the CLI):

```bash
pip install "metacooler[plotly]"
```

### Documentation

The documentation site is built with [Zensical](https://zensical.org) and
[mkdocstrings](https://mkdocstrings.github.io). To build and preview it locally:

```bash
pip install "metacooler[docs]"   # or: uv sync --extra docs
zensical serve                   # live preview at http://localhost:8000
zensical build                   # render the static site into ./site
```

## Verifying the install

```bash
metacooler --version
```
