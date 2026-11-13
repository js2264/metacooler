---
title: metacooler
hide:
  - navigation
  - toc
  - footer
---

<div class="mc-landing" markdown>

<h1 class="mc-masthead"><img src="assets/metacooler-mark.svg" alt="" width="64" height="64"><span>metacooler</span></h1>

<p class="mc-tagline">Hi-C contact maps for metagenomes.</p>

[metacooler](https://github.com/js2264/metacooler) extends the [cooler](https://github.com/open2c/cooler) HDF5 format to
handle multi-genome (metagenomic) Hi-C data. It adds a `genomes/` hierarchy on
top of cooler's `chroms/bins/pixels` model while preserving **full backward
compatibility**: any `.metacool` file opens with `cooler.Cooler()` out of the
box.
{ .mc-lede }

<div class="mc-install" markdown>

```bash
pip install metacooler
```

</div>

<div class="mc-actions" markdown>

[Why metacooler?](why.md){ .mc-button }
[Quickstart :lucide-arrow-right:](quickstart.md){ .mc-button .mc-button--filled }

</div>

<div class="mc-badges" markdown>

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue)](https://www.python.org/downloads/)
[![PyPI version](https://img.shields.io/pypi/v/metacooler.svg)](https://pypi.org/project/metacooler/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/js2264/metacooler/branch/devel/graph/badge.svg?token=XMYgQxb7KH)](https://codecov.io/gh/js2264/metacooler)

</div>

<p class="mc-kicker">Start here</p>

<div class="mc-tiles" markdown>

[<span class="mc-glyph">:lucide-rocket:</span> <strong>Quickstart</strong> <span class="mc-tile-text">From raw pairs to a genome heatmap in five commands.</span>](quickstart.md){ .mc-tile }
[<span class="mc-glyph">:lucide-lightbulb:</span> <strong>Concepts</strong> <span class="mc-tile-text">The mental model: genomes, contigs, bins, pixels.</span>](concepts.md){ .mc-tile }
[<span class="mc-glyph">:lucide-square-terminal:</span> <strong>CLI reference</strong> <span class="mc-tile-text">Every command and option.</span>](cli.md){ .mc-tile }
[<span class="mc-glyph">:lucide-braces:</span> <strong>Python API</strong> <span class="mc-tile-text">The <code>MetaCooler</code> class and functions.</span>](api/metacooler.md){ .mc-tile }

</div>

</div>
