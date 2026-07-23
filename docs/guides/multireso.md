# Multi-resolution files

metacooler supports two related operations for working at multiple resolutions:

- **`binnify`** — re-bin a single `.metacool` file to a coarser resolution,
  producing a new `.metacool` file.
- **`zoomify`** — generate a pyramid of resolutions from a base-resolution
  `.metacool`, storing all levels in a single `.mmetacool` collection file.

## Re-binning a single file (`binnify`)

```bash
metacooler binnify sample.metacool \
    --binsize 5000 \
    --output  sample_5kb.metacool
```

The new bin size must be **strictly larger** than the current bin size. Contacts
are reassigned to new bins using midpoint assignment.

```bash
# Also works at genome level (contigs concatenated within each genome)
metacooler binnify sample.metacool --binsize 10000 --level genome -o sample_10kb_genome.metacool
```

## Creating a multi-resolution collection (`zoomify`)

`zoomify` takes a base-resolution `.metacool` and creates a `.mmetacool` file
with multiple resolution levels, each doubling the bin size:

```bash
metacooler zoomify sample.metacool \
    --levels 5 \
    --output sample_multireso.mmetacool
```

This produces:

```
sample_multireso.mmetacool
└── /entries
    ├── 1000   ← base resolution (copied from sample.metacool)
    ├── 2000
    ├── 4000
    ├── 8000
    └── 16000
```

Generation stops early if fewer than 2 bins would remain at the next level.

## Working with collection files

```bash
# List all resolution levels
metacooler ls sample_multireso.mmetacool
# sample_multireso.mmetacool::entries/1000
# sample_multireso.mmetacool::entries/2000
# ...

# Inspect a specific resolution
metacooler info sample_multireso.mmetacool::entries/4000
```

Open a specific resolution in Python:

```python
import metacooler

mc = metacooler.MetaCooler("sample_multireso.mmetacool::entries/4000")
```

## Managing collection entries

```bash
# Add an existing metacool file as a new entry
metacooler stack sample_5kb.metacool sample_10kb.metacool \
    --names "5kb,10kb" \
    -o multi.mmetacool

# Copy an entry
metacooler cp multi.mmetacool::entries/5kb multi.mmetacool::entries/5kb_backup

# Rename an entry
metacooler mv multi.mmetacool::entries/5kb_backup multi.mmetacool::entries/5kb_v2

# Remove an entry
metacooler rm multi.mmetacool::entries/5kb_v2

# Cross-file copy
metacooler cp source.mmetacool::entries/rep1 target.mmetacool::entries/rep1
```

## Python API

```python
import metacooler
from metacooler.ops.binnify import binnify
from metacooler.ops.zoomify import zoomify
from metacooler.io.mmetacool import list_entries, add_entry, remove_entry

# Re-bin
binnify("sample.metacool", output="sample_5kb.metacool", binsize=5000)

# Zoomify
zoomify("sample.metacool", output="sample_multireso.mmetacool", levels=5)

# Inspect collection
entries = list_entries("sample_multireso.mmetacool")
print(entries)  # ['1000', '2000', '4000', ...]
```
