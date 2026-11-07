# Multi-resolution files

metacooler supports two related operations for working at multiple resolutions:

- **`binnify`** — re-bin a single `.metacool` file to a coarser resolution,
  producing a new `.metacool` file.
- **`zoomify`** — generate a pyramid of resolutions from a base-resolution
  `.metacool`, storing all levels in a single `.mmetacool` collection file.

## Re-binning a single file (`binnify`)

```bash
# 1 kb source → 5 kb output (5000 = 5 × 1000 ✓)
metacooler binnify sample.metacool \
    --binsize 5000 \
    --output  sample_5kb.metacool
```

Contacts are reassigned to new bins using midpoint assignment.

!!! warning "The target bin size is not arbitrary"
    For a **fixed-bin** source, `--binsize` must be a strict, *exact multiple*
    of the current bin size. From a 1 kb file, the valid targets are 2000, 3000,
    4000, 5000, … — `--binsize 3500` is now rejected:

    ```console
    $ metacooler binnify sample_1kb.metacool --binsize 3500 -o out.metacool
    Error: New binsize (3500) must be an exact multiple of the source binsize
    (1000) — 3500 % 1000 = 500. […] Nearest valid multiples: 3000 and 4000.
    ```

    Older versions accepted this. The midpoint rule then packed an oscillating
    number of source bins into each output bin (1 kb → 3.5 kb alternates 3 and 4
    source bins, a built-in ±14 % checkerboard) while the file advertised a
    uniform 3500 bp bin size. If you produced such a file, re-bin from the base
    resolution to a valid multiple.

    For a **variable-bin** source (restriction fragments, contig-level files)
    there is no exact multiple to check, so the *median* source bin width is the
    effective resolution: `--binsize` must be strictly larger than it. A target
    between the median and the widest source bin is accepted with a warning —
    those widest bins each collapse into one output bin, so the effective
    resolution is uneven across the file. Re-binning such a file to a *finer*
    size is refused outright (it used to be accepted and left almost every
    output bin empty).

```bash
# Also works at genome level (contigs ligated within each genome, longest first)
metacooler binnify sample.metacool --binsize 10000 --level genome -o sample_10kb_genome.metacool

# --chunksize bounds the peak memory of the pixel remapping pass
metacooler binnify sample.metacool --binsize 5000 --chunksize 200000 -o sample_5kb.metacool
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

Generation **stops early**: a level that would leave fewer than 2 bins, or that
would duplicate an already-written level (once the bin size exceeds the longest
contig, every contig collapses to a single bin), is skipped and reported in the
log rather than computed and thrown away.

### Choosing resolutions explicitly

`-r/--resolutions` uses cooler's own **ladder-token grammar** — deliberately
*not* the plain comma list taken by `--genomes`/`--contigs` elsewhere in the
CLI, so that `metacooler zoomify` and `cooler zoomify` stay interchangeable:

| Form | Meaning |
|------|---------|
| `1000,2000,5000` | An explicit list. |
| `1000N` | "Nice" 1-2-5 progression from 1000. |
| `1000B` | Binary ×2 progression from 1000. |
| `N` / `B` | The same, starting from the base bin size. |
| `4DN` | The 4DN standard ladder (`1000,2000,5000N`). |

```bash
metacooler zoomify sample.metacool --resolutions 1000,2000,5000,10000 -o sample_multireso.mmetacool
```

!!! warning "Every level must nest in the base bin size"
    A requested resolution that is not an exact multiple of the base is
    rejected — the pyramid is built by re-binning the base, so a level that does
    not nest in it would not be self-consistent. With a 1 kb base,
    `--resolutions 1500,3000` fails; the message names the nearest valid
    resolution. Be careful with `N`-style tokens against a non-1 kb base: with a
    4 kb base, `1000N` expands to 5000/10000/20000/50000 — not one of which is a
    multiple of 4000.

### Controlling temporary disk

Each level is written to a temporary directory, appended to the output, and
deleted immediately, so peak scratch use is one level rather than a second copy
of the whole pyramid. Point `--temp-dir` at fast local scratch on a cluster:

```bash
metacooler zoomify sample.metacool --levels 5 \
    --temp-dir /scratch/$USER \
    -o sample_multireso.mmetacool
```

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

`binnify` and `zoomify` are exported from the package root:

```python
import metacooler

# Re-bin (binsize must be an exact multiple of the source bin size)
metacooler.binnify(
    "sample.metacool", output="sample_5kb.metacool", binsize=5000, chunksize=1_000_000
)

# Zoomify
metacooler.zoomify(
    "sample.metacool",
    output="sample_multireso.mmetacool",
    levels=5,
    temp_dir="/scratch/tmp",
)

# Inspect collection
entries = metacooler.list_entries("sample_multireso.mmetacool")
print(entries)  # ['1000', '2000', '4000', ...]
```

## Viewing in HiGlass { #higlass }

HiGlass reads multi-resolution coolers laid out as `/resolutions/<binsize>`
(`.mcool`). The `.mmetacool` that `metacooler zoomify` writes stores its levels
under `/entries/<binsize>` instead, so HiGlass cannot ingest it directly. Since
every `.metacool` is also a valid cooler, cooler's own `zoomify` produces an
`.mcool` from it:

```bash
cooler zoomify -r 1000,2000,4000,8000 -o sample.mcool sample.metacool
```

The result is a plain `.mcool`: contigs are its chromosomes, and the genome
layer is not carried over.

