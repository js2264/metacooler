# Format specification

This page is the authoritative reference for the `.metacool` and `.mmetacool`
HDF5 formats. It is intended for developers integrating with metacooler, not
for end users.

## `.metacool` — single-sample file

A `.metacool` file is a valid HDF5 file that is simultaneously a valid cooler
file. Any software that reads cooler files (cooltools, HiGlass, the cooler CLI)
will open a `.metacool` file without modification.

### Top-level structure

```
file.metacool  (HDF5)
├── @attrs
│   ├── format                    = "HDF5::Cooler"         ← cooler compat
│   ├── format-url                = "https://github.com/open2c/cooler"
│   ├── format-version            = 3                      ← cooler format version
│   ├── storage-mode              = "symmetric-upper"
│   ├── metacooler-format         = "HDF5::MetaCooler"     ← identity marker
│   ├── metacooler-format-url     = "https://github.com/js2264/metacooler"
│   ├── metacooler-format-version = 1
│   ├── contig-ordering           = "genome-then-length"   ← ordering strategy
│   ├── genome-assembly           = <str>                  ← "unknown" if unset (CLI: --assembly-name)
│   ├── bin-size                  = <int | null>           ← null for variable bins
│   ├── bin-type                  = "fixed" | "variable"
│   ├── creation-date             = <ISO 8601>
│   ├── generated-by              = "metacooler-<version>"
│   ├── nchroms                   = <int>                  ← = number of contigs
│   ├── ngenomes                  = <int>
│   ├── nbins                     = <int>
│   ├── nnz                       = <int>                  ← number of non-zero pixels
│   ├── sum                       = <int>                  ← total contact count
│   └── metadata                  = <JSON string>          ← user metadata blob
│
├── /genomes    ← MAG / genome table
├── /contigs    ← contig table (canonical chromosome reference)
├── /bins       ← genomic bins
├── /pixels     ← contacts
├── /indexes    ← fast-access offset arrays
└── /chroms     → soft link to /contigs               ← cooler compat
```

### `/genomes` group

One row per genome / MAG.

| Dataset | dtype | Description |
|---------|-------|-------------|
| `name` | `S<n>` (fixed bytes) | Genome name, ASCII only (see [Names must be ASCII](#names-must-be-ascii)). `__unbinned__` is reserved for unassigned contigs. |
| `length` | `int64` | Total length in bp (sum of assigned contig lengths). |
| `n_contigs` | `int32` | Number of contigs assigned to this genome. |
| `n_bins` | `int64` | Number of bins across all contigs of this genome. |
| `qc_*` | `int64` / `float64` | Optional per-genome contact-quality columns (see [Contact-quality block](#qc-block)). Present when the file was created with `--qc` (the default) or annotated with `metacooler qc --store`. |
| *extra columns* | any | Optional metadata columns from the `--genomes` table (e.g. `phylum`, `completeness`, `contamination`). Stored as fixed-length bytes or numeric arrays. |

The row *set* is derived from the contig-to-genome assignments: one row per
distinct genome name in `/contigs`. A `--genomes` metadata table only supplies
the extra columns, left-joined on `name`; it cannot add a genome that owns no
contig, nor suppress one that does.

Row order: `__unbinned__` first (index 0), then remaining genomes sorted
alphabetically by name. This order defines the `genome_id` integers used in
`/contigs`. For genome-grouped orderings (`genome-then-length`,
`genome-then-alpha`) the contigs in `/contigs` are laid out in ascending
`genome_id` order, so `genome_offset` delimits each genome's contig block.

### `/contigs` group

One row per contig. This is the canonical chromosome reference table. The soft
link `/chroms → /contigs` makes it transparent to the cooler ecosystem.

| Dataset | dtype | Description |
|---------|-------|-------------|
| `name` | `S<n>` (fixed bytes) | Contig name, unique within the file, ASCII only (see [Names must be ASCII](#names-must-be-ascii)). |
| `length` | `int32` | Contig length in bp. |
| `genome_id` | `int32` | Zero-based index into `/genomes`. |
| *extra columns* | any | Optional metadata columns from the `--contigs` table (e.g. `gc_content`, `coverage`). ASCII only: cooler reads every column of this group. |

Row order is determined by the `--contig-order` strategy chosen at creation time
and stored in the `contig-ordering` attribute.

### `/bins` group

One row per genomic bin. Identical to the cooler `bins` table, except the
`chrom` column stores an integer contig ID (the row index into `/contigs`)
rather than a string chromosome name or categorical code.

| Dataset | dtype | Description |
|---------|-------|-------------|
| `chrom` | `int32` | Index into `/contigs`. |
| `start` | `int32` | Zero-based start position (inclusive). |
| `end` | `int32` | Zero-based end position (exclusive). |

### `/pixels` group

One row per non-zero contact, stored in upper-triangular form
(`bin1_id ≤ bin2_id`).

| Dataset | dtype | Description |
|---------|-------|-------------|
| `bin1_id` | `int64` | Row bin index (upper-triangular: always ≤ `bin2_id`). |
| `bin2_id` | `int64` | Column bin index. |
| `count` | `int32` or `int64` | Number of contacts. `int32` at base resolution; `int64` wherever counts are summed (`binnify`, `merge`), since aggregation readily exceeds the `int32` range. Readers must not assume `int32`. |

Pixels are sorted by `(bin1_id, bin2_id)` to allow the `bin1_offset` index to
work correctly.

### `/indexes` group

| Dataset | dtype | Length | Description |
|---------|-------|--------|-------------|
| `chrom_offset` | `int64` | `n_contigs + 1` | `chrom_offset[i]` is the index of the first bin belonging to contig `i`. |
| `bin1_offset` | `int64` | `n_bins + 1` | `bin1_offset[i]` is the index of the first pixel with `bin1_id == i`. Written by cooler. |
| `genome_offset` | `int64` | `n_genomes + 1` | `genome_offset[i]` is the index of the first contig belonging to genome `i`. metacooler-specific. |

`genome_offset` is written for every file, but it only *delimits* contiguous
per-genome blocks under a genome-grouped `contig-ordering` (`genome-then-length`,
`genome-then-alpha`). Under `fasta`, `alpha`, `cool` or `custom` ordering a
genome's contigs may be scattered, so readers must fall back to filtering
`/contigs/genome_id` rather than slicing with `genome_offset`.

### Identity check

A file is a valid `.metacool` if and only if:

1. It is a valid HDF5 file.
2. `f.attrs["metacooler-format"] == "HDF5::MetaCooler"`.
3. The `/genomes` group exists.
4. It passes `metacooler.schema.validate()` without raising.

```python
import metacooler
print(metacooler.is_metacooler("gut_T0.metacool"))  # True
```

---

## `.mmetacool` — multi-sample collection

A `.mmetacool` file is a container for multiple `.metacool` datasets, stored
as HDF5 groups. It is **not** itself a valid cooler file.

### Top-level structure

```
file.mmetacool  (HDF5)
├── @attrs
│   ├── format          = "HDF5::MMetaCooler"
│   └── format-version  = 1
│
└── /entries
    ├── <entry_name_1>   ← full metacooler (same layout as a .metacool root)
    ├── <entry_name_2>
    └── ...
```

Each entry under `/entries` contains a complete metacooler: its own `/genomes`,
`/contigs`, `/bins`, `/pixels`, `/indexes`, and top-level attributes. Entries
are independent; they do not share data.

### Addressing entries

Use the `::` URI syntax to address a specific entry:

```bash
metacooler info gut_T0.mmetacool::entries/rep1
```

```python
mc = metacooler.MetaCooler("gut_T0.mmetacool::entries/rep1")
```

### Multi-resolution files (zoomify output)

`metacooler zoomify` produces a `.mmetacool` file where each entry is named
after its bin size:

```
gut_T0_multireso.mmetacool
└── /entries
    ├── 1000    ← base resolution
    ├── 2000
    ├── 4000
    └── 8000
```

### Identity check

```python
print(metacooler.is_mmetacooler("gut_T0.mmetacool"))  # True
```

---

## Contig ordering strategies

The order of contigs in `/contigs` and `/bins` is fixed at creation time and
stored in the `contig-ordering` attribute. Merging two files requires identical
ordering. `reassign` is the one command that can change it after the fact: a
re-grouping rewrites the rows into the named strategy, and `--layout keep`
downgrades the attribute to `custom` rather than let the file keep claiming a
grouping its rows no longer have.

| Strategy | Description |
|----------|-------------|
| `genome-then-length` | Group by genome (alphabetical), longest contig first within each genome. **Default.** |
| `genome-then-alpha` | Group by genome (alphabetical), contigs alphabetical within each genome. |
| `fasta` | Preserve the order of contigs as they appear in the FASTA or pairs header. |
| `alpha` | All contigs sorted alphabetically, ignoring genome assignment. |
| `cool` | Preserve the exact contig order of the input cooler. `create cool` only. |
| `custom` | Row order inherited from an earlier layout and no longer described by any strategy. Written by `reassign --layout keep` when the new assignment breaks genome grouping; not a value `create` can produce. |

---

## Contact-quality block { #qc-block }

`create` records where the library's contacts landed, and `metacooler qc`
recomputes it on demand. It occupies two places in the file, both optional —
a file written with `--no-qc` has neither.

**`metadata["qc"]`**, a JSON object inside the root `metadata` attribute,
alongside the `commands` provenance list:

| Key | Meaning |
|-----|---------|
| `version` | Block schema version (currently `1`). |
| `source` | `"pairs"` — measured on read positions by `create pairs`, before binning — or `"pixels"`, measured between bin midpoints. |
| `distance_threshold` | Short/long intra-contig split, in bp (default 10 000). |
| `n_contacts` | Total contacts classified. |
| `intra_contig_short`, `intra_contig_long`, `inter_contig_intra_genome`, `inter_genome` | The four mutually exclusive categories; they sum to `n_contacts`. |
| `intra_contig`, `inter_contig`, `intra_genome` | Convenience sums of the above. |
| `trans_contig_ratio` | The "3D ratio": `inter_contig / n_contacts`. |
| `*_fraction` | Each category's share of `n_contacts`. `null` on an empty file. |
| `binsize` | Bin size at measurement time (pixel source). |
| `n_pairs_read`, `n_pairs_discarded` | Pairs seen and dropped while binning (pairs source), so the denominator the ratio does *not* use stays visible. |

**`/genomes` columns**, one row per genome:
`qc_intra_contig_short`, `qc_intra_contig_long`,
`qc_inter_contig_intra_genome`, `qc_inter_genome` (`int64`), and
`qc_intra_genome_ratio` (`float64`, `NaN` for a genome with no contact). An
inter-genome contact is credited to both genomes involved, so these columns sum
to more than the file-level totals.

Both are statements about one specific pixel table, so `binnify` and `merge`
drop them from their output instead of carrying them forward, and `merge`
ignores `qc_*` columns when it compares the inputs' `/genomes` tables for
compatibility — two libraries over one assembly legitimately differ there.

!!! note "P(s) is a companion file, not part of the format"
    `create pairs` also writes a per-MAG P(s) curve, but as a plain TSV beside
    the map (`gut_T0.metacool` → `gut_T0.ps.tsv`), never inside the HDF5. It is
    measured at base-pair resolution from the read positions, so it holds
    detail the binned file cannot represent, and it is a plotting input rather
    than something any reader of the map needs. See [`ps`](cli.md#ps).

## Cooler compatibility detail

The compatibility mechanism relies on two things written by metacooler at file
creation time:

1. **Soft link** `/chroms → /contigs`: cooler looks up `/chroms` for the
   chromosome reference table. The soft link redirects it transparently to
   `/contigs`.
2. **`format = "HDF5::Cooler"` attribute**: cooler's `is_cooler()` check passes.

The consequence is that `cooler.Cooler("gut_T0.metacool").chroms()[:]` returns
the contigs table, and all cooler range queries work at contig resolution.
The genome hierarchy is invisible to cooler.

### Names must be ASCII

Because `/contigs` *is* cooler's `/chroms`, it inherits cooler's assumption
about the text it holds: cooler decodes those byte columns with an ASCII
`.astype("U")`. A single non-ASCII byte anywhere in the group — in `name` or in
any extra metadata column — therefore yields a file that satisfies every
structural rule above and still raises `UnicodeDecodeError` in
`cooler.Cooler()`, and so in `MetaCooler()`, which opens through it. Genome
names are subject to the same rule for the opposite reason: cooler never reads
`/genomes`, but `aggregate --level genome` writes those names into a `.cool` as
chrom names, where the encode fails instead of the decode.

So, in both directions:

| Where | Rule |
|-------|------|
| `/contigs`, every string column | ASCII only |
| `/genomes/name` | ASCII only |
| `/genomes`, extra metadata columns | any UTF-8 — cooler never sees them |
| HDF5 attributes, metadata JSON | any UTF-8 |

`create` refuses non-ASCII names before it bins anything, `validate()` refuses
them on read, and `check` reports them; the limitation is upstream, in cooler,
but a file that breaks it is one metacooler can write and then not open, which
is the one thing this format promises never to do.
