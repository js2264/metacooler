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
│   ├── format              = "HDF5::Cooler"          ← cooler compat
│   ├── format-url          = "https://github.com/open2c/cooler"
│   ├── metacooler-format   = "HDF5::MetaCooler"      ← identity marker
│   ├── metacooler-format-url = "https://github.com/js2264/metacooler"
│   ├── metacooler-format-version = 1
│   ├── contig-ordering     = "genome-then-length"    ← ordering strategy
│   ├── assembly-name       = <str | null>
│   ├── bin-size            = <int | null>             ← null for variable bins
│   ├── bin-type            = "fixed" | "variable"
│   ├── creation-date       = <ISO 8601>
│   ├── generated-by        = "metacooler-<version>"
│   ├── nbins               = <int>
│   ├── nnz                 = <int>
│   └── metadata            = <JSON string>           ← user metadata blob
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
| `name` | `S<n>` (fixed bytes) | Genome name. `__unbinned__` is reserved for unassigned contigs. |
| `length` | `int64` | Total length in bp (sum of assigned contig lengths). |
| `n_contigs` | `int32` | Number of contigs assigned to this genome. |
| `n_bins` | `int64` | Number of bins across all contigs of this genome. |
| *extra columns* | any | Optional metadata columns from the `--genomes` table (e.g. `phylum`, `completeness`, `contamination`). Stored as fixed-length bytes or numeric arrays. |

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
| `name` | `S<n>` (fixed bytes) | Contig name, unique within the file. |
| `length` | `int32` | Contig length in bp. |
| `genome_id` | `int32` | Zero-based index into `/genomes`. |
| *extra columns* | any | Optional metadata columns from the `--contigs` table (e.g. `gc_content`, `coverage`). |

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
| `count` | `int32` | Number of contacts. |

Pixels are sorted by `(bin1_id, bin2_id)` to allow the `bin1_offset` index to
work correctly.

### `/indexes` group

| Dataset | dtype | Length | Description |
|---------|-------|--------|-------------|
| `chrom_offset` | `int64` | `n_contigs + 1` | `chrom_offset[i]` is the index of the first bin belonging to contig `i`. |
| `bin1_offset` | `int64` | `n_bins + 1` | `bin1_offset[i]` is the index of the first pixel with `bin1_id == i`. Written by cooler. |
| `genome_offset` | `int64` | `n_genomes + 1` | `genome_offset[i]` is the index of the first contig belonging to genome `i`. metacooler-specific. |

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
ordering.

| Strategy | Description |
|----------|-------------|
| `genome-then-length` | Group by genome (alphabetical), longest contig first within each genome. **Default.** |
| `genome-then-alpha` | Group by genome (alphabetical), contigs alphabetical within each genome. |
| `fasta` | Preserve the order of contigs as they appear in the FASTA or pairs header. |
| `alpha` | All contigs sorted alphabetically, ignoring genome assignment. |
| `custom` | Explicit order supplied via the Python API (`custom_order` argument). CLI only. |

---

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
