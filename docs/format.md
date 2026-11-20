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
│   └── metadata                  = <JSON string>          ← user metadata + reserved keys (see below)
│
├── /genomes    ← MAG / genome table
├── /contigs    ← contig table (canonical chromosome reference)
├── /bins       ← genomic bins
├── /pixels     ← contacts
├── /indexes    ← fast-access offset arrays
├── /clusters   ← optional: the clustering rounds that led to /genomes
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
| *extra columns* | any | Optional metadata columns from the `--genomes` table (e.g. `phylum`, `completeness`, `contamination`), at creation or later with `metacooler annotate`. Stored as fixed-length bytes or numeric arrays. A column called `genome` is stored as `user_genome`: `genome` is the key `dump` prints this table under. |

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
| `qc_intra_contig_short`, `qc_intra_contig_long` | `int64` / `float64` | Optional per-contig contact-quality columns (see [Contact-quality block](#qc-block)), present whenever the `metadata["qc"]` block was stored with them. |
| *extra columns* | any | Optional metadata columns from the `--contigs` table (e.g. `gc_content`, `coverage`), at creation or later with `metacooler annotate`. ASCII only: cooler reads every column of this group. Every member of the group must be a dataset; a `/` in a user column name is stored as `_`, since HDF5 would read it as a nested group. A column called `contig` or `genome` is stored as `user_contig` / `user_genome`: those are the keys `dump` prints this table under. |

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
| *weight columns* | `float64` | Optional normalization weights, one per bin — see [Normalization weights](#weights). |

#### Normalization weights { #weights }

A normalization is stored as cooler stores a balancing vector: one `float64`
column of `/bins`, one **multiplicative** weight per bin, `NaN` for a masked
bin. The normalized value of a pixel is `count × w[bin1_id] × w[bin2_id]`, which
is what `cooler.Cooler.matrix(balance=NAME)` computes. Several columns may
coexist. `metacooler normalize ice` writes one (`weight`, by default), any
other tool stores its own with `metacooler.store_weights`, and `show`, `dump`
and `network --balance NAME` read them. Pixels are never rewritten.

A column written by `normalize ice` or `store_weights` carries HDF5 attributes:

| Attribute | Description |
|-----------|-------------|
| `normalization` | The method: `ice`, `normcc` or `hiczin`, or the name the writer gave its own (`store_weights(method=...)`). |
| `divisive_weights` | Always `false`, as `cooler balance` writes it. |
| `depends_on` | What the weights were computed from, a list of: `pixels`; `genome_id`, the contig-to-genome assignment, for a weight computed per genome (or fitted on the pairs of contigs of one genome); `contigs/<column>` for each `/contigs` column read (a coverage, a site count). `normalize ice` writes `pixels`, plus `genome_id` at its default `genome` scope; `hiczin` writes `genome_id` and its coverage column. |
| `spurious_cutoff`, `spurious_rule` | Optional: a contig pair whose normalized contact is below the cutoff (`<`) or at most the cutoff (`<=`) is spurious. NormCC and HiCzin store one; `network --drop-spurious` applies it. |
| *parameters* | The method's parameters and fit statistics — `scope`, `converged`, `coefficients`, … — or the writer's own (`store_weights(parameters=...)`). |
| `metacooler_version` | The version that stored them. |

A column without these attributes — one `cooler balance` wrote — is read as
depending on the pixels alone. Names cooler reads as divisive (`KR`, `VC`,
`VC_SQRT`) are inverted on read, and refused as names for a metacooler weight.

Weights are derived data: `reassign` drops the columns whose `depends_on`
includes `genome_id` when a contig changes genome, and `annotate` those that
name a `/contigs` column it changes. Writers that rebuild the pixels or the
bins (`binnify`, `merge`, `collect`, `extract`, `flatten`) write none.

### `/pixels` group

One row per non-zero contact, stored in upper-triangular form
(`bin1_id ≤ bin2_id`).

| Dataset | dtype | Description |
|---------|-------|-------------|
| `bin1_id` | `int64` | Row bin index (upper-triangular: always ≤ `bin2_id`). |
| `bin2_id` | `int64` | Column bin index. |
| `count` | `int32` or `int64` | Number of contacts. `int32` at base resolution; `int64` wherever counts are summed (`binnify`, `merge`), since aggregation readily exceeds the `int32` range. Readers must not assume `int32`. |

Pixels are sorted by `(bin1_id, bin2_id)` to allow the `bin1_offset` index to
work correctly, and strictly so: a `(bin1_id, bin2_id)` pair appears at most
once. `validate()` and `check --deep` refuse a repeated pair, which cooler would
sum on some reads and not on others.

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

### `/clusters` group { #clusters }

Optional. The clustering rounds of the latest run of a method that groups
contigs into genomes, whatever the method: the successive steps of one binner
([metaTOR](https://github.com/koszullab/metaTOR) writes its overlapping,
recursive and final clusters), the bin sets of several binners and the bins a
refiner keeps from them, or clusters that share contigs. They are kept
without being written into `/genomes` and `/contigs`, which hold the current
assignment: `/clusters` says how a method reached one, and nothing requires
them to agree (a later `reassign` changes `/genomes`, not `/clusters`;
`reassign --round` makes a round the map's genomes). A cluster is a set of
contigs; a contig may be a member of several clusters, of one round or of
several.

```
/clusters
├── @rounds                  ← JSON list, one object per round, in order
├── table/                   ← one row per cluster, every round
├── members/                 ← one row per (cluster, contig)
└── indexes/
    ├── cluster_offset       ← members/ rows of each cluster
    ├── contig_offset        ← contig_clusters rows of each contig
    └── contig_clusters      ← members/cluster_id, ordered by contig, then cluster
```

**`@rounds`**, a JSON list of objects, one per round, in the order the tool
made them. Each has a `name` — non-empty, ASCII, unique — and any other keys:
by convention `tool`, `version`, `parameters` and `date`.

**`table/`**, one row per cluster:

| Dataset | dtype | Description |
|---------|-------|-------------|
| `round` | `int32` | Zero-based index into `@rounds`. Rows come round by round, in the order of `@rounds`. |
| `name` | `S<n>` (fixed bytes) | Cluster name, non-empty, ASCII, unique within its round. |
| `parent` | `int32` | Row of `table/` of the cluster it was split from, in an earlier round; `-1` for none. |
| `n_contigs` | `int32` | Number of member contigs. |
| `length` | `int64` | Total length of its member contigs, in bp. |
| *extra columns* | any | Optional per-cluster columns (completeness, contamination, …), named anything but `genome` (below), `parent_round` (the round of `parent`, which `MetaCooler.clusters` computes) and `cluster` (what `dump` calls `name`). |

Not stored: **`genome`**, the genome of the map whose contigs a cluster holds
exactly, which `MetaCooler.clusters` and `dump -t clusters` add after
`length`. It is computed from `/contigs/genome_id` whenever a map is read, so
it always says which clusters the map's current genomes are, whatever made
them, and follows a `reassign`; `__unbinned__` is no genome. A stored
`table/genome` is refused.

**`members/`**, one row per membership, sorted by `(cluster_id, contig_id)`,
each pair at most once:

| Dataset | dtype | Description |
|---------|-------|-------------|
| `cluster_id` | `int32` | Row of `table/`. |
| `contig_id` | `int32` | Row of `/contigs`. |

**`indexes/`** looks the memberships up both ways, each in one slice:

| Dataset | dtype | Description |
|---------|-------|-------------|
| `cluster_offset` | `int64`, `n_clusters + 1` | The members of cluster `i` are the rows `cluster_offset[i]` to `cluster_offset[i + 1]` of `members/`. |
| `contig_offset` | `int64`, `n_contigs + 1` | The clusters of contig `j` are the rows `contig_offset[j]` to `contig_offset[j + 1]` of `contig_clusters`. |
| `contig_clusters` | `int32`, one per member | `members/cluster_id`, ordered by `contig_id`, then `cluster_id`: the members transposed. Clusters are numbered round by round, so a contig's clusters come in the order of the rounds. |

Both are derived from `members/`, and `validate()` checks them against it.

The group holds one run: writing clusters (`metacooler.write_clusters`)
replaces it whole, and earlier runs are kept only as the provenance their tool
leaves in the metadata. The commands that write a map from a map keep it:
`annotate` and the `.mmetacool` entries as stored; `reassign` too, with the
contig ids renumbered when it re-groups the contigs; `binnify` by contig and
`zoomify`, whose contigs are those of their source (`binnify --level genome`
makes other contigs, and warns that it drops the group); `merge`, which keeps the first
input's and warns about an input whose clusters differ. `create` and `collect`
write maps without clusters, and `extract`, `aggregate` and `flatten` write
plain cooler files. metacooler 0.3.0 and earlier read a map with this group.
Of the maps they write from it, those they copy keep it as stored —
`annotate`, a relabelling `reassign`, the `.mmetacool` entries, and the base
level of `zoomify`, a copy of the map — and those they rewrite drop it: a
re-grouping `reassign`, `binnify`, `merge`, and the coarser levels of
`zoomify`.

### Storage filters

Every dataset is chunked and compressed with gzip (level 6) behind the shuffle
filter: cooler's defaults, which it applies to `/bins` and `/pixels`, and which
metacooler applies to the datasets it writes itself (`/genomes`, `/contigs`,
`/bins/chrom`, `indexes/chrom_offset`, `indexes/genome_offset`). HDF5 filters
are transparent to readers, so this is not part of the format's contract: a
file whose metacooler tables are stored uncompressed (as every file written
before this change is) reads and validates exactly the same.

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
| `cool` | Preserve the exact contig order of the input cooler. `create cool`, and `collect` — where every input is one genome, genomes are laid out in `genome_id` order, and each keeps its input's contig order, so the layout is genome-grouped even though `cool` is not declared as such (until a `reassign` scatters a genome — see [`reassign`](cli.md#reassign-layout)). |
| `custom` | Row order inherited from an earlier layout and no longer described by any strategy. Written by `reassign --layout keep` when the new assignment breaks genome grouping; not a value `create` can produce. |

---

## The `metadata` attribute { #metadata }

The root `metadata` attribute is a JSON object. Every key in it is user
metadata — from `--metadata` files, `metadata add`/`load` — except four that
metacooler writes itself. `--metadata` and `metadata add`/`load`/`remove` refuse
to set or remove those four; `metadata list`/`get` read them like any other key.

| Key | Written by | Content |
|-----|------------|---------|
| `commands` | every command that writes a map | List of the command lines that produced the file, oldest first. Carried forward by every derived map. |
| `qc` | `create`, `collect`, `qc --store`, `reassign`, `merge` | The contact-quality block ([below](#qc-block)). |
| `reassign` | `reassign` | List, one entry per run, oldest first: `source` (the map reassigned), `table` (the assignment table; empty when the API was given a mapping), `layout` (the `contig-ordering` written) and `n_contigs_reassigned`. |
| `annotate` | `annotate` | List, one entry per run that changed the map, oldest first: `genomes_table` and `contigs_table` (the tables given, or `null`), `replace` and `missing` (the options used), and `added`, `replaced`, `dropped` — the columns affected, each as `<table>.<column>` (e.g. `contigs.gc`). A run that changes nothing records nothing. |

## Contact-quality block { #qc-block }

`create` records where the library's contacts landed, and `metacooler qc`
recomputes it on demand. It occupies three places in the file, all optional —
a file written with `--no-qc` has none of them.

**`metadata["qc"]`**, a JSON object inside the root `metadata` attribute,
alongside the `commands` provenance list:

| Key | Meaning |
|-----|---------|
| `version` | Block schema version (currently `2`; version 1 counted a contact between two different `__unbinned__` contigs as intra-genome rather than inter-genome). |
| `source` | `"pairs"` — measured on read positions by `create pairs`, before binning — or `"pixels"`, measured between bin midpoints. |
| `distance_threshold` | Short/long intra-contig split, in bp (default 10 000). |
| `n_contacts` | Total contacts classified. |
| `intra_contig_short`, `intra_contig_long`, `inter_contig_intra_genome`, `inter_genome` | The four mutually exclusive categories; they sum to `n_contacts`. |
| `intra_contig`, `inter_contig`, `intra_genome` | Convenience sums of the above. |
| `trans_contig_ratio` | The "3D ratio": `inter_contig / n_contacts`. |
| `*_fraction` | Each category's share of `n_contacts`. `null` on an empty file. |
| `binsize` | Bin size at measurement time (pixel source). |
| `n_pairs_read`, `n_pairs_discarded` | Pairs seen and dropped while binning (pairs source), so the denominator the ratio does *not* use stays visible. Summed over every pairs file when `create pairs` reads several, and over the inputs of `merge`. |

**`/genomes` columns**, one row per genome:
`qc_intra_contig_short`, `qc_intra_contig_long`,
`qc_inter_contig_intra_genome`, `qc_inter_genome` (`int64`), and
`qc_intra_genome_ratio` (`float64`, `NaN` for a genome with no contact). An
inter-genome contact is credited to both genomes involved, so these columns sum
to more than the file-level totals.

**`/contigs` columns**, one row per contig: `qc_intra_contig_short` and
`qc_intra_contig_long` (`int64`; `float64` on a map with float counts), the
same split contig by contig. A contig's contacts with itself are credited to
its genome, so each `/genomes` column of that name is the sum of this one over
the genome's contigs, and the column totals are the block's
`intra_contig_short` and `intra_contig_long`. Written with the block by every
command that stores one (`create`, `collect`, `qc --store`, `reassign`, and
`merge` when every input has both, as numbers adding up to its block); maps
written by metacooler 0.2.0 and earlier do not have them. These two are the
only `qc_*` columns metacooler writes to `/contigs`: any other (`qc_pass`,
from a `--contigs` table) is user metadata — every writer keeps it, `annotate`
can drop it, and `merge` compares it like any metadata column.

These two categories are the only ones that do not depend on the
contig-to-genome assignment, which is why they are kept per contig:
[`reassign`](cli.md#reassign) carries them over, sums them per new genome and
re-measures only the two inter-contig categories from the pixels — exactly, at
any bin size — so a `source="pairs"` block keeps its split measured on read
positions, and every other key of the block (`source`, `distance_threshold`,
`binsize`, `n_pairs_read`, `n_pairs_discarded`) is kept as it was. A source
without these columns (or whose columns no longer add up to its block) has
its block recomputed from the pixels instead, and loses the split measured on
read positions.

All three are statements about one specific pixel table, so `binnify` drops
them from its output instead of carrying them forward. `merge` sums them
instead: a merged map's contacts are its inputs' put together, every count of
the block, of `/genomes` and of `/contigs` adds up over libraries, and every
fraction and ratio is recomputed from the sums. The result is what a single
tally of every input's contacts gives — `create pairs` given every library's
pairs at once or, for blocks measured on pixels, `qc --store` on the merged
map. That takes blocks counted alike: the same `version`, `source` and
`distance_threshold`, and the same keys, `binsize` aside (kept when every input
records the same value), with each block holding together with its columns.
Otherwise `merge` drops them, with a warning that says why (see
[`merge`](cli.md#merge-qc)). It ignores `qc_*` columns when it compares the
inputs' `/genomes` and `/contigs` tables for compatibility — two libraries over
one assembly legitimately differ there. The names are reserved: a user table
column called `qc_intra_contig_short`, say, is stored as
`user_qc_intra_contig_short`.

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
   `/contigs`. metacooler writes it *relative* (`contigs`), so it keeps
   resolving when a map is copied into an `.mmetacool` entry, and `validate()`
   requires it to resolve to the map's own `/contigs`.
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
