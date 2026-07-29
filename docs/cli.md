# CLI reference

The `metacooler` command-line tool is the primary way to build, inspect, and
convert `.metacool` files. This page documents every command and option; run any
command with `-h/--help` for the same information in your terminal.

```console
metacooler --help
metacooler <command> --help
```

## Conventions

A few conventions hold across the whole CLI:

| Convention | Meaning |
|------------|---------|
| `METACOOL_PATH` | A **positional** path to a `.metacool` file, or an `.mmetacool` entry URI (`file.mmetacool::entries/<name>`). Every command takes its input this way — there are no input *options*. |
| `-o, --output` | Output path. Required unless the command writes to stdout. |
| `-f, --force` | Overwrite the output if it already exists. Available on **every** command that writes a file, including `dump`, `network`, and `show`. |
| `-v, --verbose` | Enable debug logging. Available on every command. |
| `-q, --quiet` | Suppress informational logging; warnings and errors are still printed. Available on every command. `-v` and `-q` are mutually exclusive in effect (`-v` wins). |
| `--genomes`, `--contigs` | Multi-value filters: repeat the flag or pass a comma-separated list (`--genomes MAG_01,MAG_02`). An explicitly empty value (`--genomes ''`) is **rejected** — it used to silently disable the filter, which turned an unset shell variable into a full dump. |

!!! warning "A bare `.mmetacool` path is rejected"
    A `.mmetacool` file is a *collection*, not a single contact map. Every
    command that operates on one map (`dump`, `extract`, `info`, `binnify`,
    `show`, `metadata …`, …) requires an entry URI and refuses a bare collection
    path, listing the available entries in the error message:

    ```console
    $ metacooler dump reps.mmetacool -t genomes
    Error: Invalid value for 'METACOOL_PATH': 'reps.mmetacool' is an mmetacool
    collection; specify an entry, e.g. reps.mmetacool::entries/rep1
    ```

    Only `ls`, `tree`, and `attrs` — which describe the container itself —
    accept a bare `.mmetacool` path.

New to the tool? Start with the [Quickstart](quickstart.md); for the ideas
behind genomes, contigs, bins, and levels, see [Concepts](concepts.md).

---

## Create and edit

### `create pairs` { #create-pairs }

Create a `.metacool` file from a 4DN pairs file, a set of bins, and a
contig-to-genome assignment table. **API:** `metacooler.from_pairs`
(`metacooler.io.create.from_pairs`).

```console
metacooler create pairs --pairs P --bins B --contigs C -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--pairs` | *required* | Pairs file (plain text, `.gz`, `.bgz`, or `.zst`). |
| `--bins` | *required* | BED file of bins (no header), **or** an integer bin size in bp (e.g. `1000`) — one or the other. For an integer, contig lengths must come from the pairs `#chromsize:` header or `--fasta`. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | *required* | Contig-to-genome table; layout set by `--contigs-format`. |
| `--contigs-format` | `contigs` | `contigs`, `contig_genome`, `metator`, or a column-index spec `I,J[,noheader]`. See [Creating from pairs](guides/from-pairs.md#contig-to-genome-table). |
| `--genomes` | — | Optional genome metadata TSV (first column `genome_name`); extra columns land in `/genomes`. It only *decorates* the genome rows derived from `--contigs` — it does not decide which genomes exist. A mismatch between the two tables is reported as a warning. See [Genomes metadata table](guides/from-pairs.md#genomes-metadata-table-optional). |
| `--metadata` | — | JSON or YAML file of free-form metadata to embed. May **not** contain a `commands` key — that name is reserved for the provenance record metacooler writes itself. |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `fasta`, or `alpha`. See [Contig ordering](concepts.md#contig-ordering). |
| `--fasta` | — | Assembly FASTA, used for contig lengths when `--bins` is an integer and the pairs header has no `#chromsize:` lines. |
| `--assembly-name` | — | Name for this assembly / metagenome. Stored in the standard cooler `genome-assembly` attribute (reported as `assembly` by `info`). |
| `--drop-unassigned` | off | Drop contigs absent from the contig table instead of collecting them under `__unbinned__`. |
| `--chunksize` | `50000000` | Pairs processed per chunk. This is chunk *granularity*, not read memory (the pyarrow read block is capped separately). Under `external` sorting each chunk becomes a temporary cooler at a fixed per-call cost that scales with contig count, so smaller chunks cost disproportionately more. |
| `--sort-strategy` | `auto` | Pixel sort strategy: `auto` (`memory` unless the binned pixels exceed 100 M rows), `memory` (sort in RAM — fast for small files), or `external` (stream to cooler's external sort — scalable). |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

!!! warning "`--contig-order fasta` needs a source of order"
    `fasta` means *"preserve the order of the input assembly"*, so it now
    **requires** either `--fasta` or a pairs file carrying `#chromsize:` header
    lines. Without one of those there is no assembly order to preserve;
    metacooler used to fall back silently to alphabetical order while still
    recording `contig-ordering="fasta"` in the file, so two files with different
    layouts claimed the same ordering and `merge` would scramble them.

```bash
# Fragment-level bins from MetaTOR output
metacooler create pairs \
    --pairs contacts.pairs.gz \
    --bins  fragments.bed \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_table.tsv \
    --assembly-name gut_T0 \
    -o gut_T0.metacool

# Fixed 2 kb bins, contig lengths from FASTA
metacooler create pairs \
    --pairs contacts.pairs.gz \
    --bins  2000 --fasta assembly.fa \
    --contigs contigs2mag.tsv --contigs-format contig_genome \
    -o gut_T0_2kb.metacool
```

### `create cool` { #create-cool }

Create a `.metacool` file from an existing cooler file by adding the genome
hierarchy on top of its bins and pixels. **API:** `metacooler.from_cool`
(`metacooler.io.create.from_cool`).

```console
metacooler create cool --cool IN.cool --contigs C -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--cool` | *required* | Input `.cool` file to convert. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | *required* | Contig-to-genome table; layout set by `--contigs-format`. |
| `--contigs-format` | `contigs` | Same choices as `create pairs`. |
| `--genomes` | — | Optional genome metadata TSV. Decorates the genome rows derived from `--contigs`; same semantics and warnings as `create pairs`. |
| `--metadata` | — | JSON or YAML metadata file (no `commands` key — see `create pairs`). |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `cool` (preserve the input's order), or `alpha`. |
| `--assembly-name` | — | Assembly / metagenome name, stored as `genome-assembly`. |
| `--drop-unassigned` | off | Drop unassigned contigs. |
| `--nproc` | `1` | Parallel reader threads for pixel remapping. |
| `--sort-strategy` | `auto` | As `create pairs`, but `nnz` is known up front so the 100 M-pixel threshold is compared directly. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

- There is no `--fasta` option: contig lengths always come from the cool file's
  own bins table.
- `--contig-order cool` is the only ordering that round-trips row order exactly;
  the default regroups contigs by genome.

```bash
metacooler create cool \
    --cool existing.cool \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_quality.tsv \
    -o sample.metacool
```

### `merge` { #merge }

Sum the pixel counts of two or more metacool files that share **identical**
bins, contigs, genomes, and contig ordering (e.g. sequencing replicates of the
same assembly). Counts are summed in `int64`. **API:**
`metacooler.merge` (`metacooler.io.merge.merge`).

```console
metacooler merge A.metacool B.metacool [...] -o merged.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.metacool` path. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

```bash
metacooler merge rep1.metacool rep2.metacool -o merged.metacool
```

### `binnify` { #binnify }

Re-bin a metacool file to a **coarser** fixed bin size, writing a new metacool
file. **API:** `metacooler.binnify` (`metacooler.ops.binnify.binnify`).

```console
metacooler binnify METACOOL_PATH --binsize N -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.metacool` path. |
| `--binsize` | *required* | New bin size in bp. See the constraints below. |
| `--level` | `contig` | Tiling strategy: `contig` (tile each contig independently) or `genome` (ligate a genome's contigs by decreasing length, then tile — the output's contigs become the ligated genomes). |
| `--chunksize` | `1000000` | Source pixels read, remapped, and aggregated per chunk (peak-memory control). |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

!!! warning "`--binsize` must be an exact multiple of the source bin size"
    For a **fixed-bin** source, the target must be a strict, *exact multiple* of
    the current bin size. A non-multiple target is rejected:

    ```console
    $ metacooler binnify gut_T0_1kb.metacool --binsize 3500 -o out.metacool
    Error: New binsize (3500) must be an exact multiple of the source binsize
    (1000) — 3500 % 1000 = 500. […] Nearest valid multiples: 3000 and 4000.
    ```

    This is not pedantry. Under the midpoint rule a non-multiple target packs an
    *oscillating* number of source bins into each output bin — 1 kb → 3.5 kb
    alternates 3 and 4 source bins, a built-in ±14 % checkerboard — while the
    file advertises a uniform 3500 bp bin size. Earlier versions accepted this
    and produced a silently wrong file. From 1 kb the valid targets are 2000,
    3000, 4000, 5000, …

    For a **variable-bin** source (fragment-level or contig-level files) there is
    no exact multiple to check, so the *median* source bin width is used as the
    effective resolution: the target must be strictly larger than it. A target
    between the median and the widest source bin is accepted with a warning —
    those widest bins each collapse into a single output bin, so the effective
    resolution is uneven across the file.

- At `--level contig`, each contig is tiled independently, the last bin of each
  contig may be shorter, and the contig hierarchy is preserved.
- At `--level genome`, each genome's contigs are ligated (longest first) into one
  artificial sequence which is then tiled, so a bin may span original contig
  boundaries. The output's contigs *are* the genomes (one per genome). Genomes
  with no contigs are skipped.
- Pixels are re-assigned and re-aggregated for the new layout; integer counts are
  widened to `int64` so a coarse bin pair cannot overflow.

```bash
# 1 kb → 5 kb (5000 = 5 × 1000 ✓)
metacooler binnify gut_T0.metacool --binsize 5000 -o gut_T0_5kb.metacool

# genome-level: ligate each MAG's contigs, then tile at 10 kb
metacooler binnify gut_T0.metacool --binsize 10000 --level genome -o gut_T0_genome10kb.metacool

# smaller chunks on a memory-constrained machine
metacooler binnify gut_T0.metacool --binsize 5000 --chunksize 200000 -o gut_T0_5kb.metacool
```

### `zoomify` { #zoomify }

Build a multi-resolution `.mmetacool` collection from a fixed-bin metacool file.
By default each level doubles the bin size; `--resolutions` gives cooler-style
control. **API:** `metacooler.zoomify` (`metacooler.ops.zoomify.zoomify`).

```console
metacooler zoomify METACOOL_PATH -o OUT.mmetacool [--levels N | --resolutions R]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.mmetacool` path. |
| `--levels` | `5` | Number of bin-size doublings (used when `--resolutions` is not set). |
| `-r, --resolutions` | — | Resolution ladder in cooler's ladder-token grammar (see below). Takes precedence over `--levels`. |
| `--temp-dir` | system temp | Directory for the intermediate per-level files. Mirrors cooler's `zoomify --temp-dir`; point it at fast local scratch with room for **one** level. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

Entries are addressed as `file.mmetacool::entries/<binsize>`.

!!! note "`--resolutions` is a ladder grammar, not a plain list"
    Unlike `--genomes`/`--contigs`, which take a plain comma-separated list,
    `-r/--resolutions` takes cooler's own **ladder-token** grammar. This
    difference is deliberate — it keeps the syntax identical to `cooler zoomify`.
    Accepted forms:

    | Form | Meaning |
    |------|---------|
    | `1000,2000,5000` | An explicit list of resolutions. |
    | `1000N` | "Nice" 1-2-5 progression starting at 1000. |
    | `1000B` | Binary ×2 progression starting at 1000. |
    | `N` / `B` | The same progressions, starting from the base bin size. |
    | `4DN` | The 4DN standard ladder (`1000,2000,5000N`). |

    Tokens are comma-separated and case-insensitive, and may be mixed
    (`1000,2000,5000N`). Anything else is a `ValueError`.

!!! warning "Every level must nest in the base bin size"
    A requested resolution that is not an exact multiple of the base bin size is
    rejected — with a 1 kb base, `--resolutions 1500,3000` fails, because a
    1500 bp level does not nest in the 1 kb binning and the pyramid would not be
    self-consistent:

    ```console
    Error: […] not multiples of the base bin size (1000); such levels would not
    nest in the base binning. Nearest valid resolutions: 1500 -> 2000.
    ```

    Watch out for `1000N`-style tokens against a non-1 kb base: with a 4 kb base
    the nice ladder 5000/10000/20000/50000 contains no valid level at all.

The base resolution is always included, resolutions finer than the base are
ignored, and generation **stops early**: levels that would leave fewer than 2
bins, or that would duplicate an already-written level (once the bin size exceeds
the longest contig, every contig collapses to one bin), are skipped and reported
in the log rather than computed and discarded.

```bash
metacooler zoomify gut_T0.metacool --levels 6 -o gut_T0_multireso.mmetacool

# explicit cooler-style resolutions (all multiples of the 1 kb base)
metacooler zoomify gut_T0.metacool --resolutions 1000,2000,5000,10000 -o gut_T0_multireso.mmetacool

# keep the per-level scratch on a fast local disk
metacooler zoomify gut_T0.metacool --levels 6 --temp-dir /scratch/$USER -o gut_T0_multireso.mmetacool
```

---

## Extract and view

### `extract` { #extract }

Write a standard `.cool` file containing only the selected genomes and/or
contigs, at the original bin resolution. In single-file mode, genome and contig
selections combine. **API:** `metacooler.extract`
(`metacooler.io.extract.extract`).

```console
metacooler extract METACOOL_PATH --genomes G -o OUT.cool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.cool` path (or filename prefix with `--split`). |
| `--genomes` | — | Genome name(s) to extract. |
| `--contigs` | — | Contig name(s) to extract (single-file mode only). |
| `--split` | off | Write one `.cool` per genome (genome-only); `--output` becomes a prefix. Cannot be combined with `--contigs`. |
| `-f, --force` | off | Overwrite existing output(s), including per-genome files with `--split`. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

At least one of `--genomes`/`--contigs` is required. The output is a plain
cooler file, readable by cooltools, HiGlass, and anything else that reads
`.cool`.

```bash
# Single MAG
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool

# Multiple MAGs, one file each
metacooler extract gut_T0.metacool --genomes MAG_01,MAG_02 --split -o split
```

### `flatten` { #flatten }

Drop the genome hierarchy and write a standard `.cool` file in which every
contig is a chromosome. Bins and pixels are preserved exactly. **API:**
`metacooler.flatten` (`metacooler.io.flatten.flatten`).

```console
metacooler flatten METACOOL_PATH -o flat.cool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.cool` path. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

### `dump` { #dump }

Dump any table to TSV (stdout by default). **API:** `metacooler.dump_table`
(`metacooler.io.dump.dump_table`).

```console
metacooler dump METACOOL_PATH -t TABLE [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-t, --table` | `pixels` | `bins`, `chroms`, `contigs`, `genomes`, `pixels`, `indexes`, or `all`. |
| `--genomes` | — | Filter rows to the named genome(s). |
| `--contigs` | — | Filter rows to the named contig(s). |
| `-H, --header` | off | Include a header row. |
| `-o, --output` | stdout | Output file. |
| `--prefix` | — | Filename prefix; **required** with `-t all`, and rejected with any other table. |
| `--pixels-format` | `long` | `long` (row per pixel) or `wide` (2-D matrix). |
| `--weight` | `count` | Value column for `wide` format (pixels only). |
| `--join` | off | Annotate pixels with contig names + coordinates (pixels only). |
| `-r, --range` | — | Row range (pixels only): genome name, contig name, or `contig:start-end`. |
| `-r2, --range2` | — | Column range (pixels only); defaults to `--range`. |
| `-f, --force` | off | Overwrite `-o/--output` if it already exists (`dump` used to clobber it silently). |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

`-t contigs` adds a `genome` column next to the internal `genome_id`, and
`-t genomes`/`-t chroms` respect `--genomes`/`--contigs`. `-t all` writes
`<prefix>.genomes.tsv`, `.contigs.tsv`, `.bins.tsv`, `.pixels.tsv`, and
`.indexes.tsv`; `chroms` is omitted because it duplicates `contigs`.

!!! note "Incompatible combinations are now errors"
    Options that do not apply to the chosen table used to be accepted and
    silently ignored (exit 0, output quietly unfiltered). They are now rejected
    with a usage error:

    | Combination | Why it fails |
    |-------------|--------------|
    | `--weight` / `--join` with `-t bins/chroms/contigs/genomes/indexes` | Only the pixels table has a value column or coordinates to join. |
    | `--genomes` / `--contigs` with `-t indexes` | The indexes table has no genome or contig dimension. |
    | `-o/--output` with `-t all` | `all` writes one file per table — use `--prefix`. |
    | `--prefix` with any table other than `all` | Single-table dumps go to `-o/--output`. |
    | `--genomes ''` / `--contigs ''` | An explicitly empty filter is rejected rather than treated as "no filter". |

```bash
# Genomes table with header
metacooler dump gut_T0.metacool -t genomes -H

# Contigs of one MAG (the `genome` column holds the MAG name)
metacooler dump gut_T0.metacool -t contigs --genomes MAG_01 -H

# Pixels within a genomic range
metacooler dump gut_T0.metacool -t pixels --range contig_0_0:0-50000

# Pixels annotated with contig names + coordinates (cooler --join style)
# → contig1  start1  end1  contig2  start2  end2  count
metacooler dump gut_T0.metacool -t pixels --join -H

# Overwrite an existing dump
metacooler dump gut_T0.metacool -t pixels -o pixels.tsv --force

# Every table at once → gut_T0.{genomes,contigs,bins,pixels,indexes}.tsv
metacooler dump gut_T0.metacool -t all --prefix gut_T0
```

### `network` { #network }

Export a contact network as an edge list (TSV: `node1`, `node2`, `value`).
**API:** `metacooler.network` (`metacooler.io.network.network`).

```console
metacooler network METACOOL_PATH --level LEVEL [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--level` | `contig` | Aggregation level: `bin`, `contig`, or `genome` (singular — not `bins`/`contigs`/`genomes`). |
| `--weight` | `count` | Pixel column used as edge weight. |
| `-o, --output` | stdout | Output TSV file. |
| `--header` | off | Write a `node1  node2  value` header line. |
| `--chunksize` | `1000000` | Pixels processed per chunk (memory control). |
| `-f, --force` | off | Overwrite `-o/--output` if it already exists (`network` used to clobber it silently). |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

```bash
metacooler network gut_T0.metacool --level genome -o genome_network.tsv
```

### `aggregate` { #aggregate }

Aggregate contacts to contig or genome level and write a standard `.cool` file
(one chromosome per contig or per genome). **API:** `metacooler.aggregate`
(`metacooler.ops.aggregate.aggregate`).

```console
metacooler aggregate METACOOL_PATH --level LEVEL -o OUT.cool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.cool` path. |
| `--level` | `contig` | `contig` (one chromosome per contig) or `genome` (sum contacts per genome pair). |
| `--count-col` | `count` | Pixel value column to aggregate. |
| `--chunksize` | `1000000` | Pixels read per chunk (peak-memory control). |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

```bash
metacooler aggregate gut_T0.metacool --level genome -o genome_level.cool
```

### `show` { #show }

Render a contact-map heatmap to an image file. **API:**
`metacooler.ops.show.show`.

```console
metacooler show METACOOL_PATH -o heatmap.png [OPTIONS]
```

!!! warning "`show` changed its input and `-f` in the 0.2 CLI pass"
    `show` used to be the odd one out: it read its input from a required
    `-f/--file` **option**, so `-f` meant something different there than in
    every other command. It now takes `METACOOL_PATH` as a **positional
    argument** like every other command, and `-f` is `--force`. Old invocations
    such as `metacooler show -f gut_T0.metacool -o out.png` no longer do what
    they used to — drop the `-f`.

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output image path (`.png`, `.pdf`, `.svg`). |
| `--level` | `genome` | `bin` (finest), `contig`, or `genome` (coarsest). |
| `--genome` | — | Restrict to one genome (intra-genome view). Honoured at every level. |
| `--genomes` | — | A pair of genomes (inter-genome view). Takes exactly **two space-separated** names (`--genomes MAG_01 MAG_02`) — not the comma list used elsewhere. Honoured at every level. |
| `--transform` | `log1p` | Value transform: `log1p`, `log10`, or `none`. |
| `--cmap` | `afmhot_r` | Matplotlib colormap name. |
| `--vmax` | — | Upper colormap bound (default: 99th percentile of non-zero values). |
| `--separators / --no-separators` | on | Draw lines at genome/contig boundaries. |
| `--dpi` | `500` | Output image resolution. |
| `-r, --range` | — | Row range (implies `--level bin`): genome, contig, or `contig:start-end`. |
| `-r2, --range2` | — | Column range; defaults to `--range`. |
| `-f, --force` | off | Overwrite the output image if it already exists (`show` used to clobber it silently). |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

A dense-matrix size guard applies at **every** level, including `--level bin`
and range queries: a selection too large to materialise is refused with an error
rather than attempting a multi-hundred-gigabyte allocation.

```bash
# Genome-level heatmap (default level) — note the positional path
metacooler show gut_T0.metacool --level genome -o genome_heatmap.png

# Intra-genome view at bin resolution
metacooler show gut_T0.metacool --level bin --genome MAG_01 -o MAG_01.png

# A specific genomic range, raw counts, custom colormap, overwriting an old image
metacooler show gut_T0.metacool -r contig_0_0:0-100000 \
    --transform none --cmap Reds -f -o range.png
```

---

## Inspect

### `info` { #info }

Print core file metadata and statistics (bin size, counts, ordering, …). Prints
JSON by default; `--table` renders a rich table. Also accepts standard
`.cool`/`.mcool` files and mmetacool entries. **API:** `MetaCooler.info`.

The `assembly` field reports the standard cooler `genome-assembly` attribute
(`"unknown"` when `--assembly-name` was not given at creation time); metacooler
no longer writes a separate `assembly-name` attribute.

```bash
metacooler info gut_T0.metacool           # JSON
metacooler info gut_T0.metacool --table   # rich table
metacooler info gut_T0.mmetacool::entries/rep1
```

### `tree` { #tree }

Print the HDF5 group/dataset hierarchy of a `.metacool` or `.mmetacool` file.

```bash
metacooler tree gut_T0.metacool
```

### `metadata` { #metadata }

Read and edit the free-form `metadata` attribute. Subcommands:

```bash
metacooler metadata list   gut_T0.metacool                       # list user keys
metacooler metadata get    gut_T0.metacool -k sample             # one value
metacooler metadata add    gut_T0.metacool -k enzyme --value DpnII
metacooler metadata remove gut_T0.metacool -k enzyme
metacooler metadata load   sample_meta.yaml gut_T0.metacool      # INPUT_FILE first
```

!!! note "`metadata load` argument order"
    `load` takes the metadata file **first**, then the target:
    `metacooler metadata load <metadata.yaml> <file.metacool>`.

!!! warning "Collection roots are rejected"
    **All five** subcommands require a single map: pass an entry URI
    (`reps.mmetacool::entries/rep1`), not the collection path. `add`, `get`,
    `remove`, and `load` used to open the `.mmetacool` root and write there —
    exiting 0 while the value stayed invisible from every entry.

The `commands` key holds metacooler's own provenance record and cannot be
overwritten, either by `metadata load` or by a `--metadata` file at creation
time. It is hidden from `metadata list` but readable with
`metadata get -k commands`.

### `attrs` { #attrs }

Print the HDF5 attributes of a group as YAML (same format as `cooler attrs`). The
`GROUP_PATH` argument is optional and defaults to the file root `/`. Accepts a
bare `.mmetacool` path.

```bash
metacooler attrs gut_T0.metacool           # root attributes
metacooler attrs gut_T0.metacool /genomes  # a specific group
```

Specify the group **once**: giving both a URI group and a `GROUP_PATH`
(`metacooler attrs reps.mmetacool::entries/rep1 /genomes`) is an error. The URI
group used to be silently discarded, so the lookup hit the file root.

---

## Manage collections

### `ls` { #ls }

List the entries in a `.mmetacool` (or `.mcool`, or plain `.metacool`) file.
**API:** `metacooler.list_entries`. The printed URIs can be fed straight back
into any other command.

```bash
metacooler ls gut_T0.mmetacool
# gut_T0.mmetacool::entries/rep1
# gut_T0.mmetacool::entries/rep2

metacooler ls gut_T0.metacool
# gut_T0.metacool::/
```

### `stack` { #stack }

Combine several metacool files into one `.mmetacool` collection. **API:**
`metacooler.create_mmetacool`.

```console
metacooler stack A.metacool B.metacool -o OUT.mmetacool [--names ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.mmetacool` path. |
| `--names` | see below | Comma-separated entry names (must match the number of inputs, and must be unique — duplicates are rejected up front rather than failing mid-write). |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

Inputs may be entry URIs as well as plain files; an entry URI copies **that
entry**, and its default name is the entry name (a plain path defaults to the
file stem, `sample.metacool` → `sample`).

```bash
metacooler stack rep1.metacool rep2.metacool --names "rep1,rep2" -o combined.mmetacool
metacooler stack a.mmetacool::entries/x b.mmetacool::entries/y -o picked.mmetacool
```

### `cp` { #cp }

Copy an entry within or between `.mmetacool` files. Both endpoints are entry
URIs — a bare file path is rejected. **API:** `metacooler.copy_entry` /
`metacooler.copy_entry_cross_file`.

```bash
metacooler cp combined.mmetacool::entries/rep1 combined.mmetacool::entries/rep1_backup
metacooler cp source.mmetacool::entries/rep1 target.mmetacool::entries/rep1
```

### `mv` { #mv }

Rename an entry within a single `.mmetacool` file (an O(1) HDF5 link rename, not
a copy-and-delete). **API:** `metacooler.move_entry`.

```bash
metacooler mv combined.mmetacool::entries/rep1_backup combined.mmetacool::entries/rep1_v2
```

### `rm` { #rm }

Remove an entry from a `.mmetacool` file (asks for confirmation unless `-f`).
**API:** `metacooler.remove_entry`.

```bash
metacooler rm combined.mmetacool::entries/rep1_v2
```

!!! note "Confirmation prompts in `cp`, `mv`, and `rm`"
    - Answering **`y`** to an overwrite prompt now actually overwrites the
      destination. It previously fell through and then failed; `-f` was the only
      way through.
    - **Declining** a prompt now exits with a **non-zero** status. It used to
      return 0, so `metacooler rm … && echo ok` printed `ok` after you said no.
