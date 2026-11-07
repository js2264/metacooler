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
| `--filter` | Select by *property* rather than by name — see [Filtering with `--filter`](#filter) below. Available on `create pairs`, `create cool`, `extract`, `dump`, and `show`. |

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

!!! warning "Contig and genome names must be ASCII"
    `create` refuses a non-ASCII name before it bins anything:

    ```console
    $ metacooler create pairs --pairs r.pairs --contigs mags.tsv -o out.metacool
    Error: Refusing to write a metacool file that cooler could not open — these
    names are not ASCII:
      /contigs.name: 'contig_é'
    ```

    `/contigs` doubles as cooler's `/chroms`, and cooler decodes every column
    of it as ASCII, so an accented contig name — or an accented value in any
    extra `--contigs` metadata column — makes the finished file unopenable by
    `cooler.Cooler()`, and therefore by `metacooler` itself. Genome names are
    held to the same rule because `aggregate --level genome` writes them into a
    `.cool` as chrom names. Extra `--genomes` metadata columns are exempt:
    cooler never reads that group, so accented taxonomy is fine there. See
    [Names must be ASCII](format.md#names-must-be-ascii).

New to the tool? Start with the [Quickstart](quickstart.md); for the ideas
behind genomes, contigs, bins, and levels, see [Concepts](concepts.md).

## Filtering with `--filter` { #filter }

`--genomes`/`--contigs` select by *name*. `--filter` selects by *property* —
"every genome over 1 Mb", "every contig over 10 kb", "every MAG above 90 %
completeness". One expression per flag; repeat the flag to require several
conditions at once.

```console
--filter "<table>.<column> <op> <value>"
```

| Part | Values |
|------|--------|
| `<table>` | `genomes` or `contigs`. **Always required** — both tables have a `length` and a `name`. |
| `<column>` | Any column of that table, including the extra ones your `--contigs`/`--genomes` metadata tables added (run `dump -t genomes -H` to see them). Names with spaces, dashes, or dots need no quoting: `genomes.Present Markers >= 50`. |
| `<op>` | `>` `>=` `<` `<=` `==` `!=`, plus `=~` / `!~` for a regular-expression match. A lone `=` means `==`. |
| `<value>` | A number — plain (`1000000`), scientific (`1e6`), underscored (`1_000_000`), or with a size suffix (`10kb`, `2.5Mb`, `1Gb`, case-insensitive) — or, for `==`/`!=`/`=~`, free text. |

Two extra columns exist beyond what `dump` prints: `contigs.genome` holds the
**genome name** of each contig (`/contigs` stores only the integer
`genome_id`), and it is available during `create` too, before any id exists.

```bash
# Genomes above 1 Mb, dumped as a table
metacooler dump gut_T0.metacool -t genomes -H --filter "genomes.length > 1Mb"

# Contigs of MAG_01 longer than 50 kb — names choose the candidates,
# the filter narrows them (this is an intersection, not a union)
metacooler dump gut_T0.metacool -t contigs -H \
    --genomes MAG_01 --filter "contigs.length > 50kb"

# Several conditions at once
metacooler extract gut_T0.metacool -o big_complete.cool \
    --filter "genomes.length > 1Mb" \
    --filter "genomes.Completeness > 90" \
    --filter "contigs.length >= 10kb"

# Drop the sub-5 kb tail at creation time, so it never costs a bin
metacooler create pairs --pairs contacts.pairs.gz --bins 2000 --fasta assembly.fa \
    --contigs contig_data_final.txt --contigs-format metator \
    --filter "contigs.length >= 5kb" -o gut_T0.metacool
```

**Semantics.**

- **Repeating the flag ANDs** the conditions. There is no `or` and no
  parenthesised grammar — that is the price of letting a column be called
  `Weighted completeness` without backticks. For an `or`, run the two filters
  separately.
- **Names choose, filters narrow.** `--genomes MAG_01 --filter
  "contigs.length > 50kb"` is *MAG_01's long contigs*, not "MAG_01 plus every
  long contig in the file". With no `--genomes`/`--contigs` at all, the filter
  selects from the whole file.
- **A `genomes.*` filter selects that genome's contigs**, so it narrows the
  bins, pixels, and contigs tables too — not just the `/genomes` rows.
- **Missing values never match.** A row whose value is `NA` fails every
  comparison, including the negated `!=` and `!~`: a filter selects rows it can
  positively vouch for.
- **Text columns holding numbers still work.** MetaTOR writes `-` for an absent
  shotgun coverage, so `contigs.Shotgun_coverage > 3` coerces what it can,
  excludes the rest, and logs how many it could not read. A column where
  *nothing* is a number is rejected outright, with a pointer to `==`/`=~`.
- **Nothing selected is an error, not an empty file.** The message reports what
  each filter would have kept on its own, so a mistyped threshold is
  distinguishable from a genuinely empty intersection:

    ```console
    $ metacooler dump gut_T0.metacool -t genomes --filter "genomes.length > 1Gb"
    Error: Nothing passed --filter. On its own, each filter keeps:
      genomes.length > 1Gb  ->  0 of 49 genomes   <-- keeps nothing
    ```

- **Unknown columns are caught, with a suggestion.** `genomes.lenght > 1Mb`
  reports `Did you mean 'length'?` and lists the available columns. A malformed
  expression is rejected before the file is even opened.

!!! note "`--filter` in `create` is not the same as `--filter` in `extract`"
    In `create`, a dropped contig is never tiled into bins and its pairs are
    discarded while binning (counted in the usual drop summary), so the sub-kb
    tail of an assembly costs nothing at all — including the per-chunk
    `cooler.create` overhead, which scales with contig count. The `/genomes`
    rows, lengths, and offsets are then derived from what survives, so a genome
    left with no contigs disappears entirely.

    One column is unavailable there: `genomes.n_bins` does not exist until the
    bins are tiled, and asking for it says so rather than reporting an unknown
    column. Filter on `genomes.length` or `genomes.n_contigs` instead.

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
| `--filter` | — | Keep only the genomes/contigs matching an expression, *before* anything is binned — e.g. `--filter "contigs.length >= 5kb"`. Repeatable. See [Filtering with `--filter`](#filter). `genomes.n_bins` is not available here. |
| `--qc / --no-qc` | on | Tally contact-quality statistics while the pairs are binned and store them in `metadata.qc` plus the `/genomes` `qc_*` columns. See [`qc`](#qc). Costs a few vectorised operations per chunk, and is the only way to get an *exact* short/long split — binning throws the sub-bin distances away. |
| `--qc-threshold` | `10000` | Distance separating short- from long-range intra-contig signal. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--ps / --no-ps` | on | Accumulate a **per-MAG P(s) curve** from the same pass and write it beside the output as `<output>.ps.tsv`. Measured on read positions, so it resolves distances the binned map cannot. No plot is drawn. See [`ps`](#ps). |
| `--ps-bins-per-decade` | `10` | Log-grid resolution of the P(s) table. |
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
| `--filter` | — | Keep only the genomes/contigs matching an expression, *before* the bins are remapped. Repeatable. Same semantics and caveats as `create pairs`. |
| `--qc / --no-qc` | on | Tally contact-quality statistics off the pixel stream being copied, and store them in `metadata.qc` plus the `/genomes` `qc_*` columns. A cool file has no read positions left, so intra-contig distances are measured between bin midpoints and the block records `source="pixels"`. See [`qc`](#qc). |
| `--qc-threshold` | `10000` | Distance separating short- from long-range intra-contig signal. Accepts `bp/kb/Mb/Gb` suffixes. |
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
| `--filter` | — | Keep only the genomes/contigs matching an expression. Repeatable; narrows `--genomes`/`--contigs`, or selects from the whole file on its own. See [Filtering with `--filter`](#filter). |
| `--split` | off | Write one `.cool` per genome (genome-only); `--output` becomes a prefix (a trailing `.cool` is dropped, nothing else). Cannot be combined with `--contigs`. |
| `-f, --force` | off | Overwrite existing output(s), including per-genome files with `--split`. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

At least one of `--genomes`/`--contigs`/`--filter` is required. The output is a
plain cooler file, readable by cooltools, HiGlass, and anything else that reads
`.cool`.

With `--split`, one file is written per genome that still has a contig, so
`--filter` alone is a complete selection there too.

```bash
# Single MAG
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool

# Multiple MAGs, one file each
metacooler extract gut_T0.metacool --genomes MAG_01,MAG_02 --split -o split

# Every MAG over 1 Mb, one file each — no need to name them
metacooler extract gut_T0.metacool --filter "genomes.length > 1Mb" --split -o big
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
| `--filter` | — | Filter rows by property, e.g. `--filter "contigs.length >= 10kb"`. Repeatable; narrows `--genomes`/`--contigs`. See [Filtering with `--filter`](#filter). |
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
    | `--genomes` / `--contigs` / `--filter` with `-t indexes` | The indexes table has no genome or contig dimension. |
    | `-o/--output` with `-t all` | `all` writes one file per table — use `--prefix`. |
    | `--prefix` with any table other than `all` | Single-table dumps go to `-o/--output`. |
    | `--genomes ''` / `--contigs ''` | An explicitly empty filter is rejected rather than treated as "no filter". |

!!! note "`dump … | head` is not an error"
    Piping a dump into a program that stops reading early (`head`, `less`, a
    `grep -m`) closes the pipe under `dump`. That used to surface as a red
    `Error: [Errno 32] Broken pipe` box after the rows had already been printed,
    which read as though the dump had failed. It is now silent, as in any other
    Unix tool: nothing is written to stderr, and the exit status is `141`
    (`128 + SIGPIPE`) — the same status `yes | head` gives. A plain pipeline
    still reports the *last* command's status, so `metacooler dump … | head`
    exits `0`; only `set -o pipefail` surfaces the `141`.

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

### `ps` { #ps }

Export — and optionally plot — **P(s)**, the probability of contact as a
function of genomic distance. **API:** `metacooler.contact_probability` /
`metacooler.plot_ps` (`metacooler.ops.ps`).

```console
metacooler ps METACOOL_PATH [-o table.tsv] [--plot fig.png] [OPTIONS]
```

Only contacts whose two ends sit on the **same contig** take part — those are
the ones a genomic distance is defined for. Contacts between contigs (of the
same MAG or not) are what [`qc`](#qc) counts instead.

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | stdout | TSV destination. |
| `--plot` | — | Also render the curve here (`.png`, `.pdf`, `.svg`). |
| `--pairs` | — | Measure from this pairs file instead of the binned pixel table (plain text, `.gz`, `.bgz` — read by the same streamer `create pairs` uses). See [Two sources](#ps-sources). |
| `--by` | `all` | `all` pools every selected contig into one curve; `genome` emits one curve per MAG. |
| `--genomes` | — | Restrict to these genomes, **comma-separated or repeated**. |
| `--contigs` | — | Restrict to these contigs, same syntax. |
| `--filter` | — | Narrow the selection by property. Repeatable. See [Filtering with `--filter`](#filter). |
| `--bins-per-decade` | `10` | Log-grid resolution. |
| `--min-distance` | 1 bp with `--pairs`, else the bin size | Upper edge of the leading `[0, min)` bin, which holds separation 0. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--max-distance` | longest span | Largest separation to report. Same suffixes. |
| `--value` | `prob` | Column plotted on the y axis: `prob` or `avg_contacts`. |
| `--header / --no-header` | header on | Write a header row in the TSV. |
| `--dpi` | `300` | Output image resolution. |
| `--chunksize` | `50000000` | Pixels (or pairs, under `--pairs`) processed per chunk. |
| `-f, --force` | off | Overwrite existing outputs. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

#### The table

| Column | Meaning |
|--------|---------|
| `genome` | The MAG name under `--by genome`; `all` when the curves are pooled. |
| `dist` | Geometric centre of the log-spaced distance bin, in bp — the x coordinate to plot against. |
| `n_contacts` | Contacts observed at that separation. |
| `n_bin_pairs` | Opportunities: pairs of bins (or of positions, under `--pairs`) that *could* have contributed one. |
| `avg_contacts` | `n_contacts / n_bin_pairs`. |
| `prob` | `avg_contacts` normalised so the curve integrates to 1 over the reported range. |

`n_bin_pairs` is what makes the curve mean anything: a contig of *n* bins
offers only `n - d` pairs at separation *d*, so a raw count of long-range
contacts is starved by construction and would look like a decay even on
uniform data. Dividing it out is also what makes `--by genome` curves
comparable between a 5 Mb MAG and a 300 kb one, and `prob` is what makes two
libraries of different depth comparable. Distance bins that nothing could
populate are omitted rather than reported as a rate of zero.

#### Two sources { #ps-sources }

| | distances | `n_bin_pairs` counts | starts at |
|---|---|---|---|
| **pixels** (default) | between bin midpoints | pairs of **bins**, `n - d` per contig | one bin |
| **pairs** (`--pairs`) | between read positions | pairs of **positions**, `L - s` per contig | 1 bp |

A binned map cannot say anything below its own resolution: every sub-bin
contact has already been collapsed onto the diagonal. Since the short-range
head of P(s) is exactly where a metaHi-C library shows undigested fragments and
self-religation, `create pairs` measures the curve while it still has the read
positions, and writes it as `<output>.ps.tsv` — the same table `--pairs`
recomputes, byte for byte. Reach for `--pairs` when you need that resolution on
a file created earlier, or with a different distance grid or selection.

```bash
# The fine-grained curve create pairs already wrote, per MAG
head gut_T0.ps.tsv

# Pooled curve to stdout, from the binned map
metacooler ps gut_T0.metacool

# One curve per MAG, table + figure
metacooler ps gut_T0.metacool --by genome -o ps_by_mag.tsv --plot ps_by_mag.png

# Only the large, well-covered MAGs, on a finer distance grid
metacooler ps gut_T0.metacool --by genome \
    --filter "genomes.length > 1Mb" --bins-per-decade 20 --plot ps_big_mags.pdf

# Raw mean contacts per opportunity instead of the normalised probability
metacooler ps gut_T0.metacool --value avg_contacts --plot ps_raw.png

# Base-pair resolution from the original pairs, for one MAG, plotted
metacooler ps gut_T0.metacool --pairs contacts.pairs.gz \
    --genomes MAG_01 --plot MAG_01_ps.png -o MAG_01_ps.tsv
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
| `--genomes` | — | Genome names, **comma-separated or repeated** (`--genomes MAG_01,MAG_02`). Every name goes on **both** axes — see below. Honoured at every level. |
| `--contigs` | — | Contig names, same syntax and same rule. Bare names only; for coordinates use `-r`. Upgrades `--level genome` to `contig`. |
| `--filter` | — | Choose the axes by property instead of by name, e.g. `--filter "genomes.length > 1Mb"`. Repeatable; narrows `--genomes`/`--contigs`. A `contigs.*` filter upgrades `--level genome` to `contig`, exactly as `--contigs` does. Cannot be combined with `-r/-r2`. See [Filtering with `--filter`](#filter). |
| `--transform` | `log1p` | Value transform: `log1p`, `log10`, or `none`. |
| `--cmap` | `afmhot_r` | Matplotlib colormap name. |
| `--vmax` | — | Upper colormap bound (default: 99th percentile of the non-zero values). |
| `--separators / --no-separators` | on | Draw lines at genome/contig boundaries. |
| `--dpi` | `500` | Output image resolution. |
| `-r, --range` | — | Row range (implies `--level bin`): genome, contig, or `contig:start-end`. |
| `-r2, --range2` | — | Column range; defaults to `--range`. |
| `-f, --force` | off | Overwrite the output image if it already exists (`show` used to clobber it silently). |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

#### Two ways to choose the axes

There are exactly two, and they cannot be combined:

| Selector | Shape | Meaning |
|----------|-------|---------|
| `--genomes` / `--contigs` (± `--filter`) | **symmetric** | Every name goes on **both** axes, so the figure holds all the intra- *and* inter-block contacts between them. Uniform in the number of names: one name is that genome against itself, N names are all N against all N. Blocks appear in the order the names were given — you control the layout — with a heavier boundary line between them. |
| `-r/--range` + `-r2/--range2` | **asymmetric** | A different region per axis, at bin level. This is the only way to get a rectangular image. |

`--genomes` and `--contigs` are mutually exclusive: both put their names on
*both* axes, so combining them would have no single meaning. `--contigs` takes
bare names — a coordinate range is `-r`'s job.

`--filter` belongs to the symmetric side: it narrows whichever list you gave
(or the whole file), keeps the file's order, and is rejected alongside
`-r/-r2` for the same reason `--genomes` is.

!!! note "There is no `--genome` (singular)"
    It was removed: `--genomes MAG_01` with one name is exactly the same thing,
    and the near-identical pair of flags meant `--genome MAG_01,MAG_02` looked
    for a genome literally named `MAG_01,MAG_02`. `show` was also the only
    command in the CLI with a singular form.

#### Titles and axes

Every figure is titled with what it shows and at which resolution —
`MAG_01 — contig level`, `contig_0_0:0-100000 — bin level`, `All genomes — genome
level`. A multi-name selection lists the names, truncating past four:
`5 genomes: MAG_01, MAG_02, MAG_03, MAG_04 (+1 more) — contig level`.

Axis labels name the selection and its size (`MAG_01 (12 contigs)`,
`5 genomes (19 contigs)`), and at `--level bin` the ticks are **base-pair
coordinates**, not bin indices:

- a selection inside one contig is ticked in bp at its true position, so
  `-r contig_0_0:20000-50000` reads `20 kb … 50 kb` rather than starting at 0;
- a selection spanning several contigs has no single coordinate system, so the
  ticks name the contigs instead, one at the first bin of each (lined up with
  the separators) — or the selected genomes, when there are several.

A dense-matrix size guard applies at **every** level, including `--level bin`
and range queries: a selection too large to materialise is refused with an error
rather than attempting a multi-hundred-gigabyte allocation.

```bash
# Genome-level heatmap (default level) — note the positional path
metacooler show gut_T0.metacool --level genome -o genome_heatmap.png

# Intra-genome view at bin resolution
metacooler show gut_T0.metacool --level bin --genomes MAG_01 -o MAG_01.png

# Two genomes: both on both axes, so the image holds both intra blocks and the
# inter block between them
metacooler show gut_T0.metacool --level contig --genomes MAG_01,MAG_02 -o pair.png

# Four genomes, all against all, blocked in this order (repeats accumulate)
metacooler show gut_T0.metacool --level contig \
    --genomes MAG_01,MAG_02 --genomes MAG_03,MAG_04 -o panel.png

# A handful of contigs of interest, against each other
metacooler show gut_T0.metacool --level bin \
    --contigs contig_0_0,contig_3_7,contig_11_2 -o contigs.png

# Rectangular: a different region per axis (the only asymmetric form)
metacooler show gut_T0.metacool -r MAG_01 -r2 MAG_02 -o inter.png

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

### `qc` { #qc }

Report where a library's contacts land: the **3D ratio** and its breakdown.
**API:** `metacooler.qc` (`metacooler.ops.qc.qc`).

```console
metacooler qc METACOOL_PATH [OPTIONS]
```

Every contact falls in exactly one of four categories:

| Category | Meaning |
|----------|---------|
| `intra_contig_short` | Same contig, less than `--threshold` apart. |
| `intra_contig_long` | Same contig, at least `--threshold` apart. |
| `inter_contig_intra_genome` | Two contigs of the same MAG. |
| `inter_genome` | Two contigs of different MAGs. |

The headline `trans_contig_ratio` — the "3D ratio" — is the share of contacts
joining two *different* contigs, i.e.
`(inter_contig_intra_genome + inter_genome) / n_contacts`. The block also
reports `intra_genome` / `inter_genome` totals and a fraction for every
category. On an empty file the ratios are `null` rather than `0`.

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | stdout | Write the JSON summary here. |
| `--table` | off | Print a rich table instead of JSON. |
| `--per-genome` | — | Write the per-genome breakdown as a TSV to this path. |
| `--threshold` | `10000` | Short/long intra-contig split. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--store` | off | Write the result back into the file: `metadata.qc` plus the `/genomes` `qc_*` columns. |
| `--chunksize` | `50000000` | Pixels processed per chunk. |
| `-f, --force` | off | Overwrite `-o`/`--per-genome` if they exist. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

#### Two sources, one vocabulary

The block records which one produced it:

- **`source="pairs"`** — written by `create pairs` from the raw read positions,
  before binning. Distances are exact, and the block additionally reports
  `n_pairs_read` / `n_pairs_discarded`, so the denominator the ratio does *not*
  use is visible.
- **`source="pixels"`** — what this command computes, and what `create cool`
  stores. Intra-contig distances are gaps between bin midpoints, so the
  short/long split is only as fine as the bin size; at `binsize >= threshold`
  everything but the diagonal is "long", and the command warns. Every other
  number is identical to the pairs-based one.

`--store` on a file created with `create pairs` therefore *replaces* a more
precise block with a coarser one, and says so.

#### Per-genome columns

`--store` (and `create --qc`) add five columns to `/genomes`, readable with
`metacooler dump --table genomes`:

`qc_intra_contig_short`, `qc_intra_contig_long`,
`qc_inter_contig_intra_genome`, `qc_inter_genome`, and
`qc_intra_genome_ratio` — the share of a MAG's contacts that stay inside it.
That last one is a bin-quality score: a MAG whose contigs really belong
together sits near 1, while a chimeric or over-split one leaks contacts to its
neighbours. An inter-genome contact is credited to *both* genomes, so the
per-genome columns sum to more than the file-level total.

Because these columns describe one particular pixel table, `binnify` and
`merge` drop them (and the `metadata.qc` block) from their output rather than
carry stale numbers forward — re-run `metacooler qc --store` on the result.

```bash
# The 3D ratio and its breakdown, as JSON
metacooler qc gut_T0.metacool

# A 5 kb short/long split, as a table, stored back into the file
metacooler qc gut_T0.metacool --threshold 5kb --table --store

# Per-MAG breakdown for downstream plotting
metacooler qc gut_T0.metacool --per-genome mag_qc.tsv

# What create pairs already recorded, without recomputing anything
metacooler metadata get gut_T0.metacool -k qc
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

### `check` { #check }

Check the integrity of a `.metacool` file (or of every entry in an
`.mmetacool`). **API:** `metacooler.check` (`metacooler.schema.check`).

```console
metacooler check METACOOL_PATH [--deep]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--deep` | off | Also scan the **whole** pixel table instead of the three sampled windows, and reconcile `attrs['sum']` with the stored counts. O(nnz) — slow on large files. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

This is the reporting front end to the same validators that run on every file
open ([`validate` / `validate_mmetacooler`](api/metacooler.md)): instead of raising
on the first problem it prints one line per check, so every failure is visible in
one pass. Accepts a bare `.mmetacool` path — the per-map checks then repeat once
per entry, tagged with the entry name.

| Check | What it establishes |
|-------|---------------------|
| `readable` | The path exists, opens as HDF5, and contains the group the URI names. |
| `format` | It carries `metacooler-format` / `HDF5::MMetaCooler`. A plain `.cool` is reported as such, with the `create cool` hint. |
| `structure` | `validate()`: required groups and columns, offset lengths, pixel sort order, `genome_id` range, and every declared count against the tables. |
| `dtypes` | The 16 core columns are stored in a usable dtype **kind**. `validate` casts on read, so a `bins/start` written as float64 or a `genome_id` written as a string passes it and fails much later inside cooler. |
| `cooler-compat` | `cooler.Cooler()` opens the file — the central metacooler invariant that a `.metacool` is always also a cooler. |
| `pixels` (`--deep`) | Every pixel, not a sample: in range, upper-triangular, sorted by `(bin1_id, bin2_id)` across chunk boundaries, in the `bin1_offset` slot that claims it, and summing to `attrs['sum']`. |

Exits `0` when everything passed, `1` otherwise, listing the checks that failed.

```bash
metacooler check gut_T0.metacool
# ✓ readable: gut_T0.metacool opens as HDF5
# ✓ format: HDF5::MetaCooler
# ✓ structure: required groups, columns, indexes and counts are consistent
# ✓ dtypes: 16 core columns are stored in usable types
# ✓ cooler-compat: cooler.Cooler() opens it: 44 chroms (= contigs), 661 bins, 36249 pixels
# gut_T0.metacool: OK (5 checks passed)

metacooler check --deep gut_T0.metacool        # + a full pixel scan
metacooler check gut_T0.mmetacool              # every entry in a collection
```

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
