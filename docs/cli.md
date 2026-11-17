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
| `-f, --force` | Overwrite the output if it already exists. Available on **every** command that writes a file, including `dump`, `network`, and `show`. An output that is one of the command's own inputs is refused even with `-f`, and `create`, `reassign` and `annotate` keep the previous file until the new one has been written and validated. |
| `-v, --verbose` | Enable debug logging. Available on every command. |
| `-q, --quiet` | Suppress informational logging; warnings and errors are still printed. Available on every command. `-v` and `-q` are mutually exclusive in effect (`-v` wins). |
| `--genomes`, `--contigs` | Multi-value filters: repeat the flag or pass a comma-separated list (`--genomes MAG_01,MAG_02`). An explicitly empty value (`--genomes ''`) is **rejected** — it used to silently disable the filter, which turned an unset shell variable into a full dump. |
| `--filter` | Select by *property* rather than by name — see [Filtering with `--filter`](#filter) below. Available on `create pairs`, `create cool`, `collect`, `extract`, `dump`, and `show`. |

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
    accept a bare `.mmetacool` path, and [`annotate`](#annotate), which then
    annotates every entry of the collection.

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

## Create { #create }

Build a map: from pairs files ([`create pairs`](#create-pairs)), from a cool
file ([`create cool`](#create-cool)), or from single-genome cool files, one
genome each ([`collect`](#collect)).

### `create pairs` { #create-pairs }

Create a `.metacool` file from one or more 4DN pairs files, a set of bins, and
a contig-to-genome assignment table. **API:** `metacooler.from_pairs`
(`metacooler.io.create.from_pairs`).

```console
metacooler create pairs --pairs P [--pairs P2 ...] --bins B [--contigs C] -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--pairs` | *required* | Pairs file (plain text, `.gz`, `.bgz`, `.bz2`, or `.zst`). **Repeatable**: give one `--pairs` per Hi-C library of the same assembly and they are read back to back, in the given order, in one pass — see [Several libraries](#create-pairs-libraries). Not comma-split, since a path may contain a comma. A pair with an end on an unknown contig, with an empty position, or at or past its contig's end is dropped and counted, and a summary is logged; `metadata.qc.n_pairs_discarded` records the total. |
| `--bins` | *required* | BED file of bins (no header), **or** an integer bin size in bp (e.g. `1000`) — one or the other. For an integer, contig lengths must come from the pairs `#chromsize:` header(s) or `--fasta`. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | — | Contig-to-genome table; which columns hold the contig and the genome names is set by `--contigs-format`, and its other columns land in `/contigs`. It need not list every contig, but must name at least one contig of the assembly — see [The contig table](#create-contig-table). Leave it out to build a map with no binning yet. |
| `--contigs-format` | `contigs` | `CONTIG,GENOME[,noheader]`, each a 1-based position or a header name, or an alias for a binner's table — see [Contig and genome columns](#contigs-format) below. |
| `--genomes` | — | Optional genome metadata TSV; extra columns land in `/genomes`. It only *decorates* the genome rows derived from `--contigs` — it does not decide which genomes exist. A mismatch between the two tables is reported as a warning. See [Genomes metadata table](guides/from-pairs.md#genomes-metadata-table-optional). |
| `--genomes-format` | `genome_name` | Which column of `--genomes` holds the genome names: a 1-based position or a header name, plus `,noheader` when the table has no header line — e.g. `Name` for CheckM2's `quality_report.tsv`, `user_genome` for GTDB-Tk's summary. |
| `--metadata` | — | JSON or YAML file of free-form metadata to embed. May **not** contain a key metacooler writes itself: `commands` (the provenance record), `qc` (the contact-quality block), `reassign` or `annotate` (those commands' histories). |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `fasta`, or `alpha`. See [Contig ordering](concepts.md#contig-ordering). |
| `--fasta` | — | Assembly FASTA, used for contig lengths when `--bins` is an integer and the pairs headers have no `#chromsize:` lines. |
| `--assembly-name` | — | Name for this assembly / metagenome. Stored in the standard cooler `genome-assembly` attribute (reported as `assembly` by `info`). |
| `--drop-unassigned` | off | Drop contigs absent from the contig table instead of collecting them under `__unbinned__`. |
| `--filter` | — | Keep only the genomes/contigs matching an expression, *before* anything is binned — e.g. `--filter "contigs.length >= 5kb"`. Repeatable. See [Filtering with `--filter`](#filter). `genomes.n_bins` is not available here. |
| `--qc / --no-qc` | on | Tally contact-quality statistics while the pairs are binned and store them in `metadata.qc` plus the `qc_*` columns of `/genomes` and `/contigs`. See [`qc`](#qc). Costs a few vectorised operations per chunk, and is the only way to get an *exact* short/long split — binning throws the sub-bin distances away. The per-contig split is what lets a later [`reassign`](#reassign) keep it. |
| `--qc-threshold` | `10000` | Distance separating short- from long-range intra-contig signal. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--ps / --no-ps` | on | Accumulate a **per-MAG P(s) curve** from the same pass and write it beside the output as `<output>.ps.tsv`. Measured on read positions, so it resolves distances the binned map cannot. No plot is drawn. See [`ps`](#ps). |
| `--ps-bins-per-decade` | `10` | Log-grid resolution of the P(s) table. |
| `--chunksize` | `50000000` | Pairs processed per chunk. This is chunk *granularity*, not read memory (the pyarrow read block is capped separately). Under `external` sorting each chunk becomes a temporary cooler at a fixed per-call cost that scales with contig count, so smaller chunks cost disproportionately more. |
| `--sort-strategy` | `auto` | Pixel sort strategy: `auto` (`memory` unless the binned pixels exceed 100 M rows), `memory` (sort in RAM — fast for small files), or `external` (stream to cooler's external sort — scalable). |
| `--nproc` | `1` | Threads binning pairs while the file is read and decompressed (which stays on one thread: one gzip stream). The output, QC block, P(s) table and drop report are identical whatever the value; each thread holds one more batch of pairs. Two threads take ~15 % off a 30 M-pair run; more add little. |
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
# No binning yet: every contig in __unbinned__, to be binned with reassign later
metacooler create pairs --pairs contacts.pairs.gz --bins 2000 --fasta assembly.fa \
    -o gut_T0.metacool

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

#### The contig table { #create-contig-table }

`--contigs` assigns contigs to genomes, and need not cover the whole assembly —
nor only this assembly:

- **A contig the table does not list** is collected under `__unbinned__`, or
  dropped with `--drop-unassigned`. So is one listed with an empty or NA-like
  genome name: `NA`, `N/A`, `nan`, `ND`, `n/d`, or `-1` (SemiBin2's label for
  an unbinned contig).
- **A row naming a contig the assembly does not have** is ignored, with a
  warning naming up to 10 such rows, so one table can serve several
  assemblies.
- **No table at all** — `--contigs` left out — builds a map with **no binning
  yet**: the cooler plus the genome layer, with every contig in a single
  `__unbinned__` genome, to be binned later with [`reassign`](#reassign).
  `--genomes` and `--drop-unassigned` are refused then — there is no genome to
  describe, and every contig would be dropped.
- **A table naming no contig of the assembly at all** — or with no row — is
  **refused**, before anything is binned. Every contig would land in
  `__unbinned__`, and that is never the map you meant: it is the wrong table,
  or the right one read with the wrong `--contigs-format` — typically a leading
  ID column read as `1,2`, which takes the IDs for contig names and the contig
  names for genome names.

```console
$ metacooler create pairs ... --contigs mags.tsv --contigs-format 1,2
Error: None of the 44 contig name(s) in the contig table mags.tsv matches a
contig in the bins file fragments.bed, so no contig would be assigned to a
genome.
  table names:      '1', '2', '3', ...
  assembly contigs: 'contig_1_000', 'contig_1_001', 'contig_1_002', ...
Check that the table describes this assembly and that --contigs-format (here
'1,2') selects its contig-name and genome-name columns.
```

The table is compared with every contig of the assembly — as named by the bins
file, the FASTA or the pairs `#chromsize:` header, whichever supplies the
contig lengths — *before* `--drop-unassigned` and `--filter` narrow them. A
filter that keeps only contigs the table does not list is therefore not a
refusal: those contigs are unbinned, as asked. `create cool` applies the same
rule against the cool file's contigs, and [`reassign`](#reassign) refuses such
a table too.

#### Contig and genome columns { #contigs-format }

`--contigs-format` says which two columns of the `--contigs` table hold the
contig and the genome names, in that order:

```console
--contigs-format CONTIG,GENOME[,noheader]
```

Each column is a **1-based position** (digits only) or a **header name**
(anything else), and the two can be mixed: `2,5`, `Name,Final_bin`,
`1,Final_bin`. Add `,noheader` when the table has no header line; only
positions work then. Names are looked up wherever they are in the header, and
every other column of the table is stored in `/contigs` (named `col3`, `col4`,
… in a table without a header).

The default, `contigs`, means `contig_name,genome_name`: looked up by name, so a
table given without its format fails with its header listed instead of being
read with the wrong columns. Aliases cover the tables common binners write:

| Alias | Same as | Table |
|-------|---------|-------|
| `contigs` | `contig_name,genome_name` | metacooler's own layout (the default). |
| `contig_genome` | `1,2,noheader` | Two columns, contig then genome, no header. |
| `metator` | `Name,Final_bin` | MetaTOR's `contig_data_final.txt`. |
| `dastool` | `1,2,noheader` | DAS Tool's `*_DASTool_contigs2bin.tsv`, and the tables its `Fasta_to_Contigs2Bin.sh` builds from FASTA bins (MetaBAT2, MaxBin2, …). |
| `binette` | `1,2,noheader` | Binette's `final_contig_to_bin.tsv`. |
| `semibin` | `contig,bin` | SemiBin2's `contig_bins.tsv`. |
| `vamb` | `contigname,clustername` | VAMB's `vae_clusters_unsplit.tsv` / `vae_clusters_split.tsv`. |

- **No binning yet, with contig metadata.** Leaving `--contigs` out builds a
  map with no binning yet; to store a table's columns (depth, GC, …) in it
  meanwhile, give the table with `NA` in the genome's place — `contig_name,NA`
  for a table with a header, `1,NA,noheader` for one without. Every contig then
  lands in `__unbinned__`, as without a table. `NA` has to be written out: a
  spec that just leaves the genome column off is an error, so a forgotten
  column cannot pass for "no binning". `--genomes` and `--drop-unassigned` are
  refused with it, as without a table.
- **An alias wins** over a header name it happens to equal; give such a column
  by position.
- **CONCOCT** writes comma-separated tables, which have no alias: convert one
  first (`tr ',' '\t' < clustering_merged.csv > clustering_merged.tsv`), then
  give its columns by position — contig first, cluster second: `1,2`, or
  `1,2,noheader` if the file has no header line.

#### Several libraries { #create-pairs-libraries }

Repeat `--pairs` to build one map from several libraries of the same assembly
(replicates, sequencing runs, restriction enzymes). The files are streamed back
to back, in the given order, as if they were one concatenated file, so
everything a single pass measures covers all of them, **exactly**: the summed
contacts, the QC block and its `qc_*` columns, the `<output>.ps.tsv`
P(s) table, and the dropped-pairs report with `n_pairs_read` /
`n_pairs_discarded`. One map per library joined with [`merge`](#merge) holds the
same contacts but not the rest: `merge` drops the QC block, and no P(s) table
comes out of it.

```bash
# Three libraries of one assembly, in one pass: one map, one QC block, one P(s)
metacooler create pairs \
    --pairs gut_T0_DpnII.pairs.gz \
    --pairs gut_T0_HinfI.pairs.gz \
    --pairs gut_T0_rep2.pairs.gz \
    --bins  2000 \
    --contigs contig_data_final.txt --contigs-format metator \
    -o gut_T0.metacool
```

- **Headers.** The files' `#chromsize:` lines are united, in file order (the
  first file's contigs, then each contig new to a later file). A file without a
  header, or with only some of the contigs, is fine; one contig given two
  different lengths is an error naming both files, since the libraries were
  then aligned to different assemblies. The headers are read only where they
  are needed: for the contig lengths (integer `--bins` without `--fasta`) or
  for the order (`--contig-order fasta` without `--fasta`).
- **Refused up front:** the same file given twice (its contacts would count
  twice), and a missing file — before any file is read, rather than after the
  others have been binned.
- **Streaming.** `--chunksize`, `--sort-strategy` and `--nproc` behave as for
  one file; a chunk may span two files, and the output is the same for every
  `--nproc`.
- **Provenance.** The `commands` record lists every `--pairs`, in order.
- **Recomputing P(s).** `ps --pairs` reads one file: to recompute the table of
  a map made from several, give it their concatenation (one header, then every
  file's pairs).

### `create cool` { #create-cool }

Create a `.metacool` file from an existing cooler file by adding the genome
hierarchy on top of its bins and pixels. **API:** `metacooler.from_cool`
(`metacooler.io.create.from_cool`).

```console
metacooler create cool --cool IN.cool [--contigs C] -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--cool` | *required* | Input `.cool` file to convert. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | — | Contig-to-genome table; which columns hold the names is set by `--contigs-format`. Must name at least one contig of the cool file — see below. Leave it out to build a map with no binning yet. |
| `--contigs-format` | `contigs` | Same grammar and aliases as `create pairs` — see [Contig and genome columns](#contigs-format). |
| `--genomes` | — | Optional genome metadata TSV. Decorates the genome rows derived from `--contigs`; same semantics and warnings as `create pairs`. |
| `--genomes-format` | `genome_name` | Which column of `--genomes` holds the genome names, as in `create pairs`. |
| `--metadata` | — | JSON or YAML metadata file (no key metacooler writes itself — see `create pairs`). |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `cool` (preserve the input's order), or `alpha`. |
| `--assembly-name` | — | Assembly / metagenome name, stored as `genome-assembly`. |
| `--drop-unassigned` | off | Drop unassigned contigs. |
| `--filter` | — | Keep only the genomes/contigs matching an expression, *before* the bins are remapped. Repeatable. Same semantics and caveats as `create pairs`. |
| `--qc / --no-qc` | on | Tally contact-quality statistics off the pixel stream being copied, and store them in `metadata.qc` plus the `qc_*` columns of `/genomes` and `/contigs`. A cool file has no read positions left, so intra-contig distances are measured between bin midpoints and the block records `source="pixels"`. See [`qc`](#qc). |
| `--qc-threshold` | `10000` | Distance separating short- from long-range intra-contig signal. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--nproc` | `1` | Parallel reader threads for pixel remapping. |
| `--sort-strategy` | `auto` | As `create pairs`, but `nnz` is known up front so the 100 M-pixel threshold is compared directly. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

- There is no `--fasta` option: contig lengths always come from the cool file's
  own bins table.
- `--contig-order cool` is the only ordering that round-trips row order exactly;
  the default regroups contigs by genome. Without `--contigs`, it gives back the
  cool file's contigs, bins and pixels unchanged, plus the genome layer with one
  `__unbinned__` genome — and a later [`reassign`](#reassign) relabels such a
  map in place instead of rewriting its pixels.
- The `--contigs` table follows the rule of
  [`create pairs`](#create-contig-table), compared with the contigs of the cool
  file (its `chroms` table): unlisted contigs go to `__unbinned__`, rows naming
  a contig the cool file does not have are ignored with a warning, and a table
  naming none of its contigs is refused before a single pixel is read — before
  `--drop-unassigned` and `--filter` narrow anything.

```bash
metacooler create cool \
    --cool existing.cool \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_quality.tsv \
    -o sample.metacool
```

### `collect` { #collect }

Collect single-genome cool files — one per species, strain, or MAG — into one
`.metacool` file, **one genome per input**. No contig-to-genome table is needed:
the file a contig comes from *is* its genome. **API:** `metacooler.collect`
(`metacooler.io.collect.collect`).

```console
metacooler collect A.cool B.cool [...] -o OUT.metacool [OPTIONS]
```

Inputs are `.cool` files, or one resolution of an `.mcool` each
(`E_coli.mcool::resolutions/1000`). Every input becomes a genome named after its
file stem (`E_coli.cool` → `E_coli`) unless `--names` says otherwise, and
**every contig is renamed `<genome>__<contig>`**, so the `chromosome` or
`contig_1` that so many assemblies share cannot collide.

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.metacool` path. It may not be one of the inputs. |
| `--names` | file stems | Genome name for each input, in input order — comma-separated or repeated. Names (or, without `--names`, file stems) must be unique. |
| `--genomes` | — | Optional genome metadata TSV; its genome-name column (`--genomes-format`) is matched against the genome names, and extra columns land in `/genomes`. A genome with no row keeps empty metadata, a row matching no input is ignored, and both are reported as warnings. |
| `--genomes-format` | `genome_name` | Which column of `--genomes` holds the genome names: a 1-based position or a header name, plus `,noheader` when the table has no header line. |
| `--contigs` | — | Optional per-contig metadata TSV; extra columns land in `/contigs`. See [The contigs table](#collect-contigs-table) below. |
| `--contigs-format` | `contig_name` | Which column of `--contigs` holds the contig names: a 1-based position or a header name, plus `,noheader` when the table has no header line — e.g. `contigName` or `1,noheader`. |
| `--metadata` | — | JSON or YAML metadata file (no key metacooler writes itself — see `create pairs`). |
| `--filter` | — | Keep only the genomes/contigs matching an expression, before anything is written — e.g. `--filter "contigs.length >= 5kb"`. Repeatable. Unlike in `create`, `genomes.n_bins` is available: the bins come with the inputs. See [Filtering with `--filter`](#filter). |
| `--contig-order` | `cool` | Order of the contigs *within* each genome: `cool` keeps each input's own order, `genome-then-length` puts the longest first, `genome-then-alpha` sorts them by name. Genomes themselves are always laid out alphabetically, in `genome_id` order. |
| `--assembly-name` | — | Name for the collection, stored as `genome-assembly`. |
| `--qc / --no-qc` | on | Tally contact-quality statistics off the pixel stream being copied (`source="pixels"`, as in `create cool`). See [`qc`](#qc). |
| `--qc-threshold` | `10000` | Distance separating short- from long-range intra-contig signal. Accepts `bp/kb/Mb/Gb` suffixes. |
| `--sort-strategy` | `auto` | Only matters when `--contig-order` reorders an input's contigs: `auto` (in memory unless one such input exceeds 100 M pixels), `memory` (in RAM, one input at a time), or `external` (cooler's external sort). |
| `--chunksize` | `50000000` | Pixels read per chunk (peak-memory control), and the smallest chunk handed to cooler. |
| `-f, --force` | off | Overwrite the output if it already exists. An input is never overwritten. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

!!! warning "The inter-genome blocks are empty"
    Every input was binned on its own, usually from its own library, so no
    contact *between* two inputs was ever observed and the collected map holds
    none. It is not the map you would get by aligning one metagenomic library to
    the combined assembly — that is [`create pairs`](#create-pairs). Sequencing
    depth also differs between genomes exactly as it did between the inputs. In
    the QC block, `inter_genome` is 0 and every genome's `qc_intra_genome_ratio`
    is 1, by construction.

**All inputs must share one bin size** — or all have variable
(restriction-fragment) bins — because a metacool file has a single bin table:

```console
$ metacooler collect E_coli_1kb.cool B_subtilis_5kb.cool -o species.metacool
Error: The inputs are not binned alike, but a metacool file has a single bin
table:
  E_coli_1kb.cool: 1000 bp
  B_subtilis_5kb.cool: 5000 bp
Collect maps of one resolution — pick the same resolution from each .mcool
(file.mcool::resolutions/<binsize>), or coarsen the finer ones first ('cooler
coarsen').
```

A few inputs are refused outright, each with a message saying what to do
instead: a `.metacool` file or `.mmetacool` entry (it already holds several
genomes, which collecting would squash into one — `flatten` or
`extract --split` it first), the root of an `.mcool` (the message lists its
resolutions), and a cooler stored as a full `square` matrix.

- **No sort under the default order.** Genomes are laid out one after another,
  so each input is a single block of bins and its pixels only move within it.
  Under `--contig-order cool` every block stays sorted, and the whole collection
  is written in one streaming pass **without sorting anything**. The other
  orderings re-sort only the inputs whose contigs they actually moved, one at a
  time.
- **Round trip.** Under `cool`, `extract --genomes E_coli` gives back
  `E_coli.cool`, bin for bin and pixel for pixel; only the contig names carry
  their `E_coli__` prefix. A contig that already starts with its genome's prefix
  is not prefixed twice, so collecting such an extracted file again keeps its
  names.
- **Contig lengths** come from each input's bins, as in `create cool`.
  Balancing weights (`bins/weight`) are not carried over; re-balance the result
  with `cooler balance` if you need them.
- **Pixel counts** keep their width: `int32` inputs stay `int32`, and one
  `int64` (or float) input widens the whole column rather than being truncated.

#### The contigs table { #collect-contigs-table }

`--contigs` attaches per-contig metadata — coverage, GC content, taxonomy — to
the collected contigs. It is a TSV whose contig-name column is the one
`--contigs-format` names: a 1-based position or a header name, plus
`,noheader` when the table has no header line (its other columns are then
named `col2`, `col3`, …). By default it is the column called `contig_name`,
wherever it is. Every column other than the contig names and `genome_name` is
stored in `/contigs`.

```text
contig_name	genome_name	coverage	gc
chromosome	E_coli	120.4	0.51
plasmid1	E_coli	310.2	0.47
chromosome	B_subtilis	98.0	0.44
```

- A row may name a contig **as in its input cool** (`chromosome`) or **as
  collected** (`E_coli__chromosome`).
- The `genome_name` column is optional. It says which genome a row describes,
  which is what makes a name several inputs share — `chromosome` above —
  unambiguous. Without it, such a name is an error rather than a guess.
- A `genome_name` is **checked, never used to re-assign**: a row placing a
  contig under a genome it does not come from is an error, since `collect`
  assigns every contig to the genome of its input file. A row about a genome
  that is not being collected is simply ignored, so one master table can serve
  several collections.
- Rows matching no contig are ignored and reported; a contig described twice
  keeps its last row. A column named `name`, `length`, or `genome_id` is stored
  as `user_name`, `user_length`, `user_genome_id`.
- Values must be ASCII: cooler reads every column of `/contigs`. See
  [Names must be ASCII](format.md#names-must-be-ascii).

```bash
# One genome per cool file, named after the files
metacooler collect E_coli.cool B_subtilis.cool S_aureus.cool -o species.metacool

# One resolution from each .mcool, with explicit genome names
metacooler collect \
    E_coli.mcool::resolutions/1000 B_subtilis.mcool::resolutions/1000 \
    --names Escherichia_coli,Bacillus_subtilis \
    -o species_1kb.metacool

# Genome and contig metadata, and no contig under 5 kb
metacooler collect *.cool \
    --genomes species.tsv --contigs contig_metadata.tsv \
    --filter "contigs.length >= 5kb" \
    -o species.metacool
```

---

## Edit an existing map { #edit }

Three commands change a map that already exists, each a different part of it,
and none re-bins a single contact.

| Command | What it changes | Cost |
|---------|-----------------|------|
| [`reassign`](#reassign) | Which genome each contig belongs to — moving contigs, merging or splitting MAGs, renaming a genome, unbinning contigs — and everything that follows from it: the `/genomes` rows and their ids, the per-genome QC counts, and possibly the contig layout. | Seconds when the contigs keep their places; a rewrite of every pixel when restoring a genome-grouped layout means moving them ([why](#reassign-layout)). |
| [`annotate`](#annotate) | The metadata columns of `/genomes` and `/contigs` — completeness, taxonomy, coverage, … — never the assignment or any column derived from it. | Seconds. |
| [`metadata`](#metadata) | The free-form `metadata` attribute: sample, enzyme, sequencing run, … | Instant. |

`reassign` and `annotate` write either an edited copy (`-o`) or the map itself
(`--inplace`); `metadata` always edits the map itself.

`reassign` and `annotate` meet in one place: after a reassignment, the metadata
of the genomes whose contigs changed still describes their old contigs, and
`reassign` names them. Refresh those values with `annotate --replace`, from a
QC run on the new bins.

One option of another command also edits a map: [`qc --store`](#qc) writes the
contact-quality block and the `qc_*` columns back into the map it measured.

### `reassign` { #reassign }

Re-attribute contigs to genomes (MAGs) in a file that already exists, from a
new contig-to-genome table — a re-binning, a curation pass, a merge of two
MAGs. The contact matrix is **never recomputed**: the contigs, their lengths
and their bins are those of the source file and every pixel keeps its count.
**API:** `metacooler.reassign` (`metacooler.io.reassign.reassign`).

```console
metacooler reassign METACOOL_PATH --contigs TABLE -o OUT.metacool [OPTIONS]
metacooler reassign METACOOL_PATH --contigs TABLE --inplace [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--contigs` | *required* | New contig-to-genome table. It need not list every contig — see `--missing`. Only its two name columns are read; extra columns are ignored. |
| `--contigs-format` | `contigs` | Same grammar and aliases as `create pairs` — see [Contig and genome columns](#contigs-format). `NA` as the genome column moves every contig the table lists to `__unbinned__`. |
| `-o, --output` | *required* | Output `.metacool` path. Required unless `--inplace`. |
| `--inplace` | off | Rewrite `METACOOL_PATH` itself instead of writing a new file. Only possible when no contig has to move (see `--layout`). |
| `--genomes` | — | Genome metadata table for the *new* genomes. **Replaces** the metadata columns carried over from the source `/genomes` table rather than merging with them. To add or change a few columns without reassigning anything, use [`annotate`](#annotate). |
| `--genomes-format` | `genome_name` | Which column of `--genomes` holds the genome names, as in `create pairs`. |
| `--missing` | `keep` | What to do with a contig the table does not mention. `keep`: leave it in its current genome, so a short table is an incremental edit. `unbinned`: move it to `__unbinned__`, so the table is authoritative. Neither drops a contig. |
| `--layout` | `auto` | `auto`, `keep`, `genome-then-alpha`, or `genome-then-length` — see below. |
| `--metadata` | — | JSON or YAML metadata file (no key metacooler writes itself — see `create pairs`). |
| `--qc / --no-qc` | on | Keep the contact-quality block when the source carries one, for the new genomes — see [What a reassignment invalidates](#reassign-invalidates). With `--no-qc` the block and its `qc_*` columns are dropped and not replaced. |
| `--chunksize` | `50000000` | Pixels read per chunk: remapped on the re-ordering path, and classified when the QC block is stored. |
| `-f, --force` | off | Overwrite the output if it already exists. The previous file is replaced only once the new one has validated, and an input (the map, `--contigs`, `--genomes`, `--metadata`) is never overwritten. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

What is rewritten, and what is not:

- `/genomes` is *derived* from the assignment, so the whole table is rebuilt:
  rows appear and disappear as genomes are created and emptied, every
  `genome_id` is renumbered (`__unbinned__` first, then alphabetical), and
  `length`, `n_contigs` and `n_bins` are recomputed. So is
  `indexes/genome_offset`.
- `/contigs` keeps its rows and every extra column; only `genome_id` changes
  (and its two QC columns, with the QC block — below).
- `/bins` and `/pixels` are untouched unless the layout has to be restored
  (below).

#### The layout, and why a reassignment can cost a full rewrite { #reassign-layout }

Under a genome-grouped `contig-ordering` (`genome-then-length`,
`genome-then-alpha` — the defaults) the contigs are *physically* laid out in
ascending `genome_id` order. That is what makes `genome_offset` delimit a
genome's contig block and a genome's bins a contiguous range, which
[`show`](#show), [`ps`](#ps) and `MetaCooler.genome_extent()` all rely on. Move
one contig between two non-adjacent MAGs and the layout no longer holds —
[`check`](#check) rejects the file, and it is right to.

Restoring it means permuting the contigs, hence the bins, hence every pixel's
bin ids: an O(nnz) rewrite. `--layout` is the choice between paying that and
giving up the guarantee:

| `--layout` | What it does | Cost |
|------------|--------------|------|
| `auto` | Re-group a file whose `contig-ordering` is genome-grouped; leave the order of any other file (`fasta`, `alpha`, `cool`) alone, since it never made the promise. A file whose genomes were contiguous anyway — every [`collect`](#collect)ed map, under `cool` — is relabelled too, with a warning naming any genome the new assignment scatters. | Whichever of the two below applies. |
| `keep` | Never re-order. The file is recorded as `contig-ordering="custom"`, because it can no longer claim to be grouped. | O(n_contigs) — the file is copied and four datasets are rewritten. |
| `genome-then-alpha`, `genome-then-length` | Re-group with that strategy, whatever the source used. | O(nnz) — bins are permuted and pixels remapped and re-sorted. |

A re-grouping that turns out to reproduce the order the file already has is
detected and skipped, so merging two adjacent MAGs, renaming one, or unbinning
a contig usually costs nothing. `--inplace` is available exactly when no
re-ordering is needed, and is refused (before writing anything) when one is.

Re-grouping folds every remapped pixel into the upper triangle, which is only
right for a symmetric matrix. A map stored as `storage-mode="square"` does not
promise one, so re-grouping it is refused; `--layout keep` relabels it without
moving a pixel.

!!! warning "`--layout keep` gives up contiguous genomes"
    Under `contig-ordering="custom"` a genome's contigs may be scattered
    through `/contigs`, so `genome_offset` no longer delimits them and anything
    that needs one bin range per genome fails loudly rather than returning a
    range that includes foreign bins:

    ```console
    $ metacooler show out.metacool --genomes MAG_01 -o mag01.png
    Error: Genome 'MAG_01' is not contiguous in bin space under contig-ordering
    'custom'; a single bin range cannot represent it. […]
    ```

    Selecting by contig *membership* is unaffected: `dump --genomes`,
    `extract --genomes`, `MetaCooler.bins_for_genome()` and
    `MetaCooler.genome_bin_ids()` all still work. Re-run with
    `--layout genome-then-alpha` to get the contiguity back.

#### What a reassignment invalidates { #reassign-invalidates }

- **Half of the QC block.** Two of its four categories
  (`inter_contig_intra_genome` and `inter_genome`), and every per-genome `qc_*`
  count, are defined by the contig-to-genome assignment; the two intra-contig
  categories are not — a contact inside one contig stays inside it whichever
  genome the contig is given to. So with `--qc` (the default) the block is
  *carried*, not recomputed: the per-contig `qc_intra_contig_short` /
  `qc_intra_contig_long` columns of `/contigs` are kept and summed per new
  genome, and only the two inter-contig categories are re-measured from the
  pixels — exactly, at any bin size, since all a pixel needs for them is the
  contig and the genome of each of its bins. Every other key (`source`,
  `distance_threshold`, `binsize`, `n_pairs_read`, `n_pairs_discarded`) is kept,
  so a map made by [`create pairs`](#create-pairs) keeps the short/long split
  it measured on read positions: reassigning it gives the block that
  `create pairs` with the new table would have given. A source without the
  per-contig columns — written by metacooler 0.2.0 or earlier — has its block
  recomputed from the pixels instead, at its own distance threshold, with
  intra-contig distances measured between bin midpoints (`source="pixels"`).
  A warning names what that loses when the block was measured on read
  positions (at one bin per contig, `intra_contig_long` drops to 0); otherwise
  an INFO line says so. The new block is measured before anything is written,
  so a run stopped meanwhile — `--inplace` included — leaves the source's
  block as it was. With `--no-qc` the block is dropped and
  [`qc --store`](#qc) can add one back later.
- **Per-genome metadata columns.** Columns from the source `/genomes` table are
  carried over for the genomes that survive, and a genome whose contig set
  changed keeps values (completeness, contamination, …) that were computed for
  a different set of contigs. Pass `--genomes` with a fresh table to replace
  them. `/contigs` columns are always still true — contigs are not modified.
- **`merge` compatibility.** [`merge`](#merge) requires identical contigs,
  genomes and ordering, so a reassigned file no longer merges with un-reassigned
  replicates. Reassign every replicate with the same table, or merge first and
  reassign the result.

!!! note "`reassign` never adds or removes a contig"
    A name in the table that is not a contig of the file is ignored with a
    warning (and a table matching *no* contig at all is refused — that is a
    wrong table, or a wrong `--contigs-format` — as it is by
    [`create`](#create-contig-table)). To *drop* contigs, and the pixels that
    touch them, use [`extract`](#extract).

```bash
# Incremental: move two contigs into MAG_07, leave everything else alone
metacooler reassign gut_T0.metacool --contigs fixes.tsv -o gut_T0_curated.metacool

# A fresh MetaTOR run is authoritative: anything it does not bin becomes unbinned
metacooler reassign gut_T0.metacool \
    --contigs contig_data_final.txt --contigs-format metator \
    --missing unbinned \
    --genomes mag_quality_v2.tsv \
    -o gut_T0_v2.metacool

# Relabel a 400 GB file without rewriting its pixel table
metacooler reassign huge.metacool --contigs fixes.tsv --layout keep --inplace

# One entry of a collection; the output is a standalone map
metacooler reassign reps.mmetacool::entries/rep1 --contigs fixes.tsv -o rep1_v2.metacool
```

### `annotate` { #annotate }

Add, replace or drop the **metadata columns** of `/genomes` and `/contigs` in a
map that already exists — completeness and taxonomy from CheckM2 or GTDB-Tk,
coverage from a depth table, anything that arrives after the map was built.
Rows are matched by name; the genome assignment, the bins and the pixels are
never touched, so this costs seconds even on a large map. **API:**
`metacooler.annotate` (`metacooler.io.annotate.annotate`).

```console
metacooler annotate METACOOL_PATH [--genomes T] [--contigs T] [OPTIONS] (-o OUT | --inplace)
```

| Option | Default | Description |
|--------|---------|-------------|
| `--genomes` | — | Genomes metadata table. Its genome-name column (`--genomes-format`) is matched against `/genomes/name`; every other column is added to `/genomes`. |
| `--genomes-format` | `genome_name` | Which column of `--genomes` holds the genome names: a 1-based position or a header name, plus `,noheader` when the table has no header line — e.g. `Name` for CheckM2's `quality_report.tsv`, `user_genome` for GTDB-Tk's summary. |
| `--contigs` | — | Per-contig metadata table. Its contig-name column (`--contigs-format`) is matched against `/contigs/name`; every other column is added to `/contigs`. A `genome_name` column is **checked**, never applied: a row placing a contig under another genome than the map's is an error — that is [`reassign`](#reassign)'s job. |
| `--contigs-format` | `contig_name` | Which column of `--contigs` holds the contig names: a 1-based position or a header name, plus `,noheader` when the table has no header line — e.g. `contigName` for a MetaBAT depth table. |
| `--columns` | all | Take only these columns of the tables, as `genomes.<name>` / `contigs.<name>` (comma-separated or repeated), named as in the tables. |
| `--drop` | — | Remove these metadata columns, as `genomes.<name>` / `contigs.<name>`. |
| `--replace` | off | Allow changing the values of a column the map already has. |
| `--missing` | `keep` | For a replaced column, what the rows the table does not list get. `keep`: their current value, so a short table is a patch. `na`: an empty value, so the table is authoritative. |
| `-o, --output` | — | Write an annotated copy here. Required unless `--inplace`. |
| `--inplace` | off | Edit `METACOOL_PATH` itself. |
| `-f, --force` | off | Overwrite the output if it already exists (only once the new file validates). An input is never overwritten. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

A column gets the same name, values and type it would have had if its table
had been given to `create` (or `collect`): a column called `name`, `length`,
`n_contigs`, `n_bins` or `genome_id` is stored as `user_<name>`, and a column is
numeric when every non-empty value parses as a number. So `--filter` reads an
annotated column exactly as it reads one added at creation:

```bash
# CheckM2 completeness and contamination, straight from its report
metacooler annotate gut_T0.metacool --inplace \
    --genomes quality_report.tsv --genomes-format Name \
    --columns genomes.Completeness,genomes.Contamination

# ... then use them
metacooler extract gut_T0.metacool --split -o hq \
    --filter "genomes.Completeness >= 90" --filter "genomes.Contamination < 5"

# Per-contig depth from MetaBAT's jgi_summarize_bam_contig_depths
metacooler annotate gut_T0.metacool --inplace \
    --contigs depth.txt --contigs-format contigName --columns contigs.totalAvgDepth

# CheckM2 again, on the bins a reassignment changed: patch only those rows
metacooler annotate gut_T0_v2.metacool --inplace --replace \
    --genomes quality_report_v2.tsv --genomes-format Name

# Every resolution of a zoomified collection at once
metacooler annotate gut_T0.mmetacool --inplace --genomes gtdbtk.bac120.summary.tsv \
    --genomes-format user_genome --columns genomes.classification
```

How rows and columns are matched:

- **Names, not positions.** A table may list its rows in any order and cover
  any subset of the map. A genome or contig it does not list gets an empty
  value in a new column; a table row naming nothing in the map is ignored with
  a warning; a table that matches *nothing* is refused — it is the wrong table,
  or the format points at the wrong column. Rows are never added or removed.
- **Existing columns.** Giving a column the values it already holds is not a
  change, so re-running a command is harmless. A table that would change an
  existing column is refused unless `--replace`, with the number of rows that
  would change and an example.
- **Derived columns are off limits.** `name`, `length`, `genome_id`,
  `n_contigs`, `n_bins` and the QC columns — every `qc_*` column of `/genomes`,
  `qc_intra_contig_short` and `qc_intra_contig_long` in `/contigs` — cannot be
  dropped or replaced: the assignment changes with [`reassign`](#reassign), the
  QC columns with [`qc --store`](#qc). Any other `qc_*` column of `/contigs`
  (`qc_pass`, say) is metadata like the rest.
- **Contig names.** A row names a contig as the map stores it — or, in a map
  made by [`collect`](#collect), as in its input cool (`chromosome`, with a
  `genome_name` column to say which genome), as `collect` itself accepts. Its
  `genome_name` may use the unbinned labels `create` reads as `__unbinned__`
  (`NA`, `ND`, `-1`, …).
- **ASCII in `/contigs`.** cooler reads every `/contigs` column, so a non-ASCII
  value there is refused before anything is written; `/genomes` metadata may
  hold any text.

Where it writes:

| `METACOOL_PATH` | With `--inplace` | With `-o OUT` |
|-----------------|------------------|---------------|
| `x.metacool` | `x` itself | a copy of `x` |
| `c.mmetacool::entries/rep1` | that entry, inside `c`; the other entries are untouched | a standalone `.metacool` |
| `c.mmetacool` | **every entry** of `c` | a copy of the whole collection |

Every table is checked against every map before anything is written, and each
annotated map is validated. The run is recorded in the map's metadata: the
command line in `commands`, and the tables and the columns added, replaced and
dropped in `annotate`. A run that changes nothing writes nothing. Since
[`merge`](#merge) compares only the identity of two maps, an annotated
replicate still merges with the others.

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

Four keys belong to metacooler, which writes them itself: `commands` (the
provenance record of every command that wrote the map), `qc` (the
contact-quality block), and `reassign` and `annotate` (those commands'
histories). `metadata list` and `metadata get` show them, but no `metadata`
subcommand and no `--metadata` file can set or remove them.

---

## Combine and re-bin { #combine }

Build a new map from the pixels of existing ones: [`merge`](#merge) sums
replicates of one assembly, [`binnify`](#binnify) re-bins to a coarser bin
size, and [`zoomify`](#zoomify) builds a multi-resolution collection.

### `merge` { #merge }

Sum the pixel counts of two or more metacool files that share **identical**
bins, contigs, genomes, and contig ordering (e.g. sequencing replicates of the
same assembly). Inputs may be `.mmetacool` entry URIs. Integer counts are
summed in `int64`, float counts in `float64`. **API:**
`metacooler.merge` (`metacooler.io.merge.merge`).

"Identical" means the same contig names and lengths, the same contig-to-genome
assignment and the same genome sizes — the columns that say what the pixels
are. Metadata columns (completeness, coverage, …, whether from creation or
[`annotate`](#annotate)) and the QC columns may differ: the merged map keeps
the first input's metadata, and each difference in a metadata column is
reported as a warning naming the input and the columns.

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

!!! tip "Libraries you still have as pairs files"
    `merge` drops the inputs' QC blocks and has no P(s) table to sum. When the
    libraries' pairs files are at hand, pass them all to
    [`create pairs`](#create-pairs-libraries) instead — one `--pairs` each — for
    the same contacts with an exact QC block and P(s) table.

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
| `--chunksize` | `50000000` | Source pixels read, remapped, and aggregated per chunk (peak-memory control). |
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
- At `--level contig` the re-binning keeps every pixel's row order, so the source
  is processed one block of whole output rows at a time (located through
  `bin1_offset`) and the output is written in order, with no temporary files and
  memory proportional to one `--chunksize` block. `--level genome` reorders
  contigs, so it re-aggregates through cooler's external merge.

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
| `--temp-dir` | system temp | Directory for the intermediate per-level files. Mirrors cooler's `zoomify --temp-dir`; point it at fast local scratch with room for **one** level (with `--nproc`, for the levels being built and the levels they are coarsened from). |
| `--nproc` | `1` | Levels built at once, each in its own process; a level waits for the finer level it is coarsened from (with `-r 1000,2000,5000,10000`: 2 kb and 5 kb together, then 10 kb). Levels are appended in ladder order, so the output is identical whatever the value; memory is one level's chunk per process. |
| `-f, --force` | off | Overwrite the output if it already exists. |
| `-v, --verbose` / `-q, --quiet` | off | Debug logging / suppress INFO logging. |

Entries are addressed as `file.mmetacool::entries/<binsize>`. Each level is
re-binned from the coarsest level already written whose bin size divides it
(with `-r 1000,2000,5000,10000`: 2 kb and 5 kb from 1 kb, 10 kb from 5 kb), as
cooler's `zoomify` does; since every level nests in the base, the counts are
those of re-binning the base directly.

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
| `dist` | Geometric centre of the separations the distance bin can actually hold, in bp — whole base pairs under `--pairs`, multiples of the bin size on a fixed-bin map — and the x coordinate to plot against. The leading bin holds separation 0, which a log axis cannot place, and sits at half its upper edge. |
| `n_contacts` | Contacts observed at that separation. |
| `n_bin_pairs` | Opportunities, in pairs of **positions** (base pairs) that *could* have contributed a contact — on both paths, so the two are directly comparable. A bin pair of widths `w_i`, `w_j` holds `w_i × w_j` of them, a diagonal one `w (w + 1) / 2`. (The column name predates the unit.) |
| `avg_contacts` | `n_contacts / n_bin_pairs`: contacts per position pair. |
| `prob` | `avg_contacts` normalised so the curve integrates to 1 over the reported range. Each distance bin weighs the share of separation it actually represents, so the result does not depend on `--bins-per-decade`. |

`n_bin_pairs` is what makes the curve mean anything: a contig of *L* bp offers
only `L - s` position pairs at separation *s*, so a raw count of long-range
contacts is starved by construction and would look like a decay even on
uniform data. Dividing it out is also what makes `--by genome` curves
comparable between a 5 Mb MAG and a 300 kb one, and `prob` is what makes two
libraries of different depth comparable. Distance bins that nothing could
populate are omitted rather than reported as a rate of zero.

#### Two sources { #ps-sources }

| | distances | `n_bin_pairs` counts | starts at |
|---|---|---|---|
| **pixels** (default) | between bin midpoints | **position** pairs, summed over the bin pairs at that midpoint separation | one bin |
| **pairs** (`--pairs`) | between read positions | **position** pairs, `L - s` per contig | 1 bp |

A binned map cannot say anything below its own resolution: every sub-bin
contact has already been collapsed onto the diagonal. Since the short-range
head of P(s) is exactly where a metaHi-C library shows undigested fragments and
self-religation, `create pairs` measures the curve while it still has the read
positions, and writes it as `<output>.ps.tsv` — the same table `--pairs`
recomputes, byte for byte (for a map read from several pairs files, from their
concatenation: see [Several libraries](#create-pairs-libraries)). Reach for
`--pairs` when you need that resolution on a file created earlier, or with a
different distance grid or selection.

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

`__unbinned__` is not a MAG but the absence of one, so each unbinned contig
counts as a genome of its own: a contact between two *different* unbinned
contigs is `inter_genome`. (Before block version 2 it was counted as
intra-genome, which inflated `intra_genome_fraction` by every contact among
the unassigned contigs.)

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
| `--store` | off | Write the result back into the file: `metadata.qc` plus the `qc_*` columns of `/genomes` and `/contigs`. |
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

They come with two per-contig columns in `/contigs`, `qc_intra_contig_short`
and `qc_intra_contig_long`: the same split contig by contig, which the
`/genomes` columns of those names sum per genome. Those two categories do not
depend on the genome assignment, so [`reassign`](#reassign) carries them over
to the new genomes rather than re-measuring them (see
[What a reassignment invalidates](#reassign-invalidates)).

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
file stem, `sample.metacool` → `sample`). Every input is validated first, so a
plain `.cool` or a whole `.mmetacool` is refused, and an entry name may not
contain `/`.

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
