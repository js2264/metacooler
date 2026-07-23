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
| `METACOOL_PATH` | A positional path to a `.metacool` file, or an `.mmetacool` entry URI (`file.mmetacool::entries/<name>`). |
| `-o, --output` | Output path. Required unless the command writes to stdout. |
| `-f, --force` | Overwrite the output if it already exists. |
| `-v, --verbose` | Enable debug logging. Available on every command. |
| `--genomes`, `--contigs` | Multi-value filters: repeat the flag or pass a comma-separated list (`--genomes MAG_01,MAG_02`). |

!!! note "`show` is the exception"
    Unlike every other command, [`show`](#show) takes its input file through the
    required option `-f/--file` rather than as a positional argument.

New to the tool? Start with the [Quickstart](quickstart.md); for the ideas
behind genomes, contigs, bins, and levels, see [Concepts](concepts.md).

---

## Create and edit

### `create pairs` { #create-pairs }

Create a `.metacool` file from a 4DN pairs file, a set of bins, and a
contig-to-genome assignment table. **API:** `metacooler.from_pairs`.

```console
metacooler create pairs --pairs P --bins B --contigs C -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--pairs` | *required* | Pairs file (plain text or gzipped). |
| `--bins` | *required* | BED file of bins (no header), **or** an integer bin size in bp (e.g. `1000`). For an integer, contig lengths must come from the pairs `#chromsize:` header or `--fasta`. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | *required* | Contig-to-genome table; layout set by `--contigs-format`. |
| `--contigs-format` | `contigs` | `contigs`, `contig_genome`, `metator`, or a column-index spec `I,J[,noheader]`. See [Creating from pairs](guides/from-pairs.md#contig-to-genome-table). |
| `--genomes` | — | Optional genome metadata TSV (first column `genome_name`); extra columns land in `/genomes`. |
| `--metadata` | — | JSON or YAML file of free-form metadata to embed. |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `fasta`, or `alpha`. See [Contig ordering](concepts.md#contig-ordering). |
| `--fasta` | — | Assembly FASTA, used for contig lengths when `--bins` is an integer and the pairs header has no `#chromsize:` lines. |
| `--assembly-name` | — | Name for this assembly / metagenome, stored in file attrs. |
| `--drop-unassigned` | off | Drop contigs absent from the contig table instead of collecting them under `__unbinned__`. |
| `--chunksize` | `1000000` | Pairs processed per chunk. |
| `--sort-strategy` | `auto` | Pixel sort strategy: `auto` (choose by pairs-file size, 2 GiB threshold), `memory` (sort in RAM — fast for small files), or `external` (stream to cooler's external sort — scalable). |

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
hierarchy on top of its bins and pixels. **API:** `metacooler.from_cool`.

```console
metacooler create cool --cool IN.cool --contigs C -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--cool` | *required* | Input `.cool` file to convert. |
| `-o, --output` | *required* | Output `.metacool` path. |
| `--contigs` | *required* | Contig-to-genome table; layout set by `--contigs-format`. |
| `--contigs-format` | `contigs` | Same choices as `create pairs`. |
| `--genomes` | — | Optional genome metadata TSV. |
| `--metadata` | — | JSON or YAML metadata file. |
| `--contig-order` | `genome-then-length` | `genome-then-length`, `genome-then-alpha`, `cool` (preserve the input's order), or `alpha`. |
| `--assembly-name` | — | Assembly / metagenome name. |
| `--drop-unassigned` | off | Drop unassigned contigs. |
| `--nproc` | `1` | Parallel reader threads for pixel remapping. |
| `--sort-strategy` | `auto` | As `create pairs`, but the `auto` threshold is 100 M pixels. |

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
same assembly). **API:** `metacooler.io.merge.merge`.

```console
metacooler merge A.metacool B.metacool [...] -o merged.metacool
```

```bash
metacooler merge rep1.metacool rep2.metacool -o merged.metacool
```

### `binnify` { #binnify }

Re-bin a metacool file to a coarser fixed bin size, writing a new metacool file.
The new size must be strictly larger than the current one. **API:**
`metacooler.ops.binnify.binnify`.

```console
metacooler binnify METACOOL_PATH --binsize N -o OUT.metacool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.metacool` path. |
| `--binsize` | *required* | New bin size in bp. |
| `--level` | `contig` | Tiling strategy: `contig` (tile each contig independently) or `genome` (concatenate a genome's contigs, then tile). |

```bash
metacooler binnify gut_T0.metacool --binsize 5000 -o gut_T0_5kb.metacool
```

### `zoomify` { #zoomify }

Build a multi-resolution `.mmetacool` collection from a fixed-bin metacool file.
Each successive level doubles the bin size. **API:**
`metacooler.ops.zoomify.zoomify`.

```console
metacooler zoomify METACOOL_PATH -o OUT.mmetacool [--levels N]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.mmetacool` path. |
| `--levels` | `5` | Maximum number of resolution levels to generate. |

```bash
metacooler zoomify gut_T0.metacool --levels 6 -o gut_T0_multireso.mmetacool
```

---

## Extract and view

### `extract` { #extract }

Write a standard `.cool` file containing only the selected genomes and/or
contigs, at the original bin resolution. Genome and contig selections combine.
**API:** `metacooler.io.extract.extract`.

```console
metacooler extract METACOOL_PATH --genomes G -o OUT.cool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.cool` path (or filename prefix with `--split`). |
| `--genomes` | — | Genome name(s) to extract. |
| `--contigs` | — | Contig name(s) to extract. |
| `--split` | off | Write one `.cool` per genome; `--output` becomes a prefix. |

```bash
# Single MAG
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool

# Multiple MAGs, one file each
metacooler extract gut_T0.metacool --genomes MAG_01,MAG_02 --split -o split
```

### `flatten` { #flatten }

Drop the genome hierarchy and write a standard `.cool` file in which every
contig is a chromosome. Bins and pixels are preserved exactly. **API:**
`metacooler.io.flatten`.

```console
metacooler flatten METACOOL_PATH -o flat.cool
```

### `dump` { #dump }

Dump any table to TSV (stdout by default). **API:** `metacooler.io.dump.dump_table`.

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
| `--prefix` | — | Filename prefix; **required** with `-t all`. |
| `--pixels-format` | `long` | `long` (row per pixel) or `wide` (2-D matrix). |
| `--weight` | `count` | Value column for `wide` format. |
| `-r, --range` | — | Row range (pixels only): genome name, contig name, or `contig:start-end`. |
| `-r2, --range2` | — | Column range (pixels only); defaults to `--range`. |

```bash
# Genomes table with header
metacooler dump gut_T0.metacool -t genomes -H

# Contigs of one MAG
metacooler dump gut_T0.metacool -t contigs --genomes MAG_01 -H

# Pixels within a genomic range
metacooler dump gut_T0.metacool -t pixels --range contig_0_0:0-50000

# Every table at once → gut_T0.{genomes,contigs,bins,pixels}.tsv
metacooler dump gut_T0.metacool -t all --prefix gut_T0
```

### `network` { #network }

Export a contact network as an edge list (TSV: `node1`, `node2`, `value`).
**API:** `metacooler.io.network`.

```console
metacooler network METACOOL_PATH --level LEVEL [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--level` | `contig` | Aggregation level: `bin`, `contig`, or `genome`. |
| `--weight` | `count` | Pixel column used as edge weight. |
| `-o, --output` | stdout | Output TSV file. |
| `--header` | off | Write a `node1  node2  value` header line. |
| `--chunksize` | `1000000` | Pixels processed per chunk (memory control). |

```bash
metacooler network gut_T0.metacool --level genome -o genome_network.tsv
```

### `aggregate` { #aggregate }

Aggregate contacts to contig or genome level and write a standard `.cool` file
(one chromosome per contig or per genome). **API:**
`metacooler.ops.aggregate.aggregate`.

```console
metacooler aggregate METACOOL_PATH --level LEVEL -o OUT.cool [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output` | *required* | Output `.cool` path. |
| `--level` | `contig` | `contig` (one chromosome per contig) or `genome` (sum contacts per genome pair). |
| `--count-col` | `count` | Pixel value column to aggregate. |
| `--chunksize` | `1000000` | Pixels read per chunk (peak-memory control). |

```bash
metacooler aggregate gut_T0.metacool --level genome -o genome_level.cool
```

### `show` { #show }

Render a contact-map heatmap to an image file. **API:**
`metacooler.ops.show`.

!!! warning "Input and output flags"
    `show` reads its input from `-f/--file` (**not** a positional path) and
    always writes to `-o/--output`; both are required.

```console
metacooler show -f METACOOL_PATH -o heatmap.png [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | *required* | Input `.metacool` file or mmetacool entry URI. |
| `-o, --output` | *required* | Output image path (`.png`, `.pdf`, `.svg`). |
| `--level` | `genome` | `bin` (finest), `contig`, or `genome` (coarsest). |
| `--genome` | — | Restrict to one genome (intra-genome view). |
| `--genomes` | — | A pair of genomes (inter-genome view). |
| `--transform` | `log1p` | Value transform: `log1p`, `log10`, or `none`. |
| `--cmap` | `afmhot_r` | Matplotlib colormap name. |
| `--vmax` | — | Upper colormap bound (default: 99th percentile of non-zero values). |
| `--separators / --no-separators` | on | Draw lines at genome/contig boundaries. |
| `--dpi` | `500` | Output image resolution. |
| `-r, --range` | — | Row range (implies `--level bin`): genome, contig, or `contig:start-end`. |
| `-r2, --range2` | — | Column range; defaults to `--range`. |

```bash
# Genome-level heatmap (default level)
metacooler show -f gut_T0.metacool --level genome -o genome_heatmap.png

# Intra-genome view at bin resolution
metacooler show -f gut_T0.metacool --level bin --genome MAG_01 -o MAG_01.png

# A specific genomic range, raw counts, custom colormap
metacooler show -f gut_T0.metacool -r contig_0_0:0-100000 \
    --transform none --cmap Reds -o range.png
```

---

## Inspect

### `info` { #info }

Print core file metadata and statistics (bin size, counts, ordering, …). Prints
JSON by default; `--table` renders a rich table. Also accepts standard
`.cool`/`.mcool` files and mmetacool entries.

```bash
metacooler info gut_T0.metacool           # JSON
metacooler info gut_T0.metacool --table   # rich table
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

### `attrs` { #attrs }

Print the HDF5 attributes of a group as YAML (same format as `cooler attrs`). The
group path is optional and defaults to the file root `/`.

```bash
metacooler attrs gut_T0.metacool           # root attributes
metacooler attrs gut_T0.metacool /genomes  # a specific group
```

---

## Manage collections

### `ls` { #ls }

List the entries in a `.mmetacool` (or `.mcool`, or plain `.metacool`) file.

```bash
metacooler ls gut_T0.mmetacool
# gut_T0.mmetacool::entries/rep1
# gut_T0.mmetacool::entries/rep2
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
| `--names` | file stems | Comma-separated entry names (must match the number of inputs). |

```bash
metacooler stack rep1.metacool rep2.metacool --names "rep1,rep2" -o combined.mmetacool
```

### `cp` { #cp }

Copy an entry within or between `.mmetacool` files. Both endpoints are URIs.
**API:** `metacooler.copy_entry` / `metacooler.copy_entry_cross_file`.

```bash
metacooler cp combined.mmetacool::entries/rep1 combined.mmetacool::entries/rep1_backup
metacooler cp source.mmetacool::entries/rep1 target.mmetacool::entries/rep1
```

### `mv` { #mv }

Rename an entry within a single `.mmetacool` file. **API:**
`metacooler.move_entry`.

```bash
metacooler mv combined.mmetacool::entries/rep1_backup combined.mmetacool::entries/rep1_v2
```

### `rm` { #rm }

Remove an entry from a `.mmetacool` file (asks for confirmation unless `-f`).
**API:** `metacooler.remove_entry`.

```bash
metacooler rm combined.mmetacool::entries/rep1_v2
```
