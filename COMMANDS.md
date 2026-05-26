# metacooler CLI — Command Reference

Commands are grouped by purpose. Each section lists the CLI invocation, key
options, and the underlying API entry point.

---

## Create and edit

### `create pairs` — Create a metacool file from a pairs file

**API:** `metacooler.io.create.from_pairs`

```
metacooler create pairs \
  --pairs contacts.pairs.gz \
  --bins 1000 \
  --contigs contigs_table.tsv \
  -o output.metacool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--pairs` | required | Pairs file (plain text or gzipped). |
| `--bins` | required | BED bins file path, or an integer bin size in bp. If an integer, contig lengths must come from the pairs header or `--fasta`. |
| `--contigs` | required | Contig-to-genome assignment table. Layout depends on `--contigs-format`. |
| `--contigs-format` | `contigs` | Format of the contig table: `contigs`, `contig_genome`, `metator`, or `I,J[,noheader]`. |
| `--genomes` | — | Optional genomes metadata TSV (first column `genome_name`). Extra columns stored in `/genomes`. |
| `--contig-order` | `genome-then-length` | Contig ordering strategy: `genome-then-length`, `genome-then-alpha`, `fasta`, `alpha`, `custom`. |
| `--fasta` | — | Assembly FASTA. Required when `--bins` is an integer and the pairs header has no `#chromsize:` lines. |
| `--assembly-name` | — | Assembly/metagenome name stored in file attrs. |
| `--drop-unassigned` | off | Drop contigs absent from the contig table instead of collecting them under `__unbinned__`. |
| `--chunksize` | `10000000` | Pairs processed per chunk. |
| `-o/--output` | required | Output `.metacool` file path. |
| `--metadata` | — | JSON or YAML file with additional metadata to embed. |
| `-f/--force` | off | Overwrite output if it already exists. |
| `-v/--verbose` | off | Enable debug logging. |

Notes:
- `--bins` accepts either a BED file (variable-width bins, e.g. restriction fragments) or an integer (fixed-width tiling).
- When `--bins` is an integer and there are no `#chromsize:` lines in the pairs header, `--fasta` is required to derive contig lengths.

---

### `create cool` — Create a metacool file from an existing cool file

**API:** `metacooler.io.create.from_cool`

```
metacooler create cool \
  --cool sample.cool \
  --contigs contigs_table.tsv \
  -o output.metacool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--cool` | required | Input `.cool` file to convert. |
| `--contigs` | required | Contig-to-genome assignment table. Layout depends on `--contigs-format`. |
| `--contigs-format` | `contigs` | Format of the contig table: `contigs`, `contig_genome`, `metator`, or `I,J[,noheader]`. |
| `--genomes` | — | Optional genomes metadata TSV (first column `genome_name`). |
| `--contig-order` | `genome-then-length` | Contig ordering strategy: `genome-then-length`, `genome-then-alpha`, `cool`, `alpha`, `custom`. |
| `--assembly-name` | — | Assembly/metagenome name stored in file attrs. |
| `--drop-unassigned` | off | Drop contigs absent from the contig table. |
| `-o/--output` | required | Output `.metacool` file path. |
| `--metadata` | — | JSON or YAML file with additional metadata. |
| `-f/--force` | off | Overwrite output if it already exists. |
| `-v/--verbose` | off | Enable debug logging. |

Notes:
- `--contig-order cool` preserves the exact contig order from the input cool file's chroms table (not available for `create pairs`).
- No `--fasta` option: contig lengths are always taken from the cool file's bins.

---

### `merge` — Merge metacool files with identical axes

**API:** `metacooler.io.merge`

```
metacooler merge a.metacool b.metacool -o merged.metacool
```

Merges two or more metacool files that share the same bins/contigs/genomes by
summing pixel counts. All inputs must have identical bin axes.

---

### `binnify` — Re-bin a variable-width metacool to a fixed bin size

**API:** `metacooler.ops.binnify`

```
metacooler binnify sample.metacool --binsize 1000 -o sample_1kb.metacool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--binsize` | required | New bin size in base pairs. |
| `--level` | `contigs` | Tiling level: `contigs` (each contig tiled independently) or `genomes` (contigs concatenated per genome before tiling). |

Notes:
- Input must be a variable-width bin metacool (e.g., generated from restriction fragments).
- At `--level contigs`, each contig is tiled independently; the last bin of each contig may be shorter.
- At `--level genomes`, each genome is tiled as one continuous sequence; the last bin of each genome may be shorter.
- Pixels are re-assigned and aggregated according to the new bin layout.

---

### `zoomify` — Create a multi-resolution mmetacool file

**API:** `metacooler.ops.zoomify`

```
metacooler zoomify sample.metacool -o sample.mmetacool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--levels` | `5` | Maximum number of resolution levels to generate. |

Notes:
- Input must be a fixed-bin metacool.
- Each level doubles the bin size (e.g. 1 kb → 2 kb → 4 kb → …).
- Output is a `.mmetacool` collection containing one entry per resolution.

---

## Extract and view

### `extract` — Extract a cool file for specific genome(s) or contig(s)

**API:** `metacooler.io.extract`

```
metacooler extract sample.metacool --genomes genome_A -o genome_A.cool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--genomes` | — | Genome name(s) to extract (comma-separated or repeated). |
| `--contigs` | — | Contig name(s) to extract (comma-separated or repeated). |
| `--split` | off | Write one `.cool` per genome; `--output` is used as a filename prefix. |
| `-o/--output` | required | Output `.cool` file path (or prefix when `--split` is set). |
| `-f/--force` | off | Overwrite output if it already exists. |

Notes:
- `--genomes` and `--contigs` can be combined; the result includes all contigs belonging to the named genomes plus any individually named contigs.
- Output is a standard cooler file compatible with cooltools, HiGlass, etc.

---

### `flatten` — Flatten a metacool into a standard cool file

**API:** `metacooler.io.flatten`

```
metacooler flatten sample.metacool -o flat.cool
```

Removes the genome hierarchy: every contig becomes a "chromosome" in the output
cool file. The output is a valid cooler file.

---

### `dump` — Dump tables from a metacool file

**API:** `metacooler.io.dump`

```
metacooler dump sample.metacool -t pixels -H -o pixels.tsv
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `-t/--table` | `pixels` | Table to dump: `genomes`, `contigs`, `chroms`, `bins`, `pixels`, `indexes`, or `all`. |
| `--genomes` | — | Filter output to the named genome(s). |
| `--contigs` | — | Filter output to the named contig(s). |
| `-H/--header` | off | Include a header row. |
| `-o/--output` | stdout | Output file path. |
| `--prefix` | — | Filename prefix (required when `--table all`; `--output` is ignored). |
| `--pixels-format` | `long` | Pixel format: `long` (row-per-pixel) or `wide` (2-D matrix). |
| `--weight` | `count` | Column to use as matrix values when `--pixels-format wide`. |

Notes:
- `--table all` writes `.genomes.tsv`, `.contigs.tsv`, `.bins.tsv`, and `.pixels.tsv` under `--prefix`; the `indexes` group is not dumped.
- `--weight` is only meaningful with `--pixels-format wide`.

---

### `network` — Export a contact network as an edge list

**API:** `metacooler.io.network`

```
metacooler network sample.metacool --level contigs -o contigs_network.tsv
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--level` | `contigs` | Aggregation level: `bins`, `contigs`, or `genomes`. |
| `--weight` | `count` | Pixel table column to use as edge weight. |
| `-o/--output` | stdout | Output TSV file (no header; columns: `node1`, `node2`, `value`). |

---

### `aggregate` — Aggregate contact counts to contig or genome level

**API:** `metacooler.ops.aggregate`

```
metacooler aggregate sample.metacool --level contig -o aggregated.cool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--level` | `contig` | Aggregation level: `contig` or `genome`. |
| `--count-col` | `count` | Pixel column to aggregate. |
| `--chunksize` | `10000000` | Pixels per processing chunk. |

---

### `show` — Visualise a contact map

---

## Inspect

### `info` — Display main file attributes

```
metacooler info sample.metacool
```

Prints the core HDF5 attributes (bin size, creation date, number of bins/contigs/genomes, etc.) in a formatted table. Accepts `.metacool` files only.

---

### `tree` — Display the HDF5 hierarchy

```
metacooler tree sample.metacool
```

Prints the HDF5 group/dataset tree of the file.

---

### `attrs` — Display all HDF5 attributes

```
metacooler attrs sample.metacool
```

Prints all attributes of every group in the file. Accepts both `.metacool` and `.mmetacool` files.

---

### `metadata` — Read and edit the metadata field

```
metacooler metadata show sample.metacool
metacooler metadata set  sample.metacool --key experiment --value metaHiC
```

Subcommands to read, add, update, or remove keys in the flexible `metadata` dict stored in file attrs. Accepts `.metacool` files only.

---

## Manage collections

### `stack` — Combine metacool files into an mmetacool collection

```
metacooler stack rep1.metacool rep2.metacool -o reps.mmetacool
```

Key options:

| Option | Default | Description |
|--------|---------|-------------|
| `--names` | — | Comma-separated entry names (must match number of inputs). Defaults to the input file stems. |
| `-o/--output` | required | Output `.mmetacool` file path. |
| `-f/--force` | off | Overwrite existing output. |

---

### `ls` — List entries in an mmetacool file

```
metacooler ls reps.mmetacool
```

Prints all entry paths (e.g. `reps.mmetacool::entries/rep1`). For a plain
`.metacool` file, prints `sample.metacool::/`.

---

### `cp` — Copy an entry within or between mmetacool files

```
metacooler cp src.mmetacool::entries/rep1 dst.mmetacool::entries/rep1_copy
```

URI syntax: `file.mmetacool::entries/<entry_name>`. Cross-file copy is supported.
Use `-f/--force` to overwrite an existing destination without confirmation.

---

### `mv` — Move/rename an entry within an mmetacool file

```
metacooler mv reps.mmetacool::entries/rep1 reps.mmetacool::entries/replicate1
```

Renames an entry in place. Use `-f/--force` to overwrite without confirmation.

---

### `rm` — Remove an entry from an mmetacool file

```
metacooler rm reps.mmetacool::entries/rep1
```

Removes the named entry. Requires confirmation unless `-f/--force` is given.

---

## (m)metacool file structure

```
@attrs:
  format: HDF5::MMetaCooler
  format-version: 1
entries/
  <entry_name>/          ← full metacool layout
    @attrs: { format: HDF5::MetaCooler, … }
    genomes/
    contigs/
    bins/
    pixels/
    indexes/
    chroms → contigs     ← soft link for cooler compatibility
```
