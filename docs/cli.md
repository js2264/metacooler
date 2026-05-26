# CLI reference

## Create and edit

### `metacooler create pairs`

Create a `.metacool` file from a pairs file, bins, and contig-to-genome
assignments.

```
Usage: metacooler create pairs [OPTIONS]

Options:
  --pairs PATH              Pairs file (plain text or gzipped).  [required]
  --bins TEXT               BED file of bins (no column names), or an integer
                            bin size in bp (e.g. 1000).  [required]
  -o, --output PATH         Output .metacool file path.  [required]
  --contigs PATH            Contig-to-genome assignment table.  [required]
  --contigs-format TEXT     Format of the contig table. One of: contigs,
                            contig_genome, metator, I,J[,noheader].
                            [default: contigs]
  --genomes PATH            Optional genomes metadata table (TSV with header,
                            first column must be 'genome_name').
  --metadata PATH           Path to a JSON or YAML file with metadata.
  --contig-order CHOICE     Contig ordering strategy. Choices: genome-then-length
                            (default), genome-then-alpha, fasta, alpha.
  --fasta PATH              Assembly FASTA for contig length extraction.
  --assembly-name TEXT      Name for this assembly / metagenome.
  --drop-unassigned         Drop contigs not assigned to any genome.
  --chunksize INT           Pairs per chunk.  [default: 10000000]
  -f, --force               Overwrite output file if it exists.
  -v, --verbose             Debug logging.
  -h, --help                Show this message and exit.
```

**Examples:**

```bash
# Fragment-level bins from MetaTOR output
metacooler create pairs \
    --pairs contacts.pairs.gz \
    --bins  fragments.bed \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_table.tsv \
    --assembly-name gut_T0 \
    --output gut_T0.metacool

# Fixed 2 kb bins, contig lengths from FASTA
metacooler create pairs \
    --pairs contacts.pairs.gz \
    --bins  2000 --fasta assembly.fa \
    --contigs contigs2mag.tsv --contigs-format contig_genome \
    --output gut_T0_2kb.metacool
```

---

### `metacooler create cool`

Create a `.metacool` file from an existing cooler file.

```
Usage: metacooler create cool [OPTIONS]

Options:
  --cool PATH               Input cool file to convert.  [required]
  -o, --output PATH         Output .metacool file path.  [required]
  --contigs PATH            Contig-to-genome assignment table.  [required]
  --contigs-format TEXT     Format of the contig table.  [default: contigs]
  --genomes PATH            Optional genomes metadata table.
  --metadata PATH           Path to a JSON or YAML file with metadata.
  --contig-order CHOICE     Contig ordering strategy. Includes 'cool' to
                            preserve the input file's order.
                            [default: genome-then-length]
  --assembly-name TEXT      Name for this assembly / metagenome.
  --drop-unassigned         Drop unassigned contigs.
  -f, --force               Overwrite output.
  -v, --verbose             Debug logging.
  -h, --help                Show this message and exit.
```

---

### `metacooler merge`

Merge two or more `.metacool` files that share identical contigs and genomes
tables (e.g. technical replicates from the same assembly).

```
Usage: metacooler merge [METACOOL_PATHS]... [OPTIONS]

Options:
  -o, --output PATH    Output .metacool file path.  [required]
  -f, --force          Overwrite output.
  -v, --verbose        Debug logging.
  -h, --help           Show this message and exit.
```

```bash
metacooler merge rep1.metacool rep2.metacool -o merged.metacool
```

---

### `metacooler binnify`

Re-bin a `.metacool` file to a coarser fixed resolution.

```
Usage: metacooler binnify METACOOL_PATH [OPTIONS]

Arguments:
  METACOOL_PATH    Source .metacool file.

Options:
  -o, --output PATH     Output .metacool file path.  [required]
  --binsize INT         New bin size in bp.  [required]
  --level CHOICE        Tiling level: contig (default) or genome.
  -f, --force           Overwrite output.
  -v, --verbose         Debug logging.
  -h, --help            Show this message and exit.
```

```bash
metacooler binnify gut_T0.metacool --binsize 5000 -o gut_T0_5kb.metacool
```

---

### `metacooler zoomify`

Create a multi-resolution `.mmetacool` collection from a base-resolution
`.metacool` file.

```
Usage: metacooler zoomify METACOOL_PATH [OPTIONS]

Arguments:
  METACOOL_PATH    Source fixed-bin .metacool file.

Options:
  -o, --output PATH    Output .mmetacool file path.  [required]
  --levels INT         Maximum number of resolution levels.  [default: 5]
  -f, --force          Overwrite output.
  -v, --verbose        Debug logging.
  -h, --help           Show this message and exit.
```

```bash
metacooler zoomify gut_T0.metacool --levels 6 -o gut_T0_multireso.mmetacool
```

---

## Extract and view

### `metacooler extract`

Extract selected genomes and/or contigs as a standard `.cool` file.

```
Usage: metacooler extract METACOOL_PATH [OPTIONS]

Options:
  --genomes TEXT    Genome name(s) to extract (comma-separated or repeated).
  --contigs TEXT    Contig name(s) to extract (comma-separated or repeated).
  -o, --output PATH Output .cool file path.  [required]
  --split           Write one .cool per genome (output treated as prefix).
  -f, --force       Overwrite output.
  -v, --verbose     Debug logging.
  -h, --help        Show this message and exit.
```

```bash
# Single MAG
metacooler extract gut_T0.metacool --genomes MAG_01 -o MAG_01.cool

# Multiple MAGs, one file each
metacooler extract gut_T0.metacool --genomes MAG_01,MAG_02 --split -o split
```

---

### `metacooler flatten`

Remove the genome hierarchy and write a standard `.cool` file with contigs
as chromosomes.

```
Usage: metacooler flatten METACOOL_PATH [OPTIONS]

Options:
  -o, --output PATH    Output .cool file path.  [required]
  -f, --force          Overwrite output.
  -v, --verbose        Debug logging.
  -h, --help           Show this message and exit.
```

---

### `metacooler dump`

Dump a table to TSV (stdout by default).

```
Usage: metacooler dump METACOOL_PATH [OPTIONS]

Options:
  -t, --table CHOICE      Table to dump: bins, chroms, contigs, genomes,
                          pixels, indexes, all.  [default: pixels]
  --genomes TEXT          Filter by genome name(s).
  --contigs TEXT          Filter by contig name(s).
  -r, --range TEXT        Row range: genome name, contig name, or
                          contig:start-end. Pixels only.
  -r2, --range2 TEXT      Column range (defaults to --range). Pixels only.
  -H, --header            Include column header row.
  --pixels-format CHOICE  long (default) or wide.
  --weight TEXT           Column name for wide format values.  [default: count]
  --prefix TEXT           File prefix (required when --table all).
  -o, --output PATH       Output file (default: stdout).
  -v, --verbose           Debug logging.
  -h, --help              Show this message and exit.
```

---

### `metacooler network`

Export a network edge list (TSV, no header: `node1`, `node2`, `value`).

```
Usage: metacooler network METACOOL_PATH [OPTIONS]

Options:
  --level CHOICE    Aggregation level: bin, contig (default), genome.
  --weight TEXT     Pixel column to use as edge weight.  [default: count]
  -o, --output PATH Output file (default: stdout).
  -v, --verbose     Debug logging.
  -h, --help        Show this message and exit.
```

---

### `metacooler aggregate`

Aggregate to contig or genome level and write a standard `.cool` file.

```
Usage: metacooler aggregate METACOOL_PATH [OPTIONS]

Options:
  -o, --output PATH    Output .cool file path.  [required]
  --level CHOICE       contig (default) or genome.
  --count-col TEXT     Pixel value column to aggregate.  [default: count]
  --chunksize INT      Pixels per chunk.  [default: 10000000]
  -f, --force          Overwrite output.
  -v, --verbose        Debug logging.
  -h, --help           Show this message and exit.
```

---

### `metacooler show`

Render a contact map heatmap.

```
Usage: metacooler show METACOOL_PATH [OPTIONS]

Options:
  --level CHOICE      bin, contig, or genome (default).
  --genome TEXT       Intra-genome view at bin level.
  --genomes TEXT TEXT Inter-genome view (two genome names).
  -r, --range TEXT    Row range (implies --level bin).
  -r2, --range2 TEXT  Column range (defaults to --range).
  --transform CHOICE  log1p (default) or none.
  --cmap TEXT         Matplotlib colormap.  [default: YlOrRd]
  -o, --output PATH   Save to file (PNG/PDF/SVG).
  --dpi INT           Image DPI.  [default: 150]
  -v, --verbose       Debug logging.
  -h, --help          Show this message and exit.
```

---

## Inspect

### `metacooler info`

Display file metadata and statistics.

```bash
metacooler info gut_T0.metacool
metacooler info gut_T0.metacool --table    # rich table output
```

### `metacooler tree`

Display the HDF5 group hierarchy.

```bash
metacooler tree gut_T0.metacool
metacooler tree gut_T0.mmetacool
```

### `metacooler attrs`

Display HDF5 attributes as YAML (same format as `cooler attrs`).

```bash
metacooler attrs gut_T0.metacool
metacooler attrs gut_T0.metacool /genomes
```

### `metacooler metadata`

Manage user metadata stored in the file's `metadata` attribute.

```bash
# List all metadata
metacooler metadata list gut_T0.metacool

# Get a single key
metacooler metadata get gut_T0.metacool -k sample

# Add / update a key
metacooler metadata add gut_T0.metacool -k enzyme --value DpnII

# Load from a YAML/JSON file
metacooler metadata load sample_meta.yaml gut_T0.metacool

# Remove a key
metacooler metadata remove gut_T0.metacool -k enzyme
```

---

## Manage collections

### `metacooler ls`

List entries in a `.mmetacool` (or `.mcool`, or plain `.metacool`) file.

```bash
metacooler ls gut_T0.mmetacool
# gut_T0.mmetacool::entries/rep1
# gut_T0.mmetacool::entries/rep2
```

### `metacooler stack`

Combine multiple `.metacool` files into a single `.mmetacool` collection.

```bash
metacooler stack rep1.metacool rep2.metacool \
    --names "rep1,rep2" \
    -o combined.mmetacool
```

### `metacooler cp`

Copy an entry within or between `.mmetacool` files.

```bash
# Same-file copy
metacooler cp combined.mmetacool::entries/rep1 combined.mmetacool::entries/rep1_backup

# Cross-file copy
metacooler cp source.mmetacool::entries/rep1 target.mmetacool::entries/rep1
```

### `metacooler mv`

Rename an entry within a `.mmetacool` file (same file only).

```bash
metacooler mv combined.mmetacool::entries/rep1_backup combined.mmetacool::entries/rep1_v2
```

### `metacooler rm`

Remove an entry from a `.mmetacool` file.

```bash
metacooler rm combined.mmetacool::entries/rep1_v2
```
