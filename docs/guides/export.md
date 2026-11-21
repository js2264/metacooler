# Exporting data

metacooler provides several ways to get data out of `.metacool` files,
depending on your downstream use case.

## Dumping tables to TSV (`dump`)

Dump any table from the file to tab-separated text:

```bash
# Genome / MAG table (with header)
metacooler dump sample.metacool -t genomes -H

# Contig table, filtered to one MAG
metacooler dump sample.metacool -t contigs --genomes MAG_01 -H

# Bins table for two MAGs
metacooler dump sample.metacool -t bins --genomes MAG_01,MAG_02 -H

# Raw pixel table (long format, to file)
metacooler dump sample.metacool -t pixels -o pixels.tsv

# Pixel table in wide (matrix) format — use with caution on large files
metacooler dump sample.metacool -t pixels --pixels-format wide -o matrix.tsv

# All tables at once (writes sample.{genomes,contigs,bins,pixels,indexes}.tsv)
metacooler dump sample.metacool -t all --prefix sample

# Re-run over an existing file: dump refuses to clobber -o without --force
metacooler dump sample.metacool -t pixels -o pixels.tsv --force

# Quiet mode: drop the INFO log lines, keep warnings and errors
metacooler dump sample.metacool -t genomes -H -q
```

!!! note "Options that do not apply to the chosen table are now errors"
    `--weight`/`--join` with a non-pixels table, `--genomes`/`--contigs` with
    `-t indexes`, `-o` with `-t all`, and `--prefix` with any other table all
    used to be accepted and silently ignored. They now fail with a usage error.
    An explicitly empty filter (`--genomes ''`) is rejected too, instead of
    quietly dumping everything. See the
    [`dump` reference](../cli.md#dump).

### Range queries on pixels

```bash
# All pixels with at least one bin in MAG_01
metacooler dump sample.metacool -t pixels --range MAG_01

# All pixels between two genomic ranges
metacooler dump sample.metacool -t pixels \
    --range contig_0_0:0-50000 --range2 contig_0_1:0-30000
```

## Aggregating to a standard cool file (`aggregate`)

Collapse the per-bin contact matrix to contig or genome level and write a
standard cooler file. The output is compatible with cooltools, HiGlass, and
anything else that reads `.cool` files.

```bash
# One bin per contig
metacooler aggregate sample.metacool --level contig -o contig_level.cool

# One bin per genome
metacooler aggregate sample.metacool --level genome -o genome_level.cool
```

The output bins table has one entry per contig (or genome), spanning its full
length (`start=0, end=length`).

## Extracting genomes or contigs (`extract`)

Write a standard `.cool` file containing only the selected genomes and/or
contigs, preserving the original bin resolution:

```bash
# Single MAG
metacooler extract sample.metacool --genomes MAG_01 -o MAG_01.cool

# Multiple MAGs in one file
metacooler extract sample.metacool --genomes MAG_01,MAG_02 -o two_mags.cool

# One cool file per MAG (--split)
metacooler extract sample.metacool --genomes MAG_01,MAG_02 --split -o split

# Specific contigs
metacooler extract sample.metacool --contigs contig_0_0,contig_0_1 -o subset.cool

# Combine genome and contig filters
metacooler extract sample.metacool --genomes MAG_01 --contigs contig_3_0 -o mixed.cool
```

## Writing sequences (`fasta`)

Write the contigs of each genome as FASTA, from the assembly the map was made
from (plain or gzip-compressed). Each contig's record is copied as it is, header
and line breaks included, one `<genome>.fa` per genome:

```bash
# One file per MAG
metacooler fasta sample.metacool --fasta assembly.fa.gz -o mags

# Each MAG as one sequence, its contigs joined by 100 N; the unbinned contigs too
metacooler fasta sample.metacool --fasta assembly.fa.gz -o mags --stitch --unbinned

# The clusters of a clustering round, every record in one file
metacooler fasta sample.metacool --fasta assembly.fa.gz -o final.fa \
    --round final --single-file
```

See the [`fasta` reference](../cli.md#fasta) for the order of the records and
the names they take.

## Exporting a network edge list (`network`)

Write a three-column TSV (`node1`, `node2`, `value`) with no header, suitable
for input to network analysis tools (Gephi, NetworkX, igraph):

```bash
# Genome-level network (one edge per genome pair)
metacooler network sample.metacool --level genome -o genome_network.tsv

# Contig-level network (one edge per contig pair)
metacooler network sample.metacool --level contig -o contig_network.tsv

# Raw bin-level network
metacooler network sample.metacool --level bin -o bin_network.tsv

# Overwrite an existing edge list (network no longer clobbers silently)
metacooler network sample.metacool --level genome -o genome_network.tsv --force
```

The `--level` values are singular: `bin`, `contig`, `genome`. `--chunksize`
(default `1000000`) controls how many pixels are folded into the running
aggregate at a time.

## Flattening to a standard cool file (`flatten`)

Remove the genome hierarchy entirely. All contigs become "chromosomes" in the
output cooler. Bins and pixels are preserved exactly as-is.

```bash
metacooler flatten sample.metacool -o flat.cool
```

This is useful for passing a metacool file to tools that expect a standard
cooler but where genome-level information is not needed.

## Visualising contact maps (`show`)

Render a PNG heatmap directly from the command line:

`show` takes its input as a **positional** path (like every other command) and
writes to the required `-o/--output`. `-f` is `--force`, not `--file`: pre-0.2
invocations of the form `metacooler show -f sample.metacool …` are no longer
valid (see the [`show` reference](../cli.md#show)).

```bash
# Genome-level (default)
metacooler show sample.metacool --level genome -o genome_heatmap.png

# Contig-level
metacooler show sample.metacool --level contig -o contig_heatmap.png

# Contig-level, restricted to one MAG (--genomes is honoured at every level)
metacooler show sample.metacool --level contig --genomes MAG_01 -o MAG_01_contigs.png

# Bin-level: intra-genome view
metacooler show sample.metacool --level bin --genomes MAG_01 -o MAG_01.png

# Any number of genomes: every name goes on BOTH axes, so you get all the intra-
# and inter-genome blocks between them, in the order you listed them
metacooler show sample.metacool --level contig --genomes MAG_01,MAG_02,MAG_03 -o panel.png

# The same rule for contigs — a handful of contigs of interest against each other
metacooler show sample.metacool --level bin --contigs contig_0_0,contig_3_7 -o contigs.png

# Bin-level: genomic range. The axes are ticked in bp, at their true position
metacooler show sample.metacool -r contig_0_0:0-100000 -o range.png

# -r/-r2 is the only asymmetric form: a different region per axis
metacooler show sample.metacool -r MAG_01 -r2 MAG_02 -o inter.png

# Custom colormap and transform, overwriting an existing image
metacooler show sample.metacool --cmap Reds --transform none -f -o raw.png
```

`--genomes` and `--contigs` are **symmetric** selectors — every name lands on
both axes — and they are mutually exclusive with each other and with `-r`. There
is no singular `--genome`: a one-name `--genomes` is the same thing.

Every figure is titled with the selection and the level it shows, and the axes
are labelled with the selection and its size. At `--level bin` the ticks are
base-pair coordinates rather than bin indices — a bin index starts at 0 even when
the requested range does not — except across a multi-contig selection, where they
name the contigs (or the selected genomes) instead.

A dense-matrix size guard applies at every level, including `--level bin` and
range queries: a selection too large to materialise is refused with an error
instead of attempting the allocation.

## Checking a file (`check`)

`check` reports on the integrity of a file instead of raising on the first
problem, printing one line per check and exiting non-zero if any failed:

```bash
metacooler check sample.metacool          # readable/format/structure/dtypes/cooler-compat
metacooler check --deep sample.metacool   # + a full O(nnz) pixel scan
metacooler check reps.mmetacool           # every entry in a collection
```

See the [`check` reference](../cli.md#check) for what each check establishes.

## Python API for export

Every export operation is also a top-level function — `metacooler.extract`,
`.write_fasta`, `.flatten`, `.merge`, `.binnify`, `.zoomify`, `.aggregate`, `.network`,
`.dump_table`, `.validate`, `.validate_mmetacooler`, and `.check` are all
exported from the package root (previously only the creation and collection
helpers were).

```python
import metacooler

mc = metacooler.MetaCooler("sample.metacool")

# Same operations as the CLI, straight from Python
metacooler.extract("sample.metacool", output="MAG_01.cool", genomes=["MAG_01"])
metacooler.write_fasta("sample.metacool", "assembly.fa.gz", "mags")
metacooler.flatten("sample.metacool", output="flat.cool")
metacooler.aggregate("sample.metacool", output="genome_level.cool", level="genome")
metacooler.network("sample.metacool", out="genome_network.tsv", level="genome")
metacooler.dump_table("sample.metacool", table="genomes", out="genomes.tsv")

# Dump tables as DataFrames
print(mc.genomes)
print(mc.contigs)
print(mc.bins())
print(mc.pixels())

# Filtered bins and pixels
bins_mag1 = mc.bins_for_genome("MAG_01")
pixels_pair = mc.pixels_for_genome_pair("MAG_01", "MAG_02")

# Dense bin-level contact matrix (cooler selector syntax)
mat = mc.matrix(balance=False)[:]

# One genome's sub-matrix, via its bin-index extent
start, end = mc.genome_extent("MAG_01")
mat_mag1 = mc.matrix(balance=False)[start:end, start:end]

# Aggregate to genome level → standard cooler, then read it back as a matrix
import cooler

mc.to_cooler("genome_level.cool", level="genome")
genome_matrix = cooler.Cooler("genome_level.cool").matrix(balance=False)[:]
```

## Merging replicates (`merge`)

If you have two `.metacool` files from the same assembly (e.g. two sequencing
replicates), merge them:

```bash
metacooler merge rep1.metacool rep2.metacool -o merged.metacool
```

Both files must have **identical** contigs and genomes tables, and the same
`contig-ordering` attribute. Pixels are concatenated and re-sorted. Metadata
from all inputs is merged, and their [QC blocks](../cli.md#merge-qc) are
summed when they were counted alike.
