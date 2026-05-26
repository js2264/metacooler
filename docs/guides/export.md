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

# All tables at once (writes sample.{genomes,contigs,bins,pixels}.tsv)
metacooler dump sample.metacool -t all --prefix sample
```

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
```

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

```bash
# Genome-level (default)
metacooler show sample.metacool --level genome --out genome_heatmap.png

# Contig-level
metacooler show sample.metacool --level contig --out contig_heatmap.png

# Bin-level: intra-genome view
metacooler show sample.metacool --level bin --genome MAG_01 --out MAG_01.png

# Bin-level: inter-genome view
metacooler show sample.metacool --level bin --genomes MAG_01 MAG_02 --out inter.png

# Bin-level: genomic range
metacooler show sample.metacool --range contig_0_0:0-100000 --out range.png

# Custom colormap and transform
metacooler show sample.metacool --cmap Reds --transform none --out raw.png
```

## Python API for export

```python
import metacooler

mc = metacooler.MetaCooler("sample.metacool")

# Dump tables as DataFrames
print(mc.genomes)
print(mc.contigs)
print(mc.bins())
print(mc.pixels())

# Filtered bins and pixels
bins_mag1 = mc.bins_for_genome("MAG_01")
pixels_pair = mc.pixels_for_genome_pair("MAG_01", "MAG_02")

# Sparse contact matrix
mat_genome = mc.aggregate_matrix(level="genome", sparse=True)
mat_contig = mc.aggregate_matrix(level="contig", sparse=False)  # dense numpy array

# Export to cooler
mc.to_cooler("genome_level.cool", level="genome")

# Extent helpers (bin index ranges)
start, end = mc.genome_extent("MAG_01")
mat_slice = mc.matrix(balance=False)[start:end, start:end]
```

## Merging replicates (`merge`)

If you have two `.metacool` files from the same assembly (e.g. two sequencing
replicates), merge them:

```bash
metacooler merge rep1.metacool rep2.metacool -o merged.metacool
```

Both files must have **identical** contigs and genomes tables, and the same
`contig-ordering` attribute. Pixels are concatenated and re-sorted. Metadata
from all inputs is merged.
