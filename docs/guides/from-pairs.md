# Creating a metacool file from pairs

This guide covers all options for `metacooler create pairs` and its Python
equivalent `metacooler.from_pairs()`.

## Minimal example

```bash
metacooler create pairs \
    --pairs   contacts.pairs.gz \
    --bins    1000 \
    --contigs contigs2mag.tsv \
    --output  sample.metacool
```

`--bins 1000` generates 1 kb fixed-size bins automatically. Contig lengths
are read from the `#chromsize:` header lines in the pairs file. If they are
absent, provide `--fasta assembly.fa` instead.

## Input file formats

### Pairs file

Standard [4DN pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).
Plain text or gzipped. Required columns: `chrom1`, `pos1`, `chrom2`, `pos2`
(columns 2–5 by default).

```
#pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromsize: contig_0_0 22000
#chromsize: contig_0_1 13000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
read_0  contig_0_0  5211  contig_0_0  238  +  +
```

### Bins

Two modes:

=== "BED file"

    A three-column BED file (no header) with `chrom`, `start`, `end`. Generate
    with `cooler makebins` (fixed) or `cooler digest` (restriction fragments):

    ```bash
    cooler makebins assembly.fa.fai 2000 > bins_2kb.bed
    cooler digest assembly.fa.fai assembly.fa DpnII > fragments.bed
    ```

    ```bash
    metacooler create pairs --bins fragments.bed ...
    ```

=== "Integer bin size"

    Pass an integer (in bp) to generate fixed-size bins on the fly:

    ```bash
    metacooler create pairs --bins 1000 ...
    ```

    Contig lengths must be available from the pairs header or `--fasta`.

### Contig-to-genome table

Four formats are supported via `--contigs-format`:

=== "contigs (default)"

    TSV with header. First column must be `contig_name`, second `genome_name`.
    Additional columns are stored in the `/contigs` table.

    ```
    contig_name  genome_name  gc_content  coverage
    contig_0_0   MAG_01       0.47        71.70
    contig_0_1   MAG_01       0.59        10.00
    contig_1_0   MAG_02       0.45        87.98
    ```

    ```bash
    --contigs table.tsv --contigs-format contigs
    ```

=== "contig_genome"

    Two-column TSV, **no header**. Contig name in column 1, genome name in
    column 2.

    ```
    contig_0_0  MAG_01
    contig_0_1  MAG_01
    contig_1_0  MAG_02
    ```

    ```bash
    --contigs map.tsv --contigs-format contig_genome
    ```

=== "metator"

    MetaTOR `contig_data_final.txt` format. Expects `Name` and `Final_bin`
    columns; all other columns are ignored.

    ```bash
    --contigs contig_data_final.txt --contigs-format metator
    ```

=== "Column index (I,J)"

    Use column indices (1-based) when your table does not match the above
    layouts. Append `,noheader` if the file has no header row.

    ```bash
    # Column 1 = contig, column 5 = genome, with header
    --contigs-format "1,5"

    # Column 2 = contig, column 3 = genome, no header
    --contigs-format "2,3,noheader"
    ```

### Genomes metadata table (optional)

A TSV with header. First column must be `genome_name`. Any additional columns
are stored in the `/genomes` table.

```
genome_name  phylum         completeness  contamination  quality
MAG_01       Firmicutes     0.85          0.027          high
MAG_02       Bacteroidetes  0.74          0.011          medium
```

```bash
--genomes mag_quality.tsv
```

### File-level metadata (optional)

A JSON or YAML file with arbitrary metadata stored in the file's `metadata`
attribute:

```yaml
sample:
  name: gut_microbiome_T0
  organism: human gut metagenome
  enzyme: DpnII
sequencing:
  platform: Illumina NovaSeq 6000
  n_valid_pairs: 180000000
```

```bash
--metadata sample_meta.yaml
```

## Contig ordering

Control how contigs are ordered in the output with `--contig-order`:

```bash
# Default: longest contig first within each genome (good for visual clarity)
--contig-order genome-then-length

# Preserve order from FASTA / pairs header
--contig-order fasta

# Pure alphabetical
--contig-order alpha
```

## Handling unassigned contigs

By default, contigs absent from the contig table are assigned to a special
`__unbinned__` genome. Use `--drop-unassigned` to exclude them entirely:

```bash
metacooler create pairs ... --drop-unassigned
```

## Large files

Creation streams the pairs file in chunks, so memory stays bounded regardless of
file size. Two options tune this:

- `--chunksize` (default `1000000`) — pairs processed per chunk.
- `--sort-strategy` (default `auto`) — how pixels are sorted before writing.
  `auto` picks `memory` (sort in RAM, fastest) for files under 2 GiB and
  `external` (stream to cooler's on-disk sort) for larger ones. Force either
  explicitly if you know your machine's limits.

```bash
metacooler create pairs ... --sort-strategy external --chunksize 5000000
```

## Python API

```python
import metacooler

metacooler.from_pairs(
    uri="sample.metacool",
    pairs="contacts.pairs.gz",
    bins="fragments.bed",           # or an integer: binsize=1000
    contig_genome_map="contigs.tsv",
    contig_table_format="contigs",
    genomes_table="mag_table.tsv",  # optional
    assembly_name="gut_T0",         # optional
    contig_ordering="genome-then-length",
    drop_unassigned=False,
    metadata={"sample": "gut_T0"},  # optional
    chunksize=1_000_000,            # default; raise for very large pairs files
)
```

Pass a `dict[str, str]` directly as `contig_genome_map` to skip file parsing:

```python
cg_map = {"contig_A": "MAG_01", "contig_B": "MAG_02"}
metacooler.from_pairs(uri="out.metacool", contig_genome_map=cg_map, ...)
```
