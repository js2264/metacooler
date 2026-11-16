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
absent, provide `--fasta assembly.fa` instead. Several libraries of the same
assembly go in one map by repeating `--pairs` — see
[Several pairs files](#several-pairs-files).

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

### Several pairs files

A metagenome is often sequenced as several Hi-C libraries — replicates, runs,
or one library per restriction enzyme — each aligned to the same assembly.
Repeat `--pairs`, once per file:

```bash
metacooler create pairs \
    --pairs   lib_DpnII.pairs.gz \
    --pairs   lib_HinfI.pairs.gz \
    --bins    1000 \
    --contigs contigs2mag.tsv \
    --output  sample.metacool
```

The files are read back to back in one pass, as if concatenated, so the map,
the contact-quality block, the P(s) table and the dropped-pairs report all
cover every library exactly. That is what sets this apart from creating one map
per library and [merging](../cli.md#merge) them, which sums the same contacts
but drops the QC block and yields no P(s) table.

Their `#chromsize:` headers are united; a file may lack one, but a contig given
two different lengths by two files is an error, as is the same file given
twice. See [Several libraries](../cli.md#create-pairs-libraries) in the CLI
reference for the details.

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

`--contigs-format` says which two columns of the table hold the contig and the
genome names: `CONTIG,GENOME`, each a **1-based position** or a **header
name**, plus `,noheader` when the table has no header line. Every other column
is stored in the `/contigs` table. Aliases cover the tables of common binners
(see [Contig and genome columns](../cli.md#contigs-format) for the full list).

=== "contigs (default)"

    `contig_name,genome_name`: a TSV whose header names those two columns,
    wherever they are.

    ```
    contig_name  genome_name  gc_content  coverage
    contig_0_0   MAG_01       0.47        71.70
    contig_0_1   MAG_01       0.59        10.00
    contig_1_0   MAG_02       0.45        87.98
    ```

    ```bash
    --contigs table.tsv
    ```

=== "By header name"

    Any table, by the names of its two columns — here MetaTOR's
    `contig_data_final.txt` (the `metator` alias says the same):

    ```bash
    --contigs contig_data_final.txt --contigs-format Name,Final_bin
    ```

=== "By position"

    Column positions (1-based), for a table whose header you would rather not
    spell — or that has none:

    ```bash
    # Column 1 = contig, column 5 = genome, with header
    --contigs-format "1,5"

    # Column 2 = contig, column 3 = genome, no header
    --contigs-format "2,3,noheader"
    ```

=== "No binning yet"

    `NA` as the genome column: the table only lists contigs, and every contig
    lands in `__unbinned__`. Bin them later with
    [`metacooler reassign`](../cli.md#reassign).

    ```bash
    cut -f1 assembly.fa.fai > contigs.txt
    --contigs contigs.txt --contigs-format 1,NA,noheader
    ```

=== "Binner aliases"

    ```bash
    --contigs contig_data_final.txt     --contigs-format metator  # MetaTOR
    --contigs run_DASTool_contigs2bin.tsv --contigs-format dastool  # DAS Tool
    --contigs final_contig_to_bin.tsv   --contigs-format binette  # Binette
    --contigs contig_bins.tsv           --contigs-format semibin  # SemiBin2
    --contigs vae_clusters_unsplit.tsv  --contigs-format vamb     # VAMB
    --contigs map.tsv                   --contigs-format contig_genome  # 2 columns, no header
    ```

### Genomes metadata table (optional)

A TSV whose genome-name column is the one `--genomes-format` names (by
default the column called `genome_name`, wherever it is) — a 1-based position
or a header name, plus `,noheader` when the table has no header line. Any
additional columns are stored in the `/genomes` table. A QC tool's report works
as it is: `--genomes-format Name` for CheckM2's `quality_report.tsv`,
`--genomes-format user_genome` for GTDB-Tk's summary.

```
genome_name  phylum         completeness  contamination  quality
MAG_01       Firmicutes     0.85          0.027          high
MAG_02       Bacteroidetes  0.74          0.011          medium
```

```bash
--genomes mag_quality.tsv
```

Metadata that only exists once the map does — a QC run on the bins, coverage
from another library — can be added later with
[`metacooler annotate`](../cli.md#annotate), which reads these same tables.

!!! important "The contigs table decides which genomes exist"

    `/genomes` has one row per distinct genome name in the **contigs** table.
    `--genomes` only *decorates* those rows with extra columns, left-joined on
    the genome name — it neither adds nor removes genomes. So:

    - A genome named in the contigs table but missing from `--genomes` still
      gets a row, with `length`/`n_contigs`/`n_bins` computed from its contigs
      and empty metadata columns. This is common with MetaTOR output, where
      `Final_bin` also holds `MetaTOR_MGE_*` bins that a MAG-only genome table
      does not list.
    - A genome listed in `--genomes` with no contigs assigned to it has nothing
      to derive a row from, and is absent from the output.

    Both cases are reported as warnings naming the genomes involved, so a
    `/genomes` table larger or smaller than the input can be accounted for.

Both this table and the contigs table must have exactly as many fields on every
data row as their header has names. A header that declares more columns than
the rows carry is **rejected**: pandas would pad the trailing columns and store
every value under the name of a column to its left (turning, say, a genome's GC
content into its `Size`), and there is no way to tell from the file whether the
absent fields belong at the start or the end of each row. A trailing tab on the
header line is tolerated and ignored.

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

!!! warning "`commands` is reserved"
    A `--metadata` file may **not** contain a top-level `commands` key — that
    name holds metacooler's own provenance record (the exact command line used
    to produce each file). Supplying it is an error; the file used to overwrite
    the provenance silently.

## Contig ordering

Control how contigs are ordered in the output with `--contig-order`:

```bash
# Default: longest contig first within each genome (good for visual clarity)
--contig-order genome-then-length

# Group by genome, alphabetical within each genome
--contig-order genome-then-alpha

# Preserve order from FASTA / pairs header
--contig-order fasta

# Pure alphabetical
--contig-order alpha
```

!!! warning "`fasta` requires a FASTA or a `#chromsize:` pairs header"
    `--contig-order fasta` means "preserve the order of the input assembly", so
    there has to be one: pass `--fasta assembly.fa`, or use a pairs file whose
    header carries `#chromsize:` lines. Without either, the command now fails
    instead of quietly falling back to alphabetical order while still recording
    `contig-ordering="fasta"` in the file — a mismatch that made two files with
    different layouts claim the same ordering, which `merge` would then
    scramble.

Only `--contig-order cool` (available in [`create cool`](from-cool.md)) round-trips
the input row order exactly; the default regroups contigs by genome.

## Handling unassigned contigs

By default, contigs absent from the contig table are assigned to a special
`__unbinned__` genome, and so are contigs whose genome reads `NA`, `N/A`,
`nan`, `ND`, `n/d`, `-1` (SemiBin2's label) or is empty. Use
`--drop-unassigned` to exclude them entirely:

```bash
metacooler create pairs ... --drop-unassigned
```

A table that names *no* contig of the assembly is refused rather than sending
every contig to `__unbinned__`: that is the wrong table, or the wrong
`--contigs-format`. Rows naming contigs the assembly does not have are ignored
with a warning. See [The contig table](../cli.md#create-contig-table).

## Large files

Creation streams the pairs file in chunks, so memory stays bounded regardless of
file size. Two options tune this:

- `--chunksize` (default `50000000`) — pairs processed per chunk. Under
  `external` sorting each chunk is written as a temporary cooler, and that write
  carries a fixed cost proportional to the contig count, so small chunks are
  expensive on metagenomic assemblies.
- `--sort-strategy` (default `auto`) — how pixels are sorted before writing.
  `auto` buffers the binned pixels and picks `memory` (sort in RAM, fastest —
  the cooler is written in a single pass) unless the **deduplicated** pixel
  count exceeds 100 M, at which point it falls back to `external` (cooler's
  on-disk sort). Force either explicitly if you know your machine's limits.

!!! note "Why deduplicated, and why it matters at depth"

    The buffer is periodically re-aggregated as it fills, so the threshold is
    compared against the pixel table actually being built rather than the sum of
    per-chunk row counts. The two diverge sharply on deep libraries: Hi-C pixel
    counts *saturate* on a fixed bin grid, so sequencing 20x deeper yields well
    under 2x more distinct pixels. A 1 B-pair library on a typical metagenomic
    grid can stream a few hundred million raw rows that collapse to under 10 M
    pixels — roughly 250 MB, and easily sorted in RAM. Counting raw rows would
    have spilled that to the external sort and paid one `cooler.create` call per
    chunk for nothing.

    Peak memory for the in-memory path is roughly `3 x nnz x 24` bytes, so ~7 GB
    at the 100 M-pixel limit. To estimate `nnz` for your own data before
    committing, create from a ~5% subsample and read off the reported pixel
    count — saturation means the full run will not be far above it.

```bash
metacooler create pairs ... --sort-strategy external --chunksize 5000000
```

## Python API

```python
import metacooler

metacooler.from_pairs(
    uri="sample.metacool",
    pairs="contacts.pairs.gz",      # or a list: one pairs file per library
    bins="fragments.bed",           # or an integer: binsize=1000
    contig_genome_map="contigs.tsv",
    contig_table_format="contigs",
    genomes_table="mag_table.tsv",  # optional
    assembly_name="gut_T0",         # optional
    contig_ordering="genome-then-length",
    drop_unassigned=False,
    metadata={"sample": "gut_T0"},  # optional
    chunksize=50_000_000,           # default; lower it to cap peak memory
)
```

Pass a `dict[str, str]` directly as `contig_genome_map` to skip file parsing:

```python
cg_map = {"contig_A": "MAG_01", "contig_B": "MAG_02"}
metacooler.from_pairs(uri="out.metacool", contig_genome_map=cg_map, ...)
```
