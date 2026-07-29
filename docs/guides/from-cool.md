# Importing from an existing cooler file

If you already have a `.cool` file and want to add genome / MAG annotations,
use `metacooler create cool`. This is useful when your Hi-C processing pipeline
produces a standard cooler output that you want to enrich with metagenomic
context after the fact.

## Basic usage

```bash
metacooler create cool \
    --cool    existing.cool \
    --contigs contigs2mag.tsv \
    --output  sample.metacool
```

Bins and pixels are carried over from the input cool file. The genome hierarchy
is added on top.

## With genome metadata

```bash
metacooler create cool \
    --cool    existing.cool \
    --contigs contig_data_final.txt --contigs-format metator \
    --genomes mag_quality.tsv \
    --output  sample.metacool
```

## Contig ordering

By default, contigs are reordered to `genome-then-length`. To preserve the
original contig order from the cool file:

```bash
--contig-order cool
```

`cool` is the only ordering under which a `cool → metacool → cool` round trip
preserves row order exactly; every other ordering regroups contigs by genome.
Note there is no `fasta` ordering here (and no `--fasta` option): contig lengths
and order come from the input cool file's own bins table.

A `--metadata` file may not contain a `commands` key — it is reserved for
metacooler's provenance record.

!!! note "Pixel remapping"
    When the output contig order differs from the input cooler's, metacooler
    remaps every bin ID in the pixel table. This is streamed in chunks with an
    `auto`/`memory`/`external` sort strategy (`--sort-strategy`) and can use
    several reader threads (`--nproc`), so large files are remapped without
    loading all pixels into memory.

## Python API

```python
import metacooler

metacooler.from_cool(
    cool_uri="existing.cool",
    output="sample.metacool",
    contig_genome_map="contigs.tsv",
    contig_table_format="contigs",
    genomes_table="mag_table.tsv",   # optional
    contig_ordering="genome-then-length",
    assembly_name="gut_T0",          # optional
    drop_unassigned=False,
)
```

## Differences from `create pairs`

| Aspect | `create pairs` | `create cool` |
|--------|---------------|---------------|
| Primary input | `.pairs` file | `.cool` file |
| Binning | Done during creation | Inherited from input |
| Pixel streaming | Chunked (`--chunksize`, `--sort-strategy`) | Chunked, with parallel remapping (`--nproc`, `--sort-strategy`) |
| Bin size flexibility | Fixed or fragment | Inherited |
