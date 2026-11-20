# Collecting single-genome cool files

If you already have one `.cool` file per genome — isolates of a mock community,
reference strains, maps downloaded one species at a time — `metacooler collect`
zips them into a single `.metacool` file, **one genome per input**. There is no
contig-to-genome table to write: the file a contig comes from *is* its genome.

## Basic usage

```bash
metacooler collect E_coli.cool B_subtilis.cool S_aureus.cool -o species.metacool
```

```console
$ metacooler dump -t genomes -H species.metacool | cut -f1-4
genome      length   n_contigs  n_bins
B_subtilis  4215000  1          4215
E_coli      4705000  2          4705
S_aureus    2849000  2          2849
```

Genome rows are always alphabetical: that order defines the `genome_id` integers
of the format, whatever order the inputs were given in.

## Naming

Each genome is named after its input's file stem (`E_coli.cool` → `E_coli`).
Give explicit names with `--names`, one per input and in the same order:

```bash
metacooler collect run1/sample.cool run2/sample.cool \
    --names Escherichia_coli,Bacillus_subtilis \
    -o species.metacool
```

Every contig is then renamed `<genome>__<contig>`. Assemblies routinely reuse
the same contig names — `chromosome`, `plasmid1`, `contig_1`, `1` — but a
metacool file needs them unique, so the prefix is applied to every contig of
every input, whether or not it clashes. That keeps each genome's names the same
whichever other genomes it is collected with:

```console
$ metacooler dump -t contigs -H species.metacool | cut -f1-3
contig                  genome      length
B_subtilis__chromosome  B_subtilis  4215000
E_coli__chromosome      E_coli      4640000
E_coli__plasmid1        E_coli      65000
S_aureus__chromosome    S_aureus    2821000
S_aureus__plasmid1      S_aureus    28000
```

The new names are what every other command expects: `--contigs
E_coli__plasmid1`, `-r E_coli__chromosome:1000000-2000000`, `--filter
"contigs.contig =~ plasmid"`.

## Resolution

A metacool file has a single bin table, so **every input must share one bin
size** — or all have variable, restriction-fragment bins. Inputs binned
differently are refused with the resolution of each. From multi-resolution
files, pick the same resolution from each with cooler's URI syntax:

```bash
metacooler collect \
    E_coli.mcool::resolutions/5000 B_subtilis.mcool::resolutions/5000 \
    -o species_5kb.metacool
```

## Adding metadata

Three optional inputs decorate the collection:

- `--genomes`: one row per genome, first column `genome_name`, matched against
  the genome names; every other column (taxonomy, completeness, growth
  condition…) lands in `/genomes`.
- `--contigs`: one row per contig, first column `contig_name`; every other
  column lands in `/contigs`. A contig may be named as in its input cool
  (`chromosome`) or as collected (`E_coli__chromosome`), and a `genome_name`
  column says which genome a shared name belongs to. See
  [The contigs table](../cli.md#collect-contigs-table).
- `--metadata`: a JSON or YAML file of free-form metadata.

```bash
metacooler collect *.cool \
    --genomes species.tsv \
    --contigs contig_coverage.tsv \
    --metadata study.yaml \
    -o species.metacool
```

`--filter` works as in `create`, before anything is written, and can use any
column those tables add — and, unlike in `create`, `genomes.n_bins`:

```bash
metacooler collect *.cool --genomes species.tsv \
    --filter "genomes.completeness > 90" --filter "contigs.length >= 5kb" \
    -o species.metacool
```

## What the collected map is — and is not

Every input was binned on its own, usually from its own library, so no contact
*between* two inputs was ever observed. The collected map has **empty
inter-genome blocks**: it is a block-diagonal stack of the inputs, not the map
you would get by aligning one metagenomic library to the combined assembly (for
that, see [Creating from pairs](from-pairs.md)). Each genome also keeps its own
sequencing depth, so compare genomes with that in mind.

Within those limits, nothing is lost. Under the default `--contig-order cool`,
contigs keep each input's own order and extracting a genome gives back its
input, bin for bin and pixel for pixel — only the contig names carry the
prefix:

```bash
metacooler extract species.metacool --genomes E_coli -o E_coli_again.cool
```

That layout is also what makes `collect` fast: each input is one contiguous
block of bins, its pixels only move within that block, and under `cool` the
whole collection is streamed into place without sorting anything.
`--contig-order genome-then-length` (longest contig first, as `create` does by
default) and `genome-then-alpha` re-sort only the inputs whose contigs they
actually move, one at a time.

## Python API

```python
import metacooler

metacooler.collect(
    ["E_coli.cool", "B_subtilis.mcool::resolutions/5000"],
    "species.metacool",
    names=["E_coli", "B_subtilis"],        # optional; defaults to file stems
    genomes_table="species.tsv",           # optional
    contigs_table="contig_coverage.tsv",   # optional
    contig_ordering="cool",                # or "genome-then-length", ...
)
```

## Compared with the other creation paths

| Aspect | `create pairs` | `create cool` | `collect` |
|--------|----------------|---------------|-----------|
| Input | one `.pairs` file | one `.cool` file | one `.cool` per genome |
| Genome assignment | contig table | contig table | the input file |
| Contig names | unchanged | unchanged | `<genome>__<contig>` |
| Inter-genome contacts | as observed | as observed | none |
| Default contig order | `genome-then-length` | `genome-then-length` | `cool` (each input's own) |
