# Concepts

This page explains the **mental model** behind metacooler: what problem it
solves, how it organises data, and the vocabulary used throughout the rest of
the documentation. Read it once and the commands and API will make sense.

For the exact byte-level layout, see the [format specification](format.md).

## The problem

A [cooler](https://github.com/open2c/cooler) file stores a Hi-C contact map as
three tables — `chroms`, `bins`, and `pixels`. That model assumes a single
genome whose chromosomes are known in advance.

Metagenomic Hi-C breaks that assumption. Binning tools such as
[MetaTOR](https://github.com/koszullab/metator) and
[bin3C](https://github.com/cerebis/bin3C) reconstruct dozens or hundreds of
**metagenome-assembled genomes (MAGs)** at once, each made of many contigs. In a
plain cooler file every contig looks like an independent "chromosome", and the
grouping of contigs into genomes — the thing you actually care about — is lost.

metacooler adds that missing layer while keeping the file a valid cooler.

## The data model

metacooler layers a **genome → contig** hierarchy on top of cooler's
`bins → pixels` model. Four tables describe the data, from coarsest to finest:

```
genomes   one row per MAG            name, length, n_contigs, n_bins  (+ metadata)
   │  one genome owns many contigs
contigs   one row per contig         name, length, genome_id          (+ metadata)
   │  one contig is tiled into many bins
bins      one row per genomic bin    chrom (= contig id), start, end
   │  contacts are counted between pairs of bins
pixels    one row per contact        bin1_id, bin2_id, count
```

- **Genomes** are the MAGs. The reserved genome `__unbinned__` (always row 0)
  collects contigs that were never assigned to a MAG.
- **Contigs** are the assembled sequences. Each carries a `genome_id` pointing
  back into the genomes table — that link is the heart of the format.
- **Bins** tile each contig into fixed- or variable-width windows, exactly like
  cooler. The `chrom` column holds an integer contig index (not a name).
- **Pixels** are the contacts: for each observed bin pair, how many read pairs
  connected them. Stored upper-triangular (`bin1_id ≤ bin2_id`) and sorted.

Any table also has room for **extra metadata columns** — genome completeness or
taxonomy, contig GC content or coverage — carried straight through from your
input tables.

!!! note "Genome, MAG, bin — a word of caution"
    In metagenomics a "bin" can mean a MAG. In cooler a "bin" is a genomic
    window. **In this documentation, "bin" always means the cooler sense** (a
    window on a contig) and "genome" or "MAG" means the reconstructed organism.

## Levels: bin, contig, genome

Because the hierarchy has three tiers, most analysis commands accept a
`--level` that decides how contacts are aggregated before you look at them:

| Level | One row / cell represents | Typical use |
|-------|---------------------------|-------------|
| `bin` | a single genomic window (finest) | intra-genome structure, zoomed views |
| `contig` | one contig collapsed to a point | contig–contig association |
| `genome` | one whole MAG collapsed to a point | genome–genome interaction network |

The same idea drives [`show`](cli.md#show), [`network`](cli.md#network),
[`aggregate`](cli.md#aggregate), and [`dump`](cli.md#dump). Aggregating never
changes the file on disk — it is computed on read.

## Two file types

<div class="grid cards" markdown>

- **`.metacool`** — one contact map

    A single metagenomic Hi-C map: a cooler plus the `genomes/` hierarchy.
    This is what `create` produces and what most commands consume.

- **`.mmetacool`** — a collection

    A container holding several named metacool entries in one file — for
    replicates, samples, or resolution levels. Not itself a cooler.

</div>

Entries in a collection are addressed with a cooler-style **URI**:

```
collection.mmetacool::entries/<name>
```

Every command and the Python `MetaCooler(...)` constructor accept such a URI
wherever they accept a path. [`zoomify`](cli.md#zoomify), for example, writes an
`.mmetacool` whose entries are named after their bin size (`1000`, `2000`, …).

## Cooler compatibility

The central invariant of the format:

!!! success "Every `.metacool` file is also a valid cooler file"
    `cooler.Cooler("sample.metacool")` opens it with no conversion, and every
    cooler-based tool (cooltools, HiGlass, the cooler CLI) works at **contig
    resolution**. The genome layer is simply invisible to them.

This works because metacooler writes, at creation time:

- `format = "HDF5::Cooler"` so cooler's own `is_cooler()` check passes;
- a soft link `/chroms → /contigs`, so cooler finds its chromosome table where
  it expects it (the contigs *are* the chromosomes);
- the standard `indexes/chrom_offset`.

Alongside these it writes the metacooler-only markers
(`metacooler-format = "HDF5::MetaCooler"` and the `/genomes` group) that make it
a metacooler. The two identities coexist in one file. See
[Cooler compatibility detail](format.md#cooler-compatibility-detail).

## Contig ordering

The order of contigs in the file is chosen **once, at creation**, and recorded
in the `contig-ordering` attribute. It determines the layout of `bins` and
`pixels` and therefore how a heatmap reads visually.

| Strategy | Order of contigs |
|----------|------------------|
| `genome-then-length` | grouped by genome, longest contig first — **default** |
| `genome-then-alpha` | grouped by genome, alphabetical within each genome |
| `fasta` | as they appear in the FASTA / pairs header |
| `alpha` | all contigs alphabetical, ignoring genome |
| `cool` | preserve the input file's order (`create cool`, and the default of `collect`, where each input keeps its own order within its genome) |

Grouping by genome keeps each MAG's contigs contiguous, so a genome appears as a
solid block on the diagonal of a contact map. Because the ordering is baked in,
[`merge`](cli.md#merge) requires all inputs to share the same one.

## Where to go next

<div class="grid cards" markdown>

- :material-rocket-launch: **[Quickstart](quickstart.md)** — raw data to a heatmap in five commands
- :material-book-open-variant: **[Guides](guides/from-pairs.md)** — task-by-task walkthroughs
- :material-console: **[CLI reference](cli.md)** — every command and option
- :material-file-tree: **[Format specification](format.md)** — the on-disk layout

</div>
