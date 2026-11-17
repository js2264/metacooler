# MetaCooler

The `MetaCooler` class is the primary Python interface to a `.metacool` file.
It wraps a `cooler.Cooler` instance internally and layers genome- and
contig-aware accessors on top.

## Opening a file

```python
import metacooler

# Single-sample file
mc = metacooler.MetaCooler("sample.metacool")

# Specific entry in a collection
mc = metacooler.MetaCooler("sample.mmetacool::entries/rep1")

# Fast path for very large files: skip the structural rescan that runs on
# every construction. Only for files you know to be well-formed — a corrupt
# file then surfaces later as an opaque cooler IndexError/KeyError, or as
# silently wrong slices if only the indexes are stale.
mc = metacooler.MetaCooler("huge.metacool", validate=False)
```

Passing a bare `.mmetacool` path raises `KeyError` listing the available
entries; opening a single map always requires an entry URI.

## Reading pixels in chunks

`pixels()` reads the whole pixel table at once. On a large map, iterate over
it instead: `iter_pixels()` yields consecutive slices of `/pixels`, at most
`chunksize` rows each (50 M by default), read with cooler's pixel selector, so
memory stays at one chunk whatever the size of the map.

```python
mc = metacooler.MetaCooler("sample.metacool")
for chunk in mc.iter_pixels(chunksize=10_000_000, contig_ids=True):
    # bin1_id, bin2_id, count, then contig1_id, contig2_id
    between = chunk[chunk["contig1_id"] != chunk["contig2_id"]]
    ...
```

`contig_ids=True` adds the contig of each of the two bins, as a row number of
`mc.contigs`, looked up through `indexes/chrom_offset` without reading the bin
table. Put together, the chunks are `pixels()` exactly, index included.

::: metacooler.api.MetaCooler

::: metacooler.schema
    options:
      heading: "Schema helpers"
      toc_label: "Schema helpers"

::: metacooler.uri
    options:
      heading: "URI"
      toc_label: "URI"
