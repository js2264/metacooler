# I/O functions

Functions that read and write metacooler files. All of the entry points below
are re-exported at the top level, so `metacooler.from_pairs(...)` and
`metacooler.io.create.from_pairs(...)` are equivalent:

| Top-level name | Defined in |
|----------------|------------|
| `from_pairs`, `from_cool` | `metacooler.io.create` |
| `collect` | `metacooler.io.collect` |
| `reassign` | `metacooler.io.reassign` |
| `annotate` | `metacooler.io.annotate` |
| `extract` | `metacooler.io.extract` |
| `flatten` | `metacooler.io.flatten` |
| `merge` | `metacooler.io.merge` |
| `network` | `metacooler.io.network` |
| `dump_table` | `metacooler.io.dump` |
| `create_mmetacool`, `add_entry`, `remove_entry`, `copy_entry`, `copy_entry_cross_file`, `move_entry`, `list_entries` | `metacooler.io.mmetacool` |
| `is_metacooler`, `is_mmetacooler`, `validate`, `validate_mmetacooler` | `metacooler.schema` |
| `parse_uri` | `metacooler.uri` |

`extract`, `flatten`, `merge`, `network`, `dump_table`, `validate`, and
`validate_mmetacooler` were previously reachable only through their submodules.

Note that `network` and `dump_table` name their output parameter `out`, not
`output`.

## In-place writes { #in-place }

`annotate(inplace=True)`, `reassign(inplace=True)`, `store_weights`,
`write_clusters`, `drop_clusters`, `contact_qc(store=True)`, the
`normalize_*` functions and the entry functions edit a file through a staged
copy, under its lock taken before their first read, and replace it once the
copy validates (see [In-place writes](../format.md#in-place-writes)). A caller
that computes from a map and writes it holds the lock from before its read to
after its last write, and hands the map it read to `store_weights`:

```python
import metacooler

with metacooler.lock("sample.metacool"):  # MapLockedError if another process or thread writes it
    mc = metacooler.MetaCooler("sample.metacool")
    weights = my_weights(mc)
    metacooler.store_weights(
        "sample.metacool", "w", weights, method="tool", depends_on="pixels", source=mc
    )
    metacooler.annotate("sample.metacool", contigs_table="depth.tsv", inplace=True)
```

A `MetaCooler` opened before a replacement refuses to read on
(`metacooler.MapChangedError`): open the map again. Both errors, and
`MapLockError`, are `OSError`s.

## A failed write { #keep-partial }

`from_pairs`, `from_cool`, `collect`, `reassign` and `annotate` (the last two
without `inplace`) build their output as `<output>.partial` and move it into
place only once it has validated, so a failed run never leaves a half-written
file at the name asked for, nor touches a file already there. `from_pairs`
stages its `.ps.tsv` table the same way and moves it into place just before the
map, the previous table set aside (as `<name>.ps.tsv.previous`) until the map is
in place and put back if the map cannot be moved. Once the map is in place the
set-aside table is removed; if that fails, it is left where it is, with a
warning, and not put back: the new map and table match, and the run succeeds
unless it was interrupted. Only a process killed during those moves can leave
them half done: the previous table still set aside, or the new table beside the
previous map. Either way, `from_pairs` refuses to start while a
`<name>.ps.tsv.previous` exists, before it reads anything: check which table
matches the map, keep it as `<name>.ps.tsv`, and remove the other (a file of
your own of that name is refused alike: rename it). `reassign`
measures the QC block before it writes anything. What becomes of
`<output>.partial` is `keep_partial`'s choice:

| `keep_partial` | When the write fails |
|----------------|----------------------|
| `True` (default) | The partial file stays, for inspection, and an ERROR line says where. Interrupted with Ctrl-C (`KeyboardInterrupt`), it stays without a word. |
| `False` | The partial file is removed whatever stopped the run, `KeyboardInterrupt` included, and an INFO line says so. For a caller that stages its own output and wants the finished file or nothing. |

A run that fails before its partial file exists has nothing to keep or remove,
and says nothing about one.

```python
import metacooler

try:
    metacooler.reassign(
        "sample.metacool", "bins.tsv", "sample.binned.metacool", keep_partial=False
    )
except Exception:
    ...  # neither sample.binned.metacool nor its .partial exists
```

::: metacooler.io.create
    options:
      heading: "Create metacool files"
      toc_label: "Create metacool files"

::: metacooler.io.collect
    options:
      heading: "Collecting single-genome cool files"
      toc_label: "Collecting single-genome cool files"
      members:
        - collect
        - collected_contig_name
        - CONTIG_NAME_SEPARATOR
        - COLLECT_ORDERINGS

::: metacooler.io.reassign.reassign
    options:
      heading: "Reassigning contigs to genomes"
      toc_label: "Reassigning contigs to genomes"

::: metacooler.io.annotate.annotate
    options:
      heading: "Annotating genomes and contigs"
      toc_label: "Annotating genomes and contigs"

::: metacooler.io.merge.merge
    options:
      heading: "Merging metacool files"
      toc_label: "Merging metacool files"

::: metacooler.io.mmetacool
    options:
      heading: "Metacooler collections"
      toc_label: "Metacooler collections"

::: metacooler.io.extract
    options:
      heading: "Extracting subsets"
      toc_label: "Extracting subsets"

::: metacooler.io.dump
    options:
      heading: "Dumping tables"
      toc_label: "Dumping tables"

::: metacooler.io.flatten
    options:
      heading: "Flattening to a cooler"
      toc_label: "Flattening to a cooler"

::: metacooler.io.network
    options:
      heading: "Exporting network"
      toc_label: "Exporting network"

::: metacooler.io.contigs
    options:
      heading: "Parsing contigs tables"
      toc_label: "Parsing contigs tables"

::: metacooler.io.genomes
    options:
      heading: "Parsing genomes tables"
      toc_label: "Parsing genomes tables"
