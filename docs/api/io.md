# I/O functions

Functions that read and write metacooler files. All of the entry points below
are re-exported at the top level, so `metacooler.from_pairs(...)` and
`metacooler.io.create.from_pairs(...)` are equivalent:

| Top-level name | Defined in |
|----------------|------------|
| `from_pairs`, `from_cool` | `metacooler.io.create` |
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

::: metacooler.io.create
    options:
      heading: "Create metacool files"
      toc_label: "Create metacool files"

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
