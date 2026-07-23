# I/O functions

Functions that read and write metacooler files. The creation entry points
`from_pairs` and `from_cool`, and the collection helpers (`create_mmetacool`,
`add_entry`, `remove_entry`, `copy_entry`, `move_entry`, `list_entries`), are
re-exported at the top level — so `metacooler.from_pairs(...)` and
`metacooler.io.create.from_pairs(...)` are equivalent.

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
