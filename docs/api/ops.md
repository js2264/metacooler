# Operations

Computation on already-open metacooler files: aggregation, re-binning,
multi-resolution generation, visualisation, and range parsing. These act on
files produced by the [I/O functions](io.md).

`aggregate`, `binnify`, and `zoomify` are re-exported at the top level, so
`metacooler.binnify(...)` and `metacooler.ops.binnify.binnify(...)` are
equivalent.

!!! warning "Resolution constraints"
    `binnify` requires the target bin size to be an **exact multiple** of a
    fixed source bin size (or, for a variable-bin source, strictly larger than
    the median source bin width), and `zoomify` requires every requested
    resolution to be a multiple of the base bin size. Both raise `ValueError`
    otherwise; see the [multi-resolution guide](../guides/multireso.md).

::: metacooler.ops.aggregate
    options:
      toc_label: "Aggregation"
      heading: "Aggregation"

::: metacooler.ops.binnify
    options:
      toc_label: "Re-binning"
      heading: "Re-binning"

::: metacooler.ops.zoomify
    options:
      toc_label: "Multi-resolution"
      heading: "Multi-resolution"

::: metacooler.ops.show
    options:
      toc_label: "Visualisation"
      heading: "Visualisation"

::: metacooler.ops.range
    options:
      toc_label: "Range resolution"
      heading: "Range resolution"
