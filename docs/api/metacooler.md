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

::: metacooler.api.MetaCooler

::: metacooler.schema
    options:
      heading: "Schema helpers"
      toc_label: "Schema helpers"

::: metacooler.uri
    options:
      heading: "URI"
      toc_label: "URI"
