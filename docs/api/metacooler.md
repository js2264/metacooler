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
```

::: metacooler.api.MetaCooler

::: metacooler.schema
    options:
      heading: "Schema helpers"
      toc_label: "Schema helpers"

::: metacooler.uri
    options:
      heading: "URI"
      toc_label: "URI"
