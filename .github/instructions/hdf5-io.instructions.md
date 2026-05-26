---
description: "Use when reading or writing HDF5 files with h5py, cooler, or metacooler schema. Covers dtype handling, string encoding, schema validation."
applyTo: ["src/metacooler/io/**", "src/metacooler/schema.py", "src/metacooler/api.py"]
---
# HDF5 I/O Conventions

## String dtype handling
h5py cannot write pandas `object` or `StringDtype` arrays. Always convert:

```python
def _to_hdf5_data(arr):
    if arr.dtype == object or arr.dtype.kind in ("U", "S") or str(arr.dtype) == "string":
        str_arr = [str(x) for x in arr]
        max_len = max(len(s) for s in str_arr) if str_arr else 1
        return np.array([s.encode("utf-8") for s in str_arr], dtype=f"S{max_len}")
    return np.asarray(arr)
```

## Reading tables
Use `schema.read_table(uri, group)` which handles byte→string decoding automatically.

## Schema constants
- `METACOOLER_FORMAT = "MetaCooler::1.0"` — stored in `metacooler-format` attr
- `COOLER_FORMAT = "HDF5::Cooler"` — stored in `format` attr (cooler compat)
- `MMETACOOLER_FORMAT = "HDF5::MMetaCooler"` — stored in `format` attr for collections

## Validation
- `validate(uri)` — checks metacool structure
- `validate_mmetacool(uri)` — checks mmetacool structure
- `is_metacooler(uri)` / `is_mmetacooler(uri)` — quick format checks

## URI syntax
Use `uri.parse_uri(raw)` → `(Path, group_path)` for `file::group_path` addressing.

## Soft links
metacool files have `chroms → contigs` soft link for cooler compatibility. When iterating HDF5 groups in mmetacool files, wrap child access in try/except to handle broken soft links.

## Creating cooler files
Always use `cooler.create_cooler()` with `ordered=True, symmetric_upper=True`. After creation, add metacooler-specific groups/attrs with `h5py.File(path, "a")`.
