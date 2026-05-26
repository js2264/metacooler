---
name: metacool-format
description: "Reference for the metacool and mmetacool HDF5 file formats. Use when: understanding the data model, debugging schema issues, implementing new I/O operations, working with cooler compatibility."
user-invocable: true
disable-model-invocation: false
---
# metacool / mmetacool File Format Reference

## When to Use
- Understanding the HDF5 structure of `.metacool` files
- Debugging schema validation errors
- Implementing new I/O operations that read/write metacool files
- Understanding cooler backward compatibility

## metacool Format (single file)

### HDF5 Structure
```
/ (root)
├── @attrs:
│   ├── format: "HDF5::Cooler"              # cooler compat
│   ├── metacooler-format: "MetaCooler::1.0" # metacool identity
│   ├── format-version: 3
│   ├── bin-type: "fixed" | "variable"
│   ├── bin-size: <int> | 0
│   ├── contig-ordering: "mag-then-alpha" | ...
│   ├── nbins: <int>
│   ├── nnz: <int>
│   └── metadata: "{}"                       # JSON string
├── genomes/
│   ├── name:      S<N> array
│   ├── length:    int64 array
│   └── n_contigs: int32 array
├── contigs/
│   ├── name:      S<N> array
│   ├── length:    int32 array
│   └── genome_id: int32 array  (index into genomes)
├── bins/
│   ├── chrom:  int32 array  (index into contigs)
│   ├── start:  int32 array
│   └── end:    int32 array
├── pixels/
│   ├── bin1_id: int64 array
│   ├── bin2_id: int64 array
│   └── count:   int32 array
├── indexes/
│   ├── chrom_offset:  int64[n_contigs+1]
│   ├── bin1_offset:   int64[n_bins+1]
│   └── genome_offset: int64[n_genomes+1]  (optional)
└── chroms/ → soft link to contigs/
```

### Invariants
- `genome_id` values in contigs must be in `[0, n_genomes)`
- `chrom` values in bins must be in `[0, n_contigs)`
- Pixels are upper-triangular: `bin1_id ≤ bin2_id`
- `chrom_offset[i]` = first bin index for contig `i`
- `chroms` is always a soft link to `contigs` (cooler compatibility)

## mmetacool Format (collection)

### HDF5 Structure
```
/ (root)
├── @attrs:
│   ├── format: "HDF5::MMetaCooler"
│   └── format-version: 1
└── groups/
    ├── <entry1>/   (complete metacool)
    ├── <entry2>/   (complete metacool)
    └── ...
```

### URI Addressing
```
file.mmetacool::groups/entry1    → opens entry1 as MetaCooler
file.mmetacool::                 → root (for ls, attrs)
file.metacool                    → plain metacool (group path = "/")
```

## Python API
```python
from metacooler.schema import is_metacooler, is_mmetacooler, validate, validate_mmetacool
from metacooler.uri import parse_uri

path, group = parse_uri("file.mmetacool::groups/entry1")
mc = MetaCooler("file.mmetacool::groups/entry1")  # URI syntax supported
```

## Common Pitfalls
1. **String dtypes**: h5py cannot write pandas `object`/`StringDtype` — convert to `S<N>` bytes first
2. **Soft links in mmetacool**: entries copied from metacool files contain `chroms→contigs` soft links that may break when traversing with `h5py` — wrap in try/except
3. **Pixel ordering**: always enforce upper-triangle after any bin remapping
4. **Index consistency**: after filtering bins/pixels, rebuild `chrom_offset` and `bin1_offset`
