# Changelog

## Unreleased

### Fixed
- `_write_rebinned_metacool`: wrong `chrom_offset` caused by calling
  `np.searchsorted` on an unsorted `contig_names` array. Replaced with an
  explicit dict-based mapping.
- `_write_rebinned_metacool`: missing `validate()` call after write; corrupt
  files were not caught before being returned to the caller.

### Changed
- `ops/aggregate.py` now delegates to `MetaCooler.aggregate_matrix()` instead
  of duplicating the accumulation logic.

---

## 0.1.0

Initial alpha release.

### Added
- `.metacool` HDF5 format: `genomes/`, `contigs/`, `bins/`, `pixels/`,
  `indexes/` tables with full cooler backward compatibility via `/chroms` soft
  link.
- `.mmetacool` collection format for grouping multiple metacooler entries in a
  single file.
- `MetaCooler` Python class: genome- and contig-aware accessors, aggregation,
  range selectors, and `to_cooler()` export.
- `from_pairs()`: streaming creation from 4DN pairs files with support for
  fixed-size and fragment-level bins.
- `from_cool()`: conversion from existing cooler files.
- `merge()`: merging replicate metacool files.
- CLI commands: `create pairs`, `create cool`, `merge`, `binnify`, `zoomify`,
  `extract`, `flatten`, `dump`, `network`, `aggregate`, `show`, `info`, `tree`,
  `attrs`, `metadata`, `ls`, `stack`, `cp`, `mv`, `rm`.
- Five contig ordering strategies: `genome-then-length`, `genome-then-alpha`,
  `fasta`, `alpha`, `custom`.
- Four contig table input formats: `contigs`, `contig_genome`, `metator`,
  column-index spec.
- `--drop-unassigned` flag for excluding unassigned contigs.
- YAML/JSON metadata ingestion at creation time and via `metadata` subcommand.
- GitHub Actions CI across Python 3.10–3.13.
- Pre-commit hooks: ruff, ty, trailing whitespace, YAML/TOML/JSON validation.
