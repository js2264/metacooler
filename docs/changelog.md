# Changelog

## Unreleased

### Breaking

- **`show` takes a positional input path.** The required `-f/--file` option is
  gone; `metacooler show FILE.metacool -o out.png` replaces
  `metacooler show -f FILE.metacool -o out.png`. `-f` now means `--force`, as in
  every other command.
- **`binnify` rejects a target bin size that is not an exact multiple of the
  source bin size.** A non-multiple target produced a systematic ±14 %
  checkerboard while advertising a uniform bin size. For variable-bin sources
  the median source bin width is used as the effective resolution and a target
  at or below it is refused.
- **`zoomify` rejects a resolution ladder that is not divisible by the base bin
  size**, and stops early instead of emitting duplicate levels.
- **`dump` errors instead of silently ignoring incompatible combinations**:
  `--weight`/`--join` with a non-pixels table, `--genomes`/`--contigs` with
  `-t indexes`, `-o` with `-t all`, `--prefix` with a non-`all` table.
- **An explicitly empty `--genomes ''` / `--contigs ''` is rejected** rather
  than silently disabling the filter.
- **A bare `.mmetacool` path is rejected** by every command that needs a single
  map, with a message listing the entry URIs. `ls`, `tree`, and `attrs` still
  accept a collection.
- **`metadata add/get/remove/load` on a collection root is rejected**; they used
  to write to the root, invisibly from every entry.
- **`create pairs --contig-order fasta` requires a FASTA or a `#chromsize:`
  pairs header.** It used to fall back to alphabetical order while still
  recording `contig-ordering="fasta"`.
- **`--metadata` files may not contain a `commands` key** (reserved for the
  provenance record).
- **The assembly name is stored in the standard `genome-assembly` HDF5
  attribute**, not `assembly-name`. `info` now reports it, so `assembly` reads
  `"unknown"` where it previously read `null`.
- **Declining a `cp`/`mv`/`rm` confirmation prompt exits non-zero** (it used to
  exit 0).
- **`attrs` rejects a group given twice** (URI group *and* `GROUP_PATH`); the
  URI group used to be silently discarded.
- **A `--contigs` or `--genomes` table whose header and data rows disagree on
  the number of fields is rejected.** `pandas.read_csv` used to resolve the
  mismatch silently, and in opposite directions depending on which side was
  wider — a wider header pads the trailing columns, wider rows absorb the
  surplus leading fields into the row index — so in both cases every value was
  stored under the name of a neighbouring column. A trailing tab on the header
  line is still accepted, and no longer becomes an empty `Unnamed: N` column in
  `/genomes` or `/contigs`.

### Added

- **`check`**: an integrity report over the existing validators. Prints one line
  per check — `readable`, `format`, `structure`, `dtypes`, `cooler-compat` —
  instead of raising on the first problem, so every failure is visible in one
  pass, and repeats the per-map checks for each entry of an `.mmetacool`. Two
  checks go beyond what runs on every file open: the on-disk dtype *kind* of the
  16 core columns (`validate` casts on read, so a float64 `bins/start` used to
  pass it and fail later inside cooler), and, under `--deep`, a full O(nnz) pixel
  scan that also reconciles `attrs['sum']` with the stored counts. Exposed as
  `metacooler.check` / `metacooler.CheckResult`.
- **`show --genomes` accepts any number of genomes**, with one uniform rule:
  every name goes on both axes, so N names give all N against all N — intra and
  inter blocks alike — laid out in the order the names were given, with a
  heavier boundary line between blocks.
- **`show --contigs`**, the same selector in contig space: bare contig names,
  comma-separated or repeated, on both axes. It closes the one gap `-r` could
  not — more than one contig on the same axis — and upgrades `--level genome`
  (where a contig has no row of its own) to `--level contig`. Mutually
  exclusive with `--genomes`.
- **`show` figures are titled and their axes labelled.** The title names the
  selection and the level (`MAG_01 vs MAG_02 — contig level`), truncating a long
  genome list to `(+N more)`; axis labels name the selection and its size
  (`MAG_01 (12 contigs)`).
- **`show --level bin` ticks the axes in base pairs**, not bin indices. A bin
  index starts at 0 even when the selection does not, so `-r contig:20000-50000`
  used to read as though it began at the start of the contig. A selection
  spanning several contigs has no single coordinate system, so there the ticks
  name the contigs, one at the first bin of each.
- `-q/--quiet` on every command: suppresses INFO logging, keeps warnings and
  errors.
- `-f/--force` on `dump`, `network`, and `show`, which previously clobbered
  their output silently.
- `binnify --chunksize` to bound peak memory during pixel remapping.
- `zoomify --temp-dir` for the intermediate per-level files; peak scratch use is
  now one level rather than a second copy of the pyramid.
- `MetaCooler(..., validate=False)` fast path for very large, known-good files.
- Top-level exports: `metacooler.extract`, `.flatten`, `.merge`, `.binnify`,
  `.zoomify`, `.aggregate`, `.network`, `.dump_table`, `.validate`, and
  `.validate_mmetacooler`.
- `stack` accepts entry URIs and rejects duplicate `--names`.

### Fixed

- **`dump … | head` no longer reports a failure.** Closing the pipe early is the
  normal end of a pipeline, but the resulting `BrokenPipeError` is an `OSError`
  and so was rendered in the same red `Error: [Errno 32] Broken pipe` box as a
  bad input — *after* the rows had already been printed, which read as though the
  dump had failed. Every command now exits silently on a closed stdout with
  status `141` (`128 + SIGPIPE`), the same as `yes | head`, and no longer leaks an
  `Exception ignored in: …` from the interpreter's final flush.
- Answering `y` to a `cp`/`mv` overwrite prompt now actually overwrites; `-f`
  used to be the only way through.
- `show --genome`/`--genomes` are honoured at every `--level`, and the
  dense-matrix size guard applies at `--level bin` too.
- Counts are preserved as `int64` through `extract`, `flatten`, and `merge`, so
  re-binned sources no longer truncate.
- `.bgz` pairs files are decompressed correctly instead of being parsed as raw
  binary into an empty file; pairs dropped for an unknown contig or a missing
  position are now counted and reported.
- `mv` renames the HDF5 link in place instead of copying and deleting.
- A `--bins` BED with more than three columns (BED4/BED6) is read correctly
  instead of shifting `chrom`/`start`/`end` one column to the left. Because
  contig lengths are derived from the bins frame itself, the shifted table was
  self-consistent and creation went ahead with contigs named after coordinates.
- `create` reports when the `--contigs` and `--genomes` tables cover different
  sets of genomes: names assigned contigs but absent from the genomes table
  (kept, with derived `length`/`n_contigs`/`n_bins` and empty metadata) and
  names listed in the genomes table with no contigs (absent from `/genomes`).
  Both used to be silent, which made a `/genomes` table larger than the input
  look as though genomes had been invented.

### Documentation

- `COMMANDS.md` removed; [`docs/cli.md`](cli.md) is the single CLI reference.
- `_write_rebinned_metacool`: wrong `chrom_offset` caused by calling
  `np.searchsorted` on an unsorted `contig_names` array. Replaced with an
  explicit contig-name → index mapping.
- `_write_rebinned_metacool`: missing `validate()` call after write; corrupt
  files were not caught before being returned to the caller.

### Changed
- `ops/aggregate.py` centralises the accumulation logic in module-level
  helpers (`_aggregate_sparse` / `_aggregate_pixel_stream`) instead of
  duplicating it across call sites.

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
- Four contig ordering strategies for `create pairs` — `genome-then-length`,
  `genome-then-alpha`, `fasta`, `alpha` — plus `cool` for `create cool`.
- Four contig table input formats: `contigs`, `contig_genome`, `metator`,
  column-index spec.
- `--drop-unassigned` flag for excluding unassigned contigs.
- YAML/JSON metadata ingestion at creation time and via `metadata` subcommand.
- GitHub Actions CI across Python 3.10–3.13.
- Pre-commit hooks: ruff, ty, trailing whitespace, YAML/TOML/JSON validation.
