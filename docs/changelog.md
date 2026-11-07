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
- **`show --genomes` is a comma-separated/repeatable list, not `nargs=2`.**
  `--genomes MAG_01,MAG_02` (or `--genomes MAG_01 --genomes MAG_02`) replaces
  `--genomes MAG_01 MAG_02`, bringing it in line with `--genomes` everywhere
  else. The old space-separated form now fails with a usage error, because the
  second name would fall through to the positional `METACOOL_PATH` slot. An
  explicitly empty `--genomes ''` is rejected, as on every other command.
- **`show --genomes` is symmetric for *any* number of names**, including two:
  every name goes on **both** axes, so the figure holds all the intra- *and*
  inter-genome blocks between them. `--genomes A,B` used to mean "A on the rows,
  B on the columns" — for that, use `-r A -r2 B`, now the only asymmetric form.
- **`show --genome` (singular) is gone.** `--genomes MAG_01` with one name is
  exactly the same thing, and the near-identical pair meant
  `--genome MAG_01,MAG_02` looked for a genome literally named
  `MAG_01,MAG_02`. `show` was also the only command in the CLI with a singular
  form.
- **Contig and genome names must be ASCII.** `create` refuses a non-ASCII name
  before it bins anything, `validate()` refuses one on read, and `check`
  reports it. `/contigs` is cooler's `/chroms`, and cooler decodes those columns
  with an ASCII `.astype("U")`, so an accented contig name — or an accented
  value in any *extra* `/contigs` column, which cooler reads just the same —
  produced a file that passed every metacooler check and then raised
  `UnicodeDecodeError` in `cooler.Cooler()`, and therefore in `MetaCooler()`.
  Genome names hit the mirror image of the same assumption on the way out,
  where `aggregate --level genome` writes them into a `.cool` as chrom names
  (`UnicodeEncodeError`). Extra `/genomes` metadata columns are unaffected and
  still take any UTF-8: cooler never reads that group. The limitation is
  upstream, but writing a file metacooler itself cannot open breaks the one
  invariant the format is built on, so it is now a named refusal instead of a
  crash three commands later.
- **A `--contigs` or `--genomes` table whose header and data rows disagree on
  the number of fields is rejected.** `pandas.read_csv` used to resolve the
  mismatch silently, and in opposite directions depending on which side was
  wider — a wider header pads the trailing columns, wider rows absorb the
  surplus leading fields into the row index — so in both cases every value was
  stored under the name of a neighbouring column. A trailing tab on the header
  line is still accepted, and no longer becomes an empty `Unnamed: N` column in
  `/genomes` or `/contigs`.

- **`ps` from a binned map counts opportunities in position pairs.** Each bin
  pair of widths `w_i`, `w_j` contributes `w_i × w_j` (the diagonal
  `w (w + 1) / 2`) — the unit the `--pairs` path already used — so
  `n_bin_pairs`, `avg_contacts` and `prob` change on the pixel path and now
  mean the same thing on both. Counting bin pairs weighed the diagonal as a
  full opportunity though it holds half as many position pairs, halving the
  first point. `prob` also no longer depends on `--bins-per-decade`: it was
  normalised over the full width of every log bin, though a binned map only
  realises separations at bin-distance steps, so the pixel/pairs ratio at
  10 kb drifted from 1.07 to 1.52 as the grid got finer. `dist` now sits on
  the separations a distance bin can actually hold (`[1, 2)` is plotted at 1,
  not 1.41).
- **QC counts a contact between two different `__unbinned__` contigs as
  `inter_genome`.** The unbinned bucket is not a genome, and treating it as
  one inflated `intra_genome_fraction` by every contact among unassigned
  contigs. The `metadata.qc` block is now `version: 2`.
- **`create cool` refuses a `square` cooler.** Both triangles were folded into
  the upper one and summed, doubling every off-diagonal count.
- **`create pairs` drops pairs at or past their contig's end** instead of
  folding them into the last bin (as `cooler cload` does); they are counted
  in the discard report and in `qc.n_pairs_discarded`.
- **Zero-width bins in `--bins`, and a `--drop-unassigned` that removes every
  contig, are refused before binning** rather than failing after the whole
  pairs pass.
- **`--metadata` files may not contain a `qc` key**, which `create` overwrote.
- **`validate()` is stricter:** it refuses a repeated `(bin1_id, bin2_id)`
  pixel, a `/chroms` link that does not resolve to the map's own `/contigs`,
  and a table member that is a group rather than a column.
- **`stack`, `create_mmetacool` and `add_entry` validate their sources**, so a
  plain `.cool` or a whole `.mmetacool` is refused, and entry names containing
  `/` are rejected (also by `cp`/`mv`).
- **An output that is one of the command's inputs is refused, even with
  `-f`.** `binnify a.metacool -o a.metacool -f` deleted the input, then failed.
- **Distance options must be finite, non-negative and at most 2^53 bp**
  (`qc --threshold`, `ps --min-distance`/`--max-distance`,
  `create --qc-threshold`); `inf` and `1e30` used to escape as tracebacks.

### Added

- **`collect`**: zip single-genome cool files — one per species, strain, or
  MAG — into one `.metacool` file, **one genome per input**, with no
  contig-to-genome table to write. Inputs are `.cool` files or one resolution
  of an `.mcool` each (`E_coli.mcool::resolutions/1000`); a genome is named
  after its file stem, or by `--names`. Every contig is renamed
  `<genome>__<contig>`, so the `chromosome` or `contig_1` that many assemblies
  share cannot collide — without stacking a second prefix on a contig that
  already carries its genome's, so `collect → extract → collect` is stable.
  `--genomes` decorates the genome rows as in `create`; `--contigs` attaches
  per-contig metadata, naming a contig as in its input or as collected, with an
  optional `genome_name` column that disambiguates a shared name and is checked
  rather than used to re-assign a contig; `--metadata`, `--filter` (including
  `genomes.n_bins`, since the bins come with the inputs), `--assembly-name`,
  and `--qc` work as in `create cool`. Inputs must share one bin size, or all be
  variable-binned; a metacool, an `.mcool` root, and a `square` cooler are
  refused with a message saying what to use instead. Every input occupies one
  contiguous block of bins, so under the default `--contig-order cool` (each
  input keeps its own contig order) the pixels are streamed into place **with no
  sort at all**, and `extract --genomes X` gives back input X bin for bin and
  pixel for pixel. `genome-then-length` and `genome-then-alpha` re-sort only
  the inputs whose contigs they move, one at a time, falling back to cooler's
  external sort past 100 M pixels. The inter-genome blocks of a collected map
  are empty by construction — no contact between two inputs was ever observed —
  which the docs spell out. Exposed as `metacooler.collect`. See
  [`collect`](cli.md#collect) and the new guide,
  [Collecting single-genome cools](guides/collect.md).
- **`qc`**: contact-quality statistics — where a library's contacts actually
  land. Every contact is classified into one of four categories
  (`intra_contig_short` / `intra_contig_long`, split at `--threshold`, default
  10 kb; `inter_contig_intra_genome`; `inter_genome`), from which the headline
  **3D ratio** (`trans_contig_ratio`, the share of contacts joining two
  different contigs) and the intra-/inter-MAG breakdown follow. `create pairs`
  and `create cool` now compute the block during creation at essentially no
  cost — they are already streaming every pair or pixel — and store it in
  `metadata.qc` plus five per-genome `qc_*` columns of `/genomes`, including
  `qc_intra_genome_ratio`, a per-MAG "how self-contained is this bin?" score.
  Pass `--no-qc` to skip it. The standalone `metacooler qc` command recomputes
  the same block from any existing file's pixel table (`--table`,
  `--per-genome FILE.tsv`, `--store`). The two sources agree on every total;
  only the intra-contig short/long split differs, because binning discards the
  sub-bin distances — which is why `create pairs` measures it on raw read
  positions and records `source="pairs"`. `binnify` and `merge` drop the block
  and the columns rather than carry stale numbers onto a different pixel table,
  and `merge` ignores `qc_*` when comparing input `/genomes` tables, so two
  libraries over one assembly still merge. See [`qc`](cli.md#qc).
- **`ps`**: P(s), contact probability as a function of genomic distance, as a
  TSV (`-o`) and/or a log-log figure (`--plot`), pooled or one curve per MAG
  (`--by genome`), with the usual `--genomes`/`--contigs`/`--filter`
  selection. Columns are `genome`, `dist`, `n_contacts`, `n_bin_pairs`,
  `avg_contacts`, `prob`. Each log-spaced distance bin reports the contacts
  observed *and* the number of opportunities that could have produced one, so
  the decay is not confounded by the shrinking number of long-range pairs a
  contig offers, and a normalised `prob` column that integrates to 1 makes
  libraries of different depth comparable. Only intra-contig contacts take
  part — the ones a genomic distance is defined for. The bin-pair denominator
  is counted with a single global `searchsorted` per distance edge over a
  contig-folded key, and the position-pair one in closed form per contig, so
  neither costs a per-contig Python loop on a metagenome. See
  [`ps`](cli.md#ps).
- **P(s) at base-pair resolution, measured during creation.** A curve read off
  a binned map cannot resolve anything below one bin, which is precisely where
  a metaHi-C library shows its undigested fragments and self-religation. So
  `create pairs` now accumulates a **per-MAG** curve from the raw read
  positions during the pass it already makes over the pairs, and writes it
  beside the output as `<output>.ps.tsv` (`--no-ps` opts out,
  `--ps-bins-per-decade` sets the grid). No plot is drawn. Measured cost on
  4 M pairs: none outside run-to-run noise — one `searchsorted` and one
  `bincount` per chunk over the intra-contig rows, and a closed-form
  denominator computed once.
- **`ps --pairs`** recomputes that fine-grained curve for a file created
  earlier, with a different grid or selection. The pairs file goes through the
  same reader `create pairs` uses, so compression and column handling are
  identical by construction, and the output is byte-for-byte the table
  `create` wrote. The map still supplies the contig→genome assignment, the
  contig lengths, and the selection.
- **`--filter`**, a repeatable predicate selector on `create pairs`,
  `create cool`, `extract`, `dump`, and `show`, for choosing genomes and contigs
  by *property* rather than by name: `--filter "genomes.length > 1Mb"`,
  `--filter "contigs.length >= 10kb"`, `--filter "genomes.Completeness > 90"`.
  One comparison per flag — `> >= < <= == != =~ !~`, with `bp/kb/Mb/Gb`
  suffixes on numbers — and repeating the flag ANDs them. The `genomes.` /
  `contigs.` prefix is mandatory, since both tables have a `length` and a
  `name`. Any column of either table can be used, including the ones a
  `--contigs`/`--genomes` metadata table added, plus a derived
  `contigs.genome` holding the genome name. Names choose the candidates and
  filters narrow them, so `--genomes MAG_01 --filter "contigs.length > 50kb"`
  is an intersection, not a union. Because the column names in real MetaTOR
  tables contain spaces and dashes (`Present Markers`, `GC-content`), the first
  operator terminates the column name and no quoting or backticking is needed;
  the cost is that there is no `or`. Selecting nothing is an error that reports
  what each filter would have kept on its own, and an unknown column suggests
  the closest one. In `create` the filter runs before anything is tiled, so a
  dropped contig costs no bin, no pixel, and no `cooler.create` overhead — see
  [Filtering with `--filter`](cli.md#filter).
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

- **cooler ≥ 0.11 no longer brings back the quadratic per-contig write
  cost.** cooler moved `write_bins` into `cooler.create._writer`; metacooler's
  shims looked in the old module and silently did nothing, so on cooler 0.12
  creating an 80 k-contig map took 13.9 s instead of 2.1 s, growing about 3×
  per doubling of contigs.
- **`merge`** ran cooler's writer without those shims (65 s → 3.1 s at 150 k
  contigs), floored float counts, lost `bin-size` when every contig was a
  single bin, crashed on `.mmetacool` entry URIs, and — in the CLI — put back
  the stale `metadata.qc` block it had just removed.
- **`create --force` keeps the previous map (and its `.ps.tsv`) until the new
  one has been written and validated**; it used to delete them before any
  work, so a failed run left nothing.
- A contig table with its own `genome_id` column no longer sends every contig
  to `__unbinned__`, and the same collision in `--genomes` no longer crashes.
- A `/` in a `--contigs`/`--genomes` column name (e.g. `reads/kb`) is stored
  as `_`; it used to create a nested HDF5 group that neither cooler nor
  `MetaCooler()` could open.
- `--contigs-format I,J,noheader` with extra columns no longer crashes after
  the whole pairs pass (they are stored as `col3`, `col4`, …).
- `.zst` pairs (and `.bz2` with an integer `--bins`) are read; the header
  scans used their own openers and failed.
- An empty pixel stream no longer crashes `create` under
  `--sort-strategy external`, nor `binnify`, `zoomify` or `aggregate` on a
  zero-pixel map.
- Fragment-level (`--bins file.bed`) binning is one vectorised lookup instead
  of a loop over contigs (13.1 s → 4.7 s at 80 k contigs; identical output).
- `aggregate --level genome` refuses a genome longer than 2^31 bp instead of
  writing a negative chromosome length.
- `show` keeps float counts at contig/genome level, and
  `show --transform log10` is `log10(1 + count)` as documented (single-count
  pixels used to render as empty).
- `--filter` `!=` / `!~` no longer match a missing text value.
- `extract --split` strips only a trailing `.cool` from `-o`, so
  `-o sample.rep1` and `-o sample.rep2` no longer collide.
- `create cool` keeps the source's `genome-assembly`.
- `MetaCooler.root` is the entry path for an `.mmetacool` entry, `info["sum"]`
  keeps a fractional sum, and `MetaCooler` has a `repr`.
- A `G:start-end` range on a `binnify --level genome` file (where every contig
  is also a genome) selects the contig instead of being refused.
- `qc --store` no longer crashes on a non-object `qc` metadata value.

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

### Packaging and CI

- The MIT license text ships as `LICENSE` (PEP 639 metadata).
- Corrected lower bounds: `pyfaidx>=0.7.2.2` (older releases import
  `pkg_resources`) and `rich-click>=1.9` (older releases cannot import the
  CLI). Dropped the unused `bioframe` dependency and `plotly` extra.
- The development toolchain is a `dev` dependency group, so `uv sync` installs
  it; as an extra, the documented `uv sync --dev` installed nothing.
- Python 3.13 and 3.14 are supported and tested. A `dependency-bounds` CI job
  runs the suite at the lowest and the newest versions every dependency range
  allows, and `tests/data/cli_tests.sh` runs in CI.

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
