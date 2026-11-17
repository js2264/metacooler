# Changelog

## Unreleased

### Added

- **Per-contig QC columns.** Every command that stores the contact-quality
  block — `create pairs`, `create cool`, `collect`, `qc --store`, `reassign` —
  now also writes `qc_intra_contig_short` and `qc_intra_contig_long` to
  `/contigs` (`int64`, or `float64` on a map with float counts), from the same
  tally: read positions for `create pairs`, bin midpoints for a pixel tally.
  They are the two categories that do not depend on the genome assignment,
  counted contig by contig; the `/genomes` columns of the same names are their
  per-genome sums. The block's `version` stays `2`: no existing key changes
  meaning. Like the `/genomes` ones, the names are reserved (a user table
  column of that name is stored as `user_qc_intra_contig_short`), `binnify` and
  `merge` drop them, and `annotate` cannot drop or replace them. Python:
  `ContactQC(..., n_contigs=)`, `ContactQC.contig_table()`,
  `ContactQC.set_intra_contig()`, `metacooler.ops.qc.tally_map()` (the pixel
  pass of `qc`, with the bins' genomes as an option), and a `per_contig`
  argument to `store_qc`.
- **`annotate --unmatched warn|ignore|error`** (API: `annotate(...,
  unmatched=)`): what to do with the rows of a genomes or contigs table that
  name no genome or contig of the map. `warn`, the default, keeps the warning
  `annotate` always gave; `ignore` says nothing about them (a debug line);
  `error` refuses the tables before anything is written. A table that matches
  nothing at all is still refused in every case.
- **`keep_partial`** on every writer that stages its output — `from_pairs`,
  `from_cool`, `collect`, `reassign` and `annotate`. When the write fails,
  `keep_partial=False` removes `<output>.partial`, whatever stopped the run
  (`KeyboardInterrupt` included), and logs that at INFO, where the default,
  `True`, leaves it for inspection and logs at ERROR as before. For callers
  that stage their own output; see [A failed write](api/io.md#keep-partial).
  `from_pairs` now stages its `.ps.tsv` table with the map and moves it into
  place just before it, the previous table set aside until the map is in place
  and put back if the map cannot be moved: written after the map was in place, a
  failure there left the new map behind. Once the map is in place the set-aside
  table is removed, and one that cannot be is left, with a warning, not put
  back: the new map and table match, and the run succeeds unless it was
  interrupted. A `<name>.ps.tsv.previous` left behind, by that failure or by a
  run killed during the moves, makes `from_pairs` refuse to start, before it
  reads anything, until someone checks which table matches the map (a file
  of the user's own of that name is refused alike: rename it).
  `reassign` and `annotate` are now in the I/O API reference.
- **`MetaCooler.iter_pixels(chunksize=, contig_ids=False)`** yields the pixel
  table in consecutive chunks of at most `chunksize` rows, read with cooler's
  selector, so a caller can stream a large map without `pixels()` reading all
  of it at once. `contig_ids=True` adds `contig1_id`/`contig2_id`, looked up
  through `indexes/chrom_offset`. See
  [Reading pixels in chunks](api/metacooler.md#reading-pixels-in-chunks).

### Changed

- **`reassign` keeps the QC block instead of recomputing it from the pixels.**
  The per-contig intra-contig counts are carried over and summed per new
  genome, only the two inter-contig categories are re-measured from the pixels
  (exact at any bin size), and every other key — `source`,
  `distance_threshold`, `binsize`, `n_pairs_read`, `n_pairs_discarded` — is
  kept. So `create pairs` with a first table (or none), then `reassign` with a
  second, now stores the very block, `/genomes` and `/contigs` `qc_*` columns
  that `create pairs` with the second table gives. It used to recompute the
  block from bin midpoints: on a map at one bin per contig that turned every
  intra-contig contact into `intra_contig_short`, and it lost
  `n_pairs_read`/`n_pairs_discarded` and `source="pairs"` — on every map
  binned after creation, which is how MetaTOR works. A source without the
  per-contig columns (written by 0.2.0 or earlier), or whose columns no longer
  add up to its block, is still recomputed from the pixels, with a warning
  naming what is lost when the block was measured on read positions (an INFO
  line otherwise). `--no-qc` still drops the block, and its two `/contigs`
  columns with it. The new block is measured from the source's pixels before
  anything is written, and written with the rest of the map, so a run stopped
  meanwhile — `--inplace` included — leaves the source's block and per-contig
  counts as they were; it used to be stored after the relabelling.
- **In `/contigs`, only the two QC columns go with the QC block.** `merge`
  drops them (it copied the first input's `/contigs` verbatim), as `binnify`
  does, and `reassign` replaces or drops them with the block, whether it moves
  the contigs or not. A user column that merely starts with `qc_` (`qc_pass`,
  from a `--contigs` table) is contig metadata, kept by every writer:
  `binnify` and a re-grouping `reassign` used to drop every `qc_*` column of
  `/contigs`. `annotate --drop` now accepts it, and `merge` reports it when
  inputs differ, as for any metadata column; both treated every `qc_*` column
  of `/contigs` as metacooler's.
- **`merge` refuses the same map given twice.** `merge A A` summed every
  contact of `A` with itself and exited 0. Inputs are now identified by file
  (device and inode) and group, so an alias or a hard link of a listed input
  is refused too; two different entries of one `.mmetacool` still merge.
  `create pairs`' refusal of a pairs file given twice is keyed the same way, so
  a hard link no longer slips past it.
- **Faster genome and name selections on large maps.** On a map whose genomes
  are not stored contiguously (`--contig-order alpha` or `fasta`, or a
  `reassign --layout keep`), each of a genome's bin ranges was read with
  cooler's selector naming every contig of the map: at 200 k contigs,
  `bins_for_genome` of a 97-contig genome took 7.7 s (now 0.3 s) and
  `extract --genomes` 15 s (now 1 s). Lists of names (`dump`/`extract`
  `--contigs` and `--genomes`, `show`) scanned the contigs table once per name:
  5 000 contigs took 7 s, now 0.16 s. Output is unchanged.

### Fixed

- **`reassign` no longer crashes on a non-object `qc` metadata value** when it
  recomputes the block; it measures at the default distance threshold, as
  `qc --store` does without `--threshold`.
- **Type checkers see what every top-level name is.** `metacooler.from_pairs`,
  `metacooler.MetaCooler` and the rest of `__all__` are loaded lazily through
  the module's `__getattr__`, so pyright and ty typed each of them as `object`
  and flagged every call in client code as not callable, while mypy typed each
  as `Any` (`import-untyped`). An `if TYPE_CHECKING:` block in
  `metacooler/__init__.py` now names each one at its source for type checkers
  only: nothing more is imported at run time, so
  `import metacooler` costs what it did. The package also ships a `py.typed`
  marker, without which mypy does not read an installed metacooler at all
  (`import-untyped`).
- **`reassign` no longer claims that only provenance differs when no contig
  moves** although the output differs in more: a `--genomes` table, a
  re-grouped layout or a new `contig-ordering`, `--no-qc`, `--metadata` keys,
  or the QC block it recomputes. The message names each, and is an INFO line
  rather than a warning when the caller asked for one of them; a table that
  moves nothing, and nothing else asked for, still warns.
- **A left-over `<output>.partial` given as the input is no longer deleted.**
  Every writer that stages its output removes a stale `<output>.partial`
  before it starts, and it did so even when that file was its input:
  `reassign out.metacool.partial -o out.metacool`, rescuing a partial run,
  deleted the map before reading it (as did `annotate`, `collect`,
  `create cool` and `create pairs` given such an input). It is now refused
  before anything is read, as an output that is an input always was.
- **A `"` in a user table no longer swallows rows.** Genomes and contigs
  tables were read with CSV quote handling, so `"draft` on one row and
  `curated"` a few rows later became one field spanning every row between them:
  those rows vanished without a warning (in `create`, their contigs went to
  `__unbinned__`). Every user table, and the bins BED, is now read with quoting
  off; a quote is stored verbatim.
- **`network` no longer wraps large integer weights.** Edge weights were summed
  in the stored `int32`, so a node pair with more than 2³¹ contacts (the
  `__unbinned__` self-edge of a deep library) printed a negative weight. They
  are summed in `int64`.
- **`collect -f` can no longer delete one of its own inputs**, nor the previous
  output when the run then fails: `collect A.cool B.cool --genomes g.tsv -o
  g.tsv -f` deleted `g.tsv` before reading it. Every command now names its
  inputs to the output guard (checked by a test over every `cmd_*.py`).
- **A user's `qc_*` genome column survives.** Every `qc_*` column of `/genomes`
  is metacooler's, but only the five exact QC names were renamed on the way
  in, so a `qc_pass` from a `--genomes` table was stored as is and then dropped
  by the next `reassign`, `binnify` or `merge`. It is now stored as
  `user_qc_pass`, with the usual warning.
- **A column named `.` or `..`** (which HDF5 cannot store) is stored as `colN`
  with a warning, instead of failing after `create`'s whole pairs pass.
- **`annotate --inplace` is all or nothing.** It swapped each rebuilt table in
  as soon as it was built, so a failure on `/contigs` left `/genomes` already
  edited, no history recorded and the scratch group behind. Every table is now
  built before any is swapped in, and a failure leaves the map as it was.
- **Pairs files: upper-case suffixes and blank header lines.** The codec was
  chosen by pyarrow's case-sensitive suffix match, so `R1.PAIRS.GZ` was parsed
  as raw gzip bytes into an empty map; suffixes are now matched
  case-insensitively, and `.zstd` is read as `.zst`. A blank line inside the
  `#` header made the read fail on the first pair (and hid the `#chromsize`
  lines after it); blank lines are now part of the header.
- **One contigs table gives one column type, whichever command stores it.**
  `create` types a metadata column over every row of the table, so one `N/A` —
  even on a row naming a contig the map lacks — keeps it text; `collect` and
  `annotate` typed only the matched rows and stored the same column as
  numbers, and `reassign` re-typed a stored text column the same way. The
  table is now typed once, as read, by every command.
- **`collect` accepts an input whose every contig is a single bin** (a phage, a
  few short plasmids). cooler infers no bin size for it and calls it variable,
  but its bins are the grid of any bin size at least as long as its longest
  contig, so it joins such a collection instead of being refused.
- **Docs:** `metacooler validate` (which never existed) is `metacooler check`;
  `contact_qc`, `contact_probability` and `plot_ps` are in the API reference;
  the `reassign` and `annotate` history blocks of the `metadata` attribute are
  in the format specification; `network` and `aggregate` `--chunksize`
  defaults read 50 000 000, as they are.

## 0.2.0 - 2026-11-16

### Breaking

- **The default `--contigs-format` (`contigs`) finds `contig_name` and
  `genome_name` by name, wherever they are**, instead of requiring them to be
  columns 1 and 2; likewise `--genomes`' `genome_name` and `collect --contigs`'
  `contig_name`.
- **`--contigs-format contig_genome` keeps a table's extra columns**, as `col3`,
  `col4`, …, like every other headerless spec; it used to drop them.
- **A genome named `-1` in a contig table means `__unbinned__`**, as `NA` and
  `ND` already did: it is SemiBin2's label for unbinned contigs.
- **`merge` compares only what makes two maps the same map** — contig names and
  lengths, the contig-to-genome assignment, genome sizes. Metadata columns may
  now differ between inputs: the output keeps the first input's, and each
  difference is a warning. Integer width no longer matters either (`int32`
  lengths from `create cool` against `int64` ones).
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
- **`create pairs` and `create cool` refuse a `--contigs` table that names no
  contig of the assembly**, or has no row, before anything is binned. Every
  contig used to land in `__unbinned__`, and the command exited 0 with nothing
  logged above INFO — the signature of the wrong table, or of the wrong
  `--contigs-format` (a leading ID column read as `1,2` takes the IDs for
  contig names). The error shows a few names from each side and points at
  `--contigs-format`. A table naming contigs the assembly does not have is
  still accepted, now with a warning naming up to 10 of them. See
  [The contig table](cli.md#create-contig-table).
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

- **`metacooler annotate`** adds, replaces or drops metadata columns of
  `/genomes` and `/contigs` in an existing map (API `metacooler.annotate`):
  CheckM2, GTDB-Tk or depth tables that arrive after the map was built. Rows
  are matched by name; the assignment, bins and pixels are never touched, so it
  is cheap enough to run `--inplace`. Columns get the names, values and types
  the same table would have given `create`. A changed column needs
  `--replace` (`--missing keep|na` chooses between a patch and an authoritative
  table); derived columns cannot be dropped or replaced; a contig table's
  `genome_name` is checked, never applied; non-ASCII `/contigs` values are
  refused. A bare `.mmetacool` path annotates every entry.
- **One grammar for table columns.** `--contigs-format CONTIG,GENOME[,noheader]`
  takes 1-based positions or header names (`Name,Final_bin`, `2,5`), and
  aliases for common binners: `metator`, `dastool`, `binette`, `semibin`,
  `vamb` (besides `contigs` and `contig_genome`). `collect` gets a
  `--contigs-format` of its own, for its per-contig metadata table (see its
  section).
- **A map with no binning yet**: `--contigs` is optional on `create pairs` and
  `create cool`. Without it the map is the cooler plus the genome layer, every
  contig in `__unbinned__`, ready for `reassign` (`create cool --contig-order
  cool` then keeps the cool file's contigs, bins and pixels unchanged). To keep
  a table's other columns as contig metadata meanwhile, `--contigs-format
  CONTIG,NA` says the table has no genome column; with `reassign`, the same
  spec unbins the contigs the table lists.
- **`--genomes-format`** on `create pairs`, `create cool`, `collect`, `reassign`
  and `annotate`: which column of `--genomes` holds the genome names, so a
  CheckM2 (`Name`) or GTDB-Tk (`user_genome`) report needs no editing.
- **`create pairs` reads several pairs files in one pass** — one per Hi-C
  library of the same assembly: repeat `--pairs` (not comma-split, since a path
  may contain a comma), or pass a list as `from_pairs(pairs=[...])`; a single
  path works as before. The files are streamed back to back in the given order,
  so the map, the QC block and `/genomes` `qc_*` columns, the `.ps.tsv` P(s)
  table and the dropped-pairs report (`n_pairs_read`, `n_pairs_discarded`) are
  exactly those of the files' concatenation — which per-library maps joined by
  `merge` cannot give, since `merge` drops the QC block. Chunk folding, the sort
  strategy, `--nproc` (output identical to `--nproc 1`) and memory behaviour are
  unchanged; a chunk may span two files. The files' `#chromsize:` headers are
  united in file order; a file without one still works with `--fasta` or
  another file's header, and one contig given two different lengths is an error
  naming both files. The same file given twice, and a missing file, are refused
  before any pairs file is read. The `commands` provenance lists every
  `--pairs`.

- **`--nproc` on `create pairs` and `zoomify`** (default 1, which keeps the
  serial behaviour). `create pairs --nproc N` bins pairs batches on N threads
  while the file is read and decompressed; per-batch drop tallies, QC and P(s)
  are folded in file order, so every output is identical to `--nproc 1`
  (`_PairsDropStats.unknown_names` keeps insertion order for that). On a
  200 k-contig, 30 M-pair library: 23.9 s → 20.3 s with 2 threads, ~20 s with 4
  — one gzip stream is then the critical path. `zoomify --nproc N` builds
  independent levels in N processes, each as soon as the level it is coarsened
  from exists, and appends them in ladder order: `-r 1000,2000,5000,10000` on the
  same map goes from 24.9 s to 15.6 s with 2 processes (peak memory per
  process unchanged, so the total roughly doubles).

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
- **`reassign`**: re-attribute the contigs of an existing map to genomes from a
  new contig-to-genome table — a re-binning, a curation pass, a merge of two
  MAGs — without re-binning a pair. Every count is kept; `/genomes` is rebuilt
  from the new assignment (genomes created and emptied, every `genome_id`
  renumbered, `length`/`n_contigs`/`n_bins` recomputed). Under a genome-grouped
  `contig-ordering` a moved contig breaks the layout `genome_offset` and
  `genome_extent` rely on, so `--layout auto` re-groups such a file (bins
  permuted, pixels remapped and re-sorted, balancing weights carried) and only
  relabels any other, skipping the rewrite when the new grouping reproduces the
  old order; `--layout keep` never moves a row and records
  `contig-ordering="custom"`, a new value `validate()` does not hold to a
  grouping. A relabel copies the file and rewrites four datasets, so it can also
  run `--inplace`. `--missing keep|unbinned` makes the table an incremental edit
  or authoritative; neither drops a contig. The QC block is recomputed (or, with
  `--no-qc`, dropped), per-genome metadata is carried with a warning naming the
  genomes whose contig set changed, and a table matching no contig is refused.
  A `square` map is only ever relabelled, and a `collect`ed map's genomes,
  contiguous although `cool` does not promise it, are named in a warning if the
  new assignment scatters them. Exposed as `metacooler.reassign`. See
  [`reassign`](cli.md#reassign).
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

- **`metadata add` and `metadata remove` refuse the keys metacooler writes
  itself** — `commands`, `qc`, and the `reassign` and `annotate` histories.
  `metadata add -k commands --value '[]'` wiped the provenance record, and
  `metadata remove -k commands` deleted it, both exiting 0, although a
  `--metadata` file was already refused `commands` and `qc`. `--metadata` files
  and `metadata load` now refuse the two history keys as well.
- **`create pairs --sort-strategy external` works on large inputs with cooler
  < 0.12.** Its pixel chunks reached cooler's unordered writer unsorted, and
  cooler before 0.12 indexes each chunk as if it were sorted. Once the merge
  took more than one `mergebuf` epoch (20 M records), the output had
  out-of-order and duplicated pixels, and `validate()` rejected it at the very
  end of the run (`bin1_offset disagrees with pixels/bin1_id`), as it did on
  30 M pairs over 200 k contigs in the benchmarks. `create pairs` past its in-memory
  budget under `auto`, `create cool` with the external sort, and `reassign`'s
  regroup path streamed unsorted chunks in the same way. Every unordered write
  now sorts each chunk that is not sorted already (one vectorised check for
  those that are).
- **`MetaCooler.pixels(join=True)` names the contigs.** `chrom1` and `chrom2`
  held contig *ids* (`0`, `1`, …): `/bins/chrom` is a plain integer column in a
  metacool file, and cooler's join reports such a column as ids where
  `cooler.Cooler.bins()` converts it to names. The coordinates are now joined
  from `MetaCooler.bins()`. Pass `convert_enum=False` to get the ids.
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

- **The CLI reference gives editing its own section.** "Create and edit" is
  split into *Create* (`create`, `collect`), *Edit an existing map*
  (`reassign`, `annotate` and `metadata`, moved from *Inspect*, with a table of
  which one to use) and *Combine and re-bin* (`merge`, `binnify`, `zoomify`);
  `metacooler --help` groups the commands the same way.
- **The documentation site is restyled after doc.chromosomix.app.** The
  zensical theme keeps its layout but takes a new header, an in-column page
  footer with previous/next cards (also reached with the ← / → keys), a
  light/dark palette, framed code blocks and Starlight-style asides and tables
  (`overrides/`, `docs/stylesheets/metacooler.css`,
  `docs/javascripts/metacooler.js`). The home page is now the introduction and
  four entry tiles; its motivation, quick look, compatibility table and
  documentation map moved to a new [Why metacooler?](why.md) page, the on-disk
  layout to [Concepts](concepts.md#on-disk), and the `qc`/`ps` commands to the
  [Quickstart](quickstart.md#step-2-inspect-the-file).
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

### Performance

- **Genome and contig retrieval goes through the offset indexes.** A genome
  resolves to its contigs through `indexes/genome_offset` (or, under a
  `fasta`/`alpha`/`cool`/`custom` ordering, `/contigs/genome_id`), a contig to
  its bins through `indexes/chrom_offset`, and a bin range to its pixel rows
  through `indexes/bin1_offset`; those rows are then read with cooler's own bin
  and pixel selectors, and nothing scans the bin table or the pixels outside
  them. `pixels_for_genome_pair`, `genome_bin_ids`, `genome_extent`,
  `contig_extent`, `bins_for_genome`, `bins_for_contig`, `dump -r/-r2` (long
  and wide), `dump -t pixels --genomes/--contigs/--filter` and `extract` all use
  it. On a 200 k-contig map, `dump -t pixels -r <contig>` went from 3.0 s to
  0.8 s (cooler: 1.0 s) and `pixels_for_genome_pair(g, g)` from 119 ms to
  3.5 ms. `genome_extent` still raises for a genome that is not contiguous.
  `extract` of a `symmetric-upper` map now writes in order, without cooler's
  merge. New public `MetaCooler.genome_contig_ids()` and
  `MetaCooler.genome_bin_ranges()`.
- **Opening a file runs only the structural checks.** `MetaCooler()` used to run
  the whole of `schema.validate()`, O(n_bins), on every open (2.8 s on a
  2 M-contig map). It now checks format attributes, required groups and
  columns, column lengths, the `/chroms` link, the count attributes and the
  length and endpoints of every offset index — O(n_columns) reads.
  `schema.validate(uri, full=False)` runs that subset; `metacooler validate`,
  `schema.validate()` and every writer still run the full check. A file
  cooler cannot decode (non-ASCII names) still fails the open with the named
  `MetaCoolerFormatError`. `pixels_for_genome_pair` now returns every pixel
  column, typed, when nothing matches (it returned an untyped three-column
  frame).
- **Pixels are aggregated on one `int64` key.** Every pixel metacooler builds
  itself is summed through `io/_utils.aggregate_pixels`: `(bin1, bin2)` folds
  into `bin1 << 32 | bin2`, which is sorted (a stable, run-detecting sort when
  the inputs are already sorted frames) and reduced with `np.add.reduceat`,
  instead of a pandas `groupby` over two columns. `create pairs` bins each chunk
  by sorting its keys and counting runs, and its in-memory sort frees each
  frame as soon as it is keyed. Count dtypes are kept exactly as before.
  `aggregate` and `network` aggregate their chunks the same way, and every
  caller lets `aggregate_pixels` choose the sort: a stable sort (which merges
  existing runs) when the keys come mostly in order, as a genome-grouped remap
  does, quicksort when they are scrambled.
- **Pair contig names are dictionary-encoded.** `create pairs` and `ps --pairs`
  ask pyarrow for the two chrom columns as `dictionary<int32, string>`, so the
  names are never materialised per pair and are mapped to contig ids once per
  distinct name per batch rather than once per row. The drop report (unknown
  contigs with examples, missing positions, out-of-bounds pairs) is unchanged.
  Together with the key-based aggregation, on a 200 k-contig, 30 M-pair
  library: `create pairs` 76.0 s / 3.54 GB → 28.1 s / 2.56 GB, `--sort-strategy
  external` 78.8 s / 3.67 GB → 39.9 s / 3.50 GB; outputs identical.
- **metacooler's own tables are compressed** like cooler's (gzip level 6,
  shuffle, chunked): `/genomes`, `/contigs`, `/bins/chrom`,
  `indexes/chrom_offset`, `indexes/genome_offset` and the QC columns, through
  one `io/_utils.write_column`, in `create`, `binnify`, `merge`, `reassign` and
  QC storage. They were written uncompressed: on a 200 k-contig map they took
  11.8 of 47.5 MB (`bins/chrom` alone 3.4 MB, against 0.15 MB now), and the
  file shrinks from 48.1 to 38.9 MB. A storage change only: filters are
  transparent to readers, and existing files read as before (see
  [Storage filters](format.md#storage-filters)).
- **`binnify` and `zoomify` stream re-binned pixels in order.** At `--level
  contig` the source-to-target bin map is monotonic, so a source pixel's new
  row is the new row of its `bin1` and the upper triangle is preserved: the
  source is cut at output-row boundaries of its `bin1_offset` into blocks of
  about `--chunksize` pixels, each block is read with cooler's pixel selector,
  mapped, sorted on its pixel key and reduced, and the blocks are written in
  order (`ordered=True`) — no temporary
  coolers, no merge passes, memory proportional to one block. `--level genome`,
  which reorders contigs, and `square` sources keep the unordered path (now
  also key-aggregated). `zoomify` coarsens each level from the coarsest level
  already written whose bin size divides it, as cooler's does (10 kb from 5 kb;
  2 kb and 5 kb from 1 kb); each level's counts are those of re-binning the
  base, and its provenance still names the base. On a 200 k-contig map,
  `zoomify -r 1000,2000,5000,10000`: 82.7 s / 4.23 GB → 27.5 s / 2.62 GB, every
  level identical to before.
- **`create cool` remaps in NumPy and sorts on the pixel key.** The remap of
  each source chunk (`_stream_remapped_pixels`, shared with `reassign`'s
  regroup path) no longer copies the chunk and swaps through `.loc`; and since
  contigs are only reordered, the in-memory sort is one argsort on the pixel
  key whose reduction is a no-op — kept because a cooler written without
  cooler's duplicate check can repeat a pixel, which is summed. On a
  200 k-contig cool: 25.4 s / 2.99 GB → 16.7 s / 2.25 GB (together with the
  key-based aggregation above), output identical.
- **Per-contig overhead cut at metagenome scale.** With two million contigs,
  `create pairs` spent ~25 s on its contig and genome tables alone, before a
  single pair (measured with a pairs file of headers only). The table builders
  (`_build_contig_order`, `_build_genomes_table`, `_build_contigs_table`) now
  resolve each contig's genome once per builder and do the grouping, sorting
  and summing in NumPy (one `lexsort`, one `np.add.at`) instead of a Python
  sort per genome with a per-contig key function; the contig-table parser zips
  object arrays rather than pandas' Arrow-backed columns (4× faster) and
  normalises each distinct genome name once; `validate()` detects repeated
  names with a sort instead of a hash-based `np.unique` (4× faster on byte
  strings); the ASCII checks and `_str_to_hdf5` are vectorised for `str`
  columns (falling back to the per-value path for non-ASCII text); and
  `binnify` hands cooler contig ids rather than two million names, which
  cooler iterated as Arrow strings several times per zoom level. Outputs are
  unchanged.
- **`MetaCooler()` opens cooler lazily.** Building a `cooler.Cooler` reads
  every contig name into a dict (1.4 s at two million contigs), which the
  open, `info`, `metadata`, `binsize`, `shape`, `root` and the genome and
  contig lookups (`genome_extent`, `contig_extent`, `genome_bin_ranges`, ...)
  never need: they read the attributes with cooler's own `info()` and the
  offset indexes. It is now created on first use of anything that reads bins,
  pixels or chroms, or resolves a region (`bins()`, `pixels()`, `matrix()`,
  `extent()`, `chromnames`, `bins_for_*`, `pixels_for_genome_pair`, `dump`,
  ...), all of which go through it. A contig name is found by one vectorised
  comparison over the stored `/contigs` names (or cooler's dict, once open).
  An open still refuses a map cooler could not decode (non-ASCII bytes in
  `/contigs`) with the named error. At XL: opening 1.37 s → 0.055 s. Also,
  the pairs header scans read bytes, and `create pairs` maps each batch's
  name dictionary through an object array, halving that step.
- **The contact-QC tally sums each block with one `bincount`.** Every
  contact gets the index of its category, and one `bincount` keyed on
  `(genome, category)` sums every category of every genome at once, where five
  masked passes used to select each category's contacts first. A chunk is
  classified in blocks of 2¹⁸ contacts, whose temporaries stay in cache,
  rather than 50 M at a time (a float-weighted tally, whose sums do not
  associate, still takes each chunk in one pass, so its bits are unchanged).
  Tallying the 73 M pixels of a 2 M-contig map: 6.4 s → 2.2 s, with a few MB
  of scratch memory instead of ~3 GB per 50 M-pixel chunk (the peak of
  `create cool`); `metacooler qc` on that map: 25.8 s / 5.5 GB → 16.6 s /
  3.6 GB; outputs identical.

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
- Benchmarks against cooler 0.12.1: `benchmarks/benchmarks.qmd` times both tools
  from the command line on four reproducible synthetic metagenomes (10³ to 2 ×
  10⁶ contigs): creating a map from pairs, converting, zoomifying, querying a
  contig, a genome or a genome pair, and genome-level aggregation. The
  `benchmarks` workflow renders it for every release tag and on demand,
  attaches the page and its TSV tables to the release, publishes them on the
  documentation site at `/benchmarks/` (the docs build fetches them from the
  `benchmark-page` release), and keeps the datasets in a `benchmark-data-v<N>`
  release. The page says outright that every dataset is simulated (no
  accessions), and reports, beside wall clock: peak memory (for cooler's
  `dump | awk` pipelines, both stages added up), in the summary figure too,
  against the machine's RAM; the size on disk of the `.metacool` against the
  `.cool` with its binning tables and against one `.cool` per MAG (with what
  those hold); and 3 to 5 runs per command, with their spread drawn and
  tabulated. A run stopped at the time limit is drawn where it was stopped, as
  a lower bound, and no longer joined to the line. `cooler cload pairs` is not
  run on XL, where it has never finished within the time limit
  (`BENCH_SKIP_COOLER_CREATE`).

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
