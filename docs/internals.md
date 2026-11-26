# Internals

The rules metacooler's code keeps, and why each one holds: what a change has
to respect so that it breaks nothing the tests cannot see. `AGENTS.md`, the
instructions of the coding agents that work on metacooler, lists them in
short; this page holds them whole, for anyone changing the code.
[Contributing](contributing.md) says how to set up, test, document and
release, and [the format](format.md) what the files hold.

## Two identities in one file { #identities }

**metacooler** extends the [cooler](https://github.com/open2c/cooler) HDF5 format for **multi-genome (metagenomic) Hi-C contact maps**. It layers a `genomes/` → `contigs/` hierarchy on top of cooler's `chroms/bins/pixels` model. The central invariant: **every `.metacool` file must remain openable by `cooler.Cooler()`**. This is achieved by writing `format="HDF5::Cooler"`, a soft-link `chroms → contigs`, and standard `indexes/chrom_offset`, alongside the metacooler-specific `metacooler-format="HDF5::MetaCooler"` attribute and `genomes/` group. Do not break this dual identity when editing the format.

Two file types exist:
- `.metacool` — a single metagenomic contact map (a cooler + `genomes/` hierarchy).
- `.mmetacool` — a collection: `@format="HDF5::MMetaCooler"` with an `entries/<name>/` group per metacool. Addressed via cooler-style URIs: `file.mmetacool::entries/<name>`.

## The code { #the-code }

```text
src/metacooler/
├── api.py         # MetaCooler class — the primary Python API
├── schema.py      # format constants, dtypes, is_metacooler/is_mmetacooler, validate(), read_table()
├── uri.py         # parse_uri("file::group") → (Path, group_path); sample_name()
├── filters.py     # --filter expression parser + DataFrame predicate engine
├── normalize.py   # normalization arithmetic, pure numpy (ICE for any grouping, NB/ZINB fits)
├── qc.py          # the contact-quality tally: its four categories, pure numpy
├── ps.py          # P(s) arithmetic: distance grids and denominators
├── cli/           # one file per command: cmd_<name>.py; logging_utils.py, the logs
├── io/            # file read/write: create, collect, reassign, annotate, clusters, extract, fasta,
│                  # dump, flatten, network, merge, mmetacool, contigs, genomes, colspec, weights,
│                  # sites, staging (locks and staged writes), report
├── ops/           # computation: aggregate, binnify, zoomify, show, range, normalize, qc, ps, report
├── templates/     # the report's page (report.html.j2)
└── _report/       # chromokit's report kit: never edited (chromokit sync replaces it)
```

**Layer separation is intentional**: `io/` reads and writes files; `ops/` does computation on already-open files. Keep new code on the correct side.

**metacooler is a layer over cooler, not a second implementation of it.** Bins, pixels, chroms and regions are read through the `cooler.Cooler` object (`clr.bins()[lo:hi]`, `clr.pixels()[p0:p1]`, `clr.extent()`), and coolers are written through cooler's own writer (via `io/_utils.write_cooler`, below). What metacooler adds is the genome layer — `/genomes`, `/contigs`, `indexes/genome_offset`, the QC block — and fast genome/contig queries over the offset indexes. Do not replace a cooler read or write with raw h5py access for speed: the maintainer had that reverted across the codebase (gains under ~3×, and even an 8.8× h5py-copy `flatten`). The genome-layer and index lookups themselves may read HDF5 directly.

**`MetaCooler` (api.py)** wraps a `cooler.Cooler` and delegates standard reads (`bins`, `pixels`, `matrix`, `chroms`, `extent`) to it, adding genome/contig-aware accessors (`genomes`, `contigs`, `bins_for_genome`, `pixels_for_genome_pair`, `genome_extent`, …). The `Cooler` is built **lazily** (`_clr` is a `_cached`, see [Edits in place](#in-place)): constructing one reads every contig name into a dict (1.4 s at 2 M contigs). `MetaCooler()`, `info`/`metadata`/`binsize`/`shape` (read with cooler's own `info()` on the group), `root` and the genome/contig lookups never build it — a genome resolves to its contigs through `genome_offset` (or `/contigs/genome_id` under a non-grouped ordering), a contig to its bins through `chrom_offset`, a contig name to its row by one comparison over the stored `/contigs` names (`_find_contig`). Anything that reads bins, pixels or chroms, or resolves a region, does build it; pixel rows of a bin range are located through `bin1_offset` and read with cooler's selector (`_iter_pixel_rows`). Bin rows of a set of ranges (`_bins_in_ranges`) are read with cooler's selector too, but with `convert_enum=False`, and named once from `_chrom_dtype`: left to name them, cooler decodes the whole `/chroms` table on every read, so a genome scattered over a non-grouped map cost ranges × contigs (7.7 s for 97 contigs at 200 k). Resolve a *list* of names in one pass over its table, never one scan per name. Keep the cheap paths off `_clr`. The constructor accepts a URI: a bare path opens a `.metacool`; a `file.mmetacool::entries/<name>` URI opens one entry. Opening the root of an `.mmetacool` raises `KeyError` listing the available entries.

**On-disk layout** (`.metacool`):
```
/genomes  [name, length, n_contigs, n_bins]
/contigs  [name, length, genome_id]
/bins     [chrom(=contig id), start, end]
/pixels   [bin1_id, bin2_id, count, ...]      # sorted by (bin1_id, bin2_id)
/indexes  [chrom_offset, bin1_offset, genome_offset]
/chroms   → soft link to /contigs             # cooler compatibility
```
`schema.validate()` enforces this (required groups/columns, offset lengths = count+1, pixel sort order, genome_id range). Run it mentally when changing anything that writes the format. `validate(uri, full=False)` is its O(n_columns) structural subset, which `MetaCooler()` runs on every open (plus `_check_decodable`, so a map cooler could not decode still fails by name); every writer and `metacooler check` (`schema.check()`) run the full check.

## The rules { #rules }

### The public API is exported lazily { #lazy-api }

**Never add barrel imports to `__init__.py`.** Public symbols are exposed via a lazy `__getattr__` dispatch (one `if name == ...` branch per symbol) to keep import time low and avoid circular imports. Add new public API by extending `__all__` and adding a branch.

### Commands register themselves { #commands }

**CLI commands self-register**: each `cli/cmd_<name>.py` imports `cli` from `cli/main.py` and decorates with `@cli.command()`; `main.py` imports every `cmd_*` module at the bottom to trigger registration, and lists it in `COMMAND_GROUPS`. Adding a command = new file + one import line + one entry in the right group. A command's `help=` is one line (it is what `metacooler -h` lists); details belong in its options' help. Don't pass `cls=` to `@cli.command()`: commands are `MetaCoolerCommand` (the group's `command_class`), whose `format_options` starts every panel's help text at one column (`cli/_help_layout.py`; rich alone spreads each panel's spare width over all its columns, so no two panels lined up — `tests/test_cli_help_layout.py` checks every page). Heavy imports (`MetaCooler`, `io`/`ops` modules) go inside the command function, not at module level, so `--help` stays fast.

### URI arguments { #uri-arguments }

**URI arguments** use the `MetaCoolURI` click param type (in `cli/main.py`), which validates the file/group exists and passes the raw string through for downstream `parse_uri()`.

### metacooler's own datasets { #datasets }

**metacooler's own datasets go through `io/_utils.write_column`** (`/genomes`, `/contigs`, `bins/chrom`, `chrom_offset`, `genome_offset`, QC columns): it converts text to fixed-length bytes (`S`; h5py cannot write pandas `object`/`StringDtype` columns) and applies cooler's own storage filters (`H5_STORAGE`: gzip 6, shuffle, chunked). On read, `schema.read_table()` decodes `S`/bytes back to `str` (matching cooler).

### Summing pixels { #summing-pixels }

**Sum pixels with `io/_utils.aggregate_pixels`**, not a pandas `groupby`: `(bin1, bin2)` folds into one `int64` key (`bin1 << 32 | bin2`), sorted (`_sort_kind` picks a stable sort, which merges existing runs, when the chunks arrive sorted) and reduced with `np.add.reduceat`, keeping the count dtype. `create pairs`, `aggregate` and `network` go through it; `binnify` uses the same key helpers (`pixel_keys`, `reduce_sorted_keys`).

### `binnify --level contig` streams { #binnify-streams }

**`binnify --level contig` streams in output-row order**: the bin map is monotonic, so blocks of whole output rows (cut on `bin1_offset`, read with cooler's selector) are mapped, sorted, reduced and written `ordered=True`, with no cooler merge. `--level genome` and `square` sources keep the unordered path. `zoomify` builds each level from the coarsest level already written whose bin size divides it, but records the base map as each level's provenance.

### Normalizations are weights { #normalizations }

**A normalization is a weight column of `/bins`, never a rewrite of the pixels.** metacooler computes three itself — ICE, NormCC and HiCzin (`metacooler normalize ice|normcc|hiczin`: `cli/cmd_normalize.py` → `ops/normalize.py`, arithmetic in the top-level `normalize.py`, site counting in `io/sites.py`) — and stores any other tool's: every weight goes through `io/weights.store_weights` (public, re-exported as `metacooler.store_weights`), one `float64` column with cooler's balancing convention — multiplicative, `NaN` for a masked bin, `divisive_weights=False` — so `cooler.Cooler.matrix(balance=NAME)` reads it; `show`/`dump`/`network --balance` read it through `io/weights.read_weights`. Another tool's normalization (MetaTOR's, say) is computed by that tool and stored with `store_weights`, not re-implemented here. NormCC and HiCzin give every bin its contig's weight, so contig sums reproduce the published contig-level values. A weight must factorise (`raw × w_i × w_j`): host calling or observed/expected are not normalizations. Genome-scoped ICE is our own loop because cooler's cis/trans masks key on `bins/chrom` (a contig here); it must stay bit-identical to `cooler.balance_cooler` where the groupings coincide (`tests/test_normalize.py` oracles; cooler before 0.12 also masked the pixel between two contigs' adjacent bins under `ignore_diags`, so the single-genome oracle runs at `ignore_diags=1`, where every cooler version agrees). `fit_nb_glm` stops on the coefficients, not the deviance: Fisher scoring converges only linearly on this link, and the tests pin it to an exact Newton solution. A weight's `depends_on` attribute names what it was computed from — `pixels`, `genome_id` (the assignment), `contigs/<column>` — and `store_weights` refuses any other entry, or a column the map does not hold, since a typo would make a weight no edit can invalidate. **Any writer that changes one of those must drop the weight** (`io/weights.drop_stale_weights`, as `reassign` and `annotate` do); a writer that rebuilds `/bins` through cooler loses dataset attributes and must restore the ones it carries (`weight_attrs`/`restore_weight_attrs`). Never store a weight as `KR`, `VC` or `VC_SQRT`: cooler reads those names as divisive.

### `--nproc` { #nproc }

**`--nproc`** (`cli/_options.nproc_option`, default 1, refuses 0): `create pairs` bins batches on threads and folds their tallies in file order; `zoomify` builds independent levels in spawned processes and appends them in ladder order. Output must stay identical to `--nproc 1` (`tests/test_nproc.py`).

### Logging { #logging }

**Logging**: use `logging.getLogger("metacooler")`, never `print`. Every CLI command exposes `-v/--verbose` for debug logging. A command that writes a map or a file runs its work inside `cli/logging_utils.run_log`, after its own checks: its log, beside what it writes, gets everything from INFO up (DEBUG with `-v`, whatever `-q` does to stderr), and nothing of the log's header or last line reaches stderr. A command that only reads logs to stderr only. See [Logs](format.md#logs).

### Streaming creation { #streaming-creation }

**Streaming creation**: `io/create.from_pairs` processes pairs in chunks (`--chunksize`) with an `auto`/`memory`/`external` sort strategy so large files don't blow memory. Preserve this when touching the create path. Under `auto`, both entry points switch to cooler's external sort past `_SORT_IN_MEMORY_PIXEL_THRESHOLD` (100 M **deduplicated** pixels), but they measure it differently: `from_cool` knows `nnz` up front and compares directly, while `from_pairs` accumulates through `_buffer_pixels`, which **folds** (re-aggregates) its buffer every `max_rows // _PIXEL_FOLD_FRACTION` rows. The fold is the point: bin pairs recur across chunks and Hi-C pixel counts saturate on a fixed bin grid, so a raw per-chunk row count overestimates the real pixel table by ~an order of magnitude on a deep library and would spill to the external sort for data that fits in RAM easily. Do not "simplify" `_buffer_pixels` back into a plain append-and-count loop.

### `cooler.create`'s cost per call, and `write_cooler` { #write-cooler }

**`cooler.create` carries a fixed per-call overhead that scales with *contig count*, not pixels**, which is why the external path hands `_buffer_pixels`' folded accumulator back as one large frame rather than slicing it, and why **every** streaming default is `schema.DEFAULT_CHUNKSIZE` (50 M rows) — a single constant, deliberately not per-command, since a low value anywhere costs `cooler.create` calls. **Never call `cooler.create_cooler` directly — go through `io/_utils.write_cooler`**, which installs `cooler_writer_shims`. This matters far more than "once per file": with `ordered=False` (`binnify`, `aggregate`, `extract`) cooler's `create_from_unordered` calls `create()` *once per pixel chunk plus once per merge pass*, so the old 1 M-pixel default meant tens of calls, each paying the full per-call cost. The two terms that made this catastrophic are: `write_bins` typed `/bins/chrom` as an HDF5 enum with one member per contig — quadratic (`H5Tenum_insert` rescans every member already inserted), and *doomed*, since past ~4 k contigs the enum header exceeds HDF5's object-header limit, `create_dataset` raises `ValueError` and cooler silently retries as plain int. Measured 2.3 s at 32 k contigs, ×3.8 per doubling, extrapolating to **~2 h at 1.9 M**; `_finalize_metacool_file` rewrites that dataset as `int32` regardless, so the whole cost bought nothing. The second is `get_binsize`, a per-contig groupby loop that does not short-circuit on a uniform grid (~50 µs/contig, so ~1.6 min at 1.9 M) — metacooler substitutes `infer_binsize`, a vectorised equivalent pinned against cooler's by test. The shims restore cooler's globals unconditionally and change no byte of the output (asserted end-to-end at 20 k / 40 k / 80 k contigs). They patch whichever module holds both `write_bins` and `get_binsize` (`cooler.create._writer` from cooler 0.11, `cooler.create._create` before; see `cooler_writer_module`) — cooler 0.11 moved them, the old lookup silently installed nothing, and the quadratic cost came back on the cooler `pip` was resolving. `TestCoolerWriterShims.test_installed_cooler_is_recognised` exists so the next move fails a test instead; if it fails, extend `_COOLER_WRITER_MODULES`, never silence it. This applies to *anything* that ends in `cooler.create` — `cooler.merge_coolers` included, which `io/merge.py` therefore wraps in the shims itself. `write_cooler` also turns an empty (or all-empty-frames) stream under `ordered=False` into an ordered write, since cooler's unordered writer cannot concatenate zero chunks. It also sorts every chunk of an unordered stream that is not sorted already: cooler before 0.12 indexes each chunk as if it were, so once a merge needs a second `mergebuf` epoch (20 M records) the output interleaves them (out-of-order, duplicated pixels; `create pairs --sort-strategy external` failed `validate()` at 200 k contigs). Test-sized inputs never reach a second epoch, which is why `TestUnorderedWrites` forces a tiny `mergebuf`.

### No Python loop over contigs { #no-contig-loops }

**Never map bins per contig with a Python loop.** Both the tiling in `io/create._tile_bins` and the source→target mapping in `ops/binnify._map_old_bins_per_contig` used to iterate contigs — one materialised DataFrame per contig in the first, an `old_chrom == cid` mask per contig (hence O(n_contigs × n_bins)) in the second. Neither finished on a 1.9 M-contig metagenome. Both are now single vectorised passes built from `np.repeat` / `np.cumsum` / one global `np.searchsorted` whose key folds the contig id in above every coordinate, turning `(contig, start)` lexicographic order into plain integer order. Same shape applies to `_tile_uniform` and `_tile_per_genome_ligated`, which were already written this way.

### Python 3.10 and newer { #python }

**Python ≥3.10**: use `match`, `X | Y` unions, `from __future__ import annotations`.

### Options that take several values { #several-values }

**Comma-separated multi-value options** (`--genomes`, `--contigs`) use the shared `parse_comma_list` click callback in `cli/_options.py`, so a value can be repeated or given as a comma-joined list. An explicitly-supplied but empty value is rejected rather than silently disabling the filter.

### The command line's shared helpers { #cli-helpers }

**Shared CLI helpers live in `cli/_options.py` and `cli/_provenance.py`** — `check_output(path, force, inputs=..., defer_removal=...)` is the single "refuse to clobber the output" guard (backed by `io/_utils.check_output_path`, which the ops layer uses directly). **Always pass `inputs`**: an output that is the same file as an input is refused even under `--force` (it used to delete the input before reading it). Writers that stage their output and rename it into place (`create`) pass `defer_removal=True`, so `--force` cannot destroy the previous file before the new one exists. `load_metadata_file()` is the single ingest point for user `--metadata` files and rejects the keys metacooler writes itself (`RESERVED_METADATA_KEYS`: `commands`, `qc`, and the `reassign`/`annotate` histories); `metadata add`/`remove` refuse the same keys through `refuse_reserved_metadata_key`, while `metadata list`/`get` still read them. Do not re-inline these.

### One map, not a collection { #one-map }

**`MetaCoolURI(require_entry=True)`** is the single policy for "this command needs one map, not a collection": it rejects a bare `.mmetacool` path up front and lists the available entry URIs. Commands that legitimately operate on a whole collection (`ls`, `info`, `tree`) leave it False.

### The `count` dtype { #count-dtype }

**`count` dtype**: never hard-cast pixel counts to `int32` when copying — aggregation (`binnify`, `merge`) writes `int64` and a blind cast wraps. Use `io/_utils.resolve_count_dtype()`.

### Names are ASCII { #ascii }

**Contig and genome names must be ASCII**, enforced by `io/_utils.check_ascii_names` (rules in `schema.ascii_violations`). `/contigs` is cooler's `/chroms` and cooler decodes *every* column of it with an ASCII `.astype("U")`, so one accented byte — in `name` or in a `--contigs` metadata column — yields a file that passes every structural check and then raises `UnicodeDecodeError` in `cooler.Cooler()`, hence in `MetaCooler()`. `/genomes/name` is checked for the mirror reason: cooler never reads that group, but `aggregate --level genome` hands those names to `cooler.create_cooler` as chrom names. Extra `/genomes` columns are deliberately *not* checked — accented taxonomy is legal metadata cooler never sees. The check runs three times on purpose: early in `from_pairs`/`from_cool` (before hours of binning), inside `_finalize_metacool_file` (so no writer can emit a bad file), and in `schema.validate()` (so a file that already exists is named rather than vouched for).

### Relabelling a table's key columns { #key-columns }

**Relabel a user table's key columns with `io/_utils.rename_key_columns()`, by position — never `df.rename(columns=...)` by name.** A table carrying its own `genome_id` column otherwise gets two columns of that label and the contig→genome map silently reads the wrong one. The helper also gives headerless extras string labels (`colN`; h5py rejects integer dataset names) and replaces `/` in names (HDF5 would build a nested group cooler cannot read), logging every change.

### Reading a user's table { #user-tables }

**Never `pd.read_csv` a user-supplied text table directly** — go through `io/_utils.read_user_table()`. `read_csv` resolves a header/data field-count mismatch silently and in opposite directions: a header wider than the rows pads the *trailing* columns (every value ends up under a name to its left), while rows wider than the header absorb the surplus leading fields into the DataFrame *index* (shifting the other way). Both parse cleanly, and `keep_default_na=False` — which the readers need so a genome literally called `"NA"` survives — turns the pad into an empty string, so the frame keeps no record of how wide each line was. `read_user_table` counts fields on the raw bytes, rejects both shapes, and drops the all-empty `Unnamed: N` column a trailing tab on the header line creates.

### `--filter` { #filter }

**`--filter` is one comparison per flag, ANDed by repetition** (`filters.py`), and the `genomes.`/`contigs.` table prefix is **mandatory** — both tables have a `length` and a `name`, and inferring the target only when it happens to be unambiguous would make the same expression legal on one file and rejected on the next. There is deliberately no `or` and no parenthesised grammar: real MetaTOR tables carry columns called `Present Markers` and `GC-content`, and a grammar with boolean connectives has to quote or backtick those before it can tokenise them. Here the *leftmost* operator terminates the column name, so nothing needs quoting (and a regex value may contain `>` or `=`). Missing values never match — including under `!=`/`!~` — because a filter selects rows it can positively vouch for.

### Filters resolve to contigs { #selection }

**Filters resolve to a contig-id set** via `io/_utils.resolve_contig_selection`, the single entry point for `--genomes` + `--contigs` + `--filter` on the read side (`dump`, `extract`). Names choose the candidate set, filters *narrow* it — an intersection, not a union. `ops/show.py` has its own `_apply_filters` because it needs ordered *names* and a selection *kind* rather than ids. `create` filters via `_apply_create_filters`, which builds the throwaway `/genomes` and `/contigs` frames with the very same `_build_genomes_table`/`_build_contigs_table` the output will use, so a filter written against a `dump` of an existing file means exactly the same thing at creation time; `genomes.n_bins` is the one absentee there and gets a specific error via the `hints` argument.

### `collect` { #collect }

**`collect` (`io/collect.py`) makes each input cool one genome and renames every contig `<genome>__<contig>`** (`collected_contig_name`, which skips a contig already carrying exactly that prefix, so `collect → extract → collect` is stable). It reuses `create`'s builders — `_apply_create_filters` (with `genome_n_bins`, because the bins exist up front) and `_build_genomes_table` — with `report_divergence=False`, since it reports genomes-table mismatches itself, once, in terms of inputs. Its speed rests on one layout invariant: **every `COLLECT_ORDERINGS` value is genome-grouped**, so each input is one contiguous block of output bins and its pixels only move within it. Under `cool` every block stays sorted and the stream is written `ordered=True` with no sort; a block whose contigs were reordered is re-sorted alone. Do not add `alpha`/`fasta` there without a global sort path. Inputs must share one bin size — a mixed collection would be `variable` in name and neither resolution in practice. The one exception is an input whose every contig is a single bin: cooler calls it variable (it infers `bin-size` from the bins before each contig's last, and there are none), but its bins are exactly the grid of any bin size ≥ its longest contig, so `_check_resolutions` lets it join one.

### A reassignment { #reassign }

**A contig→genome reassignment is not a `/contigs/genome_id` overwrite.** `io/reassign.py` exists because two things follow from the assignment and have to be rebuilt with it: `/genomes` (whose row *set* and row *order* define the `genome_id` integers, so every genome is renumbered) and, under a genome-grouped `contig-ordering`, the *physical* contig layout — contigs stored in ascending `genome_id` order, which is what makes `genome_offset` delimit a genome's contig block and `genome_extent` return a real bin range. Moving a contig between two non-adjacent MAGs breaks that layout and `validate()` refuses the file. Restoring it means permuting the contigs, hence the bins, hence every pixel's bin ids, so `reassign` has two paths: **relabel** (O(n_contigs) edits; copies the file and rewrites four datasets) and **regroup** (O(nnz); the `from_cool` machinery — `_remap_bins` + `_stream_remapped_pixels` + `write_cooler`), both of which run `--inplace` through a staged copy (regroup refuses an entry of a collection there). `--layout auto` preserves whatever the source guaranteed, and a re-grouping that reproduces the existing order is detected and downgraded to the cheap path. `--layout keep` is the escape hatch and pays for it by writing `contig-ordering="custom"` — a value no writer other than `reassign` produces, and the *point* of which is that it is not in `GENOME_GROUPED_ORDERINGS`, so `validate()` stops enforcing a grouping the rows no longer have and `genome_extent` raises instead of returning a range that includes foreign bins. Two things a reassignment invalidates and must not carry over: the QC block (three of its four categories are defined by the assignment — recomputed by default, dropped under `--no-qc`) and per-genome metadata for genomes whose contig set changed (carried, but warned about by name). `--missing` chooses between an incremental edit (`keep`) and an authoritative table (`unbinned`); neither *drops* a contig — that is `extract`'s job, and doing it here would mean discarding pixels. `auto` reads only the *declared* ordering, so a `cool` source is relabelled even when its genomes are contiguous — every `collect`ed map's are — and `_kept_ordering` names any genome that scatters; the regroup path refuses a `square` map (remapping folds pixels into the upper triangle, as `from_cool` refuses for the same reason).

### Which columns hold the names { #colspec }

**Which columns of a user table hold its names is one grammar, `io/colspec.py`**: `COLUMN[,COLUMN][,noheader]`, each column a 1-based position (digits) or a header name, or an alias. Names only *find* positions; the key columns are then relabelled by position with `rename_key_columns`. Two kinds of contigs table share `--contigs-format`: an **assignment** table (`create pairs`, `create cool`, `reassign`; `contigs_format_option`) names two columns, contig then genome; a **metadata** table (`collect`, `annotate`; `contig_metadata_format_option`) names the contig column, and a second column, if given, is only *checked* as the genome. The one-column form is documented in `collect`/`annotate` only, and the two-column form (with aliases) there is accepted but deliberately undocumented. `--genomes-format` (one column) is on every command whose `--genomes` is a table. The default `contigs` is `contig_name,genome_name` **by name** — never make it positional: a positional default reads a table given without its format with the wrong columns, silently. `NA` in the genome's place (`ColumnSpec.no_genome`) says the table has no genome column — with `create`, a map with **no binning yet** that keeps the table's other columns (`reassign` unbins the listed contigs; `collect`/`annotate` read `X,NA` as `X`). Keep it explicit: a spec that merely omits the genome column must stay an error where one is required, or a forgotten column passes for "unbinned" after hours of binning. The binner aliases (`metator`, `dastool`, `binette`, `semibin`, `vamb`) were each taken from the tool's own writer; check a new one the same way. The CLI validates a spec while parsing options (`_column_spec_callback`), before any file is read.

### `create` and its contigs table { #create-contigs }

**`create` refuses a contigs table that names no contig of the assembly, or has no row** (`_check_contig_table_matches`), as `reassign` does: it runs against every contig the assembly names, before `--drop-unassigned` and `--filter` narrow them and before a pair or pixel is read, and only warns about rows naming contigs the assembly lacks. A map with **no binning yet** is asked for explicitly — by leaving `--contigs` out (`contig_genome_map=None`: the cooler plus the genome layer, nothing to check) or, to keep a table's other columns, with `NA` in `--contigs-format`'s genome place — never with an empty table. `_check_no_binning` refuses `--genomes` and `--drop-unassigned` either way, and the CLI refuses a `--contigs-format` given without `--contigs`. Genome labels in `_UNBINNED_SENTINELS` (`NA`, `ND`, `-1`, empty, …) mean `__unbinned__`.

### `annotate` { #annotate }

**`annotate` (`io/annotate.py`) adds/replaces/drops metadata columns of `/genomes` and `/contigs` in an existing map, and its invariant is `annotate(T) == create(T)`**: same parsers (`parse_genomes_table`, `collect`'s `_read_contig_metadata`/`_resolve_contig_metadata` with `command="annotate"`), same `user_` renames, same typing (`_coerce_numeric_columns`; a contigs table is typed **once, as read, over every row** — `_type_contig_metadata`, called by `create`'s `_parse_contig_map` and by `_read_contig_metadata` — never again in a builder, so a column one unmatched `N/A` keeps text is text under every command and in every entry), pinned by `TestSameMeaningAsCreate` and `TestOneTableOneType`. It plans every map first (no write until every table is checked: matches, conflicts, ASCII in `/contigs`), never touches derived columns (`name`, `length`, `genome_id`, `n_contigs`, `n_bins`, `qc_*` — those belong to `reassign` and `qc --store`), and rebuilds a changed table in a scratch group (`.annotate-scratch`) before swapping it in, copying untouched columns as they are. Identical values are not a change, so an in-place re-run writes nothing. A bare `.mmetacool` annotates every entry (the only editing command that loops; the others keep `require_entry=True`).

### Edits in place, locks and staged writes { #in-place }

**An in-place edit takes the map's lock before its first read, and writes through `io/staging.staged_inplace`** (`io/staging.lock`; `metacooler.lock` from Python): a fresh copy of the file (`copy_fresh`, object by object: it leaves behind the space earlier edits left, which HDF5 never reclaims) in a private directory beside it (`.<map>.<uuid>.<n>.d/`, `0700`, keeping an inherited set-group-id bit), edited, validated in full, given the map's owner and group (the group alone where the owner cannot be given back, with a warning), mode and extended attributes (POSIX ACLs), synced, checked unchanged since the copy began (its state is taken before `planned` is checked; a program writing the map without the lock meanwhile makes the edit fail with `MapChangedError`: never a lost update) and renamed over it (`_install`: a rename only, as its writers left it — do not copy it again there: a second copy doubled every edit's time and disk, and the fresh copy at the start already bounds the map to what one edit, or one session, replaced: `TestBoundedMaps`), so a failure of any kind (`KeyboardInterrupt` included) leaves the map as it was and the map never grows from one edit to the next. `annotate -o` and `reassign -o` (relabel) copy their source with `copy_fresh`, never byte for byte (`shutil.copy2` carried the source's dead space into every output, which grew with each update). An in-place edit through a symbolic link edits its target (lock, staging directory and rename beside the target); an output that is a link is replaced, never written through (`_staged_output`, `staged_edit(output=)`). `annotate(inplace=True)`, `reassign(inplace=True)` (relabel *and* regroup; `locked_when_inplace` takes the lock before it reads), `store_weights`, `write_clusters`, `drop_clusters`, `qc(store=True)`, `normalize_*`, `edit_metadata` (the `metadata` commands) and the `.mmetacool` entry functions do; a new in-place writer must too, never `h5py.File(path, "r+")` on the map itself, and must pass what it read (`planned=`; `store_weights(source=...)`) so a map changed since is refused. **`metacooler.staged_edit` is the batched form**: one staged copy (with `output=`, an entry copied out as a standalone map, `copy_fresh(group=)`), registered for its thread in `_SESSIONS`; `lock()` and `staged_inplace()` of that staged path from that thread write it directly (`_in_session`: no copy, no lock; a writer of a whole new map, a re-grouping, writes beside it and replaces it), a writer's own full pre-validation is skipped there (`in_session`), and the session validates what was written and renames once — not at all, in place, if the staged file is as it was copied: `_written_state` (inode, size, `mtime_ns`, `ctime_ns`; `_as_copied` sets its mtime back ten seconds after the copy, so a write in the same coarse clock tick still shows), judged by the file, never by which writers ran (a `reassign(src, output=staged, force=True)` renames a new file over it; a raw `h5py` write is a write too). The staged map has its own lock file (`_hold_staged`, naming `staged_for`), registered in `_HELD` for the session's thread (re-entered by it — an output session onto the staged map —, never taken over as a dead process's), so another process — a forked child included — or thread is refused it. **A writer of a new map at a path takes no lock: it goes through `writing_output`** (`_staged_output` — `create`, `collect`, `reassign`/`annotate -o`, `zoomify`, `stack` —, `binnify`, `aggregate`, `extract`, `flatten`, `merge`, `to_cooler`, `create_mmetacool`, `cp`'s new collection), which refuses this thread's session map or output (`refuse_a_session_output`) and a held lock (`check_output_lock`: another thread, a live process, one this host cannot check; a stale lock, or a file that names no process — another tool's lock file —, does not stop it), and removes a symbolic link at the output under force, refusing it without (`remove_an_output_link`: cooler's writer and h5py open what a link points to; a writer that renames its output into place — `_staged_output`, `writing_output(renamed=True)`, `check_output_path(defer_removal=True)` — leaves the link for its rename to replace, so a failed run leaves it); a directory, a device or a FIFO at the output, or a link to one, is refused by every writer, with or without force, and never removed (`_refuse_what_is_not_a_file` in `writing_output`, and `check_output_path`: cooler 0.12 renames its output over the path, so `-o /dev/stdout` replaced the system's); `_staged_output` checks the lock again right before its rename, and `_refuse_staging_an_input` and `check_output_path` run the same checks up front, before hours of binning. Never create a lock for a plain output: a writer killed outright in a container, or on another node, would leave a lock that no rerun elsewhere may take over (output locks were tried, and reverted for that). A new writer of an output goes through `writing_output`, with its `force`. From a session's thread an output over the map the session edits or copies out, or over its output, is refused (`refuse_a_session_output`, in `writing_output` and `check_output_path`); the staged map itself is an output its thread may write (metaTOR's `reassign(map, output=staged, force=True)`), logged as a write of anything. Each write of a session is logged with what it changed of what a weight may depend on (`staged_inplace(touches=)`: `()` for `store_weights`, the cluster writers and `edit_metadata`; the changed `contigs/<column>`s for `annotate` and `qc`; `genome_id` and the QC columns for a relabel; `None`, anything, for a re-grouping, an output written over the staged map, a raw write — `_sync_untracked` compares `_written_state`); The log records each write's entry too (`""`: the whole file), and `_changed_since` counts only the writes of the weight's own entry. `check_source_depends` refuses `store_weights(source=mc)` when `mc` holds what a logged write changed of its `depends_on`: `mc` on the staged map is judged by the log alone — the writes after `mc._opened`, a write of anything refusing it as its bins may have changed — and left out of `planned`, whose file check every write would fail (two weights from one read); `mc` on the session's source file (`session.source`, in whatever state) counts every write of the session (none is ever in it); `mc` of another file, the writes after `mc._opened`. A new in-place writer passes `touches=`, or leaves it `None`. A writer that raises with the staged map in hand fails the session (`_failing`: `_in_session`, a joined session, `MetaCooler.open("r+")` via `_session_write`), which then raises `RuntimeError` at its end even if the error was caught: a writer may still clean up cheaply (`write_clusters`' scratch group, `store_weights`' column written under a scratch name before it replaces the old), but a caught error never installs a half-written map. A planned map that is not the map written (another file or entry; `store_weights(source=)` opened before a session) must have the same bins — contig name, start, end — (`_check_planned_bins`; one current as copied is not read again unless a write of the session may have changed the bins). An output session takes the lock of the output by its own name (`_named`, `_locked(output=True)`: the lock `lock(out, follow_symlinks=False)` takes, which it re-enters), stages beside it, and takes no lock of the source, only read (its state taken before `planned` is checked, the copy checked untorn against it by `_check_copied`, which every staged copy runs; a map opened on the source as copied stays current whatever is written to it since); without `force`, an output created meanwhile is never replaced (`_install(refuse_existing=)` → `_move_no_replace`: a hard link, which no existing name takes, then the staged name removed; check-then-rename only where the file system has no hard links). Its copy's time is set back as an in-place session's is, and the output given its install's time. A writer that records a map — `reassign`'s history `source`, `binnify`'s — records `recorded_uri(uri)`: the session's map (`_Session.named`, as given to `staged_edit`; an entry's URI, or the map copied out), never the staged path. A nested in-place session on the same map joins it; other nestings, another thread, and a write of the map a session edits, copies out or replaces (`_session_on`) while it runs are refused. `zoomify` and `stack` build their output in one `_staged_output` registered with `_writing_directly`, so an appended level or entry is not a copy of the growing file. A map without write permission is refused (by mode bits too: tests run as root). The lock is re-entrant for the *thread* that holds it (metaTOR holds it for a whole run and calls the writers; another thread is refused), names its holder in `<map>.lock` (pid, host, machine_id, pid_ns, boot_id, pid_start, uuid), is taken over only from a dead process of the same host, machine and PID namespace (a lock that does not say where it runs, metaTOR's old one, is never broken), under the `<map>.lock.takeover` guard (never move or delete a lock that might be live); whoever takes a lock removes every other uuid's `.<map>.<uuid>.<n>...` beside the map (`_remove_foreign_staging`), naming each. Waits (an empty lock file, another's takeover) release the registry (`_pause`, a `Condition` on `_REGISTRY`), and a forked child gets a fresh one (`os.register_at_fork`). A forked child neither releases nor cleans up. No signal handlers. `MetaCooler` checks its file state before and after every read: `_clr` is a `_CheckedCooler` whose selectors are `_CheckedSelector`s and whose `open()` is checked (cooler's `balance_cooler` reads its chunks through it; cooler's own store stays a path: cooler < 0.12 cannot take anything else), and `MetaCooler.open()` checks on entry and exit; everything it keeps is a `_cached` (a data descriptor that checks the file before every read of the value it keeps — never `functools.cached_property`, which answers from its cache unchecked), the public tables (`genomes`, `contigs`, `cluster_rounds`, `clusters`) being properties over private caches; a write through its own `open("r+")` (locked, unstaged) refreshes the state and drops what it kept (`_forget`). A copy by path after reading (`copy_fresh` in `_materialise_copy` and `annotate -o`, `_copy_map`) is followed by `_check_unchanged()`; `merge` and `collect`, which read their inputs by path, check them at the end (`input_states`, `check_inputs_unchanged`).

### `fasta` { #fasta }

**`fasta` (`io/fasta.py`) copies records, it does not re-format them**: a contig's record is the bytes from its `>` to the next one, as pyfastx's `raw` gives them, so metaTOR's `final_clusters/*.fa` stay byte-identical when written through it; only a stitched record is built (bases alone, `spacer` N, `LINE_WIDTH` columns). The order is a contract metaTOR relies on: a caller's `groups` as given, stitched or not; a genome's or a round's contigs in map order as they are, **longest first (ties in map order) when stitched** -- a round's cluster may span genomes, and their map order follows the current assignment. A record is named and measured as the map's creator read it: `>[ \t]*(\S*)` (`> c1` is `c1`), and the map's length is either the bytes but line breaks (pyfaidx, `create --fasta`) or the bases alone. It selects on its own (`_genome_groups`, `_round_groups`) rather than through `resolve_contig_selection`: it needs groups, not a contig set, the unbinned contigs answer to `contigs.` filters only, names or genome/cluster filters that keep no genome/cluster are refused even with `unbinned`, and a `clusters.` filter (`filter_option(tables=...)`; every other command refuses the prefix by name) only exists with `--round`. The assembly is mapped in memory and scanned once for the literal `\n>` (several times faster than a multi-line `^>`); a gzip or piped one is copied, plain, into `temp_dir` first. Outputs are checked before the assembly is read, staged as short `.fasta-*.tmp` names and renamed only once every file is written. Never write an index next to the user's assembly.

### `merge` compares identity columns { #merge }

**`merge` compares identity columns only** (`_IDENTITY_COLUMNS`: contig `name`/`length`/`genome_id`, genome `name`/`length`/`n_contigs`/`n_bins`, widened before hashing). Metadata and `qc_*` columns may differ; the output keeps the first input's metadata and warns per input. Hashing every column made one `annotate` run enough to make a replicate unmergeable.

### `/genomes` derives from the contigs table { #genomes-table }

**`/genomes` rows are derived from the contigs table**, one per distinct genome name; a `--genomes` table only left-joins extra columns onto them. It can neither add a genome with no contigs nor suppress one that has them, so `_report_genome_table_divergence` warns about both directions of mismatch — without it, MetaTOR's `MetaTOR_MGE_*` bins appear in `/genomes` unexplained and look fabricated.

## Correctness, speed and the API { #more-rules }

What a review checks besides the rules above (see
[Reviewing a pull request](contributing.md#reviewing-a-pull-request)). A real
metagenome has 50,000 contigs or more and tens of millions of pairs: what a
test's toy data cannot show, these rules keep out.

### The HDF5 layer { #hdf5 }

- `chroms` is a soft link to `contigs`, `h5py.SoftLink("contigs")`: a hard
  link (`f["chroms"] = f["contigs"]`) copies the data and goes stale with the
  next write.
- A collection's entries, copied from metacool files, carry that soft link,
  which may point nowhere while a walk over the collection's groups reads it:
  the walk takes each child in a `try`/`except`.
- Pixels are sorted by `(bin1_id, bin2_id)` with `bin1_id <= bin2_id`, and bin
  ids are contiguous from 0: a remapping of the bins puts each pixel back in
  the upper triangle.
- A filter of the bins or the pixels rebuilds `chrom_offset` and
  `bin1_offset`, each one longer than what it indexes.
- A class that keeps the `mode` of an `open(mode)` opens every file with that
  mode, never a `"r"` written out.

### Speed { #speed }

- No Python loop over pixels, pairs or bins (`for i in range(len(array))`,
  `DataFrame.iterrows`): numpy, at most once per chunk, and not even once per
  contig ([No Python loop over contigs](#no-contig-loops)).
- No dense matrix but under a bound: `show` refuses one past
  `schema.DENSE_MATRIX_BIN_LIMIT` bins a side, saying how much memory it would
  take, rather than allocating it.
- A function that takes a `chunksize` reads by it.

### The Python API { #api-design }

- A parameter that is documented, and exposed on the command line, does what
  it says, or raises `NotImplementedError` for what it does not do yet.
- An upper-triangular matrix is made symmetric one way only, its diagonal
  once: `agg + agg.T - sparse.diags(agg.diagonal())`.
- A public method that opens a file is a context manager
  (`with mc.open("r+") as f:`), or is private and says that its caller closes
  what it returns.
- A deprecation warns when the deprecated thing is used, never only in a
  docstring: `FutureWarning` from Python, a warning on stderr from the command
  line, naming what replaces it, for one minor release at least before it
  goes (CONVENTIONS.md § 18). A renamed option keeps its old name as an alias.
- A hash that is not for security is `hashlib.md5(usedforsecurity=False)`:
  FIPS systems refuse the plain call.
- A function that may overwrite a file takes `force=False` and raises
  `FileExistsError` without it; its command has `-f/--force`.
- One word for one thing in `api.py`, `ops/` and `cli/`: `level="contig"`, not
  `"contigs"`.
- A temporary file is removed with whatever a library wrote beside it
  (pyfaidx's `.fai`).

### Errors { #errors }

- A user's mistake is a `ValueError`, `KeyError`, `OSError` or
  `NotImplementedError` from `io/` or `ops/`, which the command line reports
  as `Error: <message>`, without a traceback (`cli/main.py`'s `_USER_ERRORS`);
  any other exception is a bug, and keeps its traceback. Command code raises
  click's `UsageError`, `BadParameter` or `ClickException`.
- A message names the value it refuses and what is allowed: a name that is not
  found lists those there are (`Genome 'x' not found. Available: [...]`), and
  "strictly larger" says so when `>` is meant.
