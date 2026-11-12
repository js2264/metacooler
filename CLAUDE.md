# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**metacooler** extends the [cooler](https://github.com/open2c/cooler) HDF5 format for **multi-genome (metagenomic) Hi-C contact maps**. It layers a `genomes/` → `contigs/` hierarchy on top of cooler's `chroms/bins/pixels` model. The central invariant: **every `.metacool` file must remain openable by `cooler.Cooler()`**. This is achieved by writing `format="HDF5::Cooler"`, a soft-link `chroms → contigs`, and standard `indexes/chrom_offset`, alongside the metacooler-specific `metacooler-format="HDF5::MetaCooler"` attribute and `genomes/` group. Do not break this dual identity when editing the format.

Two file types exist:
- `.metacool` — a single metagenomic contact map (a cooler + `genomes/` hierarchy).
- `.mmetacool` — a collection: `@format="HDF5::MMetaCooler"` with an `entries/<name>/` group per metacool. Addressed via cooler-style URIs: `file.mmetacool::entries/<name>`.

## Commands

Uses `uv` (not pip directly) and `ty` as the type checker (**not mypy**).

```bash
uv sync                                        # install with the dev dependency group

uv run pytest tests/ -x -v                     # run tests (fail fast); --cov is on by default via pyproject
uv run pytest tests/test_create_correctness.py -x       # single test file
uv run pytest tests/test_create.py::test_name -x        # single test

uv run ruff check src/ tests/                  # lint
uv run ruff format src/ tests/                 # format (line-length 88, double quotes; no black)
uv run ty check src/ tests/                    # type check

# patch coverage: Codecov fails a PR whose changed lines are < 80% covered
uv run pytest tests/ --cov-report=xml
uv run diff-cover coverage.xml --compare-branch=origin/devel --fail-under=80

make fix          # ty check + ruff format + ruff check --fix + regenerate README.md from README.qmd
make regenerate   # rebuild tests/data/ fixtures (generate-raw-files.py + raw-to-cool.sh)
make cli-tests    # run tests/data/cli_tests.sh end-to-end CLI smoke tests

uv run metacooler --help                       # run CLI
```

**README.md is generated** from `README.qmd` via `quarto render` (see the `make fix` target). Edit `README.qmd`, never `README.md` directly — it is also excluded from ruff.

## Architecture

```
src/metacooler/
├── api.py       # MetaCooler class — the primary Python API
├── schema.py    # format constants, dtypes, is_metacooler/is_mmetacooler, validate(), read_table()
├── uri.py       # parse_uri("file::group") → (Path, group_path)
├── filters.py   # --filter expression parser + DataFrame predicate engine
├── cli/         # one file per command: cmd_<name>.py
├── io/          # file read/write: create, collect, reassign, extract, dump, flatten, network, merge, mmetacool, contigs, genomes
└── ops/         # computation: aggregate, binnify, zoomify, show, range
```

**Layer separation is intentional**: `io/` reads and writes files; `ops/` does computation on already-open files. Keep new code on the correct side.

**`MetaCooler` (api.py)** wraps a `cooler.Cooler` instance and delegates standard reads (`bins`, `pixels`, `matrix`, `chroms`) to it, adding genome/contig-aware accessors (`genomes`, `contigs`, `bins_for_genome`, `pixels_for_genome_pair`, `genome_extent`, …). The constructor accepts a URI: a bare path opens a `.metacool`; a `file.mmetacool::entries/<name>` URI opens one entry. Opening the root of an `.mmetacool` raises `KeyError` listing the available entries.

**On-disk layout** (`.metacool`):
```
/genomes  [name, length, n_contigs, n_bins]
/contigs  [name, length, genome_id]
/bins     [chrom(=contig id), start, end]
/pixels   [bin1_id, bin2_id, count, ...]      # sorted by (bin1_id, bin2_id)
/indexes  [chrom_offset, bin1_offset, genome_offset]
/chroms   → soft link to /contigs             # cooler compatibility
```
`schema.validate()` enforces this (required groups/columns, offset lengths = count+1, pixel sort order, genome_id range). Run it mentally when changing anything that writes the format.

## Conventions (non-obvious)

- **Never add barrel imports to `__init__.py`.** Public symbols are exposed via a lazy `__getattr__` dispatch (one `if name == ...` branch per symbol) to keep import time low and avoid circular imports. Add new public API by extending `__all__` and adding a branch.
- **CLI commands self-register**: each `cli/cmd_<name>.py` imports `cli` from `cli/main.py` and decorates with `@cli.command()`; `main.py` imports every `cmd_*` module at the bottom to trigger registration, and lists it in `COMMAND_GROUPS`. Adding a command = new file + one import line + one entry in the right group.
- **URI arguments** use the `MetaCoolURI` click param type (in `cli/main.py`), which validates the file/group exists and passes the raw string through for downstream `parse_uri()`.
- **HDF5 string writes**: h5py cannot write pandas `object`/`StringDtype` columns — convert to fixed-length bytes (`S`) before `create_dataset()`. On read, `schema.read_table()` decodes `S`/bytes back to `str` (matching cooler).
- **Logging**: use `logging.getLogger("metacooler")`, never `print`. Every CLI command exposes `-v/--verbose` for debug logging.
- **Streaming creation**: `io/create.from_pairs` processes pairs in chunks (`--chunksize`) with an `auto`/`memory`/`external` sort strategy so large files don't blow memory. Preserve this when touching the create path. Under `auto`, both entry points switch to cooler's external sort past `_SORT_IN_MEMORY_PIXEL_THRESHOLD` (100 M **deduplicated** pixels), but they measure it differently: `from_cool` knows `nnz` up front and compares directly, while `from_pairs` accumulates through `_buffer_pixels`, which **folds** (re-aggregates) its buffer every `max_rows // _PIXEL_FOLD_FRACTION` rows. The fold is the point: bin pairs recur across chunks and Hi-C pixel counts saturate on a fixed bin grid, so a raw per-chunk row count overestimates the real pixel table by ~an order of magnitude on a deep library and would spill to the external sort for data that fits in RAM easily. Do not "simplify" `_buffer_pixels` back into a plain append-and-count loop.
- **`cooler.create` carries a fixed per-call overhead that scales with *contig count*, not pixels**, which is why the external path hands `_buffer_pixels`' folded accumulator back as one large frame rather than slicing it, and why **every** streaming default is `schema.DEFAULT_CHUNKSIZE` (50 M rows) — a single constant, deliberately not per-command, since a low value anywhere costs `cooler.create` calls. **Never call `cooler.create_cooler` directly — go through `io/_utils.write_cooler`**, which installs `cooler_writer_shims`. This matters far more than "once per file": with `ordered=False` (`binnify`, `aggregate`, `extract`) cooler's `create_from_unordered` calls `create()` *once per pixel chunk plus once per merge pass*, so the old 1 M-pixel default meant tens of calls, each paying the full per-call cost. The two terms that made this catastrophic are: `write_bins` typed `/bins/chrom` as an HDF5 enum with one member per contig — quadratic (`H5Tenum_insert` rescans every member already inserted), and *doomed*, since past ~4 k contigs the enum header exceeds HDF5's object-header limit, `create_dataset` raises `ValueError` and cooler silently retries as plain int. Measured 2.3 s at 32 k contigs, ×3.8 per doubling, extrapolating to **~2 h at 1.9 M**; `_finalize_metacool_file` rewrites that dataset as `int32` regardless, so the whole cost bought nothing. The second is `get_binsize`, a per-contig groupby loop that does not short-circuit on a uniform grid (~50 µs/contig, so ~1.6 min at 1.9 M) — metacooler substitutes `infer_binsize`, a vectorised equivalent pinned against cooler's by test. The shims restore cooler's globals unconditionally and change no byte of the output (asserted end-to-end at 20 k / 40 k / 80 k contigs). They patch whichever module holds both `write_bins` and `get_binsize` (`cooler.create._writer` from cooler 0.11, `cooler.create._create` before; see `cooler_writer_module`) — cooler 0.11 moved them, the old lookup silently installed nothing, and the quadratic cost came back on the cooler `pip` was resolving. `TestCoolerWriterShims.test_installed_cooler_is_recognised` exists so the next move fails a test instead; if it fails, extend `_COOLER_WRITER_MODULES`, never silence it. This applies to *anything* that ends in `cooler.create` — `cooler.merge_coolers` included, which `io/merge.py` therefore wraps in the shims itself. `write_cooler` also turns an empty (or all-empty-frames) stream under `ordered=False` into an ordered write, since cooler's unordered writer cannot concatenate zero chunks. It also sorts every chunk of an unordered stream that is not sorted already: cooler before 0.12 indexes each chunk as if it were, so once a merge needs a second `mergebuf` epoch (20 M records) the output interleaves them (out-of-order, duplicated pixels; `create pairs --sort-strategy external` failed `validate()` at 200 k contigs). Test-sized inputs never reach a second epoch, which is why `TestUnorderedWrites` forces a tiny `mergebuf`.
- **Never map bins per contig with a Python loop.** Both the tiling in `io/create._tile_bins` and the source→target mapping in `ops/binnify._map_old_bins_per_contig` used to iterate contigs — one materialised DataFrame per contig in the first, an `old_chrom == cid` mask per contig (hence O(n_contigs × n_bins)) in the second. Neither finished on a 1.9 M-contig metagenome. Both are now single vectorised passes built from `np.repeat` / `np.cumsum` / one global `np.searchsorted` whose key folds the contig id in above every coordinate, turning `(contig, start)` lexicographic order into plain integer order. Same shape applies to `_tile_uniform` and `_tile_per_genome_ligated`, which were already written this way.
- **Python ≥3.10**: use `match`, `X | Y` unions, `from __future__ import annotations`.
- **Comma-separated multi-value options** (`--genomes`, `--contigs`) use the shared `parse_comma_list` click callback in `cli/_options.py`, so a value can be repeated or given as a comma-joined list. An explicitly-supplied but empty value is rejected rather than silently disabling the filter.
- **Shared CLI helpers live in `cli/_options.py` and `cli/_provenance.py`** — `check_output(path, force, inputs=..., defer_removal=...)` is the single "refuse to clobber the output" guard (backed by `io/_utils.check_output_path`, which the ops layer uses directly). **Always pass `inputs`**: an output that is the same file as an input is refused even under `--force` (it used to delete the input before reading it). Writers that stage their output and rename it into place (`create`) pass `defer_removal=True`, so `--force` cannot destroy the previous file before the new one exists. `load_metadata_file()` is the single ingest point for user `--metadata` files and rejects the reserved `commands` and `qc` keys. Do not re-inline these.
- **`MetaCoolURI(require_entry=True)`** is the single policy for "this command needs one map, not a collection": it rejects a bare `.mmetacool` path up front and lists the available entry URIs. Commands that legitimately operate on a whole collection (`ls`, `info`, `tree`) leave it False.
- **`count` dtype**: never hard-cast pixel counts to `int32` when copying — aggregation (`binnify`, `merge`) writes `int64` and a blind cast wraps. Use `io/_utils.resolve_count_dtype()`.
- **Contig and genome names must be ASCII**, enforced by `io/_utils.check_ascii_names` (rules in `schema.ascii_violations`). `/contigs` is cooler's `/chroms` and cooler decodes *every* column of it with an ASCII `.astype("U")`, so one accented byte — in `name` or in a `--contigs` metadata column — yields a file that passes every structural check and then raises `UnicodeDecodeError` in `cooler.Cooler()`, hence in `MetaCooler()`. `/genomes/name` is checked for the mirror reason: cooler never reads that group, but `aggregate --level genome` hands those names to `cooler.create_cooler` as chrom names. Extra `/genomes` columns are deliberately *not* checked — accented taxonomy is legal metadata cooler never sees. The check runs three times on purpose: early in `from_pairs`/`from_cool` (before hours of binning), inside `_finalize_metacool_file` (so no writer can emit a bad file), and in `schema.validate()` (so a file that already exists is named rather than vouched for).
- **Relabel a user table's key columns with `io/_utils.rename_key_columns()`, by position — never `df.rename(columns=...)` by name.** A table carrying its own `genome_id` column otherwise gets two columns of that label and the contig→genome map silently reads the wrong one. The helper also gives headerless extras string labels (`colN`; h5py rejects integer dataset names) and replaces `/` in names (HDF5 would build a nested group cooler cannot read), logging every change.
- **Never `pd.read_csv` a user-supplied text table directly** — go through `io/_utils.read_user_table()`. `read_csv` resolves a header/data field-count mismatch silently and in opposite directions: a header wider than the rows pads the *trailing* columns (every value ends up under a name to its left), while rows wider than the header absorb the surplus leading fields into the DataFrame *index* (shifting the other way). Both parse cleanly, and `keep_default_na=False` — which the readers need so a genome literally called `"NA"` survives — turns the pad into an empty string, so the frame keeps no record of how wide each line was. `read_user_table` counts fields on the raw bytes, rejects both shapes, and drops the all-empty `Unnamed: N` column a trailing tab on the header line creates.
- **`--filter` is one comparison per flag, ANDed by repetition** (`filters.py`), and the `genomes.`/`contigs.` table prefix is **mandatory** — both tables have a `length` and a `name`, and inferring the target only when it happens to be unambiguous would make the same expression legal on one file and rejected on the next. There is deliberately no `or` and no parenthesised grammar: real MetaTOR tables carry columns called `Present Markers` and `GC-content`, and a grammar with boolean connectives has to quote or backtick those before it can tokenise them. Here the *leftmost* operator terminates the column name, so nothing needs quoting (and a regex value may contain `>` or `=`). Missing values never match — including under `!=`/`!~` — because a filter selects rows it can positively vouch for.
- **Filters resolve to a contig-id set** via `io/_utils.resolve_contig_selection`, the single entry point for `--genomes` + `--contigs` + `--filter` on the read side (`dump`, `extract`). Names choose the candidate set, filters *narrow* it — an intersection, not a union. `ops/show.py` has its own `_apply_filters` because it needs ordered *names* and a selection *kind* rather than ids. `create` filters via `_apply_create_filters`, which builds the throwaway `/genomes` and `/contigs` frames with the very same `_build_genomes_table`/`_build_contigs_table` the output will use, so a filter written against a `dump` of an existing file means exactly the same thing at creation time; `genomes.n_bins` is the one absentee there and gets a specific error via the `hints` argument.
- **`collect` (`io/collect.py`) makes each input cool one genome and renames every contig `<genome>__<contig>`** (`collected_contig_name`, which skips a contig already carrying exactly that prefix, so `collect → extract → collect` is stable). It reuses `create`'s builders — `_apply_create_filters` (with `genome_n_bins`, because the bins exist up front) and `_build_genomes_table` — with `report_divergence=False`, since it reports genomes-table mismatches itself, once, in terms of inputs. Its speed rests on one layout invariant: **every `COLLECT_ORDERINGS` value is genome-grouped**, so each input is one contiguous block of output bins and its pixels only move within it. Under `cool` every block stays sorted and the stream is written `ordered=True` with no sort; a block whose contigs were reordered is re-sorted alone. Do not add `alpha`/`fasta` there without a global sort path. Inputs must share one bin size — a mixed collection would be `variable` in name and neither resolution in practice.
- **A contig→genome reassignment is not a `/contigs/genome_id` overwrite.** `io/reassign.py` exists because two things follow from the assignment and have to be rebuilt with it: `/genomes` (whose row *set* and row *order* define the `genome_id` integers, so every genome is renumbered) and, under a genome-grouped `contig-ordering`, the *physical* contig layout — contigs stored in ascending `genome_id` order, which is what makes `genome_offset` delimit a genome's contig block and `genome_extent` return a real bin range. Moving a contig between two non-adjacent MAGs breaks that layout and `validate()` refuses the file. Restoring it means permuting the contigs, hence the bins, hence every pixel's bin ids, so `reassign` has two paths: **relabel** (O(n_contigs); copies the file and rewrites four datasets, so it can also run `--inplace`) and **regroup** (O(nnz); the `from_cool` machinery — `_remap_bins` + `_stream_remapped_pixels` + `write_cooler`). `--layout auto` preserves whatever the source guaranteed, and a re-grouping that reproduces the existing order is detected and downgraded to the cheap path. `--layout keep` is the escape hatch and pays for it by writing `contig-ordering="custom"` — a value no writer other than `reassign` produces, and the *point* of which is that it is not in `GENOME_GROUPED_ORDERINGS`, so `validate()` stops enforcing a grouping the rows no longer have and `genome_extent` raises instead of returning a range that includes foreign bins. Two things a reassignment invalidates and must not carry over: the QC block (three of its four categories are defined by the assignment — recomputed by default, dropped under `--no-qc`) and per-genome metadata for genomes whose contig set changed (carried, but warned about by name). `--missing` chooses between an incremental edit (`keep`) and an authoritative table (`unbinned`); neither *drops* a contig — that is `extract`'s job, and doing it here would mean discarding pixels. `auto` reads only the *declared* ordering, so a `cool` source is relabelled even when its genomes are contiguous — every `collect`ed map's are — and `_kept_ordering` names any genome that scatters; the regroup path refuses a `square` map (remapping folds pixels into the upper triangle, as `from_cool` refuses for the same reason).
- **`/genomes` rows are derived from the contigs table**, one per distinct genome name; a `--genomes` table only left-joins extra columns onto them. It can neither add a genome with no contigs nor suppress one that has them, so `_report_genome_table_divergence` warns about both directions of mismatch — without it, MetaTOR's `MetaTOR_MGE_*` bins appear in `/genomes` unexplained and look fabricated.

## Testing

- Fixtures live in `tests/conftest.py`: `tiny_fasta`, `tiny_contigs_table`, `tiny_bins_bed`, `tiny_pairs`, `tiny_genomes_table`, `tiny_metacooler`, `tiny_cool`, `tiny_mmetacool`, plus variants for each contig-table format (`contigs`, `contig_genome`, `metator`, column-index) and edge cases (variable bins, reversed cool order).
- `tests/_cools.py` (`make_cool`, `ECOLI`, `BSUB`) builds tiny single-genome cools for the `collect` tests: two "species" sharing contig names, E. coli's deliberately not longest-first.
- Static test data is in `tests/data/`, regenerated with `make regenerate`.
- If you mark a known issue with `@pytest.mark.xfail`, do not remove the marker without fixing the underlying bug.
- **Every change to `src/` must keep patch coverage at or above 80%** — Codecov's `codecov/patch` status enforces it (`.github/codecov.yml`). Before pushing, run the two patch-coverage commands above and add tests for any uncovered changed line, error branches included; don't wait for Codecov to report it.
- CI's `test` job installs from `uv.lock`; the `dependency-bounds` job re-resolves at `lowest-direct` (py3.10) and `highest` (py3.14). To reproduce either locally, set `UV_RESOLUTION` for *every* uv command (`uv sync` silently re-locks when the lockfile's mode differs), and use a scratch copy so `uv.lock` is not rewritten.

## Reference

- `docs/cli.md` — **the** CLI reference (every command, options, and the API entry point each maps to). This is the single source of truth; `COMMANDS.md` was a competing, drifted copy and has been deleted.
- `docs/format.md` — the on-disk format specification.
- `benchmarks/` — metacooler vs cooler, timed from the CLI only (`benchmarks.qmd`, rendered by `.github/workflows/benchmarks.yml` on `v*` tags and on demand; see its README). The cooler it compares against is pinned as `COOLER_VERSION` in `bench.sh` (0.12.1), which the workflow installs over `uv.lock`'s and the notebook checks. Its datasets come from `simulate.py` and are deterministic per `GENERATOR_VERSION` — bump it whenever the model or a preset changes, since the workflow caches them in a `benchmark-data-v<N>` release.
- `reviewed/` — dated code-review notes from prior passes. Gitignored, so it exists only in a maintainer's local checkout.
