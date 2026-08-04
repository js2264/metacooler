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
uv sync --dev                                  # install with dev deps

uv run pytest tests/ -x -v                     # run tests (fail fast); --cov is on by default via pyproject
uv run pytest tests/test_create_correctness.py -x       # single test file
uv run pytest tests/test_create.py::test_name -x        # single test

uv run ruff check src/ tests/                  # lint
uv run ruff format src/ tests/                 # format (line-length 88, double quotes; no black)
uv run ty check src/ tests/                    # type check

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
├── io/          # file read/write: create, extract, dump, flatten, network, merge, mmetacool, contigs, genomes
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
- **Per-chunk cost is dominated by `cooler.create`'s fixed overhead**, which scales with *contig count*, not pixels (~9 s per call at 61 k contigs, measured, independent of nnz). That is why `--chunksize` defaults high (50 M pairs) and why the external path hands `_buffer_pixels`' folded accumulator back as one large frame rather than slicing it.
- **Never map bins per contig with a Python loop.** Both the tiling in `io/create._tile_bins` and the source→target mapping in `ops/binnify._map_old_bins_per_contig` used to iterate contigs — one materialised DataFrame per contig in the first, an `old_chrom == cid` mask per contig (hence O(n_contigs × n_bins)) in the second. Neither finished on a 1.9 M-contig metagenome. Both are now single vectorised passes built from `np.repeat` / `np.cumsum` / one global `np.searchsorted` whose key folds the contig id in above every coordinate, turning `(contig, start)` lexicographic order into plain integer order. Same shape applies to `_tile_uniform` and `_tile_per_genome_ligated`, which were already written this way.
- **Python ≥3.10**: use `match`, `X | Y` unions, `from __future__ import annotations`.
- **Comma-separated multi-value options** (`--genomes`, `--contigs`) use the shared `parse_comma_list` click callback in `cli/_options.py`, so a value can be repeated or given as a comma-joined list. An explicitly-supplied but empty value is rejected rather than silently disabling the filter.
- **Shared CLI helpers live in `cli/_options.py` and `cli/_provenance.py`** — `check_output(path, force)` is the single "refuse to clobber the output" guard (backed by `io/_utils.check_output_path`, which the ops layer uses directly); `load_metadata_file()` is the single ingest point for user `--metadata` files and rejects the reserved `commands` provenance key. Do not re-inline these.
- **`MetaCoolURI(require_entry=True)`** is the single policy for "this command needs one map, not a collection": it rejects a bare `.mmetacool` path up front and lists the available entry URIs. Commands that legitimately operate on a whole collection (`ls`, `info`, `tree`) leave it False.
- **`count` dtype**: never hard-cast pixel counts to `int32` when copying — aggregation (`binnify`, `merge`) writes `int64` and a blind cast wraps. Use `io/_utils.resolve_count_dtype()`.
- **Contig and genome names must be ASCII**, enforced by `io/_utils.check_ascii_names` (rules in `schema.ascii_violations`). `/contigs` is cooler's `/chroms` and cooler decodes *every* column of it with an ASCII `.astype("U")`, so one accented byte — in `name` or in a `--contigs` metadata column — yields a file that passes every structural check and then raises `UnicodeDecodeError` in `cooler.Cooler()`, hence in `MetaCooler()`. `/genomes/name` is checked for the mirror reason: cooler never reads that group, but `aggregate --level genome` hands those names to `cooler.create_cooler` as chrom names. Extra `/genomes` columns are deliberately *not* checked — accented taxonomy is legal metadata cooler never sees. The check runs three times on purpose: early in `from_pairs`/`from_cool` (before hours of binning), inside `_finalize_metacool_file` (so no writer can emit a bad file), and in `schema.validate()` (so a file that already exists is named rather than vouched for).
- **Never `pd.read_csv` a user-supplied text table directly** — go through `io/_utils.read_user_table()`. `read_csv` resolves a header/data field-count mismatch silently and in opposite directions: a header wider than the rows pads the *trailing* columns (every value ends up under a name to its left), while rows wider than the header absorb the surplus leading fields into the DataFrame *index* (shifting the other way). Both parse cleanly, and `keep_default_na=False` — which the readers need so a genome literally called `"NA"` survives — turns the pad into an empty string, so the frame keeps no record of how wide each line was. `read_user_table` counts fields on the raw bytes, rejects both shapes, and drops the all-empty `Unnamed: N` column a trailing tab on the header line creates.
- **`--filter` is one comparison per flag, ANDed by repetition** (`filters.py`), and the `genomes.`/`contigs.` table prefix is **mandatory** — both tables have a `length` and a `name`, and inferring the target only when it happens to be unambiguous would make the same expression legal on one file and rejected on the next. There is deliberately no `or` and no parenthesised grammar: real MetaTOR tables carry columns called `Present Markers` and `GC-content`, and a grammar with boolean connectives has to quote or backtick those before it can tokenise them. Here the *leftmost* operator terminates the column name, so nothing needs quoting (and a regex value may contain `>` or `=`). Missing values never match — including under `!=`/`!~` — because a filter selects rows it can positively vouch for.
- **Filters resolve to a contig-id set** via `io/_utils.resolve_contig_selection`, the single entry point for `--genomes` + `--contigs` + `--filter` on the read side (`dump`, `extract`). Names choose the candidate set, filters *narrow* it — an intersection, not a union. `ops/show.py` has its own `_apply_filters` because it needs ordered *names* and a selection *kind* rather than ids. `create` filters via `_apply_create_filters`, which builds the throwaway `/genomes` and `/contigs` frames with the very same `_build_genomes_table`/`_build_contigs_table` the output will use, so a filter written against a `dump` of an existing file means exactly the same thing at creation time; `genomes.n_bins` is the one absentee there and gets a specific error via the `hints` argument.
- **`/genomes` rows are derived from the contigs table**, one per distinct genome name; a `--genomes` table only left-joins extra columns onto them. It can neither add a genome with no contigs nor suppress one that has them, so `_report_genome_table_divergence` warns about both directions of mismatch — without it, MetaTOR's `MetaTOR_MGE_*` bins appear in `/genomes` unexplained and look fabricated.

## Testing

- Fixtures live in `tests/conftest.py`: `tiny_fasta`, `tiny_contigs_table`, `tiny_bins_bed`, `tiny_pairs`, `tiny_genomes_table`, `tiny_metacooler`, `tiny_cool`, `tiny_mmetacool`, plus variants for each contig-table format (`contigs`, `contig_genome`, `metator`, column-index) and edge cases (variable bins, reversed cool order).
- Static test data is in `tests/data/`, regenerated with `make regenerate`.
- If you mark a known issue with `@pytest.mark.xfail`, do not remove the marker without fixing the underlying bug.

## Reference

- `docs/cli.md` — **the** CLI reference (every command, options, and the API entry point each maps to). This is the single source of truth; `COMMANDS.md` was a competing, drifted copy and has been deleted.
- `docs/format.md` — the on-disk format specification.
- `reviewed/` — dated code-review notes from prior passes.
