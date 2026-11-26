# Working on metacooler

The instructions for coding agents, here and nowhere else (CONVENTIONS.md
§ 21): Copilot and Claude Code read this file. Do not add a `CLAUDE.md`, which
would switch Claude Code off reading it; a personal one holding the line
`@AGENTS.md` is ignored by git.

metacooler extends [cooler](https://github.com/open2c/cooler)'s HDF5 format to
metagenomic Hi-C contact maps: a `/genomes` → `/contigs` layer over cooler's
chroms, bins and pixels. **Every `.metacool` file stays a cooler**
(`format="HDF5::Cooler"`, `chroms` a soft link to `contigs`,
`indexes/chrom_offset`) **and a metacool** (`metacooler-format="HDF5::MetaCooler"`,
`/genomes`): never break either. A `.mmetacool` holds maps as entries,
`file.mmetacool::entries/<name>`. Read [README.md](README.md), then
[the documentation](https://js2264.github.io/metacooler/); the project's
conventions are
[koszullab/.github's CONVENTIONS.md](https://github.com/koszullab/.github/blob/main/CONVENTIONS.md).

## Commands

```bash
uv sync                                         # the environment
uv run pytest                                   # the tests, with coverage
uv run pytest tests/test_merge.py -x            # one file, to its first failure
uv run ruff check . && uv run ruff format .     # lint and format (88 columns; no black)
uv run ty check src tests                       # types (ty, not mypy)
uv run pre-commit install                       # the hooks
uv run --only-group docs zensical serve         # the documentation, at http://localhost:8000
uv run diff-cover coverage.xml --compare-branch=origin/devel --fail-under=80   # the patch's coverage
make cli-tests                                  # every command from the shell
make regenerate                                 # tests/data/, generated again
make fix                                        # ty, ruff, and README.md rendered from README.qmd
uv run bump-my-version bump minor               # a release (docs/contributing.md)
```

uv, never pip, makes the environment; hatchling builds the package;
`pyproject.toml` gives every dependency, each floor with its reason, and
`uv.lock` pins them. ruff formats with double quotes and leaves README.md out.
Run pytest through `uv run`: some tests run cooler's command line, which it
puts on `PATH`. To reproduce CI's `bounds` job, set `UV_RESOLUTION` for every
uv command, in a scratch copy (docs/contributing.md, Testing).

## Layout

```text
src/metacooler/            the package: api.py (MetaCooler), schema.py, uri.py, filters.py, qc.py, ps.py
src/metacooler/io/         what reads and writes files
src/metacooler/ops/        what computes on maps already open
src/metacooler/cli/        one cmd_<name>.py per command; main.py registers them
src/metacooler/templates/  the report's page
src/metacooler/_report/    chromokit's report kit: never edited (chromokit sync replaces it)
tests/                     pytest; tests/data/ holds the test data, which `make regenerate` writes
docs/                      the site; docs/cli.md is the command line's one reference, docs/format.md the format's
benchmarks/                metacooler against cooler, timed from the command line (its README)
```

[docs/internals.md](docs/internals.md) holds every rule below in full, and why
it holds; [docs/contributing.md](docs/contributing.md) the procedures.

## Rules

The project's (CONVENTIONS.md § 21):

- **A pull request is reviewed as a whole, and tested at its head**: the
  maintainer squash-merges it. A review says how the whole change affects the
  package: its behaviour, outputs, command line, Python API and docs. This
  applies to the brief given to an independent reviewer too.
- **A change a user can notice carries its changelog entry**, under
  `## Unreleased` in `docs/changelog.md`, ending with a link to its pull
  request.
- **Every commit is a Conventional Commit:** `type(scope): summary`, and a body
  that says why.
- **Coverage stays at 80 %**, overall and on the lines a pull request changes,
  error branches included: check the patch with diff-cover before pushing.
- **The command line is tested in-process**, with click's `CliRunner`; a
  subprocess is for third-party programs only (cooler's command line).
- **Diagnostics go through logging** (`logging.getLogger("metacooler")`), never
  `print`. A command that writes a map or a file runs its work in
  `cli/logging_utils.run_log`, which writes its log beside the output; one that
  only reads logs to stderr only.
- **Nothing a user calls is renamed** (commands, options, output files, modules,
  documentation URLs): an alias is added, which warns.
- **chromokit's files are never edited**: `overrides/`,
  `docs/stylesheets/chromokit.css`, `docs/javascripts/chromokit.js` and
  `src/metacooler/_report/`, which `chromokit.json` lists. The template's files
  take their changes from `uvx copier update`.
- **README.md is generated** from README.qmd (`make fix`): edit README.qmd.
- **A known bug is an `xfail` test** naming its issue; its marker goes only with
  the fix (`xfail_strict` fails it once fixed).

metacooler's own, each in full in [docs/internals.md](docs/internals.md):

- **The layers:** `io/` reads and writes files, `ops/` computes, `cli/` parses
  and calls them. metacooler is **a layer over cooler**: bins, pixels and chroms
  are read through `cooler.Cooler`, coolers written through
  `io/_utils.write_cooler`, never `cooler.create_cooler` directly; raw h5py only
  for the genome layer and the offset indexes ([#the-code](docs/internals.md#the-code),
  [#write-cooler](docs/internals.md#write-cooler)).
- **`MetaCooler` builds its `cooler.Cooler` lazily**: keep the cheap paths
  (`info`, the genome and contig lookups) off `_clr`, and resolve a list of
  names in one pass ([#the-code](docs/internals.md#the-code)).
- **The format** is `docs/format.md`'s, which `schema.validate()` enforces; a
  writer validates in full what it writes, `MetaCooler()` checks the structure
  on every open. A change to what any writer writes changes both.
- **The public API is lazy**: a name in `__all__`, an `if name == ...` branch in
  `__getattr__`, a line in the `TYPE_CHECKING` block; never a barrel import
  ([#lazy-api](docs/internals.md#lazy-api)).
- **Commands** register themselves (`cmd_<name>.py`, an import line and a
  `COMMAND_GROUPS` entry in `cli/main.py`), with a one-line `help=`, no `cls=`,
  their heavy imports inside the function, a `MetaCoolURI` argument
  (`require_entry=True` for one map), and `cli/_options.py`'s shared options and
  `check_output(..., inputs=...)` ([#commands](docs/internals.md#commands),
  [#cli-helpers](docs/internals.md#cli-helpers)).
- **Datasets** of the genome layer go through `io/_utils.write_column`; pixels
  are summed with `aggregate_pixels`; no Python loop over contigs, bins or
  pixels; counts keep their dtype (`resolve_count_dtype`); every streaming
  default is `schema.DEFAULT_CHUNKSIZE` ([#datasets](docs/internals.md#datasets),
  [#summing-pixels](docs/internals.md#summing-pixels),
  [#no-contig-loops](docs/internals.md#no-contig-loops),
  [#streaming-creation](docs/internals.md#streaming-creation)).
- **Names are ASCII** (contigs, genomes, `/contigs` columns), checked three
  times ([#ascii](docs/internals.md#ascii)). A user's table is read with
  `read_user_table`, its key columns relabelled by position
  (`rename_key_columns`), which columns hold its names said by `io/colspec.py`
  ([#user-tables](docs/internals.md#user-tables), [#colspec](docs/internals.md#colspec)).
- **Selections:** `--filter` is one comparison per flag, its table prefix
  required; names choose, filters narrow (`resolve_contig_selection`)
  ([#filter](docs/internals.md#filter)).
- **A normalization is a weight column of `/bins`**, stored by `store_weights`
  with its `depends_on`; a writer that changes what a weight depends on drops it
  ([#normalizations](docs/internals.md#normalizations)).
- **An edit in place** takes the map's lock before its first read and writes a
  staged copy that replaces the map only once validated (`io/staging`); never
  `h5py.File(map, "r+")`. A writer of a new output goes through
  `writing_output` ([#in-place](docs/internals.md#in-place)).
- **Each writer's own invariants** (`collect`, `reassign`, `annotate`, `create`'s
  contigs table, `binnify --level contig`'s stream, `fasta`, `merge`,
  `/genomes`, `--nproc`) are in [Internals](docs/internals.md#rules): read the
  one you change.
- **Correctness, speed, the API, errors:** the HDF5 layer's links and offsets,
  no dense matrix past its bound, `force` or `FileExistsError`, deprecations
  that warn, messages that name what they refuse
  ([#more-rules](docs/internals.md#more-rules)).
- **The report** (`metacooler report`): `ops/report.collect` computes it,
  `io/report.write_report` writes it with the kit; a key, a column or a MultiQC
  column that changes or goes raises `ops/report.FORMAT_VERSION`, and
  `tests/data/assembly_1_1kb.report.txt` with it.
- **Style:** Python 3.10 and newer (`match`, `X | Y`, `from __future__ import
  annotations`); Google-style docstrings; ruff's rules, a deviation commented
  where it is.

## Testing

- The fixtures are in `tests/conftest.py` (`tiny_metacooler`, `tiny_mmetacool`,
  `tiny_pairs`, …) and `tests/_cools.py`; the static data in `tests/data/`.
- The command line through `CliRunner` (`result.exit_code`, `result.output`); an
  error with `pytest.raises(..., match=r"...")`; a shared fixture built once per
  module. Never the network.
- How to write one: [Writing a test](docs/contributing.md#writing-a-test).

## Procedures

- **Adding a command:** [docs/contributing.md](docs/contributing.md#adding-a-command):
  its code in `io/` or `ops/`, its `cmd_<name>.py`, its registration, its log,
  its API, its tests, its section in `docs/cli.md` (`tests/test_docs.py` fails
  without it) and its changelog entry.
- **Reviewing a pull request:** [docs/contributing.md](docs/contributing.md#reviewing-a-pull-request):
  as a whole, at its head, layer by layer, each finding with its place and its
  severity.
- **Planning a change** and **debugging a map:** the sections of the same names
  in [docs/contributing.md](docs/contributing.md).
- **Documenting:** docstrings, help strings and the words to use,
  [docs/contributing.md](docs/contributing.md#docstrings-help-and-words).
- **Releasing:** [docs/contributing.md](docs/contributing.md#releasing).
- **The benchmarks:** [benchmarks/README.md](benchmarks/README.md). They time
  uv.lock's cooler (0.12.1): `bench.sh`'s `COOLER_VERSION` moves with it.
  `simulate.py`'s `GENERATOR_VERSION` moves with its model or a preset (the
  `benchmark-data-v<N>` release caches the data). The page reaches the site
  through the `benchmark-page` release, which `docs.yml` serves at every
  build, since a deployment replaces the whole site.

## Reference

- `docs/cli.md`: **the** command line's reference, every command and option,
  and the Python function each maps to; `tests/test_docs.py` keeps it true.
- `docs/format.md`: the format's specification, the report's keys and the logs.
- The site is zensical with chromokit's look (`zensical.toml`, `overrides/`);
  the version in its header is `extra.release`, which bump-my-version updates.
- `reviewed/`: dated notes of earlier reviews, ignored by git: they exist only
  in a maintainer's checkout.
