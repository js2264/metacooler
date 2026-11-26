# Contributing

Bug reports, fixes, features and corrections to these pages are all welcome,
on [GitHub](https://github.com/js2264/metacooler). metacooler follows
the project's [conventions](https://github.com/koszullab/.github/blob/main/CONVENTIONS.md).

## Setting up

The environment is [uv](https://docs.astral.sh/uv/)'s; `uv.lock` pins every
dependency and is committed.

```bash
git clone https://github.com/js2264/metacooler.git
cd metacooler
uv sync                       # the environment: metacooler, editable, and the dev group
uv run pre-commit install     # the hooks: ruff and ty, on every commit
```

## Testing

```bash
uv run pytest                                       # the tests, with coverage
uv run pytest tests/test_create_correctness.py -x   # one file, up to its first failure
uv run pytest tests/test_merge.py::test_name -x     # one test
uv run ruff check . && uv run ruff format .         # lint and format
uv run ty check src tests                           # types
make cli-tests                                      # every command from the shell (tests/data/cli_tests.sh)
```

- **Run the tests through `uv run`**: some of them run cooler's own command
  line, which `uv run` puts on `PATH` (`.venv/bin/python -m pytest` alone fails
  them with `FileNotFoundError`).
- The command line is tested in-process, with click's `CliRunner`; a
  subprocess is for third-party programs only, and a test that needs a
  program that is not installed is skipped.
- A known bug is an `xfail` test naming its issue in `reason=`.
  `xfail_strict` fails the suite once the bug is fixed, so the marker leaves
  with the fix.
- Coverage stays at 80 %, overall and on the lines a pull request changes,
  error branches included: Codecov fails a pull request below either. Before
  pushing:

    ```bash
    uv run pytest
    uv run diff-cover coverage.xml --compare-branch=origin/devel --fail-under=80
    ```

- **The fixtures** are in `tests/conftest.py`:
  - `tiny_fasta`, four contigs: `contig_A` (1000 bp), `contig_B` (600),
    `contig_C` (400) and `contig_D` (200);
  - `tiny_contigs_table`, their assignment to two MAGs, `contig_D` left
    unbinned (`__unbinned__`); `tiny_genomes_table`, the genomes' metadata,
    extra columns included;
  - `tiny_bins_bed`, fixed 100 bp bins; `tiny_pairs`, about 200 synthetic
    contacts (numpy's seed 42);
  - `tiny_metacooler`, the map `from_pairs` makes of them; `tiny_cool`, the
    cool `cooler.create_cooler` makes of the same pairs and bins;
    `tiny_mmetacool`, a collection of two entries of `tiny_metacooler`;
  - variants for each contigs-table format (`contig_genome`, `metator`, by
    column index) and edge case (variable bins, a cool whose contigs are not
    in the map's order).

  `tests/_cools.py` (`make_cool`, `ECOLI`, `BSUB`) builds the small
  single-genome cools of the `collect` tests: two "species" that share contig
  names, as real assemblies do, E. coli's deliberately not longest-first. A
  fixture that many tests share is built once per module
  (`tmp_path_factory`).
- **The test data**, in `tests/data/`, is generated: `make regenerate`
  rebuilds it (`generate-raw-files.py`, then `raw-to-cool.sh`). Two
  assemblies (`assembly_1.fa`, `assembly_2.fa`, with their `.fai` and
  `.chrom.sizes`), two libraries of each (`contacts_<n>_rep<r>.pairs.gz`),
  their binnings (`contigs_table_<n>.tsv`, `mag_table_<n>.tsv`,
  `assembly_<n>_DpnII-digested.bed`, `metadata.yaml`), and the cools, mcools
  and maps made of them at 1 kb and by restriction fragment
  (`assembly_<n>_rep<r>_{1kb,DpnII}.*`, `assembly_<n>_1kb.metacool` and
  `.mmetacool`). `assembly_1_1kb.report.txt` is the reference of `metacooler
  report`. No file there is over 5 MB, which pre-commit enforces.
- **The suite's warnings** are cooler's, on the tiny `mergebuf` that two tests
  of `TestUnorderedWrites` force on purpose; a warning from anywhere else is
  worth a look.
- **The dependency bounds.** The `test` job installs from `uv.lock`; the
  `bounds` job re-resolves without it, at `lowest-direct` on Python 3.10 and
  at `highest` on the newest Python. To reproduce one, set `UV_RESOLUTION` for
  every uv command (`uv sync` silently locks again when the lock's mode
  differs), in a scratch copy, so that `uv.lock` is not rewritten:

    ```bash
    cp -r . /tmp/metacooler-lowest && cd /tmp/metacooler-lowest
    export UV_RESOLUTION=lowest-direct
    uv lock --upgrade --python 3.10
    uv sync --python 3.10
    uv run pytest -q --no-cov
    ```

### Writing a test

- **Group tests by what they test**, in classes (`class TestReassignLayout:`),
  in `tests/test_<module>.py`; the command line's in `test_cli*.py`
  (`test_cli_mmetacool.py` for the collections' commands), and the
  collections' `info`, `tree` and `attrs` in `test_info_mmetacool.py`.
- **A task of writing tests** reads the module it tests first, for its API,
  changes the tests and not that module, and reports the tests it wrote and
  their run.
- **Each feature is tested on its happy path, its edge cases (empty inputs,
  missing files, invalid formats), its errors (`pytest.raises(ValueError,
  match=r"...")`, a raw pattern) and what it must keep working**: the
  command line (`result.exit_code`, `result.output`; an unhandled exception
  is `result.exception`), and the files it writes, opened by `MetaCooler()`
  and by `cooler.Cooler()` both.
- **Fixtures go in `tests/conftest.py`**, or in a test file when they serve
  that file alone. A test never reaches the network.
- **A fixture that many tests share is built once** (`scope="module"`,
  `tmp_path_factory`).

## Documenting

The site is built from `docs/` by [zensical](https://zensical.org), in
[chromokit](https://github.com/js2264/chromokit)'s look, and published at
<https://js2264.github.io/metacooler/>.

```bash
uv run --only-group docs zensical serve     # http://localhost:8000
```

- **Pages:** the Markdown files of `docs/`, listed in order by `nav` in
  `zensical.toml`. A new page goes in both.
- **The command line:** `docs/cli.md`, written by hand, is its one reference.
  `tests/test_docs.py` fails when it misses a command or an option, or
  documents one that does not exist.
- **The format:** `docs/format.md` specifies `.metacool` and `.mmetacool`;
  `schema.validate()` enforces it. A change to what any writer writes changes
  both.
- **The Python API:** `::: metacooler.<module>` directives in `docs/api/`,
  which mkdocstrings expands from the docstrings (Google style).
- **The README** is quarto's render of `README.qmd`: edit `README.qmd`, never
  `README.md`, then `make fix` renders it (and runs ty and ruff).
- **The benchmarks:** [`benchmarks/README.md`](https://github.com/js2264/metacooler/tree/devel/benchmarks)
  says how they run and are published, at `/benchmarks/` on this site.
- **The look:** chromokit's files (`overrides/`, `docs/stylesheets/chromokit.css`,
  `docs/javascripts/chromokit.js`, listed in `chromokit.json`) are never
  edited; the site's own styles go in `docs/stylesheets/site.css`.

    ```bash
    npx github:js2264/chromokit sync     # a newer chromokit
    npx github:js2264/chromokit check    # nothing of chromokit's edited
    ```

### Docstrings, help and words

- **A module** starts with a one-line docstring saying what it does, not what
  it is: `"""URI parsing for metacool and mmetacool files."""`.
- **A class** has a one-sentence summary, then its public responsibilities
  and its key attributes; it does not repeat what its `__init__`'s own
  docstring says of its parameters (mkdocstrings merges the two).
- **A public function or method** has a Google-style docstring: a summary
  line, then `Args:`, `Returns:` (or `Yields:`) and `Raises:` sections, which
  name the actual parameters, the units, the formats and the valid ranges, and
  when each exception is raised. A private helper has a summary line at
  least, and a property one line.
- **`src/metacooler/__init__.py`'s** module docstring says how its lazy
  `__getattr__` exports the API ([Internals](internals.md#lazy-api)).
- **A command's `help=`** is one line, an imperative sentence, which `metacooler
  -h` lists; its details go in its options' help. **Every option and argument**
  has a `help=` that says what it controls, the format or the unit of its
  value, and its default where that is not obvious. A deprecated alias says
  so: `[DEPRECATED: use --genomes]`. An `epilog=` is short, its rich markup
  for emphasis only (`[bold]`, `[cyan]`).
- **Documentation never changes behaviour**: a pass over the docstrings
  changes no algorithm and no data flow, adds no type annotation, reformats
  no code (indentation, the order of imports) and touches `src/metacooler/`
  alone, its tests only when asked.
- **A pass over the documentation**: find the modules without a docstring,
  the options and arguments without a `help=`, the functions with a summary
  alone, and what is stale (a renamed option, an old parameter's name); write
  what the code does, having read it, one precise sentence rather than two
  vague ones; apply it, then read it again. Say which files changed and what
  in each, and what was left, and why.
- **One word for one thing**, everywhere:

    | Word | What it is |
    |---|---|
    | `.metacool` file | one contact map, of one genome or of many |
    | `.mmetacool` collection | several maps in one file, its entries |
    | entry | one map in a collection, `file.mmetacool::entries/<name>` |
    | contig | a sequence of the assembly: cooler's chrom |
    | genome | a MAG, or `__unbinned__`: one or more contigs |
    | bin | an interval of a contig, of fixed or variable width |
    | pixel | the count of contacts between two bins |
    | pair | a read pair, before it is binned |
    | resolution | the bin size, in bp |
    | `genome_id`, `bin1_id`, `bin2_id` | 0-based indices into `/genomes` and `/bins` |

## Changing things

- **Every commit is a [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/):**
  `type(scope): summary`, and a body that says why.
- **A change a user can notice** gets its entry under `## Unreleased` in
  `docs/changelog.md`, ending with a link to its pull request.
- **A pull request is reviewed as a whole, tested at its head, and squash-merged.**
- **Nothing a user calls is renamed:** commands, options, output files,
  modules. An alias is added instead.

## Adding a command

1. **Where its code goes.** What reads and writes files goes in
   `src/metacooler/io/<name>.py`; what computes on an open map in
   `src/metacooler/ops/<name>.py` (see [Internals](internals.md#the-code)).
2. **The command**, `src/metacooler/cli/cmd_<name>.py`:

    ```python
    """metacooler <name> command."""

    from __future__ import annotations

    from pathlib import Path

    import rich_click as click

    from metacooler.cli._options import (
        check_output,
        force_option,
        quiet_option,
        verbose_option,
    )
    from metacooler.cli.main import MetaCoolURI, cli


    @cli.command(name="<name>", help="Do something to a map, in one line.")
    @click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
    @click.option(
        "-o",
        "--output",
        type=click.Path(path_type=Path),
        required=True,
        help="Output file.",
    )
    @force_option
    @verbose_option
    @quiet_option
    def cmd_name(metacool_path, output, force, verbose, quiet):
        """What it does, in more words."""
        from metacooler.cli.logging_utils import setup_logging
        from metacooler.io.name import do_it  # heavy imports here: --help stays fast

        setup_logging(verbose, command_name="metacooler <name>", quiet=quiet)
        check_output(output, force, inputs=[metacool_path])
        do_it(metacool_path, output)
    ```

    - its `help=` is one line, and it takes no `cls=`: the group's
      `MetaCoolerCommand` aligns its help's panels;
    - a map is a `MetaCoolURI` argument, `require_entry=True` when it needs one
      map rather than a collection; an input path is `click.Path(exists=True,
      path_type=Path)`, an output `click.Path(path_type=Path)`;
    - the shared options are `cli/_options.py`'s: `-f/--force` with
      `check_output(..., inputs=...)`, `-v`, `-q`, `--nproc`, `--genomes` and
      `--contigs` (`parse_comma_list`), `--filter`;
    - it logs through `setup_logging`, to stderr only: metacooler writes no log
      file (it is not a pipeline, an exception to CONVENTIONS.md § 11).
3. **Its registration:** one line in the import block at the end of
   `cli/main.py`, in alphabetical order, and its name in the right panel of
   `COMMAND_GROUPS` there.
4. **Its Python API**, if it has one: a name in `__all__` of
   `src/metacooler/__init__.py`, one `if name == ...` branch in its
   `__getattr__`, and a line in its `TYPE_CHECKING` block. Never a barrel
   import.
5. **Its tests**, `tests/test_<name>.py`, with `tiny_metacooler`,
   `tiny_mmetacool` and the other fixtures, the command through `CliRunner`:
   80 % of its lines covered, its errors included.
6. **Its documentation:** its section in `docs/cli.md`, with a table of its
   options (`tests/test_docs.py` fails without it), a guide if it is a task
   of its own, its API in `docs/api/`, and its changelog entry.
7. **Checked:**

    ```bash
    uv run pytest tests/test_<name>.py -x
    uv run ruff check . && uv run ruff format .
    uv run ty check src tests
    uv run pytest                                  # everything, and the coverage
    ```

## Planning a change

Before writing a feature: find the modules it touches and the patterns of the
ones like it, list the files to create and to change and the tests it needs,
and how it will be checked (the commands above).

## Reviewing a pull request

A pull request is reviewed as a whole, at its head, whatever its commits: the
maintainer squash-merges it. A review says how the whole change affects
metacooler (its behaviour, its outputs, its command line, its Python API, its
documentation), and gives each finding its file, function and lines, and a
severity that reflects what it costs: 🔴 correctness, 🟠 speed, 🟡 the API,
🔵 the documentation. It reports; it does not rewrite. A finding is specific:
"this could be cleaner" is not one, nor is a matter of taste (a name, a
comment's wording) that harms neither correctness nor use. cooler itself is
out of a review's scope: it is reported upstream, never changed here. Every
step below is taken, however small the change: faults sit in code the
change calls as well as in the code it changes.

1. **Orient:** which layers the change touches. `cli → io/ops → api → schema`:
   the lower, the wider the harm. A fault in `cli/` is rarely wrong data; one
   in `io/create.py` can be wrong data that nothing notices.
2. **Correctness:** dead code (a result computed and dropped, code that a
   later change made inert; 🔴 when it was clearly meant to be used, a check
   whose result should have raised), return values against their annotations
   (a generator `yield`s; `-> None` returns nothing), exceptions caught too
   broadly, constant expressions (`x or True`), an `x or []` in a body where
   `x` may be empty, a raw string not parsed by `parse_uri`, a URI written
   with a separator other than `::`, a path handled as a `str` rather than a
   `Path`, and the rules of [Internals](internals.md)
   (the format's invariants, `validate()`, the in-place writes, ASCII names,
   cooler's writer through `write_cooler`).
3. **Speed:** [Internals](internals.md#speed): no Python loop over pixels,
   bins or contigs, no dense matrix past its bound, `chunksize` honoured.
4. **The API:** [Internals](internals.md#api-design), and nothing a user
   calls renamed.
5. **The documentation:** docstrings against what the code does (a table's
   columns, its header), annotations as precise as the values, messages that
   name what they refuse; `docs/cli.md` and `docs/format.md` up to date.

A review is written down as `review-<YYYYMMDD-hhmm>.md`: a summary (what was
reviewed, and how many findings of each severity, in a table), then each
finding (its title; its file, function and line; its severity; the step it
comes from; two to four sentences on the fault and what it costs; the code
as it is and as it should be), then a checklist of the findings, one line
each. The maintainer keeps them in `reviewed/`.

## Debugging a map

1. Look at it: `metacooler info`, `metacooler tree` and `metacooler attrs`.
2. Check it: `metacooler check MAP`, which runs the full `validate()`.
3. Then the usual suspects: a text column written as `object` rather than
   fixed-length bytes, `chroms` not a soft link to `contigs`, a pixel with
   `bin1_id > bin2_id` or out of order, an offset index of the wrong length
   (`chrom_offset` is `n_contigs + 1` long), a `genome_id` out of range.
4. Say what caused it and what fixes it, and the test that would have caught
   it.

## Releasing

Changelog entries go under `## Unreleased` in `docs/changelog.md` as their
pull requests land on `devel`. From an up-to-date `devel` with green CI:

```bash
uv run bump-my-version bump minor     # or patch, or major
git push --follow-tags
```

bump-my-version sets the version: `__version__` in
`src/metacooler/_version.py`, which is the only place the package's version is
written, `zensical.toml`'s `release` and `CITATION.cff`'s version and date.
It turns `## Unreleased` into the version's section,
`## [X.Y.Z] - YYYY-MM-DD`, commits `chore(release): vX.Y.Z` and tags `vX.Y.Z`.

The tag starts `release.yml`: CI on the tagged commit; the sdist and the wheel
built and checked (the tag is `v` and the wheel's version, the changelog has
the version's section, the sdist ships `tests/data/`, the wheel installs and
runs); PyPI, by [trusted publishing](https://docs.pypi.org/trusted-publishers/),
with no stored token; and the GitHub release, its notes the changelog's
section, the sdist and the wheel attached. `benchmarks.yml` runs on the tag
too: it attaches the benchmarks to the release and publishes their page.
`docs.yml` redeploys the site from `devel`'s push of the release commit, not
from the tag, which the `github-pages` environment refuses.

metacooler is not on bioconda yet. Once its recipe is in bioconda-recipes,
bioconda's bot picks up each release from PyPI; when the runtime dependencies
change, the recipe's are updated there the same day.

The project's template (chromokit's) brings its changes with `uvx copier update`.
