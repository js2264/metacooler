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

- **The fixtures** are in `tests/conftest.py`: `tiny_fasta`,
  `tiny_contigs_table`, `tiny_bins_bed`, `tiny_pairs`, `tiny_genomes_table`,
  `tiny_metacooler`, `tiny_cool`, `tiny_mmetacool`, and their variants for each
  contigs-table format and edge case. `tests/_cools.py` (`make_cool`, `ECOLI`,
  `BSUB`) builds the small single-genome cools of the `collect` tests. A
  fixture that many tests share is built once per module
  (`tmp_path_factory`).
- **The test data**, in `tests/data/`, is generated: `make regenerate`
  rebuilds it (`generate-raw-files.py`, then `raw-to-cool.sh`). No file there
  is over 5 MB, which pre-commit enforces.
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

## Changing things

- **Every commit is a [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/):**
  `type(scope): summary`, and a body that says why.
- **A change a user can notice** gets its entry under `## Unreleased` in
  `docs/changelog.md`, ending with a link to its pull request.
- **A pull request is reviewed as a whole, tested at its head, and squash-merged.**
- **Nothing a user calls is renamed:** commands, options, output files,
  modules. An alias is added instead.

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
