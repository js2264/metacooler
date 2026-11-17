"""``keep_partial``: what a failed write leaves behind.

Every writer that stages its output -- built as ``<output>.partial`` and moved
into place once it validates -- takes ``keep_partial``.  ``True``, the default,
leaves the partial file for inspection and says so at ERROR; ``False`` removes
it whatever stopped the run, ``KeyboardInterrupt`` included, and says so at
INFO.  Either way the output path is never written by a failed run.
"""

from __future__ import annotations

import inspect
import logging
import os
from collections.abc import Callable
from functools import partial
from pathlib import Path

import pytest
from _cli_output import flat
from _cools import BSUB, ECOLI, make_cool

import metacooler.io.annotate as annotate_mod
import metacooler.io.collect as collect_mod
import metacooler.io.create as create_mod
import metacooler.io.reassign as reassign_mod
from metacooler.io.annotate import annotate
from metacooler.io.collect import collect
from metacooler.io.create import _staged_output, from_cool, from_pairs
from metacooler.io.reassign import reassign
from metacooler.qc import ContactQC
from metacooler.schema import Layout, validate


def _partial(path: Path) -> Path:
    return path.with_name(path.name + ".partial")


def _levels(caplog: pytest.LogCaptureFixture, text: str) -> set[str]:
    return {r.levelname for r in caplog.records if text in r.getMessage()}


# ---------------------------------------------------------------------------
# The staging context itself
# ---------------------------------------------------------------------------


class TestStagedOutput:
    def test_success_moves_the_file_into_place(self, tmp_path: Path):
        final = tmp_path / "out.metacool"
        with _staged_output(final, keep_partial=False) as staged:
            staged.write_bytes(b"done")
        assert final.read_bytes() == b"done"
        assert not _partial(final).exists()

    @pytest.mark.parametrize("error", [ValueError("boom"), KeyboardInterrupt()])
    def test_without_keep_partial_a_failure_removes_it(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture, error: BaseException
    ):
        caplog.set_level(logging.INFO, logger="metacooler")
        final = tmp_path / "out.metacool"
        with (
            pytest.raises(type(error)),
            _staged_output(final, keep_partial=False) as staged,
        ):
            staged.write_bytes(b"half")
            raise error
        assert not _partial(final).exists()
        assert not final.exists()
        assert _levels(caplog, "was removed") == {"INFO"}
        assert not _levels(caplog, "for inspection")

    def test_by_default_an_error_leaves_it_and_says_so(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ):
        """Guard: the default is the behaviour every writer always had."""
        final = tmp_path / "out.metacool"
        with pytest.raises(ValueError, match="boom"), _staged_output(final) as staged:
            staged.write_bytes(b"half")
            raise ValueError("boom")
        assert _partial(final).read_bytes() == b"half"
        assert not final.exists()
        assert _levels(caplog, "for inspection") == {"ERROR"}

    def test_by_default_an_interrupt_leaves_it_quietly(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ):
        """Guard: Ctrl-C was never reported as a failed creation."""
        final = tmp_path / "out.metacool"
        with pytest.raises(KeyboardInterrupt), _staged_output(final) as staged:
            staged.write_bytes(b"half")
            raise KeyboardInterrupt
        assert _partial(final).exists()
        assert not caplog.records

    @pytest.mark.parametrize("keep", [True, False])
    def test_a_failure_before_anything_was_written_says_nothing(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture, keep: bool
    ):
        """No partial file to keep or remove, so no line claiming either."""
        caplog.set_level(logging.INFO, logger="metacooler")
        final = tmp_path / "out.metacool"
        with pytest.raises(ValueError), _staged_output(final, keep_partial=keep):
            raise ValueError("before the first byte")
        assert not _partial(final).exists()
        assert not caplog.records


# ---------------------------------------------------------------------------
# Every public writer that stages its output
# ---------------------------------------------------------------------------


@pytest.fixture
def table(tmp_path: Path) -> Path:
    path = tmp_path / "moves.tsv"
    path.write_text("contig_name\tgenome_name\ncontig_C\tMAG_1\n")
    return path


@pytest.fixture
def genomes(tmp_path: Path) -> Path:
    path = tmp_path / "genomes.tsv"
    path.write_text("genome_name\tcompleteness\nMAG_1\t90\nMAG_2\t80\n")
    return path


@pytest.fixture
def writers(
    tmp_path: Path,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    tiny_cool: Path,
    tiny_metacooler: Path,
    table: Path,
    genomes: Path,
) -> dict[str, tuple[object, Callable[[Path, bool], None]]]:
    """Each writer, as the module it validates from and ``write(output, keep)``."""
    ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
    bsub = make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)
    return {
        "from_pairs": (
            create_mod,
            lambda out, keep: from_pairs(
                uri=out,
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                ps=False,
                keep_partial=keep,
            ),
        ),
        "from_cool": (
            create_mod,
            lambda out, keep: from_cool(
                tiny_cool,
                out,
                contig_genome_map=tiny_contigs_table,
                keep_partial=keep,
            ),
        ),
        "collect": (
            collect_mod,
            lambda out, keep: collect([ecoli, bsub], out, keep_partial=keep),
        ),
        "reassign": (
            reassign_mod,
            lambda out, keep: reassign(tiny_metacooler, table, out, keep_partial=keep),
        ),
        "annotate": (
            annotate_mod,
            lambda out, keep: annotate(
                tiny_metacooler, out, genomes_table=genomes, keep_partial=keep
            ),
        ),
    }


WRITERS = ["from_pairs", "from_cool", "collect", "reassign", "annotate"]


def _fail_validation(monkeypatch: pytest.MonkeyPatch, module: object, error) -> None:
    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(module, "validate", fail)


@pytest.mark.parametrize("writer", WRITERS)
@pytest.mark.parametrize("error", [ValueError("invalid"), KeyboardInterrupt()])
def test_without_keep_partial_a_failed_write_leaves_nothing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    writers,
    writer: str,
    error: BaseException,
):
    caplog.set_level(logging.INFO, logger="metacooler")
    module, write = writers[writer]
    _fail_validation(monkeypatch, module, error)
    out = tmp_path / "out.metacool"

    with pytest.raises(type(error)):
        write(out, False)

    assert not _partial(out).exists()
    assert not out.exists()
    assert _levels(caplog, "was removed") == {"INFO"}
    assert not _levels(caplog, "for inspection")


@pytest.mark.parametrize(
    "function", [from_pairs, from_cool, collect, reassign, annotate]
)
def test_every_writer_keeps_the_partial_file_by_default(function):
    default = inspect.signature(function).parameters["keep_partial"].default
    assert default is True


@pytest.mark.parametrize("writer", WRITERS)
def test_keeping_the_partial_file_is_what_a_failed_write_always_did(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    writers,
    writer: str,
):
    """``keep_partial=True`` leaves ``<output>.partial`` and says so, as before."""
    module, write = writers[writer]
    _fail_validation(monkeypatch, module, ValueError("invalid"))
    out = tmp_path / "out.metacool"

    with pytest.raises(ValueError, match="invalid"):
        write(out, True)

    assert _partial(out).exists()
    assert not out.exists()
    assert _levels(caplog, "for inspection") == {"ERROR"}


@pytest.mark.parametrize("writer", WRITERS)
def test_a_successful_write_is_the_same_either_way(
    tmp_path: Path, writers, writer: str
):
    """Guard: ``keep_partial`` only matters once something has failed."""
    _, write = writers[writer]
    out = tmp_path / "out.metacool"
    write(out, False)
    assert out.exists()
    assert not _partial(out).exists()


# ---------------------------------------------------------------------------
# Failures past the map's own writing: nothing is moved into place before them
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("keep", [True, False])
@pytest.mark.parametrize("error", [ValueError("no P(s)"), KeyboardInterrupt()])
def test_a_failed_ps_table_leaves_the_previous_map_and_table(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    keep: bool,
    error: BaseException,
):
    """The P(s) table is staged with the map, not written after it is in place."""
    out = tmp_path / "out.metacool"
    ps = tmp_path / "out.ps.tsv"
    out.write_bytes(b"previous map")
    ps.write_bytes(b"previous table")

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(create_mod, "write_ps_table", fail)
    with pytest.raises(type(error)):
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            force=True,
            keep_partial=keep,
        )

    assert out.read_bytes() == b"previous map"
    assert ps.read_bytes() == b"previous table"
    assert _partial(out).exists() is keep
    assert not _partial(ps).exists()


@pytest.mark.parametrize("keep", [True, False])
@pytest.mark.parametrize("error", [OSError("cannot move it"), KeyboardInterrupt()])
def test_a_map_that_cannot_be_moved_into_place_puts_the_table_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    keep: bool,
    error: BaseException,
):
    """The table is moved first; when the map cannot follow, that is undone.

    Left in place, the new table would sit beside the previous map, and with
    ``keep_partial=False`` that map's partial file removed: a table describing
    a map that exists nowhere.
    """
    out = tmp_path / "out.metacool"
    ps = tmp_path / "out.ps.tsv"
    out.write_bytes(b"previous map")
    ps.write_bytes(b"previous table")
    replace = os.replace

    def fail_for_the_map(source, destination):
        if Path(destination) == out:
            raise error
        replace(source, destination)

    monkeypatch.setattr(os, "replace", fail_for_the_map)
    with pytest.raises(type(error)):
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            force=True,
            keep_partial=keep,
        )

    assert out.read_bytes() == b"previous map"
    assert ps.read_bytes() == b"previous table"
    assert not list(tmp_path.glob("*.previous"))
    assert _partial(out).exists() is keep
    assert _partial(ps).exists() is keep


def test_the_previous_table_is_removed_once_the_map_is_in_place(
    tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, tiny_contigs_table: Path
):
    """Set aside while the map is moved, the previous table does not outlive it."""
    fresh, again = tmp_path / "fresh", tmp_path / "again"
    fresh.mkdir()
    again.mkdir()
    (again / "out.metacool").write_bytes(b"previous map")
    (again / "out.ps.tsv").write_bytes(b"previous table")

    for directory in (fresh, again):
        from_pairs(
            uri=directory / "out.metacool",
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            force=True,
        )

    assert sorted(p.name for p in again.iterdir()) == ["out.metacool", "out.ps.tsv"]
    ps = (again / "out.ps.tsv").read_bytes()
    assert ps == (fresh / "out.ps.tsv").read_bytes()
    validate(again / "out.metacool")


@pytest.mark.parametrize("ps", [True, False], ids=["ps", "no-ps"])
@pytest.mark.parametrize("table", [True, False], ids=["with-table", "no-table"])
def test_a_table_left_set_aside_is_refused_before_anything_is_read(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    table: bool,
    ps: bool,
):
    """Which table matches the map is unknown; setting one aside would lose one.

    A run killed between its moves leaves ``<name>.ps.tsv.previous`` beside the
    new table (its map still the previous one) or alone.
    """
    run = tmp_path / "run"
    run.mkdir()
    out = run / "out.metacool"
    previous = run / "out.ps.tsv.previous"
    out.write_bytes(b"previous map")
    if table:
        (run / "out.ps.tsv").write_bytes(b"new table")
    previous.write_bytes(b"previous table")
    before = {p.name: p.read_bytes() for p in run.iterdir()}

    def never(*args, **kwargs):
        raise AssertionError("read before the refusal")

    monkeypatch.setattr(create_mod, "_parse_contig_map", never)
    monkeypatch.setattr(create_mod, "_stream_pairs", never)
    with pytest.raises(FileExistsError, match="killed while it moved") as refused:
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            force=True,
            ps=ps,
        )

    assert str(previous) in str(refused.value)
    assert "If so, check which table matches" in str(refused.value)
    assert "If the file is one of yours, rename it." in str(refused.value)
    assert {p.name: p.read_bytes() for p in run.iterdir()} == before


def test_the_cli_refuses_a_table_left_set_aside(
    tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, tiny_contigs_table: Path
):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out = tmp_path / "out.metacool"
    previous = tmp_path / "out.ps.tsv.previous"
    previous.write_bytes(b"previous table")
    args = ["create", "pairs", "--pairs", str(tiny_pairs), "--bins", str(tiny_bins_bed)]
    args += ["--contigs", str(tiny_contigs_table), "-o", str(out)]

    result = CliRunner().invoke(cli, args)

    assert result.exit_code != 0
    assert "killed while it moved" in flat(result.output)
    assert previous.read_bytes() == b"previous table"
    assert not out.exists()


@pytest.mark.parametrize("error", [OSError("cannot remove it"), KeyboardInterrupt()])
def test_a_set_aside_table_that_cannot_be_removed_is_left_not_put_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    error: BaseException,
):
    """Once the map is in place the new map and table match: nothing is undone.

    The table left set aside then stops the next run until someone removes it.
    """
    fresh, run = tmp_path / "fresh", tmp_path / "run"
    fresh.mkdir()
    run.mkdir()
    out, ps = run / "out.metacool", run / "out.ps.tsv"
    previous = run / "out.ps.tsv.previous"
    out.write_bytes(b"previous map")
    ps.write_bytes(b"previous table")

    def create(uri: Path) -> None:
        from_pairs(
            uri=uri,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            force=True,
        )

    create(fresh / "out.metacool")
    unlink = Path.unlink

    def fail_for_the_set_aside_table(self: Path, *args, **kwargs) -> None:
        if self == previous:
            raise error
        unlink(self, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", fail_for_the_set_aside_table)
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        if isinstance(error, Exception):
            create(out)  # the output is complete: the run succeeds
        else:
            with pytest.raises(KeyboardInterrupt):
                create(out)
    monkeypatch.undo()

    assert ps.read_bytes() == (fresh / "out.ps.tsv").read_bytes()
    validate(out)
    assert previous.read_bytes() == b"previous table"
    assert sorted(p.name for p in run.iterdir()) == [
        "out.metacool",
        "out.ps.tsv",
        "out.ps.tsv.previous",
    ]
    assert _levels(caplog, str(previous)) == {"WARNING"}
    with pytest.raises(FileExistsError, match="killed while it moved"):
        create(out)


@pytest.mark.parametrize("keep", [True, False])
@pytest.mark.parametrize("layout", ["keep", "genome-then-alpha"])
def test_a_failed_qc_pass_leaves_the_previous_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tiny_metacooler: Path,
    keep: bool,
    layout: Layout,
):
    """reassign measures the QC block before it writes anything.

    So a failure in that pass, relabelling or re-grouping, leaves an existing
    output as it was and no partial file, whatever *keep_partial* says.
    """
    moves = tmp_path / "moves.tsv"
    # contig_A joins MAG_2: re-grouped, the contigs move.
    moves.write_text("contig_name\tgenome_name\ncontig_A\tMAG_2\n")
    out = tmp_path / "out.metacool"
    out.write_bytes(b"previous output")

    def fail(*args, **kwargs):
        raise RuntimeError("the QC pass failed")

    monkeypatch.setattr(ContactQC, "tally_pixels", fail)
    with pytest.raises(RuntimeError, match="QC pass"):
        reassign(
            tiny_metacooler, moves, out, layout=layout, force=True, keep_partial=keep
        )

    assert out.read_bytes() == b"previous output"
    assert not _partial(out).exists()


# ---------------------------------------------------------------------------
# A partial file given as the input: rescued, never removed as stale
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("command", ["reassign", "annotate"])
def test_a_partial_file_given_as_input_is_never_removed(
    tmp_path: Path, tiny_metacooler: Path, table: Path, genomes: Path, command: str
):
    """``<output>.partial`` used to be removed as stale before it was read."""
    out = tmp_path / "out.metacool"
    rescued = _partial(out)
    rescued.write_bytes(tiny_metacooler.read_bytes())
    before = rescued.read_bytes()

    with pytest.raises(ValueError, match="same file as the input"):
        if command == "reassign":
            reassign(rescued, table, out)
        else:
            annotate(rescued, out, genomes_table=genomes)

    assert rescued.read_bytes() == before
    assert not out.exists()


@pytest.mark.parametrize("command", ["reassign", "annotate"])
def test_the_cli_refuses_a_partial_file_as_input(
    tmp_path: Path, tiny_metacooler: Path, table: Path, genomes: Path, command: str
):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out = tmp_path / "out.metacool"
    rescued = _partial(out)
    rescued.write_bytes(tiny_metacooler.read_bytes())
    before = rescued.read_bytes()
    option = ["--contigs", str(table)] if command == "reassign" else []
    option += ["--genomes", str(genomes)] if command == "annotate" else []

    result = CliRunner().invoke(cli, [command, str(rescued), *option, "-o", str(out)])

    assert result.exit_code != 0
    assert rescued.read_bytes() == before
    assert not out.exists()


@pytest.mark.parametrize("writer", WRITERS)
def test_a_partial_file_given_as_input_is_refused_before_any_work(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    tiny_cool: Path,
    tiny_metacooler: Path,
    table: Path,
    genomes: Path,
    writer: str,
):
    """Where the output path is checked, before a pair, a bin or a pixel is read.

    ``_staged_output`` refuses it too, but only once a run gets there: after
    every pair is binned, after reassign's whole QC pass.
    """
    out = tmp_path / "out.metacool"
    rescued = _partial(out)

    def never(*args, **kwargs):
        raise AssertionError("read before the refusal")

    if writer == "from_pairs":
        rescued.write_bytes(tiny_pairs.read_bytes())
        monkeypatch.setattr(create_mod, "_stream_pairs", never)
        run = partial(
            from_pairs,
            uri=out,
            bins=tiny_bins_bed,
            pairs=rescued,
            contig_genome_map=tiny_contigs_table,
        )
    elif writer == "from_cool":
        rescued.write_bytes(tiny_cool.read_bytes())
        monkeypatch.setattr(create_mod, "_stream_remapped_pixels", never)
        run = partial(from_cool, rescued, out, contig_genome_map=tiny_contigs_table)
    elif writer == "collect":
        rescued.write_bytes(make_cool(tmp_path / "E.cool", ECOLI, seed=1).read_bytes())
        other = make_cool(tmp_path / "B.cool", BSUB, seed=2)
        monkeypatch.setattr(collect_mod, "_read_input", never)
        run = partial(collect, [rescued, other], out, names=["E", "B"])
    elif writer == "reassign":
        rescued.write_bytes(tiny_metacooler.read_bytes())
        monkeypatch.setattr(ContactQC, "tally_pixels", never)
        run = partial(reassign, rescued, table, out)
    else:
        rescued.write_bytes(tiny_metacooler.read_bytes())
        monkeypatch.setattr(annotate_mod, "_plan_map", never)
        run = partial(annotate, rescued, out, genomes_table=genomes)
    before = rescued.read_bytes()

    with pytest.raises(ValueError, match="same file as the input"):
        run()

    assert rescued.read_bytes() == before
    assert not out.exists()
