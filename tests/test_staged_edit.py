"""metacooler.staged_edit: several edits of a map through one staged copy.

Every writer given the session's staged map writes it directly -- no copy, no
validation, no lock of its own -- and the session validates it once and
renames it once, over the map or onto an output.
"""

from __future__ import annotations

import os
from pathlib import Path

import cooler
import h5py
import pytest
from test_staging import (
    WRITERS,
    _in_thread,
    _left_beside,
    _rounds,
    _store,
    _touch,
    _walk,
    holder,  # noqa: F401 - a fixture
)

import metacooler
from metacooler.api import MetaCooler
from metacooler.io import annotate as annotate_module
from metacooler.io import clusters as clusters_module
from metacooler.io import mmetacool as mmetacool_module
from metacooler.io import reassign as reassign_module
from metacooler.io import staging
from metacooler.io.clusters import write_clusters
from metacooler.io.reassign import reassign
from metacooler.io.staging import (
    MapChangedError,
    MapLockedError,
    holds_lock,
    in_session,
    lock,
    lock_path,
    staged_edit,
)
from metacooler.schema import validate, validate_mmetacooler


def _counted(
    monkeypatch: pytest.MonkeyPatch, modules, name: str, keep=lambda kwargs: True
) -> list:
    """Count the calls of *name* in each of *modules* (its own references),
    those *keep* accepts."""
    calls: list = []
    for module in modules:
        real = getattr(module, name, None)
        if real is None:
            continue

        def counted(*args, _real=real, **kwargs):
            if keep(kwargs):
                calls.append((args, kwargs))
            return _real(*args, **kwargs)

        monkeypatch.setattr(module, name, counted)
    return calls


def _full_validations(monkeypatch: pytest.MonkeyPatch) -> list:
    """Every full validation, wherever a writer calls it from."""
    import metacooler.schema as schema

    return _counted(
        monkeypatch,
        [schema, clusters_module, annotate_module, reassign_module, mmetacool_module],
        "validate",
        keep=lambda kwargs: kwargs.get("full", True),
    )


# ---------------------------------------------------------------------------
# One copy, one validation
# ---------------------------------------------------------------------------
class TestOneCopy:
    #: Writers that edit a map without re-grouping its contigs.
    EDITS = [
        "annotate",
        "store_weights",
        "write_clusters",
        "reassign, relabel",
        "qc --store",
        "metadata",
        "normalize",
        "drop_clusters",
    ]

    def test_several_writers_make_one_copy_and_one_validation(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        copies = _counted(monkeypatch, [staging], "copy_fresh")
        validations = _full_validations(monkeypatch)
        inode = tiny_metacooler.stat().st_ino
        with staged_edit(tiny_metacooler) as staged:
            assert in_session(staged) and holds_lock(staged)
            for writer in self.EDITS:
                WRITERS[writer](staged, tmp_path)
            # Nothing has replaced the map yet.
            assert tiny_metacooler.stat().st_ino == inode
            assert (len(copies), len(validations)) == (1, 0)
        assert (len(copies), len(validations)) == (1, 1)
        assert tiny_metacooler.stat().st_ino != inode
        assert not in_session(staged)
        mc = MetaCooler(tiny_metacooler)
        assert {"w", "ice"} <= set(mc.bins()[:].columns)
        assert "coverage" in mc.contigs.columns
        assert mc.metadata.get("sample") == "S1"
        assert "qc" in mc.metadata
        assert mc.cluster_rounds == []  # written, then dropped
        validate(tiny_metacooler, full=True)
        assert _left_beside(tiny_metacooler) == []

    def test_a_re_grouping_in_a_session(self, tiny_metacooler: Path):
        """A writer of a whole new map writes it beside the staged one, which it
        then replaces: later writers of the session see the new map."""
        with staged_edit(tiny_metacooler) as staged:
            reassign(staged, {"contig_A": "MAG_2"}, inplace=True)
            _store(staged)
            names = list(MetaCooler(staged).contigs["name"])
        mc = MetaCooler(tiny_metacooler)
        assert list(mc.contigs["name"]) == names
        assert "w" in mc.bins()[:].columns
        validate(tiny_metacooler, full=True)
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.parametrize("error", [RuntimeError, KeyboardInterrupt, SystemExit])
    def test_a_failure_mid_session_leaves_the_map_as_it_was(
        self, tiny_metacooler: Path, tmp_path: Path, error: type[BaseException]
    ):
        before = tiny_metacooler.read_bytes()
        with pytest.raises(error), staged_edit(tiny_metacooler) as staged:
            _store(staged)
            WRITERS["write_clusters"](staged, tmp_path)
            raise error("stop")
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_a_session_that_fails_validation_writes_nothing(
        self, tiny_metacooler: Path
    ):
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises(metacooler.schema.MetaCoolerFormatError),
            staged_edit(tiny_metacooler) as staged,
            h5py.File(staged, "r+") as f,
        ):
            del f["genomes"]
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_an_entry_in_place(
        self, tiny_mmetacool: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """An entry is edited in a staged copy of its collection, and that entry
        alone validated."""
        validations = _full_validations(monkeypatch)
        uri = f"{tiny_mmetacool}::entries/entry1"
        with staged_edit(uri) as staged:
            assert staged.endswith("::entries/entry1")
            _store(staged)
        assert [call[0][0] for call in validations] == [
            f"{staged.split('::')[0]}::entries/entry1"
        ]
        assert "w" in MetaCooler(uri).bins()[:].columns
        assert "w" not in MetaCooler(f"{tiny_mmetacool}::entries/entry2").bins()[:]
        validate_mmetacooler(tiny_mmetacool)

    def test_a_whole_collection_in_place(self, tiny_mmetacool: Path, tiny_metacooler):
        with staged_edit(tiny_mmetacool) as staged:
            metacooler.add_entry(Path(staged), "entry3", tiny_metacooler)
            metacooler.remove_entry(Path(staged), "entry2")
        assert metacooler.list_entries(tiny_mmetacool) == ["entry1", "entry3"]
        validate_mmetacooler(tiny_mmetacool)


# ---------------------------------------------------------------------------
# To an output
# ---------------------------------------------------------------------------
class TestOutput:
    def test_an_entry_copied_out_is_the_entry(self, tiny_mmetacool: Path, tmp_path):
        """An entry of a collection becomes a standalone map, as stored: an h5py
        walk of the output is a walk of the entry."""
        out = tmp_path / "entry1.metacool"
        with staged_edit(f"{tiny_mmetacool}::entries/entry1", output=out) as staged:
            assert "::" not in staged
        with h5py.File(out, "r") as f, h5py.File(tiny_mmetacool, "r") as g:
            assert _walk(f) == _walk(g["entries/entry1"], "/entries/entry1")
        validate(out, full=True)
        assert len(cooler.Cooler(str(out)).chroms()[:]) > 0
        assert _left_beside(out) == []

    def test_edits_go_to_the_output_only(self, tiny_metacooler: Path, tmp_path):
        before = tiny_metacooler.read_bytes()
        out = tmp_path / "edited.metacool"
        with staged_edit(tiny_metacooler, output=out) as staged:
            _store(staged)
            write_clusters(staged, *_rounds(staged))
        assert tiny_metacooler.read_bytes() == before
        mc = MetaCooler(out)
        assert "w" in mc.bins()[:].columns and mc.cluster_rounds
        assert out.stat().st_mode & 0o777 == 0o666 & ~staging._umask()

    def test_a_whole_collection_copied_out(self, tiny_mmetacool: Path, tmp_path):
        out = tmp_path / "copy.mmetacool"
        with staged_edit(tiny_mmetacool, output=out) as staged:
            _store(f"{staged}::entries/entry2")
        validate_mmetacooler(out)
        assert "w" in MetaCooler(f"{out}::entries/entry2").bins()[:].columns

    def test_an_existing_output_needs_force(self, tiny_metacooler: Path, tmp_path):
        out = tmp_path / "out.metacool"
        out.write_text("a previous output")
        with pytest.raises(FileExistsError), staged_edit(tiny_metacooler, output=out):
            pass
        assert out.read_text() == "a previous output"
        with staged_edit(tiny_metacooler, output=out, force=True) as staged:
            _store(staged)
        validate(out, full=True)

    def test_an_output_created_meanwhile_needs_force(
        self, tiny_metacooler: Path, tmp_path
    ):
        out = tmp_path / "out.metacool"
        with (
            pytest.raises(FileExistsError, match="created meanwhile"),
            staged_edit(tiny_metacooler, output=out),
        ):
            out.write_text("someone else's")
        assert out.read_text() == "someone else's"

    def test_an_output_link_is_replaced_not_followed(
        self, tiny_metacooler: Path, tmp_path
    ):
        target = tmp_path / "target.metacool"
        target.write_text("not to be written through")
        link = tmp_path / "link.metacool"
        link.symlink_to(target)
        with staged_edit(tiny_metacooler, output=link, force=True):
            pass
        assert not link.is_symlink()
        assert target.read_text() == "not to be written through"
        validate(link, full=True)

    def test_the_output_is_not_the_map(self, tiny_metacooler: Path, tmp_path):
        with (
            pytest.raises(ValueError, match="the file of the map"),
            staged_edit(tiny_metacooler, output=tiny_metacooler, force=True),
        ):
            pass
        with (
            pytest.raises(ValueError, match="names the group"),
            staged_edit(tiny_metacooler, output=f"{tmp_path}/x.mmetacool::entries/x"),
        ):
            pass

    def test_a_source_torn_while_copied_is_refused(
        self, tiny_metacooler: Path, tmp_path, monkeypatch: pytest.MonkeyPatch
    ):
        real = staging.copy_fresh

        def torn(source, target, **kwargs):
            real(source, target, **kwargs)
            _touch(Path(source))

        monkeypatch.setattr(staging, "copy_fresh", torn)
        out = tmp_path / "out.metacool"
        with (
            pytest.raises(MapChangedError, match="while it was copied"),
            staged_edit(tiny_metacooler, output=out),
        ):
            pass
        assert not out.exists() and _left_beside(out) == []


# ---------------------------------------------------------------------------
# The session's guarantees
# ---------------------------------------------------------------------------
class TestGuarantees:
    def test_a_write_by_another_program_meanwhile_is_not_lost(
        self, tiny_metacooler: Path
    ):
        """A program that takes no lock writes the map while the session runs:
        the session is refused, and the map keeps that program's write."""
        with (
            pytest.raises(MapChangedError, match="another program"),
            staged_edit(tiny_metacooler) as staged,
        ):
            _store(staged)
            with h5py.File(tiny_metacooler, "r+") as f:
                f.attrs["by-cooler-balance"] = 1
        with h5py.File(tiny_metacooler, "r") as f:
            assert f.attrs["by-cooler-balance"] == 1
            assert "w" not in f["bins"]
        assert _left_beside(tiny_metacooler) == []

    def test_a_changed_plan_is_refused(self, tiny_metacooler: Path):
        planned = MetaCooler(tiny_metacooler)
        _store(tiny_metacooler)
        with (
            pytest.raises(MapChangedError),
            staged_edit(tiny_metacooler, planned=planned),
        ):
            pass

    def test_another_thread_is_refused_the_staged_map(self, tiny_metacooler: Path):
        with staged_edit(tiny_metacooler) as staged:
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(lambda: _store(staged))
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(lambda: lock(staged).__enter__())
            seen: list[bool] = []
            _in_thread(lambda: seen.append(in_session(staged)))
            assert seen == [False]
            with lock(staged) as lock_file:  # this thread's: at no cost
                assert lock_file == lock_path(tiny_metacooler)
        validate(tiny_metacooler, full=True)

    def test_the_map_itself_is_refused_while_its_session_runs(
        self, tiny_metacooler: Path
    ):
        with staged_edit(tiny_metacooler) as staged:
            with pytest.raises(ValueError, match="write its staged map"):
                _store(tiny_metacooler)
            with (
                pytest.raises(ValueError, match="write its staged map"),
                MetaCooler(tiny_metacooler).open("r+"),
            ):
                pass
            _store(staged)
        assert "w" in MetaCooler(tiny_metacooler).bins()[:].columns

    @pytest.mark.usefixtures("holder")
    def test_another_process_is_refused(self, tiny_metacooler: Path):
        with (
            pytest.raises(MapLockedError, match="test holder"),
            staged_edit(tiny_metacooler),
        ):
            pass

    def test_a_read_only_map_is_refused(self, tiny_metacooler: Path):
        os.chmod(tiny_metacooler, 0o444)
        with pytest.raises(PermissionError), staged_edit(tiny_metacooler):
            pass

    def test_a_map_that_does_not_exist(self, tmp_path: Path):
        with (
            pytest.raises(FileNotFoundError, match="No map to edit"),
            staged_edit(tmp_path / "none.metacool"),
        ):
            pass

    def test_inside_a_held_lock(self, tiny_metacooler: Path):
        """A caller holding the map's lock for a whole run opens sessions in it."""
        with lock(tiny_metacooler, command="metator pipeline"):
            for name in ("a", "b"):
                with staged_edit(tiny_metacooler) as staged:
                    _store(staged, name=name)
            assert lock_path(tiny_metacooler).exists()
        assert {"a", "b"} <= set(MetaCooler(tiny_metacooler).bins()[:].columns)
        assert _left_beside(tiny_metacooler) == []


# ---------------------------------------------------------------------------
# Nested sessions: joined in place, refused otherwise
# ---------------------------------------------------------------------------
class TestNested:
    def test_a_nested_session_in_place_joins(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        copies = _counted(monkeypatch, [staging], "copy_fresh")
        with staged_edit(tiny_metacooler) as outer:
            with staged_edit(tiny_metacooler) as inner:
                assert inner == outer
                _store(inner, name="a")
            with staged_edit(outer) as again:
                assert again == outer
                _store(again, name="b")
            # Nothing written yet: the outer session writes.
            assert "a" not in MetaCooler(tiny_metacooler).bins()[:].columns
        assert len(copies) == 1
        assert {"a", "b"} <= set(MetaCooler(tiny_metacooler).bins()[:].columns)

    def test_other_nestings_are_refused(self, tiny_metacooler: Path, tmp_path: Path):
        with staged_edit(tiny_metacooler) as staged:
            with (
                pytest.raises(ValueError, match="only joins"),
                staged_edit(tiny_metacooler, output=tmp_path / "a.metacool"),
            ):
                pass
            with (
                pytest.raises(ValueError, match="only joins"),
                staged_edit(staged, output=tmp_path / "b.metacool"),
            ):
                pass
        with (
            staged_edit(tiny_metacooler, output=tmp_path / "c.metacool"),
            pytest.raises(ValueError, match="only joins"),
            staged_edit(tiny_metacooler),
        ):
            pass


def test_exported():
    assert metacooler.staged_edit is staged_edit


class TestCorners:
    def test_a_failed_re_grouping_in_a_session(self, tiny_metacooler: Path):
        """A writer of a whole new map that fails leaves the session's staged
        map as it was, and nothing beside it."""
        with staged_edit(tiny_metacooler) as staged:
            _store(staged)
            with (
                pytest.raises(RuntimeError),
                staging.staged_inplace(staged, copy=False) as fresh,
            ):
                Path(fresh).write_text("half a map")
                raise RuntimeError("stop")
            assert not Path(fresh).exists()
            assert "w" in MetaCooler(staged).bins()[:].columns
        validate(tiny_metacooler, full=True)

    def test_a_session_on_a_whole_collection_validates_it_all(
        self, tiny_mmetacool: Path, monkeypatch: pytest.MonkeyPatch
    ):
        import metacooler.schema as schema

        whole = _counted(monkeypatch, [schema, staging], "validate_mmetacooler")
        with staged_edit(tiny_mmetacool):
            pass
        assert len(whole) == 1

    def test_an_absolute_link_into_an_entry_copied_out(
        self, tiny_mmetacool: Path, tmp_path: Path
    ):
        with h5py.File(tiny_mmetacool, "r+") as f:
            del f["entries/entry1/chroms"]
            f["entries/entry1/chroms"] = h5py.SoftLink("/entries/entry1/contigs")
        out = tmp_path / "entry1.metacool"
        with staged_edit(f"{tiny_mmetacool}::entries/entry1", output=out):
            pass
        with h5py.File(out, "r") as f:
            assert f.get("chroms", getlink=True).path == "/contigs"
        validate(out, full=True)
