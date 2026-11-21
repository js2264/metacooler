"""metacooler.staged_edit: several edits of a map through one staged copy.

Every writer given the session's staged map writes it directly -- no copy, no
validation, no lock of its own -- and the session validates it once and
renames it once, over the map or onto an output.
"""

from __future__ import annotations

import contextlib
import errno
import json
import os
import shutil
import signal
from pathlib import Path

import cooler
import h5py
import numpy as np
import pytest
from _cli_output import flat
from _cools import BSUB, ECOLI, make_cool
from test_staging import (
    DEAD_RUN,
    SPAWN,
    WRITERS,
    _annotate,
    _in_thread,
    _left_beside,
    _rounds,
    _store,
    _touch,
    _walk,
    _write_lock,
    holder,  # noqa: F401 - a fixture
)

import metacooler
from metacooler.api import MetaCooler
from metacooler.io import annotate as annotate_module
from metacooler.io import clusters as clusters_module
from metacooler.io import mmetacool as mmetacool_module
from metacooler.io import reassign as reassign_module
from metacooler.io import staging
from metacooler.io import weights as weights_module
from metacooler.io.clusters import write_clusters
from metacooler.io.reassign import reassign
from metacooler.io.staging import (
    MapChangedError,
    MapLockedError,
    MapLockError,
    edit_metadata,
    holds_lock,
    in_session,
    lock,
    lock_path,
    staged_edit,
    staged_inplace,
)
from metacooler.io.weights import store_weights
from metacooler.schema import MetaCoolerFormatError, validate, validate_mmetacooler


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
                assert lock_file == lock_path(staged)
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
        with staged_edit(tiny_mmetacool) as staged, h5py.File(staged, "r+") as f:
            f.attrs["note"] = "written"
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


# ---------------------------------------------------------------------------
# The third review: what a session lets through
# ---------------------------------------------------------------------------
def _bin_keys(path: str | Path) -> np.ndarray:
    """A number naming each bin of a map: its contig's last letter, its start."""
    bins = MetaCooler(path).bins()[:]
    letter = bins["chrom"].astype(str).str[-1].map(ord).to_numpy()
    return (letter * 1000 + bins["start"].to_numpy() // 100).astype(float)


def _put_weights(path: str | Path, weights: np.ndarray, source: MetaCooler) -> None:
    store_weights(path, "w", weights, method="test", depends_on="pixels", source=source)


class TestSourceBins:
    """A map read for an edit that is not the map written must have its bins."""

    def test_a_source_read_before_a_regrouping_is_refused(self, tiny_metacooler: Path):
        with lock(tiny_metacooler):  # metaTOR holds it for its run
            source = MetaCooler(tiny_metacooler)
            weights = _bin_keys(tiny_metacooler)
            with staged_edit(tiny_metacooler) as staged:
                reassign(staged, {"contig_A": "MAG_2"}, inplace=True)  # re-groups
                with pytest.raises(MapChangedError, match="bins of the map written"):
                    _put_weights(staged, weights, source)
                # Computed from the staged map: its bins, written.
                _put_weights(staged, _bin_keys(staged), MetaCooler(staged))
        bins = MetaCooler(tiny_metacooler).bins()[:]
        assert list(MetaCooler(tiny_metacooler).contigs["name"]) == [
            "contig_D",
            "contig_B",
            "contig_A",
            "contig_C",
        ]
        np.testing.assert_array_equal(bins["w"], _bin_keys(tiny_metacooler))

    def test_a_source_whose_bins_hold_is_written(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        source = MetaCooler(tiny_metacooler)
        weights = _bin_keys(tiny_metacooler)
        with staged_edit(tiny_metacooler) as staged:
            _annotate(staged, tmp_path)
            reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
            _put_weights(staged, weights, source)
        bins = MetaCooler(tiny_metacooler).bins()[:]
        np.testing.assert_array_equal(bins["w"], _bin_keys(tiny_metacooler))

    def test_a_joined_session_planned_before_a_regrouping_is_refused(
        self, tiny_metacooler: Path
    ):
        planned = MetaCooler(tiny_metacooler)
        with staged_edit(tiny_metacooler) as staged:
            reassign(staged, {"contig_A": "MAG_2"}, inplace=True)
            with (
                pytest.raises(MapChangedError, match="bins of the map written"),
                staged_edit(tiny_metacooler, planned=planned),
            ):
                pass
            with staged_edit(tiny_metacooler, planned=MetaCooler(staged)) as again:
                assert again == staged
        validate(tiny_metacooler, full=True)

    def test_another_file_with_other_bins_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        regrouped = tmp_path / "regrouped.metacool"
        reassign(tiny_metacooler, {"contig_A": "MAG_2"}, output=regrouped)
        source = MetaCooler(tiny_metacooler)
        with pytest.raises(MapChangedError, match="bins of the map written"):
            _put_weights(regrouped, _bin_keys(tiny_metacooler), source)
        copy = tmp_path / "copy.metacool"
        shutil.copy(tiny_metacooler, copy)
        _put_weights(copy, _bin_keys(tiny_metacooler), source)
        np.testing.assert_array_equal(MetaCooler(copy).bins()[:]["w"], _bin_keys(copy))

    def test_another_entry_with_other_bins_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        from metacooler.io.mmetacool import create_mmetacool

        regrouped = tmp_path / "regrouped.metacool"
        reassign(tiny_metacooler, {"contig_A": "MAG_2"}, output=regrouped)
        collection = tmp_path / "both.mmetacool"
        create_mmetacool(
            collection, {"a": tiny_metacooler, "b": regrouped, "c": tiny_metacooler}
        )
        source = MetaCooler(f"{collection}::entries/a")
        weights = _bin_keys(source.uri)
        with pytest.raises(MapChangedError, match="bins of the map written"):
            _put_weights(f"{collection}::entries/b", weights, source)
        _put_weights(f"{collection}::entries/c", weights, source)
        written = MetaCooler(f"{collection}::entries/c").bins()[:]["w"]
        np.testing.assert_array_equal(written, weights)

    def test_an_entry_copied_out_is_compared_with_its_root(
        self, tiny_mmetacool: Path, tmp_path: Path
    ):
        entry = f"{tiny_mmetacool}::entries/entry1"
        source = MetaCooler(entry)
        out = tmp_path / "entry1.metacool"
        with staged_edit(entry, output=out) as staged:
            with staged_edit(staged, planned=source):
                pass
            _put_weights(staged, _bin_keys(entry), source)
        np.testing.assert_array_equal(MetaCooler(out).bins()[:]["w"], _bin_keys(out))

    def test_a_new_map_written_from_another_is_not_compared(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """metaTOR's pattern: a whole new map (``copy=False``) written from a
        map read for it has whatever bins its writer gives it."""
        source = MetaCooler(tiny_metacooler)
        out = tmp_path / "new.metacool"
        for _ in range(2):  # no map there, then a map of other bins
            with staged_inplace(out, copy=False, planned=source) as staged:
                reassign(tiny_metacooler, {"contig_A": "MAG_2"}, output=staged)
            validate(out, full=True)


class TestWritesOfTheMapOfASession:
    """A write of the map a session copies out, or of its output, is refused."""

    def test_the_map_copied_out(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "edited.metacool"
        before = tiny_metacooler.read_bytes()
        writes = {
            "store_weights": lambda: _store(tiny_metacooler),
            "annotate": lambda: _annotate(tiny_metacooler, tmp_path),
            "metadata": lambda: edit_metadata(
                tiny_metacooler, lambda meta: meta.update({"sample": "S9"})
            ),
            "open r+": lambda: MetaCooler(tiny_metacooler).open("r+").__enter__(),
        }
        with staged_edit(tiny_metacooler, output=out) as staged:
            for name, write in writes.items():
                with pytest.raises(ValueError, match="being copied to") as caught:
                    write()
                assert staged in str(caught.value), name
            _store(staged)
        assert tiny_metacooler.read_bytes() == before
        assert "w" in MetaCooler(out).bins()[:].columns
        assert _left_beside(tiny_metacooler) == []

    def test_the_output_it_replaces(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "edited.metacool"
        shutil.copy(tiny_metacooler, out)
        with staged_edit(tiny_metacooler, output=out, force=True) as staged:
            with pytest.raises(ValueError, match="replaced by the output"):
                _store(out, "x")
            with pytest.raises(ValueError, match="only joins"):
                staged_edit(out).__enter__()
            _store(staged)
        columns = MetaCooler(out).bins()[:].columns
        assert "w" in columns and "x" not in columns
        assert _left_beside(out) == []


@pytest.mark.skipif(not hasattr(os, "fork"), reason="needs os.fork")
@pytest.mark.filterwarnings(
    "ignore:This process .* is multi-threaded:DeprecationWarning"
)
class TestTheStagedMapsLock:
    """The staged map has a lock file of its own, naming the session."""

    def test_named_and_released(self, tiny_metacooler: Path):
        with staged_edit(tiny_metacooler) as staged:
            own = lock_path(staged)
            holder = json.loads(own.read_text())
            assert holder["pid"] == os.getpid()
            assert holder["staged_for"] == str(tiny_metacooler.resolve())
            with lock(staged) as lock_file:
                assert lock_file == own
        assert not own.exists()
        assert _left_beside(tiny_metacooler) == []

    def test_a_forked_child_is_refused(self, tiny_metacooler: Path):
        readable, writable = os.pipe()
        with staged_edit(tiny_metacooler) as staged:
            inode = Path(staged).stat().st_ino
            pid = os.fork()
            if pid == 0:  # the child: report what a write did, and leave
                try:
                    _store(staged, "child")
                    said = "written"
                except BaseException as error:  # noqa: BLE001 - reported
                    said = f"{type(error).__name__}: {error}"
                os.write(writable, said.encode()[:4000])
                os._exit(0)
            os.waitpid(pid, 0)
            os.close(writable)
            said = os.read(readable, 4096).decode()
            os.close(readable)
            assert said.startswith("MapLockedError"), said
            assert "staged map of an edit session of another process" in said
            assert Path(staged).stat().st_ino == inode
            _store(staged, "parent")
        columns = MetaCooler(tiny_metacooler).bins()[:].columns
        assert "parent" in columns and "child" not in columns

    def test_another_process_is_refused(self, tiny_metacooler: Path):
        from _staging_workers import try_store

        with staged_edit(tiny_metacooler) as staged, SPAWN.Pool(1) as pool:
            said = pool.apply(try_store, (staged,))
        assert said.startswith("MapLockedError"), said
        assert f"writes it to {tiny_metacooler.resolve()}" in said
        assert "w" not in MetaCooler(tiny_metacooler).bins()[:].columns


@pytest.mark.filterwarnings(
    "ignore:This process .* is multi-threaded:DeprecationWarning"
)
class TestOutputSessionLock:
    """An output session holds the lock of its output, not of the map."""

    def test_the_lock_of_the_output(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "out" / "edited.metacool"
        out.parent.mkdir()
        beside_the_map = sorted(os.listdir(tiny_metacooler.parent))

        def lock_briefly(path: Path) -> None:
            with lock(path):
                pass

        with staged_edit(tiny_metacooler, output=out) as staged:
            holder = json.loads(lock_path(out).read_text())
            assert holder["pid"] == os.getpid()
            assert Path(staged).parent == out.with_name(
                f".edited.metacool.{holder['uuid']}.1.d"
            )
            assert sorted(os.listdir(tiny_metacooler.parent)) == beside_the_map
            _in_thread(lambda: lock_briefly(tiny_metacooler))  # the map is free
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(lambda: lock_briefly(out))
            _store(staged)
        assert "w" in MetaCooler(out).bins()[:].columns
        assert _left_beside(out) == []
        assert sorted(os.listdir(tiny_metacooler.parent)) == beside_the_map

    @pytest.mark.usefixtures("holder")
    def test_a_holder_of_the_map_does_not_block_it(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "edited.metacool"
        with staged_edit(tiny_metacooler, output=out) as staged:
            _store(staged)
        assert "w" in MetaCooler(out).bins()[:].columns

    def test_a_writer_of_the_output_is_refused_meanwhile(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "edited.metacool"
        shutil.copy(tiny_metacooler, out)
        with staged_edit(tiny_metacooler, output=out, force=True) as staged:
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(lambda: _store(out, "x"))
            assert Path(staged).exists()
            _store(staged)
        columns = MetaCooler(out).bins()[:].columns
        assert "w" in columns and "x" not in columns
        assert _left_beside(out) == []

    @pytest.mark.skipif(not hasattr(os, "fork"), reason="needs os.fork")
    def test_what_a_killed_session_left_is_removed_by_the_next(
        self, tiny_metacooler: Path, tmp_path: Path, caplog
    ):
        out = tmp_path / "out" / "edited.metacool"
        out.parent.mkdir()
        pid = os.fork()
        if pid == 0:  # the child: killed inside its session
            with staged_edit(tiny_metacooler, output=out):
                os.kill(os.getpid(), signal.SIGKILL)
        os.waitpid(pid, 0)
        left = [n for n in os.listdir(out.parent) if n.startswith(".edited.")]
        assert lock_path(out).exists() and len(left) == 1
        assert not lock_path(tiny_metacooler).exists()
        with (
            caplog.at_level("WARNING", logger="metacooler"),
            staged_edit(tiny_metacooler, output=out) as staged,
        ):
            _store(staged)
        assert "Took over" in caplog.text
        assert f"Removed {left[0]}, left beside {out}" in caplog.text
        assert _left_beside(out) == []
        validate(out, full=True)

    def test_a_refusal_names_what_is_staged_beside_the_output(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "out" / "edited.metacool"
        out.parent.mkdir()
        _write_lock(out, 4321, host="node12", run=DEAD_RUN)
        with (
            pytest.raises(MapLockedError) as caught,
            staged_edit(tiny_metacooler, output=out),
        ):
            pass
        message = str(caught.value)
        assert str(lock_path(out)) in message
        assert str(out.parent / f".edited.metacool.{DEAD_RUN}.*.d/") in message

    def test_a_missing_output_directory_is_named(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "missing" / "edited.metacool"
        with (
            pytest.raises(FileNotFoundError) as caught,
            staged_edit(tiny_metacooler, output=out),
        ):
            pass
        assert str(caught.value).count(str(out)) == 1
        assert ".edited.metacool." not in str(caught.value)

    def test_an_unwritable_output_directory_is_named(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        out = tmp_path / "edited.metacool"
        real = staging._create

        def refused(path: Path, record: dict) -> bool:
            if path == lock_path(out):
                raise PermissionError(errno.EACCES, "Permission denied", str(path))
            return real(path, record)

        monkeypatch.setattr(staging, "_create", refused)
        with (
            pytest.raises(MapLockError) as caught,
            staged_edit(tiny_metacooler, output=out),
        ):
            pass
        assert str(caught.value).startswith(f"{out} cannot be written")
        assert "in-place" not in str(caught.value)

    @pytest.mark.skipif(
        not hasattr(os, "geteuid") or os.geteuid() == 0,
        reason="root writes a directory without write permission",
    )
    def test_a_map_in_a_read_only_directory(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "out" / "edited.metacool"
        out.parent.mkdir()
        directory = tiny_metacooler.parent
        mode = directory.stat().st_mode
        directory.chmod(0o555)
        try:
            with staged_edit(tiny_metacooler, output=out) as staged:
                _store(staged)
        finally:
            directory.chmod(mode)
        assert "w" in MetaCooler(out).bins()[:].columns


class TestFailedWrites:
    """A write that fails with the staged map in hand fails the session."""

    @staticmethod
    def _disk_full(*args, **kwargs):
        raise OSError(errno.ENOSPC, "No space left on device")

    def test_write_clusters(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises(RuntimeError, match="its error was caught") as caught,
            staged_edit(tiny_metacooler) as staged,
        ):
            _store(staged)
            with monkeypatch.context() as patch:
                patch.setattr(clusters_module, "_check_clusters", self._disk_full)
                with pytest.raises(OSError, match="No space"):
                    write_clusters(staged, *_rounds(staged))
            with h5py.File(staged, "r") as f:
                assert ".clusters-scratch" not in f
        assert isinstance(caught.value.__cause__, OSError)
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_store_weights_keeps_the_column_it_replaces(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        _store(tiny_metacooler, "w", 1.0)
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises(RuntimeError, match="its error was caught"),
            staged_edit(tiny_metacooler) as staged,
        ):
            with monkeypatch.context() as patch:
                patch.setattr(weights_module, "write_column", self._disk_full)
                with pytest.raises(OSError, match="No space"):
                    _store(staged, "w", 2.0)
            with h5py.File(staged, "r") as f:
                assert set(f["bins"]) >= {"w"}
                assert not any(name.startswith(".") for name in f["bins"])
                assert (f["bins/w"][()] == 1.0).all()
        assert tiny_metacooler.read_bytes() == before

    def test_a_write_of_its_own_and_a_joined_session(self, tiny_metacooler: Path):
        before = tiny_metacooler.read_bytes()
        # The caught error leaves each ``with`` but the outer two.
        with (
            pytest.raises(RuntimeError, match="its error was caught"),
            staged_edit(tiny_metacooler) as staged,
            contextlib.suppress(KeyError),
            MetaCooler(staged).open("r+") as grp,
        ):
            grp.attrs["half"] = "written"
            raise KeyError("stop")
        assert tiny_metacooler.read_bytes() == before
        with (
            pytest.raises(RuntimeError, match="its error was caught"),
            staged_edit(tiny_metacooler),
            contextlib.suppress(KeyError),
            staged_edit(tiny_metacooler),
        ):
            raise KeyError("stop")
        assert tiny_metacooler.read_bytes() == before

    def test_an_input_refused_before_the_write_fails_nothing(
        self, tiny_metacooler: Path
    ):
        with staged_edit(tiny_metacooler) as staged:
            with pytest.raises(ValueError):
                store_weights(
                    staged, "w", np.ones(3), method="test", depends_on="pixels"
                )
            with pytest.raises(ZeroDivisionError):
                edit_metadata(staged, lambda meta: meta.update({"x": 1 / 0}))
            _store(staged)
        assert "w" in MetaCooler(tiny_metacooler).bins()[:].columns


class TestThirdReviewCorners:
    def test_a_session_that_writes_nothing_leaves_the_map(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        validations = _full_validations(monkeypatch)
        reader = MetaCooler(tiny_metacooler)
        inode = tiny_metacooler.stat().st_ino
        with staged_edit(tiny_metacooler) as staged:
            assert len(MetaCooler(staged).contigs)  # a read is no write
        assert tiny_metacooler.stat().st_ino == inode
        assert validations == []
        assert len(reader.pixels())  # still the map it opened
        assert _left_beside(tiny_metacooler) == []

    def test_a_copy_torn_by_another_program_is_refused_at_once(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        real = staging.copy_fresh

        def torn(source, target, **kwargs):
            real(source, target, **kwargs)
            _touch(Path(source))  # written meanwhile, by a program without lock

        monkeypatch.setattr(staging, "copy_fresh", torn)
        ran: list[bool] = []
        with (
            pytest.raises(MapChangedError, match="while it was copied"),
            staged_edit(tiny_metacooler),
        ):
            ran.append(True)
        assert ran == []
        with pytest.raises(MapChangedError, match="while it was copied"):
            _store(tiny_metacooler)
        assert _left_beside(tiny_metacooler) == []


# ---------------------------------------------------------------------------
# Every install writes through a fresh file; nothing records a staged path
# ---------------------------------------------------------------------------
def _fresh_size(path: Path, tmp_path: Path) -> int:
    """The size of a fresh copy of the file *path*."""
    fresh = tmp_path / "fresh-copy.h5"
    staging.copy_fresh(path, fresh)
    try:
        return fresh.stat().st_size
    finally:
        fresh.unlink()


def _dead_space(path: str | Path, megabytes: int = 4) -> None:
    """Leave *megabytes* of space no object holds in the file *path*, as an
    edit that replaced a column did."""
    rng = np.random.default_rng(0)
    with h5py.File(path, "r+") as f:
        f.create_dataset("replaced", data=rng.random(megabytes * 131_072))
        # Written after it, so that the space it leaves is not at the end of
        # the file, which HDF5 would give back.
        bins = f["bins"]
        f.copy(bins["start"], bins, name="start.moved")
        del bins["start"]
        bins.move("start.moved", "start")
        del f["replaced"]
    assert Path(path).stat().st_size > megabytes * 1_000_000


def _writes(staged: str, tmp_path: Path) -> None:
    """A metaTOR step: annotate, weights, rounds, a relabelling."""
    _annotate(staged, tmp_path)
    _store(staged)
    write_clusters(staged, *_rounds(staged))
    reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")


def _metators_session(path: Path, tmp_path: Path) -> None:
    """metaTOR's step: the map reassigned into the staged map, then annotated,
    weighted and given its rounds."""
    with staged_edit(path, planned=MetaCooler(path)) as staged:
        reassign(path, _assignment(tmp_path), output=staged, force=True)
        _annotate(staged, tmp_path)
        _store(staged)
        write_clusters(staged, *_rounds(staged))


class TestBoundedMaps:
    """A map edited again and again does not grow: every staged map starts as
    a fresh copy, so a map holds at most the space one edit, or one session,
    replaced -- never what earlier ones did."""

    def test_sessions_in_a_row(self, tiny_metacooler: Path, tmp_path: Path):
        _dead_space(tiny_metacooler)
        dead = []
        for _ in range(6):
            _metators_session(tiny_metacooler, tmp_path)
            size = tiny_metacooler.stat().st_size
            dead.append(size - _fresh_size(tiny_metacooler, tmp_path))
        # One session's rewritten tables, measured from a fresh copy.
        compact = tmp_path / "compact.metacool"
        staging.copy_fresh(tiny_metacooler, compact)
        _metators_session(compact, tmp_path)
        one_session = compact.stat().st_size - _fresh_size(compact, tmp_path)
        assert dead[0] < 1_000_000  # the 4 MB left before are gone
        assert max(dead) <= one_session + 8_192, (dead, one_session)
        assert max(dead[1:]) - min(dead[1:]) <= 8_192, dead  # no growth

    def test_the_installed_map_is_the_staged_map(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        with staged_edit(tiny_metacooler) as staged:
            _writes(staged, tmp_path)
            with h5py.File(staged, "r") as f:
                staged_objects = _walk(f)
        with h5py.File(tiny_metacooler, "r") as f:
            assert _walk(f) == staged_objects

    def test_a_single_writer(self, tiny_metacooler: Path, tmp_path: Path):
        _dead_space(tiny_metacooler)
        sizes = []
        for _ in range(3):
            _store(tiny_metacooler)
            sizes.append(tiny_metacooler.stat().st_size)
        assert sizes[0] < 1_000_000
        # Once the column exists, each write replaces it: the same size.
        assert len(set(sizes[1:])) == 1, sizes

    def test_an_output_session(self, tiny_metacooler: Path, tmp_path: Path):
        _dead_space(tiny_metacooler)
        out = tmp_path / "out.metacool"
        with staged_edit(tiny_metacooler, output=out) as staged:
            _writes(staged, tmp_path)
        assert out.stat().st_size < 1_000_000

    def test_a_relabelled_copy_carries_no_space_of_its_source(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        _dead_space(tiny_metacooler)
        out = tmp_path / "relabelled.metacool"
        reassign(tiny_metacooler, {"contig_D": "MAG_2"}, output=out, layout="keep")
        annotated = tmp_path / "annotated.metacool"
        table = tmp_path / "depth.tsv"
        table.write_text("contig_name\tdepth\ncontig_A\t3.5\n")
        metacooler.annotate(tiny_metacooler, contigs_table=table, output=annotated)
        # metaTOR's pattern: a new map written from another.
        new = tmp_path / "new.metacool"
        with staged_inplace(new, copy=False) as staged:
            reassign(tiny_metacooler, {"contig_D": "MAG_2"}, output=staged)
        for written in (out, annotated, new):
            assert written.stat().st_size < 1_000_000


def _recorded(path: str | Path) -> list[str]:
    """Every attribute and every text value of the file *path*, as text."""
    found: list[str] = []

    def visit(name: str, node: h5py.HLObject) -> None:
        found.extend(f"{name}@{key}={node.attrs[key]!r}" for key in node.attrs)
        if isinstance(node, h5py.Dataset) and node.dtype.kind in "SOU":
            found.extend(repr(value) for value in np.asarray(node[()]).ravel())

    with h5py.File(path, "r") as f:
        found.extend(f"/@{key}={f.attrs[key]!r}" for key in f.attrs)
        f.visititems(visit)
    return found


class TestRecordedPaths:
    """A writer in a session names the session's map, never its staged map."""

    def test_every_writer_in_place(self, tiny_metacooler: Path, tmp_path: Path):
        from metacooler.ops.binnify import binnify

        coarse = tmp_path / "coarse.metacool"
        with staged_edit(tiny_metacooler) as staged:
            for write in WRITERS.values():
                write(staged, tmp_path)
            binnify(staged, coarse, 200)
            workdir = Path(staged).parent.name
        for path in (tiny_metacooler, coarse):
            assert not [text for text in _recorded(path) if workdir in text]
        metadata = MetaCooler(tiny_metacooler).metadata
        sources = [entry["source"] for entry in metadata["reassign"]]
        assert sources == [str(tiny_metacooler)] * 2
        assert MetaCooler(coarse).metadata["binnify"][-1]["source"] == str(
            tiny_metacooler
        )

    def test_an_entry_and_an_output(self, tiny_mmetacool: Path, tmp_path: Path):
        entry = f"{tiny_mmetacool}::entries/entry1"
        with staged_edit(entry) as staged:
            reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        out = tmp_path / "entry2.metacool"
        other = f"{tiny_mmetacool}::entries/entry2"
        with staged_edit(other, output=out) as staged:
            reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        assert MetaCooler(entry).metadata["reassign"][-1]["source"] == entry
        assert MetaCooler(out).metadata["reassign"][-1]["source"] == other

    def test_as_the_session_was_given(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.chdir(tiny_metacooler.parent)
        with staged_edit(tiny_metacooler.name) as staged:
            assert staging.recorded_uri(staged) == tiny_metacooler.name
            assert staging.recorded_uri(f"{staged}::/") == tiny_metacooler.name
            reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        history = MetaCooler(tiny_metacooler).metadata["reassign"]
        assert history[-1]["source"] == tiny_metacooler.name
        assert staging.recorded_uri(tiny_metacooler) == str(tiny_metacooler)


# ---------------------------------------------------------------------------
# The lock of a path itself; an output that is the map; what counts as written
# ---------------------------------------------------------------------------
class TestTheLockOfAPathItself:
    """``lock(path, follow_symlinks=False)``: the lock of an output that is a
    symbolic link, which is replaced, not followed."""

    @pytest.mark.usefixtures("holder")
    def test_a_link_to_a_map_another_process_holds(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        before = tiny_metacooler.read_bytes()
        out = tmp_path / "out" / "map.metacool"
        out.parent.mkdir()
        out.symlink_to(tiny_metacooler)
        with pytest.raises(MapLockedError, match="test holder"), lock(out):
            pass  # the target's lock
        with lock(out, follow_symlinks=False) as lock_file:
            assert lock_file == lock_path(out, follow_symlinks=False)
            assert lock_file == out.with_name("map.metacool.lock")
            with staged_edit(tiny_metacooler, output=out, force=True) as staged:
                _store(staged)
            assert lock_file.exists()  # the outer block's, still
        assert not out.is_symlink()
        assert "w" in MetaCooler(out).bins()[:].columns
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(out) == []

    def test_a_dangling_link(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "out" / "map.metacool"
        out.parent.mkdir()
        out.symlink_to(tmp_path / "gone" / "map.metacool")
        with pytest.raises(MapLockError), lock(out):
            pass  # beside a target whose directory is gone
        with (
            lock(out, follow_symlinks=False),
            staged_edit(tiny_metacooler, output=out, force=True) as staged,
        ):
            _store(staged)
        assert not out.is_symlink()
        validate(out, full=True)
        assert _left_beside(out) == []

    def test_re_entered_by_the_session(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "edited.metacool"
        with lock(out, follow_symlinks=False) as held:
            run = json.loads(held.read_text())["uuid"]
            with staged_edit(tiny_metacooler, output=out) as staged:
                assert Path(staged).parent.name == f".edited.metacool.{run}.1.d"
                assert holds_lock(out, follow_symlinks=False)
                _store(staged)
            assert held.exists()
        assert not held.exists()
        assert "w" in MetaCooler(out).bins()[:].columns


class TestAnOutputThatIsTheMap:
    def test_a_link_to_the_map_is_replaced(self, tiny_metacooler: Path, tmp_path: Path):
        before = tiny_metacooler.read_bytes()
        out = tmp_path / "link.metacool"
        out.symlink_to(tiny_metacooler)
        with staged_edit(tiny_metacooler, output=out, force=True) as staged:
            _store(staged)
        assert not out.is_symlink()
        assert "w" in MetaCooler(out).bins()[:].columns
        assert tiny_metacooler.read_bytes() == before

    def test_the_maps_own_file_is_refused(self, tiny_metacooler: Path, tmp_path: Path):
        hard = tmp_path / "hard.metacool"
        os.link(tiny_metacooler, hard)
        directory = tmp_path / "directory-link"
        directory.symlink_to(tiny_metacooler.parent)
        for output in (
            tiny_metacooler,
            f"{tiny_metacooler.parent}/./{tiny_metacooler.name}",
            directory / tiny_metacooler.name,
            hard,
        ):
            with (
                pytest.raises(ValueError, match="is the file of the map edited"),
                staged_edit(tiny_metacooler, output=output, force=True),
            ):
                pass
        assert _left_beside(tiny_metacooler) == []


def _assignment(tmp_path: Path) -> Path:
    table = tmp_path / "bins.tsv"
    table.write_text(
        "contig_name\tgenome_name\ncontig_A\tMAG_1\ncontig_B\tMAG_1\n"
        "contig_C\tMAG_2\ncontig_D\tMAG_2\n"
    )
    return table


class TestWhatIsWritten:
    """A session in place writes its map back when the staged file changed --
    however it changed -- and only then."""

    def test_a_reassignment_into_the_staged_map(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        with staged_edit(tiny_metacooler) as staged:
            reassign(tiny_metacooler, _assignment(tmp_path), output=staged, force=True)
        mc = MetaCooler(tiny_metacooler)
        genome = dict(zip(mc.contigs["name"], mc.contigs["genome_id"], strict=True))
        assert mc.genomes["name"][genome["contig_D"]] == "MAG_2"
        assert mc.metadata["reassign"][-1]["source"] == str(tiny_metacooler)

    def test_a_raw_write_as_soon_as_the_copy_is_made(self, tiny_metacooler: Path):
        with staged_edit(tiny_metacooler) as staged, h5py.File(staged, "r+") as f:
            f.attrs["note"] = "written"
        with h5py.File(tiny_metacooler, "r") as f:
            assert f.attrs["note"] == "written"

    def test_a_writer_with_nothing_to_write(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        _annotate(tiny_metacooler, tmp_path)
        inode = tiny_metacooler.stat().st_ino
        with staged_edit(tiny_metacooler) as staged:
            _annotate(staged, tmp_path)  # the values it holds already
        assert tiny_metacooler.stat().st_ino == inode

    def test_a_staged_map_deleted_writes_nothing(self, tiny_metacooler: Path):
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises((OSError, MetaCoolerFormatError)),
            staged_edit(tiny_metacooler) as staged,
        ):
            os.unlink(staged)
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.parametrize("where", ["in place", "to an output", "an entry out"])
    def test_metators_session(
        self, tiny_metacooler: Path, tiny_mmetacool: Path, tmp_path: Path, where
    ):
        src = (
            f"{tiny_mmetacool}::entries/entry1"
            if where == "an entry out"
            else str(tiny_metacooler)
        )
        output = None if where == "in place" else tmp_path / "out.metacool"
        mc = MetaCooler(src)
        before = mc.store.read_bytes()
        with staged_edit(src, output=output, planned=mc) as staged:
            reassign(src, _assignment(tmp_path), output=staged, force=True)
            _annotate(staged, tmp_path)
            _store(staged)
            write_clusters(staged, *_rounds(staged))
        result = output or src
        validate(result, full=True)
        written = MetaCooler(result)
        assert written.metadata["reassign"][-1]["source"] == src
        assert "coverage" in written.contigs.columns
        assert "w" in written.bins()[:].columns
        assert [r["name"] for r in written.cluster_rounds] == ["r1"]
        if output is not None:  # the map was only read
            assert mc.store.read_bytes() == before


# ---------------------------------------------------------------------------
# The fourth review: the staged map's lock held, outputs locked, what a map
# read for a weight holds
# ---------------------------------------------------------------------------
def _other_map(tiny_metacooler: Path, tmp_path: Path, name: str = "other") -> Path:
    """A copy of the map, in a directory of its own."""
    other = tmp_path / name / tiny_metacooler.name
    other.parent.mkdir()
    shutil.copy(tiny_metacooler, other)
    return other


class TestTheStagedMapsLockIsHeld:
    """The staged map's lock is held by its session's thread, as any lock of
    this process: never taken for one a process that is gone left."""

    def test_an_output_session_onto_it_re_enters_it(
        self, tiny_metacooler: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ):
        from _staging_workers import try_store

        other = _other_map(tiny_metacooler, tmp_path)
        with staged_edit(tiny_metacooler) as staged:
            with (
                caplog.at_level("WARNING", logger="metacooler"),
                staged_edit(other, output=staged, force=True) as inner,
            ):
                edit_metadata(inner, lambda meta: meta.update({"by": "inner"}))
            assert "Took over" not in caplog.text
            assert lock_path(staged).exists()
            with SPAWN.Pool(1) as pool:
                said = pool.apply(try_store, (staged,))
            assert said.startswith("MapLockedError"), said
        written = MetaCooler(tiny_metacooler)
        assert written.metadata["by"] == "inner"
        assert "w" not in written.bins()[:].columns
        assert _left_beside(tiny_metacooler) == []

    def test_another_thread_is_refused_it(self, tiny_metacooler: Path, tmp_path: Path):
        other = _other_map(tiny_metacooler, tmp_path)

        def onto(staged: str) -> None:
            with staged_edit(other, output=staged, force=True) as inner:
                edit_metadata(inner, lambda meta: meta.update({"by": "other thread"}))

        with staged_edit(tiny_metacooler) as staged:
            inode = Path(staged).stat().st_ino
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(lambda: onto(staged))
            assert Path(staged).stat().st_ino == inode
            assert lock_path(staged).exists()
            _store(staged)
        written = MetaCooler(tiny_metacooler)
        assert "by" not in written.metadata
        assert "w" in written.bins()[:].columns


class TestTheOutputsOfASession:
    """An output written where an edit session of this thread would replace
    it, or over the map it copies, is refused; its staged map is written as
    any output, by its thread alone."""

    def test_the_output_of_an_output_session(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "out.metacool"
        other = _other_map(tiny_metacooler, tmp_path)
        table = _assignment(tmp_path)
        depth = tmp_path / "depth.tsv"
        depth.write_text("contig_name\tdepth\ncontig_A\t3.5\n")
        mc = MetaCooler(tiny_metacooler)
        with staged_edit(tiny_metacooler, output=out, force=True, planned=mc) as staged:
            with pytest.raises(ValueError, match="is the output of an edit session"):
                reassign(tiny_metacooler, table, output=out, force=True)
            with pytest.raises(ValueError, match="is the output of an edit session"):
                metacooler.annotate(
                    tiny_metacooler, out, contigs_table=depth, force=True
                )
            with pytest.raises(ValueError, match="is being copied to"):
                reassign(other, table, output=tiny_metacooler, force=True)
            assert not out.exists()
            # metaTOR's step: the map reassigned into the staged map.
            reassign(tiny_metacooler, table, output=staged, force=True)
            _store(staged)
        written = MetaCooler(out)
        assert written.metadata["reassign"][-1]["source"] == str(tiny_metacooler)
        assert "w" in written.bins()[:].columns
        assert "reassign" not in MetaCooler(tiny_metacooler).metadata
        assert _left_beside(out) == []

    def test_a_map_edited_in_place(self, tiny_metacooler: Path, tmp_path: Path):
        other = _other_map(tiny_metacooler, tmp_path)
        table = _assignment(tmp_path)
        with staged_edit(tiny_metacooler) as staged:
            _store(staged, "outer")
            with (
                pytest.raises(ValueError, match="being edited in an edit session"),
                staged_edit(other, output=tiny_metacooler, force=True),
            ):
                pass
            with pytest.raises(ValueError, match="being edited in an edit session"):
                reassign(other, table, output=tiny_metacooler, force=True)
        assert "outer" in MetaCooler(tiny_metacooler).bins()[:].columns
        assert _left_beside(tiny_metacooler) == []

    def test_its_staged_map_from_another_thread(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        table = _assignment(tmp_path)
        with staged_edit(tiny_metacooler) as staged:
            _annotate(staged, tmp_path)
            with pytest.raises(MapLockedError, match="another thread"):
                _in_thread(
                    lambda: reassign(tiny_metacooler, table, output=staged, force=True)
                )
        written = MetaCooler(tiny_metacooler)
        assert "coverage" in written.contigs.columns
        assert "reassign" not in written.metadata

    def test_its_staged_map_from_another_process(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        from _staging_workers import try_reassign

        table = _assignment(tmp_path)
        with staged_edit(tiny_metacooler) as staged, SPAWN.Pool(1) as pool:
            _annotate(staged, tmp_path)
            said = pool.apply(try_reassign, (str(tiny_metacooler), str(table), staged))
        assert said.startswith("MapLockedError"), said
        written = MetaCooler(tiny_metacooler)
        assert "coverage" in written.contigs.columns
        assert "reassign" not in written.metadata


class _CliFailed(Exception):
    """A metacooler command that failed: its output, flattened."""


def _cli(*args: str | Path) -> None:
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(cli, [str(arg) for arg in args])
    if result.exit_code:
        raise _CliFailed(flat(result.output))


def _from_pairs(out: Path, request: pytest.FixtureRequest) -> None:
    from metacooler.io.create import from_pairs

    from_pairs(
        uri=out,
        bins=request.getfixturevalue("tiny_bins_bed"),
        pairs=request.getfixturevalue("tiny_pairs"),
        contig_genome_map=request.getfixturevalue("tiny_contigs_table"),
        contig_table_format="contigs",
        fasta=request.getfixturevalue("tiny_fasta"),
        chunksize=1000,
    )


def _collect(out: Path) -> None:
    from metacooler.io.collect import collect

    cools = out.parent.parent / "cools"
    cools.mkdir()
    collect(
        [
            make_cool(cools / "E_coli.cool", ECOLI, seed=1),
            make_cool(cools / "B_subtilis.cool", BSUB, seed=2),
        ],
        out,
    )


def _merge(src: Path, out: Path) -> None:
    from metacooler.io.merge import merge

    replicate = out.parent.parent / "replicate.metacool"
    shutil.copy(src, replicate)
    merge([src, replicate], out)


def _an_output_session(src: Path, out: Path) -> None:
    with staged_edit(src, output=out):
        pass


def _depth(out: Path) -> Path:
    table = out.parent.parent / "depth.tsv"
    table.write_text("contig_name\tdepth\ncontig_A\t3.5\n")
    return table


#: Every writer of a new map at a path: the output's name, and a call writing
#: it from the map *src*.
OUTPUT_WRITERS = {
    "create pairs": ("out.metacool", lambda src, out, rq: _from_pairs(out, rq)),
    "create cool": (
        "out.metacool",
        lambda src, out, rq: metacooler.from_cool(rq.getfixturevalue("tiny_cool"), out),
    ),
    "collect": ("out.metacool", lambda src, out, rq: _collect(out)),
    "reassign, relabelled": (
        "out.metacool",
        lambda src, out, rq: reassign(
            src, {"contig_D": "MAG_2"}, output=out, layout="keep"
        ),
    ),
    "reassign, re-grouped": (
        "out.metacool",
        lambda src, out, rq: reassign(src, {"contig_A": "MAG_2"}, output=out),
    ),
    "annotate": (
        "out.metacool",
        lambda src, out, rq: metacooler.annotate(src, out, contigs_table=_depth(out)),
    ),
    "binnify": ("out.metacool", lambda src, out, rq: metacooler.binnify(src, out, 200)),
    "zoomify": (
        "out.mmetacool",
        lambda src, out, rq: metacooler.zoomify(src, out, levels=2),
    ),
    "aggregate": (
        "out.cool",
        lambda src, out, rq: metacooler.aggregate(src, out, level="genome"),
    ),
    "extract": (
        "out.cool",
        lambda src, out, rq: metacooler.extract(src, genomes=["MAG_1"], output=out),
    ),
    "flatten": ("out.cool", lambda src, out, rq: metacooler.flatten(src, out)),
    "merge": ("out.metacool", lambda src, out, rq: _merge(src, out)),
    "a new collection": (
        "out.mmetacool",
        lambda src, out, rq: mmetacool_module.create_mmetacool(out, {"a": src}),
    ),
    "stack": (
        "out.mmetacool",
        lambda src, out, rq: _cli("stack", src, src, "--names", "a,b", "-o", out),
    ),
    "cp, a new collection": (
        "out.mmetacool",
        lambda src, out, rq: _cli(
            "cp",
            f"{rq.getfixturevalue('tiny_mmetacool')}::entries/entry1",
            f"{out}::entries/entry1",
        ),
    ),
    "an output session": (
        "out.metacool",
        lambda src, out, rq: _an_output_session(src, out),
    ),
}


@pytest.mark.parametrize("writer", list(OUTPUT_WRITERS))
def test_an_output_another_process_writes_is_refused(
    writer: str, tiny_metacooler: Path, tmp_path: Path, request: pytest.FixtureRequest
):
    """Every writer of a new map holds the lock of its output -- the path
    itself -- while it writes it: a writer of the same path in another process
    refuses it, and nothing is written there."""
    from _staging_workers import hold

    name, write = OUTPUT_WRITERS[writer]
    out = tmp_path / "out" / name
    out.parent.mkdir()
    ready, release = SPAWN.Event(), SPAWN.Event()
    process = SPAWN.Process(
        target=hold, args=(str(out), ready, release), kwargs={"follow_symlinks": False}
    )
    process.start()
    try:
        assert ready.wait(60), "the holder did not take the lock"
        with pytest.raises((MapLockedError, _CliFailed), match="test holder"):
            write(tiny_metacooler, out, request)
        assert sorted(p.name for p in out.parent.iterdir()) == [f"{name}.lock"]
    finally:
        release.set()
        process.join(60)


class TestAMapCopiedOut:
    """An output session does not lock the map it copies: what planned was
    read from must be what is copied, and what was read from the map as it
    was copied holds, whatever is written to the map since."""

    def test_a_map_written_before_its_copy_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        out = tmp_path / "out.metacool"
        planned = MetaCooler(tiny_metacooler)
        real = staging._hold_staged

        def hold_while_another_program_writes(session, command):
            # Once planned is checked, before the map is copied.
            with h5py.File(tiny_metacooler, "r+") as f:
                f.attrs["by"] = "another program"
            real(session, command)

        monkeypatch.setattr(staging, "_hold_staged", hold_while_another_program_writes)
        with (
            pytest.raises(MapChangedError),
            staged_edit(tiny_metacooler, output=out, planned=planned),
        ):
            pass
        assert not out.exists()
        assert _left_beside(out) == []

    def test_a_map_written_after_its_copy(self, tiny_metacooler: Path, tmp_path: Path):
        out = tmp_path / "out.metacool"
        source = MetaCooler(tiny_metacooler)
        weights = _bin_keys(tiny_metacooler)
        with staged_edit(tiny_metacooler, output=out, planned=source) as staged:
            with h5py.File(tiny_metacooler, "r+") as f:  # not locked: anyone may
                f.attrs["by"] = "another program"
            with staged_edit(staged, planned=source) as again:
                assert again == staged
            _put_weights(staged, weights, source)
        np.testing.assert_array_equal(MetaCooler(out).bins()[:]["w"], weights)
        with h5py.File(out, "r") as f:
            assert "by" not in f.attrs

    def test_a_map_read_after_a_regrouping_has_its_old_bins(
        self, tiny_metacooler: Path
    ):
        """In place, the map is not written before the session ends: read
        after a re-grouping of the staged map, it still has the bins the
        session started from."""
        with staged_edit(tiny_metacooler) as staged:
            reassign(staged, {"contig_A": "MAG_2"}, inplace=True)  # re-groups
            after = MetaCooler(tiny_metacooler)
            with pytest.raises(MapChangedError, match="bins of the map written"):
                _put_weights(staged, _bin_keys(tiny_metacooler), after)
        validate(tiny_metacooler, full=True)


def _per_genome(mc: MetaCooler) -> np.ndarray:
    """1 / the number of bins of each bin's genome: a weight that depends on
    the assignment."""
    genome_of = dict(zip(mc.contigs["name"], mc.contigs["genome_id"], strict=True))
    genome = mc.bins()[:]["chrom"].astype(str).map(genome_of).to_numpy()
    return 1.0 / np.bincount(genome)[genome]


def _by_contig(mc: MetaCooler, column: str) -> np.ndarray:
    """The value of a ``/contigs`` column at each bin."""
    value = dict(zip(mc.contigs["name"], mc.contigs[column], strict=True))
    return mc.bins()[:]["chrom"].astype(str).map(value).to_numpy(dtype=float)


def _store_from(
    path: str | Path, name: str, weights, depends_on: list[str], source: MetaCooler
) -> None:
    store_weights(
        path, name, weights, method="test", depends_on=depends_on, source=source
    )


class TestWhatAWeightDependsOn:
    """In a session, a weight computed from a map that holds the map as it was
    before a write that changed what the weight depends on is refused; one
    whose dependencies no write changed, or read from the staged map, is
    written."""

    def test_a_relabelling(self, tiny_metacooler: Path):
        with lock(tiny_metacooler):  # metaTOR holds it for its run
            before = MetaCooler(tiny_metacooler)
            weights = _per_genome(before)
            with staged_edit(tiny_metacooler) as staged:
                reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
                with pytest.raises(MapChangedError, match="changed genome_id"):
                    _store_from(staged, "pg", weights, ["genome_id"], before)
                # Read after the relabelling, the map as it was is no newer.
                after = MetaCooler(tiny_metacooler)
                with pytest.raises(MapChangedError, match="changed genome_id"):
                    _store_from(staged, "pg", _per_genome(after), ["genome_id"], after)
                _store_from(staged, "px", weights, ["pixels"], before)
                current = MetaCooler(staged)
                _store_from(staged, "pg", _per_genome(current), ["genome_id"], current)
        written = MetaCooler(tiny_metacooler)
        bins = written.bins()[:]
        np.testing.assert_array_equal(bins["pg"], _per_genome(written))
        assert "px" in bins.columns

    def test_an_annotation(self, tiny_metacooler: Path, tmp_path: Path):
        _annotate(tiny_metacooler, tmp_path, 2.5)
        before = MetaCooler(tiny_metacooler)
        weights = 1.0 / _by_contig(before, "coverage")
        with staged_edit(tiny_metacooler) as staged:
            _annotate(staged, tmp_path, 3.5)  # replaces the column
            with pytest.raises(MapChangedError, match="changed contigs/coverage"):
                _store_from(
                    staged, "c", weights, ["pixels", "contigs/coverage"], before
                )
            write_clusters(staged, *_rounds(staged))  # changes no dependency
            _store_from(staged, "p", weights, ["pixels"], before)
        assert "p" in MetaCooler(tiny_metacooler).bins()[:].columns

    def test_another_file(self, tiny_metacooler: Path, tmp_path: Path):
        """A map of another file is current for what no write changed after
        it was opened."""
        copy = tmp_path / "copy.metacool"
        shutil.copy(tiny_metacooler, copy)
        early = MetaCooler(copy)
        with staged_edit(tiny_metacooler) as staged:
            reassign(staged, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
            with pytest.raises(MapChangedError, match="changed genome_id"):
                _store_from(staged, "g", _per_genome(early), ["genome_id"], early)
            late = MetaCooler(copy)  # a copy made since, say
            _store_from(staged, "g", _per_genome(late), ["genome_id"], late)

    def test_metators_pattern(self, tiny_metacooler: Path, tmp_path: Path):
        """metaTOR's step: the map reassigned into the staged map, annotated,
        then weighted from the staged map as it then is."""
        with staged_edit(
            tiny_metacooler, planned=MetaCooler(tiny_metacooler)
        ) as staged:
            reassign(tiny_metacooler, _assignment(tmp_path), output=staged, force=True)
            _annotate(staged, tmp_path, 3.5)
            mc = MetaCooler(staged)
            store_weights(
                staged,
                "hit",
                1.0 / _by_contig(mc, "coverage"),
                method="test",
                depends_on=["pixels", "contigs/coverage", "genome_id"],
                source=mc,
            )
            write_clusters(staged, *_rounds(staged))
        written = MetaCooler(tiny_metacooler)
        np.testing.assert_allclose(
            written.bins()[:]["hit"], 1.0 / _by_contig(written, "coverage")
        )
        assert [r["name"] for r in written.cluster_rounds] == ["r1"]

    def test_a_raw_write_changes_anything(self, tiny_metacooler: Path):
        before = MetaCooler(tiny_metacooler)
        with staged_edit(tiny_metacooler) as staged:
            with h5py.File(staged, "r+") as f:
                f.attrs["note"] = "raw"
            with pytest.raises(MapChangedError, match="changed pixels"):
                _put_weights(staged, _bin_keys(tiny_metacooler), before)


class TestInstalls:
    def test_an_output_created_while_it_is_installed(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Without force, an output another program creates at any time before
        the rename is kept: checked right before it."""
        out = tmp_path / "out.metacool"
        real = staging._fsync

        def fsync_while_another_program_writes(path: Path) -> None:
            real(path)
            if not path.is_dir() and not out.exists():
                out.write_text("another program's result\n")

        monkeypatch.setattr(staging, "_fsync", fsync_while_another_program_writes)
        with (
            pytest.raises(FileExistsError, match="created meanwhile"),
            staged_edit(tiny_metacooler, output=out),
        ):
            pass
        assert out.read_text() == "another program's result\n"
        assert _left_beside(out) == []

    def test_a_write_in_the_tick_of_the_copy(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Where file times are coarse -- here, ten seconds, and no size or
        change time -- a write in the same tick as the copy is still seen: the
        staged map's modification time was set back when it was copied."""

        def coarse(path: Path):
            try:
                info = os.stat(path)
            except FileNotFoundError:
                return None
            return info.st_ino, info.st_mtime_ns // 10_000_000_000

        monkeypatch.setattr(staging, "_written_state", coarse)
        with staged_edit(tiny_metacooler) as staged, h5py.File(staged, "r+") as f:
            f.attrs["note"] = "written"
        with h5py.File(tiny_metacooler, "r") as f:
            assert f.attrs["note"] == "written"


def test_names_are_found_with_one_check(
    tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
):
    """Contig names resolved through the cooler check the map once, not once
    per name."""
    import metacooler.api as api_module

    mc = MetaCooler(tiny_metacooler)
    names = [str(name) for name in mc.contigs["name"]] * 50
    assert mc._clr is not None  # built
    checks = _counted(monkeypatch, [api_module], "file_state")
    ids = mc._contig_ids_of(names)
    assert len(checks) == 1
    assert list(ids[:4]) == [0, 1, 2, 3]


@pytest.mark.parametrize("writer", list(OUTPUT_WRITERS))
def test_an_output_in_a_missing_directory_is_named(
    writer: str, tiny_metacooler: Path, tmp_path: Path, request: pytest.FixtureRequest
):
    name, write = OUTPUT_WRITERS[writer]
    out = tmp_path / "out" / "missing" / name
    out.parent.parent.mkdir()
    with pytest.raises((FileNotFoundError, _CliFailed), match="does not exist"):
        write(tiny_metacooler, out, request)
    assert not out.parent.exists()
