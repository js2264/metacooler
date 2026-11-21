"""Locked, staged in-place writes, and reads that span a replacement.

metaTOR 2.0 updates a user's only copy of a map step by step, so every writer
it calls edits a map in place through a staged copy, under the map's lock
(``metacooler.io.staging``).  These tests are metaTOR #64's guarantees of its
own staged writes, ported: the lock (holders, stale locks, the races of taking
one), the staged copy (failures, interrupts, the mode, symbolic links, a map
that does not grow), and a read that spans a replacement.
"""

from __future__ import annotations

import json
import os
import socket
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import metacooler
from metacooler.api import MetaCooler
from metacooler.io import staging
from metacooler.io.annotate import annotate
from metacooler.io.clusters import drop_clusters, write_clusters
from metacooler.io.reassign import reassign
from metacooler.io.staging import (
    MapChangedError,
    MapLockedError,
    MapLockError,
    holds_lock,
    lock,
    lock_path,
    map_output,
    staged_inplace,
)
from metacooler.io.weights import store_weights
from metacooler.schema import MetaCoolerFormatError, validate

#: A pid no process of this host has (above any pid_max).
DEAD_PID = 999_999_999


def _write_lock(
    path: Path, pid: int, *, host: str | None = None, run: str = "deadrun"
) -> dict:
    """A lock file of *path*, as a holder of *pid* on *host* wrote it."""
    holder = {
        "pid": pid,
        "host": host or socket.gethostname(),
        "start": "2026-11-20T10:00:00+0000",
        "command": "metator pipeline map.metacool",
        "uuid": run,
    }
    lock_path(path).write_text(json.dumps(holder))
    return holder


def _left_beside(path: Path) -> list[str]:
    """What a writer left beside the map: staged files, partial files, locks."""
    return sorted(
        p.name
        for p in path.parent.iterdir()
        if p.name.startswith(".") or p.name.endswith(".lock")
    )


def _n_bins(path: str | Path) -> int:
    return int(MetaCooler(path).info["nbins"])


def _store(path: str | Path, name: str = "w", value: float = 1.0) -> None:
    store_weights(
        path,
        name,
        np.full(_n_bins(path), value),
        method="test",
        depends_on="pixels",
        force=True,
    )


def _rounds(path: str | Path) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    contig = str(MetaCooler(path).contigs["name"].iloc[0])
    clusters = pd.DataFrame({"round": ["r1"], "name": ["c1"]})
    members = pd.DataFrame({"round": ["r1"], "cluster": ["c1"], "contig": [contig]})
    return ["r1"], clusters, members


def _annotate(path: str | Path, tmp_path: Path, value: float = 2.5) -> None:
    table = tmp_path / f"coverage_{value}.tsv"
    entry = f"{path}::entries/entry1" if str(path).endswith(".mmetacool") else path
    names = MetaCooler(entry).contigs["name"].astype(str)
    table.write_text(
        "contig_name\tcoverage\n" + "".join(f"{n}\t{value}\n" for n in names)
    )
    annotate(path, contigs_table=table, inplace=True, replace=True)


#: Each in-place writer, as a call on a map path.
WRITERS = {
    "annotate": lambda path, tmp: _annotate(path, tmp),
    "store_weights": lambda path, tmp: _store(path),
    "write_clusters": lambda path, tmp: write_clusters(path, *_rounds(path)),
    "drop_clusters": lambda path, tmp: (
        write_clusters(path, *_rounds(path)),
        drop_clusters(path),
    ),
    "reassign, relabel": lambda path, tmp: reassign(
        path, {"contig_D": "MAG_2"}, inplace=True, layout="keep"
    ),
    "reassign, regroup": lambda path, tmp: reassign(
        path, {"contig_A": "MAG_2"}, inplace=True
    ),
}


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------
class TestLock:
    def test_held_for_the_block_and_named(self, tiny_metacooler: Path):
        with lock(tiny_metacooler) as lock_file:
            assert lock_file == tiny_metacooler.with_name("tiny.metacool.lock")
            holder = json.loads(lock_file.read_text())
            assert set(holder) == {"pid", "host", "start", "command", "uuid"}
            assert holder["pid"] == os.getpid()
            assert holder["host"] == socket.gethostname()
            assert holds_lock(tiny_metacooler)
        assert not lock_file.exists()
        assert not holds_lock(tiny_metacooler)

    def test_a_uri_locks_its_file(self, tiny_metacooler: Path, tiny_mmetacool: Path):
        with lock(f"{tiny_metacooler}::/") as lock_file:
            assert lock_file == lock_path(tiny_metacooler)
        with lock(f"{tiny_mmetacool}::entries/entry1") as lock_file:
            assert lock_file == tiny_mmetacool.with_name("tiny.mmetacool.lock")
            assert holds_lock(f"{tiny_mmetacool}::entries/entry2")

    def test_next_to_a_symbolic_link_not_its_target(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """The link is what an in-place write replaces, so it is what is
        locked."""
        link = tmp_path / "link.metacool"
        link.symlink_to(tiny_metacooler)
        with lock(link) as lock_file:
            assert lock_file == tmp_path / "link.metacool.lock"

    def test_re_entrant_for_its_holder(self, tiny_metacooler: Path):
        with lock(tiny_metacooler) as outer:
            mine = outer.read_text()
            with lock(tiny_metacooler, command="ignored") as inner:
                assert inner == outer and inner.read_text() == mine
            assert outer.exists()
        assert not outer.exists()

    def test_the_command_is_recorded(self, tiny_metacooler: Path):
        with lock(tiny_metacooler, command="metator partition map") as lock_file:
            assert json.loads(lock_file.read_text())["command"] == (
                "metator partition map"
            )

    def test_a_live_holder_is_refused_by_name(self, tiny_metacooler: Path):
        holder = _write_lock(tiny_metacooler, os.getpid(), run="otherrun")
        with pytest.raises(MapLockedError) as caught, lock(tiny_metacooler):
            pass
        message = str(caught.value)
        assert message.startswith(f"{tiny_metacooler} is being written by another")
        assert f"pid {os.getpid()} on {socket.gethostname()}" in message
        assert "since 2026-11-20T10:00:00+0000: metator pipeline map.metacool" in (
            message
        )
        assert f"delete {lock_path(tiny_metacooler)}" in message
        assert caught.value.holder == holder
        assert caught.value.lock_file == lock_path(tiny_metacooler)
        assert json.loads(lock_path(tiny_metacooler).read_text()) == holder

    def test_another_hosts_lock_is_never_broken(self, tiny_metacooler: Path):
        holder = _write_lock(tiny_metacooler, DEAD_PID, host="elsewhere")
        with (
            pytest.raises(MapLockedError, match="another host") as caught,
            lock(tiny_metacooler),
        ):
            pass
        assert f"delete {lock_path(tiny_metacooler)}" in str(caught.value)
        assert json.loads(lock_path(tiny_metacooler).read_text()) == holder

    @pytest.mark.parametrize("text", ["{not json", "[1]", ""])
    def test_a_lock_that_says_nothing_is_refused(
        self, tiny_metacooler: Path, text: str
    ):
        lock_path(tiny_metacooler).write_text(text)
        with (
            pytest.raises(MapLockedError, match="does not say which process"),
            lock(tiny_metacooler),
        ):
            pass
        assert lock_path(tiny_metacooler).read_text() == text

    def test_a_dead_holders_lock_is_taken_over(
        self, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
    ):
        """A lock of a process of this host that is gone (killed outright) is
        taken over, and the files it staged are removed: only its own."""
        _write_lock(tiny_metacooler, DEAD_PID, run="deadrun")
        mine = tiny_metacooler.parent
        staged = mine / ".tiny.metacool.deadrun.1.staged"
        partial = mine / ".tiny.metacool.deadrun.1.staged.__writing__"
        other = mine / ".tiny.metacool.otherrun.1.staged"
        for path in (staged, partial, other):
            path.write_text("x")
        with (
            caplog.at_level("WARNING", logger="metacooler"),
            lock(tiny_metacooler) as lock_file,
        ):
            assert json.loads(lock_file.read_text())["pid"] == os.getpid()
        assert f"Taking over {lock_path(tiny_metacooler)}" in caplog.text
        assert f"pid {DEAD_PID}" in caplog.text
        assert "Removed .tiny.metacool.deadrun.1.staged," in caplog.text
        assert not staged.exists() and not partial.exists() and other.exists()
        assert not lock_path(tiny_metacooler).exists()

    def test_only_its_holder_removes_it(self, tiny_metacooler: Path):
        """A lock another process took during this one (this one's deleted by
        hand) is that process's: it stays."""
        with lock(tiny_metacooler):
            other = _write_lock(tiny_metacooler, 1, host="elsewhere", run="other")
        assert json.loads(lock_path(tiny_metacooler).read_text()) == other

    def test_a_forked_child_holds_none_of_its_parents_locks(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        with lock(tiny_metacooler):
            monkeypatch.setattr(staging.os, "getpid", lambda: os.getppid())
            assert not holds_lock(tiny_metacooler)
            with pytest.raises(MapLockedError), lock(tiny_metacooler):
                pass
            monkeypatch.undo()
        assert not lock_path(tiny_metacooler).exists()

    def test_a_stale_lock_is_taken_once(self, tiny_metacooler: Path):
        """Two processes that both found a lock stale: the first moves it aside
        and takes the lock; the second moves the first's aside, sees that it is
        not the stale one, puts it back, and is refused."""
        stale = _write_lock(tiny_metacooler, DEAD_PID)
        host = socket.gethostname()
        first = {"pid": 1, "host": host, "start": "now", "uuid": "first"}
        second = {"pid": 2, "host": host, "start": "now", "uuid": "second"}
        lock_file = lock_path(tiny_metacooler)
        assert staging._take_stale_lock(lock_file, stale, first)
        assert not staging._take_stale_lock(lock_file, stale, second)
        assert json.loads(lock_file.read_text()) == first
        assert [p.name for p in lock_file.parent.iterdir() if ".lock." in p.name] == []
        # Gone in the meantime (taken and released by another process): taken.
        lock_file.unlink()
        assert staging._take_stale_lock(lock_file, stale, second)
        assert json.loads(lock_file.read_text()) == second

    def test_a_stale_lock_taken_by_a_third_process(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """The case _take_stale_lock leaves open: a second process moves aside
        the lock a first one took, and a third creates the lock meanwhile. The
        second cannot put the first's back: it leaves the third's, and is
        refused."""
        stale = _write_lock(tiny_metacooler, DEAD_PID)
        host = socket.gethostname()
        runs = {
            name: {"pid": i, "host": host, "start": "now", "uuid": name}
            for i, name in enumerate(["first", "second", "third"], start=1)
        }
        lock_file = lock_path(tiny_metacooler)
        assert staging._take_stale_lock(lock_file, stale, runs["first"])
        rename = os.rename

        def then_a_third(source, target):
            rename(source, target)
            assert staging._create_lock(lock_file, runs["third"])

        monkeypatch.setattr(staging.os, "rename", then_a_third)
        assert not staging._take_stale_lock(lock_file, stale, runs["second"])
        assert json.loads(lock_file.read_text()) == runs["third"]
        assert [p.name for p in lock_file.parent.iterdir() if ".lock." in p.name] == []

    def test_a_stale_lock_gone_meanwhile_is_taken(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A stale lock another process moved aside first: the lock is to
        take."""
        stale = _write_lock(tiny_metacooler, DEAD_PID)
        lock_file = lock_path(tiny_metacooler)
        lock_file.unlink()
        mine = {"pid": os.getpid(), "host": socket.gethostname(), "uuid": "mine"}
        assert staging._take_stale_lock(lock_file, stale, mine)
        assert json.loads(lock_file.read_text()) == mine

    def test_a_lock_that_cannot_be_created(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A map whose directory is not writable cannot be locked, nor written
        in place: refused before any work, saying why."""
        opened = os.open
        target = str(lock_path(tiny_metacooler))

        def read_only(path, *args, **kwargs):
            if str(path) == target:
                raise PermissionError(13, "Permission denied", path)
            return opened(path, *args, **kwargs)

        monkeypatch.setattr(staging.os, "open", read_only)
        with pytest.raises(MapLockError) as caught, lock(tiny_metacooler):
            pass
        assert not isinstance(caught.value, MapLockedError)
        assert str(caught.value).startswith(
            f"{target} cannot be created (Permission denied)"
        )

    def test_the_errors(self):
        assert issubclass(MapLockedError, MapLockError)
        assert issubclass(MapLockError, OSError)
        assert issubclass(MapChangedError, OSError)
        assert metacooler.lock is lock
        assert metacooler.MapChangedError is MapChangedError
        assert str(MapLockError("why", lock_file="x.lock")) == "why"
        assert str(MapChangedError("why")) == "why"

    def test_pid_checks(self, monkeypatch: pytest.MonkeyPatch):
        assert staging._pid_is_running(os.getpid())
        assert not staging._pid_is_running(DEAD_PID)
        assert not staging._pid_is_running(2**70)

        def another_users(pid, signal):
            raise PermissionError(1, "Operation not permitted")

        # A process of another user, which this one may not signal: alive.
        monkeypatch.setattr(staging.os, "kill", another_users)
        assert staging._pid_is_running(1)

    def test_what_cannot_be_removed_is_left(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        (tmp_path / ".map.run.1.staged").write_text("x")
        (tmp_path / ".map.run.1.fresh").write_text("x")
        removed = os.remove

        def refuse_one(path):
            if str(path).endswith(".staged"):
                raise PermissionError(13, "Permission denied", str(path))
            removed(path)

        monkeypatch.setattr(staging.os, "remove", refuse_one)
        assert staging._remove_named(tmp_path, ".map.run.") == [".map.run.1.fresh"]
        assert staging._remove_named(tmp_path / "missing", ".map.") == []


# ---------------------------------------------------------------------------
# Staged in-place writes
# ---------------------------------------------------------------------------
class TestStagedWrites:
    @pytest.mark.parametrize("writer", list(WRITERS))
    def test_every_writer_replaces_the_map_through_a_staged_copy(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str
    ):
        before = tiny_metacooler.stat().st_ino
        WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.stat().st_ino != before
        validate(tiny_metacooler, full=True)
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.parametrize("writer", list(WRITERS))
    @pytest.mark.parametrize("mode", [0o600, 0o640, 0o444])
    def test_the_map_keeps_its_mode(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str, mode: int
    ):
        tiny_metacooler.chmod(mode)
        WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.stat().st_mode & 0o777 == mode

    @pytest.mark.parametrize("writer", list(WRITERS))
    def test_a_symbolic_link_is_replaced_not_its_target(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str
    ):
        link = tmp_path / "link.metacool"
        link.symlink_to(tiny_metacooler)
        target = tiny_metacooler.read_bytes()
        WRITERS[writer](link, tmp_path)
        assert not link.is_symlink()
        assert tiny_metacooler.read_bytes() == target
        validate(link, full=True)

    @pytest.mark.parametrize(
        "stop", [RuntimeError("boom"), KeyboardInterrupt(), SystemExit(1)]
    )
    def test_a_failure_leaves_the_map_as_it_was(
        self, tiny_metacooler: Path, stop: BaseException
    ):
        """Whatever stops the write, the staged file and every partial file
        named after it (cooler's ``.__writing__``) are removed, and the map is
        as it was."""
        before = tiny_metacooler.read_bytes()
        with pytest.raises(type(stop)), staged_inplace(tiny_metacooler) as staged:
            Path(f"{staged}.__writing__").write_text("partial")
            with h5py.File(staged, "r+") as f:
                f.attrs["edited"] = True
            raise stop
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.parametrize("stop", [KeyboardInterrupt(), SystemExit(1)])
    def test_an_interrupt_as_the_map_is_replaced(
        self,
        tiny_metacooler: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        stop: BaseException,
    ):
        before = tiny_metacooler.read_bytes()

        def interrupted(*args):
            raise stop

        monkeypatch.setattr(staging.os, "replace", interrupted)
        with pytest.raises(type(stop)):
            _store(tiny_metacooler)
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_a_staged_map_that_fails_validation_never_replaces_the_map(
        self, tiny_metacooler: Path
    ):
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises(MetaCoolerFormatError),
            staged_inplace(tiny_metacooler) as staged,
            h5py.File(staged, "r+") as f,
        ):
            del f["genomes"]
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_from_nothing(self, tiny_metacooler: Path):
        """A writer that writes a whole new map starts from no file."""
        with staged_inplace(tiny_metacooler, copy=False) as staged:
            assert not staged.exists()
            staging.copy_fresh(tiny_metacooler, staged)
        validate(tiny_metacooler, full=True)

    def test_the_mode_of_a_new_file(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A path where no map is yet gets the mode of any new file, from the
        umask, not the private mode of a temporary file."""
        path = tmp_path / "new.metacool"
        previous = os.umask(0o027)
        try:
            with staged_inplace(path, copy=False) as staged:
                staging.copy_fresh(tiny_metacooler, staged)
        finally:
            os.umask(previous)
        assert path.stat().st_mode & 0o777 == 0o640

    @pytest.mark.parametrize("writer", ["annotate", "store_weights", "write_clusters"])
    def test_the_map_does_not_grow(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str
    ):
        """Updated again and again, the map keeps its size: an edit in place
        left the space of what it replaced in the file."""
        sizes = []
        for i in range(5):
            if writer == "annotate":
                _annotate(tiny_metacooler, tmp_path, value=float(i))
            elif writer == "store_weights":
                _store(tiny_metacooler, value=float(i + 1))
            else:
                WRITERS[writer](tiny_metacooler, tmp_path)
            sizes.append(tiny_metacooler.stat().st_size)
        assert max(sizes[1:]) - min(sizes[1:]) < 4_000, sizes

    def test_a_fresh_copy_is_the_map_as_stored(self, tiny_metacooler: Path, tmp_path):
        copy = tmp_path / "copy.metacool"
        staging.copy_fresh(tiny_metacooler, copy)

        def walk(path):
            seen: dict[str, object] = {}
            with h5py.File(path, "r") as f:

                def visit(name, obj):
                    link = f.get(name, getlink=True)
                    seen[name] = (
                        type(link).__name__,
                        dict(obj.attrs) if not isinstance(link, h5py.SoftLink) else {},
                        obj[()].tolist() if isinstance(obj, h5py.Dataset) else None,
                    )

                f.visititems(visit)
                seen["/"] = dict(f.attrs)
            return seen

        assert repr(walk(copy)) == repr(walk(tiny_metacooler))
        with h5py.File(copy, "r") as f:
            assert isinstance(f.get("chroms", getlink=True), h5py.SoftLink)

    def test_a_writer_takes_the_lock_unless_its_caller_holds_it(
        self, tiny_metacooler: Path
    ):
        with lock(tiny_metacooler) as lock_file:
            mine = lock_file.read_text()
            _store(tiny_metacooler)
            assert lock_file.read_text() == mine
        assert not lock_file.exists()
        assert "w" in MetaCooler(tiny_metacooler).bins()[:].columns

    @pytest.mark.parametrize("writer", list(WRITERS))
    def test_a_writer_is_refused_a_map_another_process_writes(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str
    ):
        _write_lock(tiny_metacooler, os.getpid(), run="otherrun")
        before = tiny_metacooler.read_bytes()
        with pytest.raises(MapLockedError):
            WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.read_bytes() == before

    def test_the_next_writer_removes_what_a_killed_one_left(
        self, tiny_metacooler: Path
    ):
        _write_lock(tiny_metacooler, DEAD_PID, run="deadrun")
        left = tiny_metacooler.with_name(".tiny.metacool.deadrun.1.staged")
        left.write_text("half a map")
        _store(tiny_metacooler)
        assert _left_beside(tiny_metacooler) == []

    def test_an_entry_of_a_collection(self, tiny_mmetacool: Path):
        """An entry is edited through a staged copy of its collection, the
        other entries copied as they are."""
        _store(f"{tiny_mmetacool}::entries/entry1")
        bins1 = MetaCooler(f"{tiny_mmetacool}::entries/entry1").bins()[:]
        bins2 = MetaCooler(f"{tiny_mmetacool}::entries/entry2").bins()[:]
        assert "w" in bins1.columns and "w" not in bins2.columns
        for entry in ("entry1", "entry2"):
            validate(f"{tiny_mmetacool}::entries/{entry}", full=True)
        assert _left_beside(tiny_mmetacool) == []

    def test_a_whole_collection_validates_every_entry(self, tiny_mmetacool: Path):
        before = tiny_mmetacool.read_bytes()
        with (
            pytest.raises(MetaCoolerFormatError),
            staged_inplace(tiny_mmetacool) as staged,
            h5py.File(staged, "r+") as f,
        ):
            del f["entries/entry2/genomes"]
        assert tiny_mmetacool.read_bytes() == before

    def test_commands_stored_with_the_weights(self, tiny_metacooler: Path):
        """The commands append to the provenance, stored as text or as
        bytes."""
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["metadata"] = np.bytes_(json.dumps({"commands": ["create"]}))
        store_weights(
            tiny_metacooler,
            "w",
            np.ones(_n_bins(tiny_metacooler)),
            method="test",
            depends_on="pixels",
            commands=["normalize"],
        )
        assert MetaCooler(tiny_metacooler).metadata["commands"] == [
            "create",
            "normalize",
        ]

    def test_annotating_a_whole_collection_in_place(
        self, tiny_mmetacool: Path, tmp_path: Path
    ):
        before = tiny_mmetacool.stat().st_ino
        _annotate(tiny_mmetacool, tmp_path)
        assert tiny_mmetacool.stat().st_ino != before
        for entry in ("entry1", "entry2"):
            contigs = MetaCooler(f"{tiny_mmetacool}::entries/{entry}").contigs
            assert (contigs["coverage"] == 2.5).all()

    def test_an_unchanged_annotation_writes_nothing(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        _annotate(tiny_metacooler, tmp_path)
        before = tiny_metacooler.stat().st_ino
        _annotate(tiny_metacooler, tmp_path)
        assert tiny_metacooler.stat().st_ino == before

    def test_dropping_no_clusters_writes_nothing(self, tiny_metacooler: Path):
        before = tiny_metacooler.stat().st_ino
        assert not drop_clusters(tiny_metacooler)
        assert tiny_metacooler.stat().st_ino == before

    def test_regrouping_an_entry_in_place_is_refused(self, tiny_mmetacool: Path):
        before = tiny_mmetacool.read_bytes()
        with pytest.raises(ValueError, match="cannot re-group an entry"):
            reassign(
                f"{tiny_mmetacool}::entries/entry1", {"contig_A": "MAG_2"}, inplace=True
            )
        assert tiny_mmetacool.read_bytes() == before

    def test_relabelling_an_entry_in_place(self, tiny_mmetacool: Path):
        entry = f"{tiny_mmetacool}::entries/entry1"
        reassign(entry, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        genome = MetaCooler(entry).contigs.set_index("name")["genome_id"]
        names = MetaCooler(entry).genomes["name"]
        assert names[genome["contig_D"]] == "MAG_2"
        validate(f"{tiny_mmetacool}::entries/entry2", full=True)

    def test_regrouping_in_place_gives_the_map_of_a_copy(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        copy = tmp_path / "copy.metacool"
        reassign(tiny_metacooler, {"contig_A": "MAG_2"}, copy)
        reassign(tiny_metacooler, {"contig_A": "MAG_2"}, inplace=True)
        in_place, written = MetaCooler(tiny_metacooler), MetaCooler(copy)
        pd.testing.assert_frame_equal(in_place.contigs, written.contigs)
        pd.testing.assert_frame_equal(in_place.genomes, written.genomes)
        pd.testing.assert_frame_equal(in_place.pixels()[:], written.pixels()[:])


class TestOutputs:
    @pytest.mark.parametrize("output", ["out.mmetacool::entries/x", "out.metacool::x"])
    def test_an_output_naming_a_group_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path, output: str
    ):
        target = f"{tmp_path}/{output}"
        with pytest.raises(ValueError, match="an output is a whole file"):
            reassign(tiny_metacooler, {"contig_D": "MAG_2"}, target)
        with pytest.raises(ValueError, match="an output is a whole file"):
            annotate(tiny_metacooler, target, drop=["contigs.nothing"])
        assert sorted(p.name for p in tmp_path.iterdir() if "out" in p.name) == []

    def test_an_output_given_as_its_file_root(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "out.metacool"
        reassign(tiny_metacooler, {"contig_D": "MAG_2"}, f"{out}::/")
        validate(out, full=True)

    def test_map_output(self, tmp_path: Path):
        assert map_output("a.metacool") == Path("a.metacool")
        assert map_output("a.metacool::/") == Path("a.metacool")
        assert map_output("a.metacool::") == Path("a.metacool")


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
class TestReads:
    def test_a_read_spanning_a_replacement_is_refused(self, tiny_metacooler: Path):
        """metacooler reads the pixels a chunk at a time, opening the file for
        each: a replacement between two chunks is refused, not mixed in."""
        mc = MetaCooler(tiny_metacooler)
        chunks = mc.iter_pixels(chunksize=2)
        next(chunks)
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError, match="replaced or written after"):
            next(chunks)

    @pytest.mark.parametrize(
        "read",
        [
            lambda mc: mc.pixels()[:],
            lambda mc: mc.bins()[:],
            lambda mc: mc.matrix(balance=False)[:],
            lambda mc: mc.info,
            lambda mc: mc.contigs,
            lambda mc: mc.genomes,
        ],
    )
    def test_every_read_of_a_replaced_map_is_refused(self, tiny_metacooler: Path, read):
        mc = MetaCooler(tiny_metacooler)
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError):
            read(mc)
        read(MetaCooler(tiny_metacooler))  # opened again: read

    def test_a_cooler_built_before_the_replacement(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        mc.pixels()[:1]
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError):
            mc.pixels()[:]

    def test_a_write_through_the_maps_own_open_is_not_a_change(
        self, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        with mc.open("r+") as grp:
            grp.attrs["note"] = "edited"
        assert len(mc.pixels()[:]) > 0
        with mc.open() as grp:
            assert grp.attrs["note"] == "edited"

    def test_a_write_of_another_process_is_a_change(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["note"] = "x" * 10_000
        with pytest.raises(MapChangedError):
            mc.pixels()[:]

    def test_the_checked_path(self, tiny_metacooler: Path):
        calls = []
        path = staging.StateCheckedPath(str(tiny_metacooler), lambda: calls.append(1))
        assert os.fspath(path) == str(tiny_metacooler) and calls == [1]
        assert str(path) == str(tiny_metacooler)
        assert repr(path) == repr(str(tiny_metacooler))
        assert staging.file_state(tiny_metacooler.with_name("nothing")) is None
