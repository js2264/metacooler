"""Locked, staged in-place writes, and reads that span a replacement.

metaTOR 2.0 updates a user's only copy of a map step by step, so every writer
edits a map in place through a staged copy, under the map's lock taken before
its first read (``metacooler.io.staging``).  These tests are metaTOR #64's
guarantees of its own staged writes, ported -- the lock (holders, stale locks
and the races of taking one over), the staged copy (failures, interrupts, the
mode, links, a map that does not grow), and a read across a replacement -- and
those of the independent review of the first version: the takeover of a stale
lock by many processes at once, threads, every writer behind the lock, a map
changed between the read and the write, a private staging directory, a forked
child, and what a re-grouping would drop.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import pickle
import socket
import threading
import time
from pathlib import Path
from typing import Any, cast

import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

import metacooler
from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io import staging
from metacooler.io.annotate import annotate
from metacooler.io.clusters import drop_clusters, write_clusters
from metacooler.io.mmetacool import (
    add_entry,
    copy_entry,
    copy_entry_cross_file,
    move_entry,
    remove_entry,
)
from metacooler.io.reassign import reassign
from metacooler.io.staging import (
    MapChangedError,
    MapLockedError,
    MapLockError,
    edit_metadata,
    holds_lock,
    lock,
    lock_path,
    map_output,
    staged_inplace,
)
from metacooler.io.weights import store_weights
from metacooler.ops.normalize import normalize_ice
from metacooler.ops.qc import qc
from metacooler.schema import MetaCoolerFormatError, validate

#: A pid no process of this host has (above any pid_max).
DEAD_PID = 999_999_999
SPAWN = multiprocessing.get_context("spawn")


def _write_lock(
    path: Path, pid: int, *, host: str | None = None, run: str = "deadrun", **extra
) -> dict:
    """A lock file of *path*, as a holder of *pid* on *host* wrote it."""
    holder = {
        "pid": pid,
        "host": host or socket.gethostname(),
        "start": "2026-11-20T10:00:00+0000",
        "command": "metator pipeline map.metacool",
        "uuid": run,
        **extra,
    }
    lock_path(path).write_text(json.dumps(holder))
    return holder


def _left_beside(path: Path) -> list[str]:
    """What a writer left beside the map: staging directories, partial files,
    locks and takeover guards."""
    return sorted(
        p.name
        for p in path.parent.iterdir()
        if p.name.startswith(".") or ".lock" in p.name
    )


def _touch(path: Path) -> None:
    """Change the file's modification time, as a write would (h5py could not
    open it for writing while another handle has it open)."""
    info = os.stat(path)
    os.utime(path, ns=(info.st_atime_ns, info.st_mtime_ns + 1_000_000_000))


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


#: Each in-place writer of a map, as a call on its path.
WRITERS = {
    "annotate": lambda path, tmp: _annotate(path, tmp),
    "store_weights": lambda path, tmp: _store(path),
    "normalize": lambda path, tmp: normalize_ice(path, name="ice", force=True),
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
    "qc --store": lambda path, tmp: qc(path, store=True),
    "metadata": lambda path, tmp: edit_metadata(
        path, lambda meta: meta.update({"sample": "S1"})
    ),
}

#: Each writer of a collection, as a call on its path.
COLLECTION_WRITERS = {
    "add_entry": lambda mm, src, tmp: add_entry(mm, "entry3", src),
    "remove_entry": lambda mm, src, tmp: remove_entry(mm, "entry2"),
    "copy_entry": lambda mm, src, tmp: copy_entry(mm, "entry1", "entry3"),
    "move_entry": lambda mm, src, tmp: move_entry(mm, "entry2", "entry3"),
    "copy_entry_cross_file": lambda mm, src, tmp: copy_entry_cross_file(
        _second_collection(tmp, src), "entry1", mm, "entry3"
    ),
}


def _second_collection(tmp: Path, src: Path) -> Path:
    from metacooler.io.mmetacool import create_mmetacool

    other = tmp / "other.mmetacool"
    if not other.exists():
        create_mmetacool(other, {"entry1": src})
    return other


@pytest.fixture
def holder(tiny_metacooler: Path):
    """Another process holding the map's lock, for the test's time."""
    from _staging_workers import hold

    ready, release = SPAWN.Event(), SPAWN.Event()
    process = SPAWN.Process(target=hold, args=(str(tiny_metacooler), ready, release))
    process.start()
    assert ready.wait(60), "the holder did not take the lock"
    yield process
    release.set()
    process.join(60)


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------
class TestLock:
    def test_held_for_the_block_and_named(self, tiny_metacooler: Path):
        with lock(tiny_metacooler) as lock_file:
            assert lock_file == tiny_metacooler.with_name("tiny.metacool.lock")
            holder = json.loads(lock_file.read_text())
            assert {"pid", "host", "start", "command", "uuid", "boot_id"} <= set(holder)
            assert holder["pid"] == os.getpid()
            assert holder["host"] == socket.gethostname()
            assert holder["pid_start"] == staging._pid_start(os.getpid())
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

    def test_one_lock_whatever_the_name_of_its_directory(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A map reached through a linked directory, or by a relative path, has
        the lock of its directory's real path; a chdir in the block changes
        nothing."""
        linked = tmp_path / "linked"
        linked.symlink_to(tiny_metacooler.parent)
        assert lock_path(linked / "tiny.metacool") == lock_path(tiny_metacooler)
        monkeypatch.chdir(tiny_metacooler.parent)
        with lock("tiny.metacool") as lock_file:
            assert lock_file == lock_path(tiny_metacooler)
            monkeypatch.chdir(tmp_path.parent)
            assert holds_lock(tiny_metacooler)
            with pytest.raises(MapLockedError):
                # Another thread's view: the same lock.
                _in_thread(lambda: lock(linked / "tiny.metacool").__enter__())
        assert not lock_path(tiny_metacooler).exists()

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

    def test_a_live_holder_is_refused_by_name(self, tiny_metacooler: Path, holder):
        with pytest.raises(MapLockedError) as caught, lock(tiny_metacooler):
            pass
        message = str(caught.value)
        assert message.startswith(f"{tiny_metacooler} is being written by another")
        assert f"pid {holder.pid} on {socket.gethostname()}" in message
        assert "test holder" in message
        assert f"delete {lock_path(tiny_metacooler)}" in message
        assert caught.value.holder["pid"] == holder.pid
        assert caught.value.lock_file == lock_path(tiny_metacooler)
        assert lock_path(tiny_metacooler).exists()

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
        self, tiny_metacooler: Path, text: str, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setattr(staging, "_EMPTY_WAIT", 0.05)
        lock_path(tiny_metacooler).write_text(text)
        with (
            pytest.raises(MapLockedError, match="does not say which process"),
            lock(tiny_metacooler),
        ):
            pass
        assert lock_path(tiny_metacooler).read_text() == text

    def test_a_lock_being_written_is_waited_for(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A lock seen empty may be one its creator is writing: read again."""
        lock_file = lock_path(tiny_metacooler)
        lock_file.write_text("")
        live = {"pid": os.getppid(), "host": socket.gethostname(), "uuid": "live"}

        def written_later():
            time.sleep(0.05)
            lock_file.write_text(json.dumps(live))

        thread = threading.Thread(target=written_later)
        thread.start()
        with (
            pytest.raises(MapLockedError, match="being written by another process"),
            lock(tiny_metacooler),
        ):
            pass
        thread.join()

    def test_a_dead_holders_lock_is_taken_over(
        self, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
    ):
        """A lock of a process of this host that is gone (killed outright) is
        taken over, and what it staged is removed: its own only."""
        _write_lock(tiny_metacooler, DEAD_PID, run="deadrun")
        mine = tiny_metacooler.parent
        staged = mine / ".tiny.metacool.deadrun.1.d"
        staged.mkdir()
        (staged / "tiny.metacool").write_text("half a map")
        partial = mine / ".tiny.metacool.deadrun.1.staged"
        other = mine / ".tiny.metacool.otherrun.1.d"
        partial.write_text("x")
        other.mkdir()
        with (
            caplog.at_level("WARNING", logger="metacooler"),
            lock(tiny_metacooler) as lock_file,
        ):
            assert json.loads(lock_file.read_text())["pid"] == os.getpid()
        assert f"Took over {lock_path(tiny_metacooler)}" in caplog.text
        assert f"pid {DEAD_PID}" in caplog.text
        assert "Removed .tiny.metacool.deadrun.1.d," in caplog.text
        assert not staged.exists() and not partial.exists() and other.exists()
        assert not lock_path(tiny_metacooler).exists()

    @pytest.mark.parametrize("why", ["own pid", "earlier boot", "reused pid"])
    def test_what_is_stale(self, tiny_metacooler: Path, why: str):
        """A lock naming this process but none of its locks (an earlier
        process of its pid), one of an earlier boot, and one whose pid a later
        process reuses, are stale: taken over."""
        if why == "own pid":
            _write_lock(tiny_metacooler, os.getpid(), run="foreign")
        elif why == "earlier boot":
            _write_lock(tiny_metacooler, os.getppid(), boot_id="an-earlier-boot")
        else:
            start = staging._pid_start(os.getppid())
            if start is None:
                pytest.skip("no /proc here")
            _write_lock(tiny_metacooler, os.getppid(), pid_start=start + 1)
        with lock(tiny_metacooler) as lock_file:
            assert json.loads(lock_file.read_text())["uuid"] != "foreign"

    def test_a_live_pid_with_its_start_time_is_live(self, tiny_metacooler: Path):
        _write_lock(
            tiny_metacooler,
            os.getppid(),
            boot_id=staging._boot_id(),
            pid_start=staging._pid_start(os.getppid()),
        )
        with pytest.raises(MapLockedError), lock(tiny_metacooler):
            pass

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

    def test_a_lock_file_is_never_left_empty(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        def fail(*args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr(staging.json, "dump", fail)
        with pytest.raises(MapLockError, match="No space left"), lock(tiny_metacooler):
            pass
        assert not lock_path(tiny_metacooler).exists()

    def test_the_errors(self):
        assert issubclass(MapLockedError, MapLockError)
        assert issubclass(MapLockError, OSError)
        assert issubclass(MapChangedError, OSError)
        assert metacooler.lock is lock
        assert metacooler.MapChangedError is MapChangedError
        assert str(MapLockError("why", lock_file="x.lock")) == "why"
        assert str(MapChangedError("why")) == "why"

    @pytest.mark.parametrize("cls", [MapLockError, MapLockedError])
    def test_the_errors_survive_pickling(self, cls):
        """A lock error raised in a worker process reaches its parent."""
        error = cls("why", lock_file="/x/m.lock", holder={"pid": 3})
        back = pickle.loads(pickle.dumps(error))
        assert type(back) is cls and str(back) == "why"
        assert back.lock_file == Path("/x/m.lock") and back.holder == {"pid": 3}

    def test_pid_checks(self, monkeypatch: pytest.MonkeyPatch):
        assert staging._pid_is_running(os.getpid())
        assert not staging._pid_is_running(DEAD_PID)
        assert not staging._pid_is_running(2**70)
        assert staging._pid_start(DEAD_PID) is None

        def another_users(pid, signal):
            raise PermissionError(1, "Operation not permitted")

        # A process of another user, which this one may not signal: alive.
        monkeypatch.setattr(staging.os, "kill", another_users)
        assert staging._pid_is_running(1)

    def test_what_proc_cannot_tell(self, tmp_path: Path, monkeypatch):
        """Without /proc, or with a stat line it cannot parse, there is no boot
        id nor start time: the pid alone is checked."""
        assert staging._read(tmp_path) == {}  # a directory: unreadable
        staging._fsync(tmp_path / "missing")  # nothing to sync: no error
        real = Path.read_text

        def garbled(self, *args, **kwargs):
            if str(self).startswith("/proc/"):
                return "12 (name) S"
            return real(self, *args, **kwargs)

        monkeypatch.setattr(Path, "read_text", garbled)
        assert staging._pid_start(os.getpid()) is None
        assert staging._boot_id() == "12 (name) S"

        def no_proc(self, *args, **kwargs):
            if str(self).startswith("/proc/"):
                raise OSError("no /proc")
            return real(self, *args, **kwargs)

        monkeypatch.setattr(Path, "read_text", no_proc)
        assert staging._boot_id() is None and staging._pid_start(os.getpid()) is None

    def test_a_lock_gone_while_read_is_contended_again(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A lock or a guard released between a failed create and its read:
        contend again."""
        _write_lock(tiny_metacooler, DEAD_PID, run="dead")
        guard = Path(f"{lock_path(tiny_metacooler)}.takeover")
        guard.write_text("{}")
        real = staging._read_settled
        gone = {"lock": False, "guard": False}

        def released_once(path):
            for what, name in (("guard", ".takeover"), ("lock", ".lock")):
                if str(path).endswith(name) and not gone[what]:
                    gone[what] = True
                    if what == "guard":
                        guard.unlink()
                    return staging._GONE
            return real(path)

        monkeypatch.setattr(staging, "_read_settled", released_once)
        with lock(tiny_metacooler) as lock_file:
            assert json.loads(lock_file.read_text())["pid"] == os.getpid()
        assert gone == {"lock": True, "guard": True}

    def test_no_takeover_on_windows(self, monkeypatch: pytest.MonkeyPatch):
        """os.kill(pid, 0) terminates a process on Windows: no lock is judged
        stale there."""
        holder = {"pid": DEAD_PID, "host": socket.gethostname(), "uuid": "x"}
        assert staging._is_stale(holder)
        monkeypatch.setattr(staging.os, "name", "nt")
        assert not staging._is_stale(holder)

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


def _in_thread(call):
    """Run *call* in another thread, raising what it raised."""
    errors: list[BaseException] = []

    def run():
        try:
            call()
        except BaseException as error:  # noqa: BLE001 - re-raised below
            errors.append(error)

    thread = threading.Thread(target=run, name="other-thread")
    thread.start()
    thread.join()
    if errors:
        raise errors[0]


# ---------------------------------------------------------------------------
# Taking over a stale lock, with many processes at once
# ---------------------------------------------------------------------------
class TestTakeover:
    def test_eight_processes_one_holder_every_round(self, tiny_metacooler: Path):
        """Eight spawned processes, released together by a barrier, contend for
        a lock whose holder is gone, round after round: one exactly holds it
        each round (the first version let two hold it at once with three
        contenders), and nothing is left once it is released."""
        from _staging_workers import contend

        n, rounds = 8, 10
        barrier, results = SPAWN.Barrier(n + 1), SPAWN.Queue()
        processes = [
            SPAWN.Process(
                target=contend, args=(str(tiny_metacooler), rounds, barrier, results)
            )
            for _ in range(n)
        ]
        for process in processes:
            process.start()
        try:
            for r in range(rounds):
                _write_lock(tiny_metacooler, DEAD_PID, run=f"dead{r}")
                barrier.wait(60)  # contend
                barrier.wait(60)  # every contender has tried
                barrier.wait(60)  # the holder has released
                assert _left_beside(tiny_metacooler) == []
        finally:
            for process in processes:
                process.join(60)
        outcomes = [results.get(timeout=10) for _ in range(n * rounds)]
        for r in range(rounds):
            held = [pid for rr, pid, what in outcomes if rr == r and what == "held"]
            assert len(held) == 1, (r, outcomes)

    def test_the_stale_lock_only_is_removed(self, tiny_metacooler: Path):
        """Under the guard, the lock is removed only if it is still the stale
        one: one another process has created meanwhile is left."""
        stale = _write_lock(tiny_metacooler, DEAD_PID, run="dead")
        live = _write_lock(tiny_metacooler, os.getppid(), run="live")
        staging._take_over(lock_path(tiny_metacooler), tiny_metacooler, stale)
        assert json.loads(lock_path(tiny_metacooler).read_text()) == live
        lock_path(tiny_metacooler).unlink()
        # Gone meanwhile: nothing to remove.
        staging._take_over(lock_path(tiny_metacooler), tiny_metacooler, stale)
        assert not lock_path(tiny_metacooler).exists()
        assert _left_beside(tiny_metacooler) == []

    def test_a_guard_left_by_a_dead_process_is_refused_by_name(
        self, tiny_metacooler: Path
    ):
        stale = _write_lock(tiny_metacooler, DEAD_PID, run="dead")
        guard = Path(f"{lock_path(tiny_metacooler)}.takeover")
        guard.write_text(
            json.dumps({"pid": DEAD_PID, "host": socket.gethostname(), "uuid": "g"})
        )
        with (
            pytest.raises(MapLockedError, match="stopped while it took over") as caught,
            lock(tiny_metacooler),
        ):
            pass
        assert f"delete {guard}" in str(caught.value)
        assert (
            guard.exists()
            and json.loads(lock_path(tiny_metacooler).read_text()) == stale
        )

    def test_a_live_takeover_is_waited_for(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A guard held by a live process: wait for it, then contend; one that
        stays too long is refused, the map's holder named."""
        _write_lock(tiny_metacooler, DEAD_PID, run="dead")
        guard = Path(f"{lock_path(tiny_metacooler)}.takeover")
        guard.write_text(
            json.dumps({"pid": os.getppid(), "host": socket.gethostname(), "uuid": "g"})
        )
        monkeypatch.setattr(staging, "_TAKEOVER_WAIT", 0.1)
        with (
            pytest.raises(MapLockedError, match="has been taking over"),
            lock(tiny_metacooler),
        ):
            pass

        def released_later():
            time.sleep(0.05)
            guard.unlink()

        monkeypatch.setattr(staging, "_TAKEOVER_WAIT", 10)
        thread = threading.Thread(target=released_later)
        thread.start()
        with lock(tiny_metacooler) as lock_file:
            assert json.loads(lock_file.read_text())["pid"] == os.getpid()
        thread.join()

    def test_an_empty_guard_is_refused(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setattr(staging, "_EMPTY_WAIT", 0.05)
        _write_lock(tiny_metacooler, DEAD_PID, run="dead")
        Path(f"{lock_path(tiny_metacooler)}.takeover").write_text("")
        with (
            pytest.raises(MapLockedError, match="stopped while it took over"),
            lock(tiny_metacooler),
        ):
            pass


# ---------------------------------------------------------------------------
# Threads
# ---------------------------------------------------------------------------
class TestThreads:
    def test_another_thread_is_refused_by_name(self, tiny_metacooler: Path):
        holding, done = threading.Event(), threading.Event()

        def hold():
            with lock(tiny_metacooler):
                holding.set()
                done.wait(30)

        thread = threading.Thread(target=hold, name="holder-thread")
        thread.start()
        try:
            assert holding.wait(30)
            with (
                pytest.raises(
                    MapLockedError,
                    match=r"another thread of this process \(holder-thread\)",
                ),
                lock(tiny_metacooler),
            ):
                pass
            before = tiny_metacooler.read_bytes()
            with pytest.raises(MapLockedError):
                _store(tiny_metacooler)
            assert tiny_metacooler.read_bytes() == before
        finally:
            done.set()
            thread.join()
        with lock(tiny_metacooler):
            pass

    def test_two_threads_one_holder(self, tiny_metacooler: Path):
        barrier, release = threading.Barrier(2), threading.Event()
        outcomes: list[str] = []

        def contend():
            barrier.wait()
            try:
                with lock(tiny_metacooler):
                    outcomes.append("held")
                    release.wait(5)
            except MapLockedError:
                outcomes.append("refused")
                release.set()

        threads = [threading.Thread(target=contend) for _ in range(2)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        assert sorted(outcomes) == ["held", "refused"]


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
    @pytest.mark.parametrize("mode", [0o600, 0o640, 0o664])
    def test_the_map_keeps_its_mode(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str, mode: int
    ):
        tiny_metacooler.chmod(mode)
        WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.stat().st_mode & 0o777 == mode

    @pytest.mark.parametrize("writer", list(WRITERS))
    def test_a_read_only_map_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str
    ):
        """Replacing a 0444 map would succeed in a writable directory, where an
        edit in place failed: it is refused, as it was."""
        tiny_metacooler.chmod(0o444)
        before = tiny_metacooler.read_bytes()
        with pytest.raises(PermissionError, match="read-only"):
            WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.skipif(os.geteuid() != 0, reason="only root may give a file away")
    def test_the_map_keeps_its_owner(self, tiny_metacooler: Path):
        os.chown(tiny_metacooler, 1234, 4321)
        _store(tiny_metacooler)
        info = tiny_metacooler.stat()
        assert (info.st_uid, info.st_gid) == (1234, 4321)

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

    def test_the_staged_map_is_private(self, tiny_metacooler: Path):
        """The staged copy of a private map is never readable by others, while
        it is written or after a kill: its directory is 0700."""
        tiny_metacooler.chmod(0o600)
        with staged_inplace(tiny_metacooler) as staged:
            assert staged.parent.stat().st_mode & 0o777 == 0o700
            assert staged.parent.name.startswith(".tiny.metacool.")
            assert staged.parent.parent == tiny_metacooler.parent
        assert tiny_metacooler.stat().st_mode & 0o777 == 0o600

    @pytest.mark.parametrize(
        "stop", [RuntimeError("boom"), KeyboardInterrupt(), SystemExit(1)]
    )
    def test_a_failure_leaves_the_map_as_it_was(
        self, tiny_metacooler: Path, stop: BaseException
    ):
        """Whatever stops the write, the staging directory and what cooler
        wrote in it are removed, and the map is as it was."""
        before = tiny_metacooler.read_bytes()
        with pytest.raises(type(stop)), staged_inplace(tiny_metacooler) as staged:
            Path(f"{staged}.__writing__").write_text("partial")
            (staged.parent / "tmpabc.multi.cool").write_text("partial")
            with h5py.File(staged, "r+") as f:
                f.attrs["edited"] = True
            raise stop
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    @pytest.mark.parametrize("stop", [KeyboardInterrupt(), SystemExit(1)])
    def test_an_interrupt_as_the_map_is_replaced(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch, stop
    ):
        before = tiny_metacooler.read_bytes()

        def interrupted(*args):
            raise stop

        monkeypatch.setattr(staging.os, "replace", interrupted)
        with pytest.raises(type(stop)):
            _store(tiny_metacooler)
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_an_interrupted_regroup_leaves_nothing(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """cooler's temporary files of a re-grouping are in the staging
        directory, removed with it."""
        from metacooler.io import reassign as reassign_module

        def interrupted(*args, **kwargs):
            raise KeyboardInterrupt

        monkeypatch.setattr(reassign_module, "_finalize_metacool_file", interrupted)
        before = tiny_metacooler.read_bytes()
        with pytest.raises(KeyboardInterrupt):
            reassign(tiny_metacooler, {"contig_A": "MAG_2"}, inplace=True)
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

    def test_the_mode_of_a_new_file(self, tiny_metacooler: Path, tmp_path: Path):
        """A path where no map is yet gets the mode of any new file, from the
        umask, not the private mode of the staging directory."""
        path = tmp_path / "new.metacool"
        previous = os.umask(0o027)
        try:
            assert staging._umask() == 0o027
            with staged_inplace(path, copy=False) as staged:
                staging.copy_fresh(tiny_metacooler, staged)
        finally:
            os.umask(previous)
        assert path.stat().st_mode & 0o777 == 0o640

    def test_the_umask_without_proc(self, monkeypatch: pytest.MonkeyPatch):
        def no_proc(self, *args, **kwargs):
            raise OSError("no /proc")

        previous = os.umask(0o027)
        try:
            monkeypatch.setattr(Path, "read_text", no_proc)
            assert staging._umask() == 0o027
        finally:
            os.umask(previous)

    def test_the_staged_map_is_synced(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        synced = []
        real = os.fsync
        monkeypatch.setattr(
            staging.os, "fsync", lambda fd: (synced.append(fd), real(fd))
        )
        _store(tiny_metacooler)
        # The lock file, the staged map, and the directory.
        assert len(synced) >= 3

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
    def test_a_writer_is_refused_a_map_another_process_holds(
        self, tiny_metacooler: Path, tmp_path: Path, writer: str, holder
    ):
        """Another process holds the lock (a spawned one): every writer is
        refused before it reads the map, which is left as it was."""
        before = tiny_metacooler.read_bytes()
        with pytest.raises(MapLockedError, match="test holder"):
            WRITERS[writer](tiny_metacooler, tmp_path)
        assert tiny_metacooler.read_bytes() == before

    def test_the_cli_writers_are_refused_too(self, tiny_metacooler: Path, holder):
        for argv in (
            ["metadata", "add", str(tiny_metacooler), "-k", "sample", "--value", "S1"],
            ["qc", str(tiny_metacooler), "--store", "-q"],
        ):
            result = CliRunner().invoke(cli, argv)
            assert result.exit_code != 0
            assert "being written by another process" in flat(result.output)

    def test_open_for_writing_takes_the_lock(self, tiny_metacooler: Path, holder):
        """MetaCooler.open in a write mode edits the file itself (unstaged),
        under the lock."""
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(MapLockedError), mc.open("r+"):
            pass

    def test_open_for_writing_refuses_a_changed_map(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError), mc.open("r+"):
            pass
        with MetaCooler(tiny_metacooler).open("r+") as grp:
            grp.attrs["note"] = "x"

    def test_the_next_writer_removes_what_a_killed_one_left(
        self, tiny_metacooler: Path
    ):
        _write_lock(tiny_metacooler, DEAD_PID, run="deadrun")
        left = tiny_metacooler.with_name(".tiny.metacool.deadrun.1.d")
        left.mkdir()
        (left / "tiny.metacool").write_text("half a map")
        (left / "tmp8x2.multi.cool").write_text("cooler's")
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

    @pytest.mark.parametrize("writer", list(COLLECTION_WRITERS))
    def test_every_collection_writer_is_staged(
        self, tiny_mmetacool: Path, tiny_metacooler: Path, tmp_path: Path, writer
    ):
        before = tiny_mmetacool.stat().st_ino
        COLLECTION_WRITERS[writer](tiny_mmetacool, tiny_metacooler, tmp_path)
        assert tiny_mmetacool.stat().st_ino != before
        metacooler.validate_mmetacooler(tiny_mmetacool)
        assert _left_beside(tiny_mmetacool) == []

    @pytest.mark.parametrize("writer", list(COLLECTION_WRITERS))
    def test_every_collection_writer_takes_the_lock(
        self, tiny_mmetacool: Path, tiny_metacooler: Path, tmp_path: Path, writer
    ):
        from _staging_workers import hold

        ready, release = SPAWN.Event(), SPAWN.Event()
        process = SPAWN.Process(target=hold, args=(str(tiny_mmetacool), ready, release))
        process.start()
        try:
            assert ready.wait(60)
            before = tiny_mmetacool.read_bytes()
            with pytest.raises(MapLockedError):
                COLLECTION_WRITERS[writer](tiny_mmetacool, tiny_metacooler, tmp_path)
            assert tiny_mmetacool.read_bytes() == before
        finally:
            release.set()
            process.join(60)

    def test_a_collection_source_changed_while_copied(
        self,
        tiny_mmetacool: Path,
        tiny_metacooler: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        from metacooler.io import mmetacool

        copy = mmetacool._copy_map

        def then_changed(*args):
            copy(*args)
            _touch(tiny_metacooler)

        monkeypatch.setattr(mmetacool, "_copy_map", then_changed)
        before = tiny_mmetacool.read_bytes()
        with pytest.raises(MapChangedError):
            add_entry(tiny_mmetacool, "entry3", tiny_metacooler)
        assert tiny_mmetacool.read_bytes() == before

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

    def test_metadata_stored_as_bytes(self, tiny_metacooler: Path):
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["metadata"] = np.bytes_(json.dumps({"sample": "S0"}))
        edit_metadata(tiny_metacooler, lambda meta: meta.update({"sample": "S1"}))
        assert MetaCooler(tiny_metacooler).metadata["sample"] == "S1"

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

    def test_regrouping_an_entry_in_place_is_refused(self, tiny_mmetacool: Path):
        before = tiny_mmetacool.read_bytes()
        with pytest.raises(ValueError, match="cannot re-group an entry"):
            reassign(
                f"{tiny_mmetacool}::entries/entry1", {"contig_A": "MAG_2"}, inplace=True
            )
        assert tiny_mmetacool.read_bytes() == before

    def test_the_cli_refusal_names_its_options(self, tiny_mmetacool: Path, tmp_path):
        table = tmp_path / "t.tsv"
        table.write_text("contig_name\tgenome_name\ncontig_A\tMAG_2\n")
        result = CliRunner().invoke(
            cli,
            [
                "reassign",
                f"{tiny_mmetacool}::entries/entry1",
                "--contigs",
                str(table),
                "--inplace",
            ],
        )
        assert result.exit_code != 0
        assert "-o" in flat(result.output) and "--layout keep" in flat(result.output)

    def test_relabelling_an_entry_in_place(self, tiny_mmetacool: Path):
        entry = f"{tiny_mmetacool}::entries/entry1"
        reassign(entry, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        genome = MetaCooler(entry).contigs.set_index("name")["genome_id"]
        names = MetaCooler(entry).genomes["name"]
        assert names[genome["contig_D"]] == "MAG_2"
        validate(f"{tiny_mmetacool}::entries/entry2", full=True)


# ---------------------------------------------------------------------------
# What is read is what is written
# ---------------------------------------------------------------------------
class TestPlanned:
    def test_a_map_changed_since_it_was_read_is_refused(self, tiny_metacooler: Path):
        """A writer hands the map it read to the staged write, which refuses
        to copy a file that changed since (another process's unlocked edit,
        here)."""
        mc = MetaCooler(tiny_metacooler)
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["changed"] = "x" * 10_000
        before = tiny_metacooler.read_bytes()
        with (
            pytest.raises(MapChangedError),
            staged_inplace(tiny_metacooler, planned=mc),
        ):
            pass
        assert tiny_metacooler.read_bytes() == before
        assert _left_beside(tiny_metacooler) == []

    def test_weights_computed_from_a_changed_map(self, tiny_metacooler: Path):
        source = MetaCooler(tiny_metacooler)
        n = _n_bins(tiny_metacooler)
        _store(tiny_metacooler, name="other")
        with pytest.raises(MapChangedError):
            store_weights(
                tiny_metacooler,
                "w",
                np.ones(n),
                method="test",
                depends_on="pixels",
                source=source,
            )

    def test_normalize_holds_the_lock_while_it_computes(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """The weights are stored under the lock taken before the map was
        read: another writer is refused meanwhile (it was a misaligned map)."""
        from metacooler.ops import normalize as normalize_module

        compute = normalize_module._ice_by_genome
        seen = []

        def while_computing(mc, **kwargs):
            seen.append(holds_lock(tiny_metacooler))
            with pytest.raises(MapLockedError):
                _in_thread(lambda: _store(tiny_metacooler, name="other"))
            return compute(mc, **kwargs)

        monkeypatch.setattr(normalize_module, "_ice_by_genome", while_computing)
        normalize_ice(tiny_metacooler)
        assert seen == [True]
        assert not holds_lock(tiny_metacooler)

    @pytest.mark.parametrize("writer", ["reassign, relabel", "reassign, regroup"])
    def test_reassign_holds_the_lock_from_its_first_read(
        self, tiny_metacooler: Path, tmp_path, writer, monkeypatch: pytest.MonkeyPatch
    ):
        from metacooler.io import reassign as reassign_module

        build = reassign_module._build_metadata
        seen = []

        def while_planning(*args, **kwargs):
            seen.append(holds_lock(tiny_metacooler))
            return build(*args, **kwargs)

        monkeypatch.setattr(reassign_module, "_build_metadata", while_planning)
        WRITERS[writer](tiny_metacooler, tmp_path)
        assert seen == [True]

    @pytest.mark.parametrize("copier", ["annotate", "reassign"])
    def test_a_copy_by_path_of_a_changed_map_is_refused(
        self, tiny_metacooler: Path, tmp_path: Path, copier, monkeypatch
    ):
        """annotate and reassign with an output copy the map by path after
        reading it: a map changed meanwhile is refused."""
        import shutil as shutil_module

        copy2 = shutil_module.copy2

        def changed_first(source, target, *args, **kwargs):
            with h5py.File(source, "r+") as f:
                f.attrs["changed"] = "x" * 10_000
            return copy2(source, target, *args, **kwargs)

        table = tmp_path / "coverage.tsv"
        names = MetaCooler(tiny_metacooler).contigs["name"].astype(str)
        table.write_text(
            "contig_name\tcoverage\n" + "".join(f"{n}\t1\n" for n in names)
        )
        monkeypatch.setattr(shutil_module, "copy2", changed_first)
        out = tmp_path / "out.metacool"
        with pytest.raises(MapChangedError):
            if copier == "annotate":
                annotate(tiny_metacooler, out, contigs_table=table)
            else:
                reassign(tiny_metacooler, {"contig_D": "MAG_2"}, out, layout="keep")
        assert not out.exists()


# ---------------------------------------------------------------------------
# What a re-grouping does not carry
# ---------------------------------------------------------------------------
class TestRegroupCarries:
    @staticmethod
    def _add_extras(path: Path, extra: str) -> str:
        """Give the map one thing a re-grouping does not write; return how the
        refusal names it."""
        with h5py.File(path, "r+") as f:
            if extra == "group":
                f.create_group("metator").create_dataset("scores", data=[1, 2])
                return "/metator"
            if extra == "root attribute":
                f.attrs["sample"] = "S1"
                return "the root attribute 'sample'"
            if extra == "pixel column":
                n = f["pixels/count"].shape[0]
                f["pixels"].create_dataset("balanced_by_tool", data=np.ones(n))
                return "/pixels/balanced_by_tool"
            if extra == "table attribute":
                f["contigs"].attrs["provenance"] = "assembly v2"
                return "the attributes of /contigs (provenance)"
            if extra == "column attribute":
                f["contigs/length"].attrs["unit"] = "bp"
                return "the attributes of /contigs/length (unit)"
            n = f["indexes/chrom_offset"].shape[0]
            f["indexes"].create_dataset("extra_offset", data=np.zeros(n))
            return "/indexes/extra_offset"

    EXTRAS = [
        "group",
        "root attribute",
        "pixel column",
        "table attribute",
        "column attribute",
        "index",
    ]

    @pytest.mark.parametrize("extra", EXTRAS)
    def test_an_in_place_regroup_that_would_drop_is_refused_by_name(
        self, tiny_metacooler: Path, extra: str
    ):
        named = self._add_extras(tiny_metacooler, extra)
        before = tiny_metacooler.read_bytes()
        with pytest.raises(ValueError, match="would re-group its contigs") as caught:
            reassign(tiny_metacooler, {"contig_A": "MAG_2"}, inplace=True)
        assert named in str(caught.value)
        assert "--layout keep" in str(caught.value) and "-o" in str(caught.value)
        assert tiny_metacooler.read_bytes() == before

    @pytest.mark.parametrize("extra", EXTRAS)
    def test_a_regroup_to_a_new_file_warns_of_what_it_drops(
        self, tiny_metacooler: Path, tmp_path: Path, extra, caplog
    ):
        named = self._add_extras(tiny_metacooler, extra)
        with caplog.at_level("WARNING", logger="metacooler"):
            reassign(tiny_metacooler, {"contig_A": "MAG_2"}, tmp_path / "out.metacool")
        assert named in caplog.text

    def test_a_relabel_in_place_keeps_everything(self, tiny_metacooler: Path):
        for extra in self.EXTRAS:
            self._add_extras(tiny_metacooler, extra)
        reassign(tiny_metacooler, {"contig_D": "MAG_2"}, inplace=True, layout="keep")
        with h5py.File(tiny_metacooler, "r") as f:
            assert "metator" in f and f.attrs["sample"] == "S1"
            assert "balanced_by_tool" in f["pixels"]

    def test_a_map_metacooler_made_carries_everything(self, tiny_metacooler: Path):
        """What metacooler writes itself -- weights, rounds, QC, metadata -- a
        re-grouping carries: nothing to refuse."""
        from metacooler.io.reassign import _not_carried_by_regroup

        normalize_ice(tiny_metacooler, scope="contig", name="ice")
        write_clusters(tiny_metacooler, *_rounds(tiny_metacooler))
        qc(tiny_metacooler, store=True)
        edit_metadata(tiny_metacooler, lambda meta: meta.update({"sample": "S1"}))
        assert _not_carried_by_regroup(MetaCooler(tiny_metacooler)) == []


# ---------------------------------------------------------------------------
# A forked child
# ---------------------------------------------------------------------------
class _ChildExit(Exception):
    pass


@pytest.mark.skipif(not hasattr(os, "fork"), reason="needs os.fork")
@pytest.mark.filterwarnings(
    "ignore:This process .* is multi-threaded:DeprecationWarning"
)
class TestFork:
    def test_a_child_leaving_by_an_error_touches_nothing(self, tiny_metacooler: Path):
        """A child forked inside an edit that leaves it by an exception removes
        neither its parent's staging directory nor its lock."""
        try:
            with staged_inplace(tiny_metacooler) as staged:
                pid = os.fork()
                if pid == 0:
                    raise _ChildExit
                os.waitpid(pid, 0)
                assert staged.parent.is_dir() and staged.exists()
                assert lock_path(tiny_metacooler).exists()
        except _ChildExit:
            os._exit(0)
        validate(tiny_metacooler, full=True)
        assert _left_beside(tiny_metacooler) == []

    def test_a_child_leaving_normally_replaces_nothing(self, tiny_metacooler: Path):
        before = tiny_metacooler.stat().st_ino
        with staged_inplace(tiny_metacooler) as staged:
            with h5py.File(staged, "r+") as f:
                f.attrs["edited"] = True
            pid = os.fork()
            if pid == 0:
                # Leaves the block normally, as its parent will.
                pass
            else:
                os.waitpid(pid, 0)
                assert tiny_metacooler.stat().st_ino == before
                assert staged.exists()
        if pid == 0:
            os._exit(0)
        with h5py.File(tiny_metacooler, "r") as f:
            assert f.attrs["edited"]


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
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

    def test_map_output(self):
        assert map_output("a.metacool") == Path("a.metacool")
        assert map_output("a.metacool::/") == Path("a.metacool")
        assert map_output("a.metacool::") == Path("a.metacool")
        assert map_output("a.metacool:: / ") == Path("a.metacool")
        with pytest.raises(ValueError, match="names the group 'entries/x'"):
            map_output("a.mmetacool:: entries/x ")


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

    def test_a_write_during_a_read_is_refused(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Checked after the read too: a file written while a slice was read
        is refused, not returned."""
        mc = MetaCooler(tiny_metacooler)
        pixels = cast(Any, mc._clr.pixels())
        selector = pixels._selector
        real = type(selector).__getitem__

        def written_meanwhile(self, key):
            out = real(self, key)
            with h5py.File(tiny_metacooler, "r+") as f:
                f.attrs["changed"] = "x" * 10_000
            return out

        monkeypatch.setattr(type(selector), "__getitem__", written_meanwhile)
        with pytest.raises(MapChangedError):
            pixels[0:2]

    @pytest.mark.parametrize(
        "read",
        [
            lambda mc: mc.pixels()[:],
            lambda mc: mc.bins()[:],
            lambda mc: mc.matrix(balance=False)[:],
            lambda mc: mc.matrix(balance=False).fetch("contig_A"),
            lambda mc: mc.chroms()[:],
            lambda mc: mc.extent("contig_A"),
            lambda mc: mc.offset("contig_A"),
            lambda mc: mc.info,
            lambda mc: mc._clr.info,
            lambda mc: mc._clr.pixels().columns,
            lambda mc: mc._clr.pixels().keys(),
            lambda mc: "count" in mc._clr.pixels(),
            lambda mc: mc._clr.pixels()["count"][0:2],
            lambda mc: mc._clr._load_dset("indexes/chrom_offset"),
            lambda mc: mc.contigs,
            lambda mc: mc.genomes,
        ],
    )
    def test_every_read_of_a_replaced_map_is_refused(self, tiny_metacooler: Path, read):
        mc = MetaCooler(tiny_metacooler)
        mc._clr  # noqa: B018 - built before the replacement
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError):
            read(mc)
        read(MetaCooler(tiny_metacooler))  # opened again: read

    def test_a_selector_is_still_one(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        pixels = cast(Any, mc._clr.pixels())
        assert len(pixels) == int(mc.info["nnz"])
        assert pixels.shape == pixels._selector.shape
        assert "RangeSelector1D" in repr(pixels)
        assert pixels["count"][0:2].tolist() == mc.pixels()[0:2]["count"].tolist()

    def test_a_cooler_built_before_the_replacement(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        mc.pixels()[:1]
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError):
            mc.pixels()[:]

    def test_a_cooler_built_after_the_replacement(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        _store(tiny_metacooler)
        with pytest.raises(MapChangedError):
            mc._clr  # noqa: B018

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

    def test_a_read_open_checks_when_it_ends(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(MapChangedError), mc.open() as grp:
            grp.attrs.get("format")
            _touch(tiny_metacooler)

    def test_file_state(self, tiny_metacooler: Path):
        assert staging.file_state(tiny_metacooler.with_name("nothing")) is None
        assert staging.file_state(tiny_metacooler) is not None
