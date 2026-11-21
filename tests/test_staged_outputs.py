"""Writers of a new map at an output path: what they check, and what they leave.

A writer of an output -- ``create``, ``collect``, ``reassign -o``, ``binnify``,
``flatten``, ``merge``, a new collection, ... -- takes no lock of its own: it
refuses an output whose lock another thread, a live process, or a process
this host cannot check holds, before it writes and again right before its
final rename; a writer killed outright leaves nothing that stops the next one.
A symbolic link at an output is replaced, never written through.
"""

from __future__ import annotations

import errno
import os
import shutil
import threading
from pathlib import Path

import h5py
import numpy as np
import pytest
from test_staged_edit import _cli, _CliFailed
from test_staging import (
    DEAD_RUN,
    SPAWN,
    _store,
    _write_lock,
    holder,  # noqa: F401 - a fixture
)

import metacooler
from metacooler.api import MetaCooler
from metacooler.io import staging
from metacooler.io.mmetacool import create_mmetacool, list_entries
from metacooler.io.staging import (
    MapLockedError,
    MapLockError,
    lock,
    lock_path,
    staged_edit,
)
from metacooler.schema import validate


def _output_lock(out: Path) -> Path:
    return lock_path(out, follow_symlinks=False)


def _kill_while_writing(src: str, out: str) -> None:
    """binnify *src* into *out*, the process killed outright in the middle of
    the write (the OOM killer, a scheduler's time limit)."""
    import signal

    import metacooler.ops.binnify as binnify_module

    def killed(*args, **kwargs):
        os.kill(os.getpid(), signal.SIGKILL)

    setattr(binnify_module, "_binnify", killed)  # noqa: B010 - a test double
    binnify_module.binnify(src, out, 200, force=True)


# ---------------------------------------------------------------------------
# A writer of an output takes no lock
# ---------------------------------------------------------------------------
class TestNoLockOfItsOwn:
    def test_a_writer_creates_no_lock(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        import metacooler.ops.binnify as binnify_module

        out = tmp_path / "out.metacool"
        seen = []
        real = binnify_module._binnify

        def watched(*args, **kwargs):
            seen.append(_output_lock(out).exists())
            return real(*args, **kwargs)

        monkeypatch.setattr(binnify_module, "_binnify", watched)
        metacooler.binnify(tiny_metacooler, out, 200)
        assert seen == [False]
        validate(out, full=True)

    def test_a_lock_this_thread_holds_is_its_own(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """A run that holds its output's lock (metaTOR's --output-map) writes
        it with every writer; the lock is still its own after."""
        out = tmp_path / "out.metacool"
        with lock(out, follow_symlinks=False) as held:
            metacooler.binnify(tiny_metacooler, out, 200)
            metacooler.flatten(tiny_metacooler, tmp_path / "flat.cool")
            assert held.exists()
        assert not held.exists()
        validate(out, full=True)

    @pytest.mark.skipif(not hasattr(os, "fork"), reason="needs os.fork")
    @pytest.mark.filterwarnings(
        "ignore:This process .* is multi-threaded:DeprecationWarning"
    )
    def test_a_writer_killed_outright_leaves_nothing_in_the_way(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """Killed in a container, or on another node, it leaves no lock that
        a rerun there could not take over: the rerun is not refused."""
        out = tmp_path / "out.metacool"
        process = SPAWN.Process(
            target=_kill_while_writing, args=(str(tiny_metacooler), str(out))
        )
        process.start()
        process.join(60)
        assert process.exitcode == -9
        assert not _output_lock(out).exists()
        assert sorted(p.name for p in tmp_path.iterdir() if "out" in p.name) == []
        metacooler.binnify(tiny_metacooler, out, 200, force=True)
        validate(out, full=True)

    def test_a_held_lock_refuses_before_the_work(
        self,
        tmp_path: Path,
        tiny_metacooler: Path,
        request: pytest.FixtureRequest,
        monkeypatch: pytest.MonkeyPatch,
    ):
        """``create pairs`` refuses an output whose lock is held before it bins
        a single pair, not once every pair is binned and sorted."""
        from test_staged_edit import _from_pairs

        import metacooler.io.create as create_module

        out = tmp_path / "out" / "out.metacool"
        out.parent.mkdir()
        _write_lock(out, 4321, host="node12", run=DEAD_RUN)
        binned = []
        real = create_module._buffer_pixels

        def counted(*args, **kwargs):
            binned.append(1)
            return real(*args, **kwargs)

        monkeypatch.setattr(create_module, "_buffer_pixels", counted)
        with pytest.raises(MapLockedError, match="this host cannot check"):
            _from_pairs(out, request)
        assert binned == []
        assert not out.exists()

    def test_a_lock_taken_meanwhile_stops_the_rename(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Checked again right before the final rename: an output whose lock
        another thread took while it was written is not replaced."""
        import metacooler.io.create as create_module

        out = tmp_path / "out.metacool"
        taken, release = threading.Event(), threading.Event()

        def hold() -> None:
            with lock(out, follow_symlinks=False):
                taken.set()
                release.wait(60)

        thread = threading.Thread(target=hold, name="holder")
        real = create_module._move_into_place

        def taken_meanwhile(*args, **kwargs):
            raise AssertionError("renamed over an output whose lock is held")

        real_check = staging.check_output_lock
        checks: list = []

        def check(path):
            checks.append(path)
            if Path(path) == out and len(checks) > 1 and not thread.is_alive():
                thread.start()  # the write is done: another thread locks it
                assert taken.wait(60)
            return real_check(path)

        monkeypatch.setattr(staging, "check_output_lock", check)
        monkeypatch.setattr(create_module, "_move_into_place", taken_meanwhile)
        try:
            with pytest.raises(MapLockedError, match="another thread"):
                metacooler.reassign(
                    tiny_metacooler, {"contig_D": "MAG_2"}, output=out, layout="keep"
                )
        finally:
            release.set()
            if thread.is_alive() or taken.is_set():
                thread.join(60)
        del real
        assert not out.exists()
        assert not Path(f"{out}.partial").exists()


# ---------------------------------------------------------------------------
# Another tool's lock file, and a lock whose holder is gone
# ---------------------------------------------------------------------------
class TestOtherLockFiles:
    @pytest.mark.parametrize(
        "content",
        ["", "12345\n", '{"locked": true}', "not json at all"],
        ids=["empty, as filelock leaves it", "a pid", "no pid", "text"],
    )
    def test_another_tools_lock_file_does_not_stop_an_output(
        self, tiny_metacooler: Path, tmp_path: Path, content: str
    ):
        out = tmp_path / "out.metacool"
        _output_lock(out).write_text(content)
        metacooler.binnify(tiny_metacooler, out, 200)
        validate(out, full=True)
        assert _output_lock(out).read_text() == content  # left as it is

    def test_it_still_stops_an_edit_in_place(self, tiny_metacooler: Path):
        """In-place writers and edit sessions refuse a lock file that says
        nothing, as they did."""
        lock_path(tiny_metacooler).write_text("")
        with pytest.raises(MapLockedError, match="does not say which process"):
            _store(tiny_metacooler)

    def test_a_lock_whose_holder_is_gone(self, tiny_metacooler: Path, tmp_path: Path):
        """A lock of this host whose holder is gone does not stop an output,
        nor is it taken over: the next lock() does that."""
        out = tmp_path / "out.metacool"
        dead = _write_lock(out, 2**22 + 12345)
        metacooler.binnify(tiny_metacooler, out, 200)
        validate(out, full=True)
        assert staging._read(_output_lock(out)) == dead

    def test_what_cannot_be_listed_is_not_named(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        def unreadable(path):
            raise PermissionError(errno.EACCES, "Permission denied", str(path))

        monkeypatch.setattr(staging.os, "listdir", unreadable)
        assert staging._staging_of(tmp_path / "out.metacool", {"uuid": DEAD_RUN}) == ""

    def test_a_lock_of_another_host_names_what_to_delete(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "out.metacool"
        _write_lock(out, 4321, host="node12", run=DEAD_RUN)
        with pytest.raises(MapLockedError) as caught:
            metacooler.binnify(tiny_metacooler, out, 200)
        message = str(caught.value)
        assert f"delete {_output_lock(out)}." in message  # nothing else is there
        assert ".*.d/" not in message and ".partial" not in message


# ---------------------------------------------------------------------------
# A symbolic link at an output is replaced, never written through
# ---------------------------------------------------------------------------
#: The writers that write their output in place, with cooler or h5py, which
#: follow a link: (the output's name, a call writing it).
LINKED_WRITERS = {
    "extract": (
        "out.cool",
        lambda src, out: metacooler.extract(src, genomes=["MAG_1"], output=out),
    ),
    "merge": ("out.metacool", lambda src, out: _merge_two(src, out)),
    "flatten": ("out.cool", lambda src, out: metacooler.flatten(src, out)),
    "aggregate": (
        "out.cool",
        lambda src, out: metacooler.aggregate(src, out, level="genome"),
    ),
    "to_cooler": ("out.cool", lambda src, out: MetaCooler(src).to_cooler(out)),
    "a new collection": (
        "out.mmetacool",
        lambda src, out: create_mmetacool(out, {"a": src}, force=True),
    ),
    "binnify": (
        "out.metacool",
        lambda src, out: metacooler.binnify(src, out, 200, force=True),
    ),
    "extract --split": ("out_MAG_1.cool", lambda src, out: _split(src, out)),
    "reassign -o, renamed into place": (
        "out.metacool",
        lambda src, out: metacooler.reassign(
            src, {"contig_D": "MAG_2"}, output=out, layout="keep", force=True
        ),
    ),
}


def _merge_two(src: Path, out: Path) -> None:
    replicate = out.parent.parent / "replicate.metacool"
    shutil.copy(src, replicate)
    metacooler.merge([src, replicate], out)


def _split(src: Path, out: Path) -> None:
    metacooler.extract(
        src, genomes=["MAG_1"], output=out.parent / "out", split=True, force=True
    )


@pytest.mark.parametrize("dangling", [False, True], ids=["to a file", "dangling"])
@pytest.mark.parametrize("writer", list(LINKED_WRITERS))
def test_an_output_link_is_replaced(
    writer: str, dangling: bool, tiny_metacooler: Path, tmp_path: Path
):
    name, write = LINKED_WRITERS[writer]
    out = tmp_path / "out" / name
    out.parent.mkdir()
    target = tmp_path / "elsewhere" / "target.cool"
    target.parent.mkdir()
    if not dangling:
        target.write_text("another file, which the link points to")
    out.symlink_to(target)
    write(tiny_metacooler, out)
    assert not out.is_symlink()
    with h5py.File(out, "r"):
        pass
    if dangling:
        assert not target.exists()
    else:
        assert target.read_text() == "another file, which the link points to"


class TestOutputLinks:
    def test_without_force_a_link_is_an_existing_output(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        gone = tmp_path / "gone.metacool"
        link = tmp_path / "link.metacool"
        link.symlink_to(gone)  # dangling: still an existing output
        with pytest.raises(FileExistsError, match="symbolic link"):
            metacooler.binnify(tiny_metacooler, link, 200)
        with pytest.raises(FileExistsError, match="already exists"):
            create_mmetacool(link, {"a": tiny_metacooler})
        split = tmp_path / "split_MAG_1.cool"
        split.symlink_to(gone)
        with pytest.raises(FileExistsError, match="symbolic link"):
            metacooler.extract(
                tiny_metacooler,
                genomes=["MAG_1"],
                output=tmp_path / "split",
                split=True,
            )
        assert link.is_symlink() and split.is_symlink() and not gone.exists()

    def test_a_failed_run_leaves_a_link_it_would_have_replaced(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A writer that renames its output into place replaces a link there
        with its rename: a run that fails leaves the link as it was, as it
        leaves a previous output."""
        import metacooler.io.create as create_module
        from metacooler.io._utils import check_output_path

        target = tmp_path / "dated.metacool"
        shutil.copy(tiny_metacooler, target)
        link = tmp_path / "latest.metacool"
        link.symlink_to(target)
        check_output_path(link, True, defer_removal=True)
        assert link.is_symlink()
        with pytest.raises(FileExistsError, match="symbolic link"):
            check_output_path(link, False, defer_removal=True)

        def fails(*args, **kwargs):
            raise RuntimeError("the run failed")

        monkeypatch.setattr(create_module, "_move_into_place", fails)
        with pytest.raises(RuntimeError, match="the run failed"):
            metacooler.reassign(
                tiny_metacooler, {"contig_D": "MAG_2"}, output=link, force=True
            )
        assert link.is_symlink() and link.resolve() == target
        assert target.read_bytes() == tiny_metacooler.read_bytes()

    def test_a_link_to_a_directory_is_a_directory(self, tmp_path: Path):
        from metacooler.io._utils import check_output_path

        directory = tmp_path / "results"
        directory.mkdir()
        link = tmp_path / "out.metacool"
        link.symlink_to(directory)
        for force in (False, True):
            with pytest.raises(IsADirectoryError, match="is a directory"):
                check_output_path(link, force)
        with pytest.raises(IsADirectoryError, match="is a directory"):
            staging.remove_an_output_link(link, force=True)
        assert link.is_symlink() and directory.is_dir()

    @pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs os.mkfifo")
    def test_a_special_file_is_never_removed(self, tmp_path: Path):
        """``-o /dev/stdout --force`` removed the system's /dev/stdout (a link
        to a descriptor): a device, a FIFO, or a link to one, is never
        removed, with or without force."""
        from metacooler.io._utils import check_output_path

        fifo = tmp_path / "fifo"
        os.mkfifo(fifo)
        link = tmp_path / "stdout"
        link.symlink_to(fifo)
        for path in (fifo, link):
            for force in (False, True):
                with pytest.raises(FileExistsError, match="regular file"):
                    check_output_path(path, force)
        assert link.is_symlink() and fifo.exists()

    @pytest.mark.parametrize("cmd", [["flatten"], ["binnify", "--binsize", "200"]])
    def test_the_cli_with_a_dangling_link(
        self, tiny_metacooler: Path, tmp_path: Path, cmd: list[str]
    ):
        link = tmp_path / "out.cool"
        target = tmp_path / "elsewhere" / "target.cool"
        target.parent.mkdir()
        link.symlink_to(target)
        with pytest.raises(_CliFailed, match="already exists"):
            _cli(*cmd, tiny_metacooler, "-o", link)
        assert link.is_symlink() and not target.exists()
        _cli(*cmd, tiny_metacooler, "-o", link, "--force", "-q")
        assert not link.is_symlink() and not target.exists()
        with h5py.File(link, "r"):
            pass

    def test_a_link_to_the_map_of_a_session(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """An output through a link onto the map an edit session of this
        thread edits replaces the link; the map, and the session, are left as
        they are."""
        other = tmp_path / "other" / "tiny.metacool"
        other.parent.mkdir()
        shutil.copy(tiny_metacooler, other)
        link = tmp_path / "latest.metacool"
        link.symlink_to(tiny_metacooler)
        before = tiny_metacooler.read_bytes()
        with staged_edit(tiny_metacooler) as staged:
            _store(staged)
            metacooler.extract(other, genomes=["MAG_2"], output=link)
            assert tiny_metacooler.read_bytes() == before
        assert not link.is_symlink()
        assert "w" in MetaCooler(tiny_metacooler).bins()[:].columns
        validate(tiny_metacooler, full=True)

    @pytest.mark.usefixtures("holder")
    def test_a_link_to_a_map_another_process_holds(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """merge through a link onto a map whose lock another process holds
        replaces the link: the locked map is not written."""
        link = tmp_path / "latest.metacool"
        link.symlink_to(tiny_metacooler)
        before = tiny_metacooler.read_bytes()
        first = tmp_path / "first" / "replicate.metacool"
        second = tmp_path / "second" / "replicate.metacool"
        for replicate in (first, second):
            replicate.parent.mkdir()
            shutil.copy(tiny_metacooler, replicate)
        metacooler.merge([first, second], link)
        assert tiny_metacooler.read_bytes() == before
        assert not link.is_symlink()
        validate(link, full=True)


# ---------------------------------------------------------------------------
# A new collection: created exclusively
# ---------------------------------------------------------------------------
class TestNewCollections:
    def test_cp_into_one_created_meanwhile(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Another cp creates the collection, an entry copied into it, between
        this cp's check and its creation: that entry is kept."""
        src = tmp_path / "src.mmetacool"
        create_mmetacool(src, {"s1": tiny_metacooler, "s2": tiny_metacooler})
        dst = tmp_path / "dst.mmetacool"
        real = staging.writing_output

        def another_cp_meanwhile(path, **kwargs):
            if not dst.exists():
                create_mmetacool(dst, {"first": tiny_metacooler})
            return real(path, **kwargs)

        monkeypatch.setattr(staging, "writing_output", another_cp_meanwhile)
        _cli("cp", f"{src}::entries/s2", f"{dst}::entries/second", "-q")
        assert sorted(list_entries(dst)) == ["first", "second"]

    def test_create_mmetacool_never_truncates_without_force(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        out = tmp_path / "new.mmetacool"
        import metacooler.io.mmetacool as mmetacool_module

        real = mmetacool_module.writing_output
        calls: list = []

        def created_meanwhile(path, **kwargs):
            calls.append(path)
            if len(calls) == 1:  # another writer, between the check and h5py
                create_mmetacool(out, {"first": tiny_metacooler})
            return real(path, **kwargs)

        monkeypatch.setattr(mmetacool_module, "writing_output", created_meanwhile)
        with pytest.raises(FileExistsError):
            create_mmetacool(out, {"second": tiny_metacooler})
        assert list_entries(out) == ["first"]

    def test_cp_into_a_new_collection(self, tiny_metacooler: Path, tmp_path: Path):
        src = tmp_path / "src.mmetacool"
        create_mmetacool(src, {"s1": tiny_metacooler})
        dst = tmp_path / "new.mmetacool"
        _cli("cp", f"{src}::entries/s1", f"{dst}::entries/x", "-q")
        assert list_entries(dst) == ["x"]

    def test_force_onto_a_writable_file_in_a_read_only_directory(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """No lock file is created beside an output: one that exists and may be
        written is written, whatever its directory allows."""
        out = tmp_path / "shared" / "results.mmetacool"
        out.parent.mkdir()
        out.write_text("an earlier result")
        real = os.open

        def no_new_file(path, flags, *args, **kwargs):
            if (
                str(path).startswith(str(out.parent))
                and flags & os.O_CREAT
                and not (os.path.exists(path))
            ):
                raise PermissionError(errno.EACCES, "Permission denied", str(path))
            return real(path, flags, *args, **kwargs)

        monkeypatch.setattr(os, "open", no_new_file)
        monkeypatch.setattr(staging.os, "open", no_new_file)
        create_mmetacool(out, {"a": tiny_metacooler}, force=True)
        assert list_entries(out) == ["a"]
        assert not _output_lock(out).exists()


# ---------------------------------------------------------------------------
# Nits
# ---------------------------------------------------------------------------
class TestNits:
    def test_a_refused_binnify_fails_no_session(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        """binnify's output check runs before its write begins: a refusal its
        caller catches leaves the session as it was."""
        n = MetaCooler(tiny_metacooler).shape[0]
        with staged_edit(tiny_metacooler) as staged:
            with pytest.raises(FileExistsError):
                metacooler.binnify(tiny_metacooler, staged, 200)
            metacooler.store_weights(
                staged, "w", np.ones(n), method="test", depends_on="pixels"
            )
        assert "w" in MetaCooler(tiny_metacooler).bins()[:].columns

    def test_the_cli_prints_no_errno(self, tiny_metacooler: Path, tmp_path: Path):
        with pytest.raises(_CliFailed) as caught:
            _cli("flatten", tiny_metacooler, "-o", tmp_path / "missing" / "out.cool")
        message = str(caught.value)
        assert "does not exist" in message and "[Errno" not in message

    def test_a_lock_in_a_directory_that_may_not_be_written(
        self, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        def refused(path: Path, record: dict) -> bool:
            raise PermissionError(errno.EACCES, "Permission denied", str(path))

        monkeypatch.setattr(staging, "_create", refused)
        for error in (PermissionError, MapLockError):
            with pytest.raises(error), lock(tiny_metacooler):
                pass

    def test_a_lock_in_a_directory_that_does_not_exist(self, tmp_path: Path):
        for error in (FileNotFoundError, MapLockError):
            with pytest.raises(error), lock(tmp_path / "gone" / "map.metacool"):
                pass

    def test_a_new_file_is_never_replaced(self, tmp_path: Path):
        staged, target = tmp_path / "staged", tmp_path / "target"
        staged.write_text("new")
        target.write_text("another program's")
        with pytest.raises(FileExistsError, match="created meanwhile"):
            staging._move_no_replace(staged, target)
        assert target.read_text() == "another program's"
        assert staged.read_text() == "new"
        target.unlink()
        staging._move_no_replace(staged, target)
        assert target.read_text() == "new" and not staged.exists()

    def test_another_error_of_the_link_is_raised(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        def broken(*args, **kwargs):
            raise OSError(errno.EIO, "Input/output error")

        monkeypatch.setattr(staging.os, "link", broken)
        staged = tmp_path / "staged"
        staged.write_text("new")
        with pytest.raises(OSError, match="Input/output error"):
            staging._move_no_replace(staged, tmp_path / "target")
        assert staged.exists() and not (tmp_path / "target").exists()

    def test_where_there_are_no_hard_links(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        def no_links(*args, **kwargs):
            raise OSError(errno.EPERM, "Operation not permitted")

        monkeypatch.setattr(staging.os, "link", no_links)
        staged, target = tmp_path / "staged", tmp_path / "target"
        staged.write_text("new")
        staging._move_no_replace(staged, target)
        assert target.read_text() == "new"
        staged.write_text("newer")
        with pytest.raises(FileExistsError):
            staging._move_no_replace(staged, target)
        assert target.read_text() == "new"

    def test_an_output_sessions_copy_is_set_back(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Where file times are coarse, a raw write as soon as an output
        session's copy is made is still seen as a write: what a weight depends
        on may have changed."""

        def coarse(path: Path):
            try:
                info = os.stat(path)
            except FileNotFoundError:
                return None
            return info.st_ino, info.st_mtime_ns // 10_000_000_000

        monkeypatch.setattr(staging, "_written_state", coarse)
        out = tmp_path / "out.metacool"
        source = MetaCooler(tiny_metacooler)
        n = source.shape[0]
        with staged_edit(tiny_metacooler, output=out) as staged:
            with h5py.File(staged, "r+") as f:
                f.attrs["note"] = "raw"
            with pytest.raises(staging.MapChangedError, match="changed pixels"):
                metacooler.store_weights(
                    staged,
                    "w",
                    np.ones(n),
                    method="t",
                    depends_on="pixels",
                    source=source,
                )
        import time

        assert abs(out.stat().st_mtime - time.time()) < 5  # its install's time


@pytest.mark.usefixtures("holder")
def test_a_lock_another_process_holds_refuses_an_output_there(
    tiny_metacooler: Path, tmp_path: Path
):
    """A map whose lock another process holds (metaTOR holds it for a run) is
    no output for any writer, as it is no map for an in-place one."""
    src = tmp_path / "other" / "tiny.metacool"
    src.parent.mkdir()
    shutil.copy(tiny_metacooler, src)
    before = tiny_metacooler.read_bytes()
    with pytest.raises(MapLockedError, match="test holder"):
        metacooler.binnify(src, tiny_metacooler, 200, force=True)
    with pytest.raises(MapLockedError, match="test holder"):
        metacooler.flatten(src, tiny_metacooler)
    assert tiny_metacooler.read_bytes() == before
