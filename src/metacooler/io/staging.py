"""Locked, staged in-place writes of a map, and the file state reads check.

A writer that edits a map in place -- ``annotate(inplace=True)``,
``reassign(inplace=True)``, :func:`~metacooler.io.weights.store_weights` (and
``normalize``), :func:`~metacooler.io.clusters.write_clusters`,
:func:`~metacooler.io.clusters.drop_clusters`, ``qc --store``, the ``metadata``
edits and the ``.mmetacool`` entry commands -- never writes the map itself.
Under the map's lock (:func:`lock`), taken before its first read, it writes a
*staged* copy in a private directory beside the map (:func:`staged_inplace`):
a fresh copy of the map, edited, validated in full, given the map's mode and
owner, and only then renamed over the map.  What that guarantees:

- **The map is whole, always.**  On any failure -- an error, a validation that
  fails, ``KeyboardInterrupt``, ``SystemExit`` -- the staging directory is
  removed and the map stays as it was.  A process killed outright (``kill
  -9``, or a scheduler's ``SIGTERM`` at a time limit: metacooler installs no
  signal handler) leaves its lock and its staging directory, both named after
  its lock's ``uuid``: the next writer that takes the stale lock removes them.
- **The map does not grow.**  HDF5 never reclaims the space of what an edit
  replaces, so a map edited in place grew with every edit.  The staged map is a
  fresh copy -- every object copied as stored -- so it carries at most what its
  own edit replaced.
- **One writer at a time.**  The lock, ``<map>.lock``, names its holder.  A
  live holder is refused by name; the lock of a process of this host that is
  gone is taken over, under a guard file, so that one process exactly takes
  it; a lock of another host is never broken.  Only its holder removes it.  It
  is re-entrant for the thread that holds it, so a caller that holds a map's
  lock for a whole run (metaTOR does) calls the writers, which take it unless
  their caller holds it; another thread is refused.
- **What was read is what is replaced.**  A writer that read the map before
  it took the lock (to compute weights, say) hands its
  :class:`metacooler.MetaCooler` to :func:`staged_inplace`, which refuses to
  copy a file that changed since it was read.
- **Private stays private.**  The staging directory is ``0700``; the staged
  map gets the map's mode, and its owner where the process may set it.  A
  symbolic link is replaced, not its target.  A map without write permission
  is refused, as editing it in place was.

A read is refused when the map is replaced or written while it is read
(:class:`MapChangedError`, raised by :class:`metacooler.MetaCooler`, which
checks the file before and after every read): metacooler reads the pixels a
chunk at a time, so a read that spanned a replacement would mix the bins of
one map with the pixels of the other.
"""

from __future__ import annotations

import contextlib
import dataclasses
import errno
import functools
import json
import logging
import os
import shutil
import socket
import stat
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, ParamSpec, TypeVar

import h5py

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterable, Sequence

    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

#: Suffix of a map's lock file, next to the map (next to a symbolic link, not
#: its target: the link is what an in-place write replaces).
LOCK_SUFFIX = ".lock"
#: Suffix of the guard a process holds while it takes over a stale lock.
TAKEOVER_SUFFIX = ".takeover"

#: How long a lock or a guard file may be seen empty (its creator is writing
#: it) before it is refused as one that says nothing.
_EMPTY_WAIT = 2.0
#: How long a process waits on another's takeover of a stale lock.
_TAKEOVER_WAIT = 10.0
_POLL = 0.005


class MapLockError(OSError):
    """The lock of a map cannot be taken.

    Attributes:
        lock_file: The lock file (or the takeover guard it names).
        holder: What the lock file says of its holder (``pid``, ``host``,
            ``start``, ``command``, ``uuid``, ...): empty when it cannot be
            read, or when the lock could not be created at all.
    """

    def __init__(
        self, message: str, *, lock_file: str | Path, holder: dict | None = None
    ) -> None:
        super().__init__(message)
        self.lock_file = Path(lock_file)
        self.holder = dict(holder or {})

    def __str__(self) -> str:
        return str(self.args[0])

    def __reduce__(self):
        return (
            _rebuild_lock_error,
            (type(self), str(self), str(self.lock_file), self.holder),
        )


class MapLockedError(MapLockError):
    """Another process, or another thread, holds the lock of a map; or the lock
    file does not say who does."""


class MapChangedError(OSError):
    """A map was replaced or written while it was read, or since it was read
    for an edit."""

    def __str__(self) -> str:
        return str(self.args[0])


def _rebuild_lock_error(
    cls: type[MapLockError], message: str, lock_file: str, holder: dict
) -> MapLockError:
    """Unpickle a :class:`MapLockError` (its constructor takes keywords)."""
    return cls(message, lock_file=lock_file, holder=holder)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
def map_file(uri: str | Path) -> Path:
    """The file of a map, as given: the path itself, or a URI's file part."""
    from metacooler.uri import parse_uri

    file_path, _ = parse_uri(str(uri))
    return Path(file_path)


def _canonical(path: str | Path) -> Path:
    """*path* absolute, its directory resolved (symbolic links followed) and its
    own name kept: a link is what an in-place write replaces, and two names of
    one directory share one lock."""
    absolute = Path(os.path.abspath(path))
    return Path(os.path.realpath(absolute.parent)) / absolute.name


def lock_path(uri: str | Path) -> Path:
    """The lock file of the map *uri*: ``<file>.lock``, next to its file.

    A URI locks its file: every entry of a ``.mmetacool`` shares the lock of the
    collection, and ``file::/`` is the file itself.  The path is absolute, its
    directory resolved, so the lock is the same whatever the working directory
    or the name of the directory the map is reached through.
    """
    path = _canonical(map_file(uri))
    return path.with_name(path.name + LOCK_SUFFIX)


def map_output(output: str | Path) -> Path:
    """The file an output path names: *output*, or the file of ``file::/``.

    An output is a whole file: a URI naming a group inside a file is refused,
    as writing it would make a file named after the URI.  The group part is
    read as :func:`metacooler.uri.parse_uri` reads it, spaces stripped.

    Raises:
        ValueError: If *output* names a group inside a file.
    """
    text = str(output)
    if "::" not in text:
        return Path(output)
    file_text, _, group = text.partition("::")
    if group.strip().strip("/"):
        raise ValueError(
            f"{text} names the group {group.strip()!r} inside {file_text}, but an "
            "output is a whole file: give a path (to edit an entry of a "
            "collection, edit it in place)."
        )
    return Path(file_text.strip())


# ---------------------------------------------------------------------------
# Who holds a lock
# ---------------------------------------------------------------------------
def _boot_id() -> str | None:
    """This host's boot id (Linux), which a reboot changes."""
    with contextlib.suppress(OSError):
        return Path("/proc/sys/kernel/random/boot_id").read_text().strip() or None
    return None


def _pid_start(pid: int) -> int | None:
    """When the process *pid* started, in clock ticks since boot (Linux), which
    tells a process from a later one that reuses its pid."""
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
    except (OSError, ValueError):
        return None
    # The command name, in parentheses, may hold spaces: count from its end.
    fields = text[text.rfind(")") + 2 :].split()
    try:
        return int(fields[19])
    except (IndexError, ValueError):
        return None


def _pid_is_running(pid: int) -> bool:
    """Whether the process *pid* runs on this host.

    Never call on Windows: there ``os.kill(pid, 0)`` terminates the process.
    """
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Another user's process, which this one may not signal: alive.
        return True
    except (OverflowError, ValueError):
        return False
    return True


@dataclasses.dataclass
class _Held:
    """A lock this process holds."""

    holder: dict
    pid: int
    thread: int
    thread_name: str
    depth: int = 1
    #: Edits staged under the lock so far, which number the next one.
    staged: int = 0


#: The locks this process holds, by lock file (absolute, canonical).
_HELD: dict[str, _Held] = {}
#: Taken while a thread acquires, enters or releases a lock: so that the
#: threads of one process see each other's locks, and a lock file naming this
#: process is either one of _HELD or left by an earlier process of its pid.
_REGISTRY = threading.RLock()


def _holding(lock_file: Path) -> _Held | None:
    """This process's hold of *lock_file*, if it has one."""
    held = _HELD.get(str(lock_file))
    # A forked child inherits the table, not the lock.
    return held if held is not None and held.pid == os.getpid() else None


def holds_lock(uri: str | Path) -> bool:
    """Whether this thread holds the lock of the map *uri*."""
    held = _holding(lock_path(uri))
    return held is not None and held.thread == threading.get_ident()


def _is_stale(holder: dict) -> bool:
    """Whether a lock's holder is a process of this host that is gone.

    A holder of another host is never stale: its process cannot be checked.
    Nor is any holder on Windows, where checking a process kills it.  On this
    host: a lock of an earlier boot; one naming this process's pid but none of
    its locks (an earlier process of that pid); a pid that runs no process;
    and a pid that a later process reuses (its start time differs) are stale.
    """
    pid = holder.get("pid")
    if (
        os.name == "nt"
        or not isinstance(pid, int)
        or isinstance(pid, bool)
        or holder.get("host") != socket.gethostname()
    ):
        return False
    boot = holder.get("boot_id")
    if boot and _boot_id() not in (None, boot):
        return True
    if pid == os.getpid():
        return not any(
            held.holder.get("uuid") == holder.get("uuid") for held in _HELD.values()
        )
    if not _pid_is_running(pid):
        return True
    start, now = holder.get("pid_start"), _pid_start(pid)
    return start is not None and now is not None and start != now


def _describe(holder: dict) -> str:
    """A lock's holder, as messages name it."""
    command = holder.get("command")
    return (
        f"pid {holder.get('pid')} on {holder.get('host')}, since "
        f"{holder.get('start')}" + (f": {command}" if command else "")
    )


def _refusal(path: Path, lock_file: Path, holder: dict) -> str:
    """Why the lock of *path* is not taken, and what to do."""
    if not holder:
        return (
            f"{lock_file} does not say which process holds {path}: if none is "
            f"writing it, delete {lock_file}."
        )
    if holder.get("host") != socket.gethostname():
        return (
            f"{path} is locked by a process of another host ({_describe(holder)}), "
            f"which this host cannot check: if it is no longer running, delete "
            f"{lock_file}."
        )
    return (
        f"{path} is being written by another process ({_describe(holder)}): wait "
        f"for it to end, or, if no such process is running, delete {lock_file}."
    )


def _record(command: str | None) -> dict:
    """What a lock (or a takeover guard) says of this process."""
    pid = os.getpid()
    return {
        "pid": pid,
        "host": socket.gethostname(),
        "boot_id": _boot_id(),
        "pid_start": _pid_start(pid),
        "start": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command": command if command is not None else _default_command(),
        "uuid": uuid.uuid4().hex,
    }


def _default_command() -> str:
    """The command line of this process, its program by name."""
    argv = list(sys.argv) or ["python"]
    return " ".join([Path(argv[0]).name, *argv[1:]])


def _create(path: Path, record: dict) -> bool:
    """Create *path* holding *record*; ``False`` when it exists.

    Created exclusively (``O_EXCL``), written and synced; a file that cannot
    be written is removed, never left empty.  Any other failure (a directory
    that is missing, or not writable) is raised as it is.
    """
    try:
        handle = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    try:
        with os.fdopen(handle, "w") as stream:
            json.dump(record, stream)
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(path)
        raise
    return True


_GONE = object()


def _read(path: Path) -> dict | object:
    """What *path* holds: a dict (``{}`` when it is empty or holds no JSON
    object), or :data:`_GONE` when there is no file."""
    try:
        with open(path) as stream:
            text = stream.read()
    except FileNotFoundError:
        return _GONE
    except OSError:
        return {}
    try:
        found = json.loads(text)
    except ValueError:
        return {}
    return found if isinstance(found, dict) else {}


def _read_settled(path: Path) -> dict | object:
    """What *path* holds, waiting a moment for a file seen empty, which its
    creator may be writing."""
    deadline = time.monotonic() + _EMPTY_WAIT
    while True:
        found = _read(path)
        if found is _GONE or found or time.monotonic() > deadline:
            return found
        time.sleep(_POLL)


def _staged_prefix(path: Path, run: str) -> str:
    """The start of the name of everything a holder of uuid *run* stages
    beside the map *path*."""
    return f".{path.name}.{run}."


def _remove_named(directory: Path, prefix: str) -> list[str]:
    """Remove what *directory* holds under names that start with *prefix*:
    what a writer staged, and the temporary files cooler wrote beside it."""
    try:
        names = os.listdir(directory)
    except OSError:
        return []
    removed = []
    for name in names:
        if not name.startswith(prefix):
            continue
        target = directory / name
        try:
            if target.is_dir() and not target.is_symlink():
                shutil.rmtree(target)
            else:
                os.remove(target)
        except OSError:
            # What cannot be removed is left: the error being raised, if any,
            # is the one to report.
            continue
        removed.append(name)
    return removed


def _take_over(lock_file: Path, path: Path, stale: dict) -> None:
    """Remove a stale lock (*stale*, what it says), under the takeover guard.

    One process at a time holds the guard, ``<map>.lock.takeover``, created
    exclusively.  Under it, the lock is read again and removed only if it is
    still the stale one: a lock another process has meanwhile created is never
    touched.  The process then contends for the lock as every other does, by
    an exclusive create, so one exactly holds it.  What the dead holder staged
    is removed with its lock.  A guard held by a live process is waited for;
    one whose holder is gone, or that says nothing, is refused by name: only
    a person knows that no process is taking the lock over.
    """
    guard = lock_file.with_name(lock_file.name + TAKEOVER_SUFFIX)
    mine = _record(None)
    deadline = time.monotonic() + _TAKEOVER_WAIT
    while not _create(guard, mine):
        found = _read_settled(guard)
        if found is _GONE:
            continue
        assert isinstance(found, dict)
        if not found or _is_stale(found):
            raise MapLockedError(
                f"{guard} was left by a process that stopped while it took over "
                f"the stale lock of {path}"
                + (f" ({_describe(found)})" if found else "")
                + f": if no process is writing {path}, delete {guard}.",
                lock_file=guard,
                holder=found,
            )
        if time.monotonic() > deadline:
            raise MapLockedError(
                f"{path}: another process ({_describe(found)}) has been taking "
                f"over its stale lock for {_TAKEOVER_WAIT:.0f} s; try again, or, "
                f"if no process is, delete {guard}.",
                lock_file=guard,
                holder=found,
            )
        time.sleep(_POLL)
    try:
        current = _read(lock_file)
        if current == stale:
            os.unlink(lock_file)
            logger.warning(
                "Took over %s: the process that held it (%s) is gone.",
                lock_file,
                _describe(stale),
            )
            run = stale.get("uuid")
            if isinstance(run, str) and run.isalnum():
                for name in _remove_named(lock_file.parent, _staged_prefix(path, run)):
                    logger.warning(
                        "Removed %s, which that process left beside %s.", name, path
                    )
    finally:
        if os.getpid() == mine["pid"]:
            with contextlib.suppress(FileNotFoundError):
                os.unlink(guard)


def _acquire(lock_file: Path, path: Path, mine: dict) -> None:
    """Create the lock file for *mine*, taking over a stale one.

    Raises:
        MapLockedError: If a live process holds the lock, the lock of another
            host or one that says nothing is in the way, or a takeover cannot
            be completed.
    """
    while True:
        if _create(lock_file, mine):
            return
        found = _read_settled(lock_file)
        if found is _GONE:
            continue
        assert isinstance(found, dict)
        if not found or not _is_stale(found):
            raise MapLockedError(
                _refusal(path, lock_file, found), lock_file=lock_file, holder=found
            )
        _take_over(lock_file, path, found)


@contextlib.contextmanager
def lock(uri: str | Path, *, command: str | None = None) -> Generator[Path, None, None]:
    """Hold the lock of a map for the time of the ``with`` block.

    The lock file, ``<map>.lock`` (see :func:`lock_path`), is created next to
    the map, exclusively, with what says who holds it: ``pid``, ``host``,
    ``boot_id`` and ``pid_start`` (to tell a process from a later one of the
    same pid), ``start``, ``command``, and a ``uuid`` of its own, which the
    release compares, so that only its holder removes it.  Take it before the
    first read of an edit, so that what is read is what is written.

    It is re-entrant for the thread that holds it: a nested ``with lock(...)``
    of the same map holds it on, and only the outermost one releases it.  Every
    in-place writer takes it (see :func:`staged_inplace`), so a caller that
    holds it for a whole run keeps other writers out from its first read to its
    last write, and calls the writers meanwhile.  Another thread of the process
    is refused, as another process is.

    A lock whose holder is a process of this host that is gone (killed with
    ``SIGKILL``, or ``SIGTERM``: metacooler installs no signal handler) is taken
    over, with a warning, by one process exactly, and what that process staged
    beside the map is removed.  A lock of another host is never broken: the
    error names the file to delete.

    Args:
        uri: A map, as a path or a URI: a URI locks its file.
        command: The command the lock names (default: this process's command
            line).

    Yields:
        The lock file.

    Raises:
        MapLockedError: If another process or thread holds the lock, or the
            lock file does not say which (it holds no JSON object).
        MapLockError: If the lock file cannot be created (the directory of the
            map is missing, or not writable).
    """
    path = map_file(uri)
    lock_file = lock_path(uri)
    key = str(lock_file)
    with _REGISTRY:
        held = _holding(lock_file)
        if held is not None:
            if held.thread != threading.get_ident():
                raise MapLockedError(
                    f"{path} is being written by another thread of this process "
                    f"({held.thread_name}): a map's lock is held by one thread at "
                    "a time.",
                    lock_file=lock_file,
                    holder=held.holder,
                )
            held.depth += 1
        else:
            mine = _record(command)
            try:
                _acquire(lock_file, path, mine)
            except MapLockedError:
                raise
            except OSError as error:
                raise MapLockError(
                    f"{lock_file} cannot be created ({error.strerror or error}): "
                    f"an in-place write of {path} locks it there, so the "
                    "directory of the map must exist and be writable. Write a "
                    "copy of the map instead.",
                    lock_file=lock_file,
                ) from error
            current = threading.current_thread()
            held = _HELD[key] = _Held(
                holder=mine,
                pid=os.getpid(),
                thread=threading.get_ident(),
                thread_name=current.name,
            )
    try:
        yield lock_file
    finally:
        with _REGISTRY:
            held.depth -= 1
            # A forked child neither releases its parent's lock nor removes
            # its file; nor does a nested block.
            if held.depth == 0 and os.getpid() == held.pid:
                _HELD.pop(key, None)
                found = _read(lock_file)
                # This process's lock only: not one that another process took
                # after it was deleted by hand.
                if isinstance(found, dict) and found.get("uuid") == held.holder["uuid"]:
                    with contextlib.suppress(FileNotFoundError):
                        os.unlink(lock_file)


P = ParamSpec("P")
R = TypeVar("R")


def locked_when_inplace(writer: Callable[P, R]) -> Callable[P, R]:
    """Run a writer under the lock of its map (its first argument) when it is
    called with ``inplace=True``, from before its first read."""

    @functools.wraps(writer)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
        if not kwargs.get("inplace"):
            return writer(*args, **kwargs)
        uri = args[0] if args else kwargs["uri"]
        with lock(map_file(str(uri))):
            return writer(*args, **kwargs)

    return wrapped


# ---------------------------------------------------------------------------
# Staged in-place writes
# ---------------------------------------------------------------------------
def _umask() -> int:
    """The umask of the process.

    Read from ``/proc/self/status`` where Linux has it, which changes nothing;
    elsewhere it can only be read by setting it, which is not thread-safe.
    """
    with contextlib.suppress(OSError, ValueError):
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Umask:"):
                return int(line.split()[1], 8)
    mask = os.umask(0o022)
    os.umask(mask)
    return mask


def copy_fresh(source: str | Path, target: str | Path) -> None:
    """Copy the HDF5 file *source* into a new file *target*, object by object.

    Every object of the root group is copied as stored -- dtype, chunks,
    filters, attributes, the columns of a table in their order -- soft links
    (``/chroms``) as links, then the root's attributes with their dtypes.  The
    copy holds none of the space that edits left in *source*.  *target* is
    opened with ``"w"``: a file already there, created private, keeps its mode.
    """
    with h5py.File(source, "r") as src, h5py.File(target, "w") as dst:
        for key in src:
            link = src.get(key, getlink=True)
            if isinstance(link, h5py.SoftLink):
                dst[key] = h5py.SoftLink(link.path)
            else:
                src.copy(src[key], dst, name=key)
        for key in src.attrs:
            dst.attrs.create(key, src.attrs[key], dtype=src.attrs.get_id(key).dtype)


def check_unchanged(planned: MetaCooler | Iterable[MetaCooler] | None) -> None:
    """Refuse an edit of a map that changed since it was read for it.

    Raises:
        MapChangedError: If the file of a map in *planned* was replaced or
            written since that map was opened.
    """
    from metacooler.api import MetaCooler

    if planned is None:
        return
    maps = [planned] if isinstance(planned, MetaCooler) else list(planned)
    for mc in maps:
        mc._check_unchanged()


def refuse_read_only(path: Path) -> None:
    """Refuse to replace a map that has no write permission.

    Replacing it would succeed in a writable directory, where editing it in
    place failed; a map made read-only is not to be edited.  Judged by its mode
    bits too, since a process run as root may write any file.

    Raises:
        PermissionError: If the map is not writable.
    """
    try:
        mode = os.stat(path).st_mode
    except FileNotFoundError:
        return
    if not (mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)) or not os.access(
        path, os.W_OK
    ):
        raise PermissionError(
            errno.EACCES,
            f"{path} is read-only (mode {stat.filemode(mode)}): an in-place edit "
            "would replace it. Make it writable (chmod u+w) to edit it, or write "
            "a copy",
        )


def _fsync(path: Path) -> None:
    """Flush *path* (a file, or a directory's entries) to disk, where the
    system allows it."""
    flags = (
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) if path.is_dir() else os.O_RDONLY
    )
    try:
        handle = os.open(path, flags)
    except OSError:
        return
    try:
        with contextlib.suppress(OSError):
            os.fsync(handle)
    finally:
        os.close(handle)


def _owned_like(path: Path, staged: Path) -> None:
    """Give *staged* the owner and group of *path*, where this process may."""
    with contextlib.suppress(OSError, AttributeError):
        info = os.stat(path)
        os.chown(staged, info.st_uid, info.st_gid)


@contextlib.contextmanager
def staged_inplace(
    uri: str | Path,
    *,
    copy: bool = True,
    validate_groups: Sequence[str] | None = None,
    planned: MetaCooler | Iterable[MetaCooler] | None = None,
) -> Generator[Path, None, None]:
    """Edit the file of a map through a staged copy, which then replaces it.

    Under the map's lock (taken unless this thread holds it), the maps in
    *planned* -- those the edit was planned on -- are checked unchanged, and a
    map without write permission is refused.  The staged map is then made in a
    private directory beside the map, ``.<map>.<uuid>.<n>.d/`` (mode ``0700``,
    named after the lock's uuid): by default a fresh copy of the map's file
    (:func:`copy_fresh`), to edit in place; with ``copy=False``, a path where
    nothing is yet, for a writer that writes a whole new map (cooler's
    temporary files land in that directory too).  Once the ``with`` block ends,
    the staged map is validated in full, given the map's mode and owner (a new
    file's mode, from the umask, where there is no map yet), synced, and
    renamed over the map's path: a symbolic link there is replaced, not its
    target.  On any failure, ``KeyboardInterrupt`` and ``SystemExit`` included,
    the staging directory is removed and the map is left as it was -- by the
    process that made it only, not by a forked child.

    Args:
        uri: The map: a ``.metacool`` path, an ``.mmetacool`` path or one of
            its entries (the whole file is staged either way).
        copy: Start from a copy of the map (``False``: from nothing).
        validate_groups: The groups of the staged file to validate in full
            (default: the entry *uri* names; every entry, with
            :func:`metacooler.schema.validate_mmetacooler`, for a bare
            ``.mmetacool``; the map otherwise).
        planned: The map(s) the edit was planned on, as opened for it: the
            edit is refused if the file changed since.

    Yields:
        The staged map's path.

    Raises:
        MapChangedError: If a map of *planned* changed since it was opened.
        PermissionError: If the map has no write permission.
        MapLockedError: If another process or thread holds the lock.
    """
    from metacooler.schema import is_mmetacooler, validate, validate_mmetacooler
    from metacooler.uri import parse_uri

    file_text, group = parse_uri(str(uri))
    # Absolute, so that a chdir in the block changes nothing.
    path = Path(os.path.abspath(file_text))
    with lock(path):
        held = _holding(lock_path(path))
        assert held is not None  # taken just above
        check_unchanged(planned)
        refuse_read_only(path)
        held.staged += 1
        prefix = f"{_staged_prefix(path, held.holder['uuid'])}{held.staged}."
        workdir = path.with_name(prefix + "d")
        os.mkdir(workdir, 0o700)
        os.chmod(workdir, 0o700)
        staged = workdir / path.name
        try:
            if copy:
                copy_fresh(path, staged)
            yield staged
            if os.getpid() != held.pid:
                # A forked child that left the block: the edit is its parent's.
                return
            if validate_groups is not None:
                for name in validate_groups:
                    validate(f"{staged}::{name}", full=True)
            elif group.strip("/"):
                validate(f"{staged}::{group}", full=True)
            elif is_mmetacooler(staged):
                validate_mmetacooler(staged)
            else:
                validate(staged, full=True)
            if path.exists():
                shutil.copymode(path, staged)
                _owned_like(path, staged)
            else:
                os.chmod(staged, 0o666 & ~_umask())
            _fsync(staged)
            os.replace(staged, path)
            _fsync(path.parent)
        except BaseException:
            if os.getpid() == held.pid:
                shutil.rmtree(workdir, ignore_errors=True)
            raise
        # What cooler left in the directory (its temporary files).
        shutil.rmtree(workdir, ignore_errors=True)


def edit_metadata(uri: str | Path, change: Callable[[dict], None]) -> None:
    """Change the ``metadata`` attribute of a map, in place.

    Under the map's lock, through a staged copy (see :func:`staged_inplace`):
    *change* gets the parsed metadata and edits it; an exception it raises
    leaves the map as it was.

    Args:
        uri: A ``.metacool`` path or an ``.mmetacool`` entry URI.
        change: Edits the metadata dict in place.
    """
    from metacooler.uri import parse_uri

    _, group = parse_uri(str(uri))
    key = group.strip("/")
    with staged_inplace(uri) as staged, h5py.File(staged, "r+") as f:
        grp = f[key] if key else f
        raw = grp.attrs.get("metadata", "{}")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        metadata = json.loads(raw) if isinstance(raw, str) else dict(raw)
        change(metadata)
        grp.attrs["metadata"] = json.dumps(metadata)


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
def file_state(path: str | Path) -> tuple[int, int, int, int] | None:
    """What changes when the file at *path* is replaced (another inode) or
    written (its size, its time of modification): ``None`` when there is no
    file there."""
    try:
        info = os.stat(path)
    except FileNotFoundError:
        return None
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns
