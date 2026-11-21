"""Locked, staged in-place writes of a map, edit sessions, and the file state
reads check.

A writer that edits a map in place -- ``annotate(inplace=True)``,
``reassign(inplace=True)``, :func:`~metacooler.io.weights.store_weights` (and
``normalize``), :func:`~metacooler.io.clusters.write_clusters`,
:func:`~metacooler.io.clusters.drop_clusters`, ``qc --store``, the ``metadata``
edits and the ``.mmetacool`` entry functions -- never writes the map itself.
Under the map's lock (:func:`lock`), taken before its first read, it writes a
*staged* copy in a private directory beside the map (:func:`staged_inplace`):
a fresh copy of the map, edited, validated in full, given the map's mode,
owner, group and extended attributes, and only then renamed over the map.
:func:`staged_edit` opens one such copy for several edits: every writer given
its path writes into it directly, and it is validated and renamed once, over
the map or onto an output path.  What that guarantees:

- **The map is whole, always.**  On any failure -- an error, a validation that
  fails, ``KeyboardInterrupt``, ``SystemExit`` -- the staging directory is
  removed and the map stays as it was.  A process killed outright (``kill
  -9``, or a scheduler's ``SIGTERM`` at a time limit: metacooler installs no
  signal handler) leaves its lock and its staging directory, both named after
  its lock's ``uuid``.  Whoever takes the map's lock next -- by a takeover, or
  once the lock was deleted by hand -- removes every staging directory of the
  map that is not its own, naming each: whoever left it no longer holds the
  lock.
- **The map does not grow.**  HDF5 never reclaims the space of what an edit
  replaces, so a map edited in place grew with every edit.  The staged map is a
  fresh copy -- every object copied as stored -- so it carries at most what its
  own edit replaced.
- **One writer at a time.**  The lock, ``<map>.lock``, names its holder.  A
  live holder is refused by name; the lock of a process of this host that is
  gone is taken over, under a guard file, so that one process exactly takes
  it; a lock of another host -- or of another machine or PID namespace under
  the same host name -- is never broken.  Only its holder removes it.  It is
  re-entrant for the thread that holds it, so a caller that holds a map's lock
  for a whole run (metaTOR does) calls the writers, which take it unless their
  caller holds it; another thread is refused.
- **What was read is what is replaced.**  A writer that read the map before
  it took the lock (to compute weights, say) hands its
  :class:`metacooler.MetaCooler` to :func:`staged_inplace`, which refuses to
  copy a file that changed since it was read.  The map is checked once more
  just before the rename: a program that writes it without metacooler's lock
  (``cooler balance``, say) while its copy is edited makes the edit fail
  (:class:`MapChangedError`) rather than lose that change -- and a copy torn by
  such a write is refused alike.  Only a write in the instant between that
  check and the rename goes unseen.
- **Private stays private.**  The staging directory is ``0700`` (keeping the
  set-group-id bit it inherits, so that the staged map takes its directory's
  group); the staged map gets the map's owner and group where the process may
  set them -- a warning names what it may not -- its mode, and its extended
  attributes, POSIX ACLs included.  A map without write permission is
  refused, as editing it in place was.
- **Symbolic links.**  An in-place edit through a symbolic link edits its
  target: the lock, the staging directory and the replacement are beside the
  target, and the link keeps pointing at the edited map.  An output path that
  is a link is replaced, never written through.

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
import re
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

#: Suffix of a map's lock file, next to the map (next to the target of a
#: symbolic link: the target is what an in-place edit replaces).
LOCK_SUFFIX = ".lock"
#: Suffix of the guard a process holds while it takes over a stale lock.
TAKEOVER_SUFFIX = ".takeover"

#: How long a lock or a guard file may be seen empty (its creator is writing
#: it) before it is refused as one that says nothing.
_EMPTY_WAIT = 2.0
#: How long a process waits on another's takeover of a stale lock.
_TAKEOVER_WAIT = 10.0
_POLL = 0.005
#: A lock's uuid, as it names what its holder stages.
_UUID = re.compile(r"[0-9a-f]{32}")


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


def _resolved(path: str | Path) -> Path:
    """*path* absolute, every symbolic link followed, its own name's too: the
    file an in-place edit of *path* edits."""
    return Path(os.path.realpath(os.path.abspath(path)))


def lock_path(uri: str | Path) -> Path:
    """The lock file of the map *uri*: ``<file>.lock``, next to its file.

    A URI locks its file: every entry of a ``.mmetacool`` shares the lock of the
    collection, and ``file::/`` is the file itself.  The path is absolute and
    resolved -- a symbolic link to the map is followed, as is any in its
    directory -- so the lock is the same whatever the working directory, the
    name of the directory, or the link the map is reached through.
    """
    path = _resolved(map_file(uri))
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


def _machine_id() -> str | None:
    """This machine's id (systemd's, or D-Bus's), where it has one: two
    machines -- or containers -- may share a host name."""
    for name in ("/etc/machine-id", "/var/lib/dbus/machine-id"):
        with contextlib.suppress(OSError):
            text = Path(name).read_text().strip()
            if text:
                return text
    return None


def _pid_namespace() -> str | None:
    """The PID namespace of this process (Linux): two processes see each
    other's pids only within one, whatever host name they share."""
    with contextlib.suppress(OSError, AttributeError):
        return os.readlink("/proc/self/ns/pid")
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


#: The locks this process holds, by lock file (absolute, resolved).
_HELD: dict[str, _Held] = {}
#: Taken while a thread acquires, enters or releases a lock: so that the
#: threads of one process see each other's locks, and a lock file naming this
#: process is either one of _HELD or left by an earlier process of its pid.  A
#: thread that waits (on a lock file being written, on another's takeover)
#: releases it meanwhile (:func:`_pause`).
_REGISTRY = threading.RLock()
_WAKE = threading.Condition(_REGISTRY)


def _after_fork_in_child() -> None:
    """Give a forked child a free registry: a thread of its parent that held it
    does not exist in the child, and would hold it for ever."""
    global _REGISTRY, _WAKE
    _REGISTRY = threading.RLock()
    _WAKE = threading.Condition(_REGISTRY)


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_fork_in_child)


def _pause(seconds: float) -> None:
    """Wait *seconds*, releasing the registry if this thread holds it, so that
    other threads lock other maps meanwhile."""
    try:
        _WAKE.wait(seconds)
    except RuntimeError:  # the registry is not held
        time.sleep(seconds)


def _holding(lock_file: Path) -> _Held | None:
    """This process's hold of *lock_file*, if it has one."""
    held = _HELD.get(str(lock_file))
    # A forked child inherits the table, not the lock.
    return held if held is not None and held.pid == os.getpid() else None


def holds_lock(uri: str | Path) -> bool:
    """Whether this thread holds the lock of the map *uri* (or writes it through
    an edit session: see :func:`staged_edit`)."""
    session = _session_of(map_file(uri))
    if session is not None:
        return session.thread == threading.get_ident()
    held = _holding(lock_path(uri))
    return held is not None and held.thread == threading.get_ident()


def _same_host(holder: dict) -> bool:
    """Whether a lock's holder runs where this process can check it: the same
    host name, machine id and PID namespace (one that a lock does not record
    differs from one this process has)."""
    return (
        holder.get("host") == socket.gethostname()
        and holder.get("machine_id") == _machine_id()
        and holder.get("pid_ns") == _pid_namespace()
    )


def _is_stale(holder: dict) -> bool:
    """Whether a lock's holder is a process of this host that is gone.

    A holder this process cannot check is never stale: one of another host,
    machine or PID namespace (see :func:`_same_host`), or any holder on
    Windows, where checking a process kills it.  On this host: a lock of an
    earlier boot; one naming this process's pid but none of its locks (an
    earlier process of that pid); a pid that runs no process; and a pid that a
    later process reuses (its start time differs) are stale.
    """
    pid = holder.get("pid")
    if (
        os.name == "nt"
        or not isinstance(pid, int)
        or isinstance(pid, bool)
        or not _same_host(holder)
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


def _staging_of(path: Path, holder: dict) -> str:
    """What the holder of a lock stages beside the map *path*, as a pattern to
    delete, or ``""`` when its lock names no uuid."""
    run = holder.get("uuid")
    if isinstance(run, str) and _UUID.fullmatch(run):
        return f" and {path.parent / _staged_prefix(path, run)}*.d/ (its staged copies)"
    return ""


def _refusal(path: Path, lock_file: Path, holder: dict) -> str:
    """Why the lock of *path* is not taken, and what to do."""
    if not holder:
        return (
            f"{lock_file} does not say which process holds {path}: if none is "
            f"writing it, delete {lock_file}."
        )
    if not _same_host(holder):
        return (
            f"{path} is locked by a process this host cannot check "
            f"({_describe(holder)}: another host, machine or PID namespace): if "
            f"it is no longer running, delete {lock_file}"
            f"{_staging_of(path, holder)}."
        )
    return (
        f"{path} is being written by another process ({_describe(holder)}): wait "
        f"for it to end, or, if no such process is running, delete {lock_file}"
        f"{_staging_of(path, holder)}."
    )


def _record(command: str | None) -> dict:
    """What a lock (or a takeover guard) says of this process."""
    pid = os.getpid()
    return {
        "pid": pid,
        "host": socket.gethostname(),
        "machine_id": _machine_id(),
        "pid_ns": _pid_namespace(),
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
        _pause(_POLL)


def _staged_prefix(path: Path, run: str) -> str:
    """The start of the name of everything a holder of uuid *run* stages
    beside the map *path*."""
    return f".{path.name}.{run}."


def _remove(target: Path) -> None:
    """Remove a file, or a directory and what it holds."""
    if target.is_dir() and not target.is_symlink():
        shutil.rmtree(target)
    else:
        os.remove(target)


def _remove_named(directory: Path, prefix: str) -> list[str]:
    """Remove what *directory* holds under names that start with *prefix*,
    returning their names; what cannot be removed is named in a warning."""
    try:
        names = os.listdir(directory)
    except OSError:
        return []
    removed = []
    for name in names:
        if not name.startswith(prefix):
            continue
        try:
            _remove(directory / name)
        except OSError as error:
            logger.warning(
                "Could not remove %s (%s): delete it by hand.",
                directory / name,
                error.strerror or error,
            )
            continue
        removed.append(name)
    return removed


def _remove_foreign_staging(path: Path, mine: str) -> None:
    """Remove what other holders of the lock of the map *path* staged beside it
    -- every ``.<map>.<uuid>.<n>...`` of another uuid than *mine* --, naming
    each: whoever left it no longer holds the lock (it was killed, and its lock
    taken over or deleted by hand)."""
    pattern = re.compile(rf"\.{re.escape(path.name)}\.([0-9a-f]{{32}})\.\d+\.")
    try:
        names = os.listdir(path.parent)
    except OSError:
        return
    for name in sorted(names):
        found = pattern.match(name)
        if found is None or found.group(1) == mine:
            continue
        target = path.parent / name
        try:
            _remove(target)
        except OSError as error:
            logger.warning(
                "Could not remove %s, left beside %s by a writer that no longer "
                "holds its lock (%s): delete it by hand.",
                target,
                path,
                error.strerror or error,
            )
            continue
        logger.warning(
            "Removed %s, left beside %s by a writer that no longer holds its lock.",
            name,
            path,
        )


def _take_over(lock_file: Path, path: Path, stale: dict) -> None:
    """Remove a stale lock (*stale*, what it says), under the takeover guard.

    One process at a time holds the guard, ``<map>.lock.takeover``, created
    exclusively.  Under it, the lock is read again and removed only if it is
    still the stale one: a lock another process has meanwhile created is never
    touched.  The process then contends for the lock as every other does, by
    an exclusive create, so one exactly holds it (and removes what the dead
    holder staged: :func:`lock`).  A guard held by a live process is waited
    for; one whose holder is gone, or that says nothing, is refused by name:
    only a person knows that no process is taking the lock over.
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
        _pause(_POLL)
    try:
        current = _read(lock_file)
        if current == stale:
            os.unlink(lock_file)
            logger.warning(
                "Took over %s: the process that held it (%s) is gone.",
                lock_file,
                _describe(stale),
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
        if found.get("pid") == os.getpid():
            for held in _HELD.values():
                if held.holder.get("uuid") == found.get("uuid"):
                    raise MapLockedError(
                        f"{path} is being written by another thread of this "
                        f"process ({held.thread_name}): a map's lock is held by "
                        "one thread at a time.",
                        lock_file=lock_file,
                        holder=found,
                    )
        if not found or not _is_stale(found):
            raise MapLockedError(
                _refusal(path, lock_file, found), lock_file=lock_file, holder=found
            )
        _take_over(lock_file, path, found)


@contextlib.contextmanager
def lock(uri: str | Path, *, command: str | None = None) -> Generator[Path, None, None]:
    """Hold the lock of a map for the time of the ``with`` block.

    The lock file, ``<map>.lock`` (see :func:`lock_path`), is created next to
    the map (next to the target of a symbolic link), exclusively, with what
    says who holds it: ``pid``, ``host``, ``machine_id`` and ``pid_ns`` (where
    the process runs), ``boot_id`` and ``pid_start`` (to tell a process from a
    later one of the same pid), ``start``, ``command``, and a ``uuid`` of its
    own, which the release compares, so that only its holder removes it.  Take
    it before the first read of an edit, so that what is read is what is
    written.  Once it is taken, what other holders staged beside the map is
    removed (each named in a warning): they no longer hold it.

    It is re-entrant for the thread that holds it: a nested ``with lock(...)``
    of the same map holds it on, and only the outermost one releases it.  Every
    in-place writer takes it (see :func:`staged_inplace`), so a caller that
    holds it for a whole run keeps other writers out from its first read to its
    last write, and calls the writers meanwhile.  Another thread of the process
    is refused, as another process is.  The staged map of an edit session
    (:func:`staged_edit`) is held by its session: the session's thread takes
    its "lock" at no cost, and any other thread is refused.

    A lock whose holder is a process of this host that is gone (killed with
    ``SIGKILL``, or ``SIGTERM``: metacooler installs no signal handler) is taken
    over, with a warning, by one process exactly.  A lock of another host,
    machine or PID namespace is never broken: the error names the files to
    delete.

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
    session = _session_of(map_file(uri))
    if session is not None:
        _check_session_thread(session, map_file(uri))
        yield session.lock_file
        return
    path = _resolved(map_file(uri))
    lock_file = path.with_name(path.name + LOCK_SUFFIX)
    key = str(lock_file)
    fresh = False
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
            fresh = True
    try:
        if fresh:
            _remove_foreign_staging(path, held.holder["uuid"])
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
# Edit sessions: one staged copy for several writes
# ---------------------------------------------------------------------------
@dataclasses.dataclass
class _Session:
    """A staged map that the writers of one thread write directly."""

    #: The staged file (resolved).
    staged: Path
    #: Where a writer writes a whole new map before it takes the staged one's
    #: place (a re-grouping).
    workdir: Path
    #: The lock file that keeps the session's map ours.
    lock_file: Path
    pid: int
    thread: int
    thread_name: str
    #: The map file a session edits in place or copies out (resolved), or
    #: ``None``.
    source: Path | None = None
    #: Whether the staged map replaces *source* (``False``: an output).
    in_place: bool = False
    #: The group of the map the session was opened on (``""``: the file).
    group: str = ""
    #: The groups the writers wrote (``""``: the whole file), validated when
    #: the session ends.
    touched: set[str] = dataclasses.field(default_factory=set)
    rewrites: int = 0


#: The edit sessions of this process, by staged file (resolved).
_SESSIONS: dict[str, _Session] = {}


def _session_of(path: str | Path) -> _Session | None:
    """The session whose staged file *path* is, if this process has one."""
    session = _SESSIONS.get(os.path.realpath(os.path.abspath(path)))
    return session if session is not None and session.pid == os.getpid() else None


def _session_on(source: Path) -> _Session | None:
    """This thread's session of the map file *source* (resolved)."""
    me = threading.get_ident()
    for session in list(_SESSIONS.values()):
        if (
            session.source == source
            and session.thread == me
            and session.pid == os.getpid()
        ):
            return session
    return None


def in_session(uri: str | Path) -> bool:
    """Whether *uri* is the staged map of an edit session of this thread, which
    its writers write directly (and validate once, when the session ends)."""
    session = _session_of(map_file(uri))
    return session is not None and session.thread == threading.get_ident()


def _check_session_thread(session: _Session, path: str | Path) -> None:
    """Refuse a write of a session's staged map from another thread."""
    if session.thread != threading.get_ident():
        raise MapLockedError(
            f"{path} is the staged map of an edit session of another thread of "
            f"this process ({session.thread_name}): only that thread writes it.",
            lock_file=session.lock_file,
        )


@contextlib.contextmanager
def _registered(session: _Session) -> Generator[_Session, None, None]:
    """Let *session*'s thread write its staged map directly while the block
    runs."""
    key = os.path.realpath(session.staged)
    with _REGISTRY:
        _SESSIONS[key] = session
    try:
        yield session
    finally:
        with _REGISTRY:
            if _SESSIONS.get(key) is session:
                del _SESSIONS[key]


@contextlib.contextmanager
def _writing_directly(path: str | Path) -> Generator[Path, None, None]:
    """Let this thread's writers write *path* directly -- unstaged, unlocked,
    unvalidated -- while the block runs: for a new file that only its writer
    knows of, which stages it as a whole (``zoomify`` and ``stack`` build their
    output so, and validate it once)."""
    staged = _resolved(path)
    session = _Session(
        staged=staged,
        workdir=staged.parent,
        lock_file=staged.with_name(staged.name + LOCK_SUFFIX),
        pid=os.getpid(),
        thread=threading.get_ident(),
        thread_name=threading.current_thread().name,
    )
    with _registered(session):
        yield Path(path)


@contextlib.contextmanager
def _in_session(
    session: _Session,
    path: str | Path,
    group: str,
    *,
    copy: bool,
    validate_groups: Sequence[str] | None,
    planned: MetaCooler | Iterable[MetaCooler] | None,
) -> Generator[Path, None, None]:
    """What :func:`staged_inplace` yields to a writer of a session's staged map:
    the map itself, or, to a writer of a whole new map (*copy* false), a path
    beside it whose file then takes its place."""
    _check_session_thread(session, path)
    check_unchanged(planned)
    if validate_groups is None:
        session.touched.add(group.strip("/"))
    else:
        session.touched.update(name.strip("/") for name in validate_groups)
    if copy:
        yield session.staged
        return
    session.rewrites += 1
    fresh = session.workdir / f".rewrite.{session.rewrites}.{session.staged.name}"
    try:
        yield fresh
        if os.getpid() == session.pid:
            os.replace(fresh, session.staged)
    except BaseException:
        if os.getpid() == session.pid:
            with contextlib.suppress(FileNotFoundError):
                os.unlink(fresh)
        raise


def _refuse_the_map_of_a_session(path: str | Path) -> None:
    """Refuse a write of a map that this thread edits in a session: the
    session replaces it when it ends, and would lose that write."""
    session = _session_on(_resolved(path))
    if session is not None and session.in_place:
        staged = _session_uri(session.staged, session.group)
        raise ValueError(
            f"{path} is being edited in an edit session of this thread: write its "
            f"staged map, {staged}, which replaces the map when the session "
            "ends."
        )


def _session_uri(staged: Path, group: str) -> str:
    """The staged map a session yields: its file, or a URI into it."""
    key = group.strip("/")
    return f"{staged}::{key}" if key else str(staged)


def _validate_session(session: _Session) -> None:
    """Validate in full what a session's writers wrote."""
    from metacooler.schema import is_mmetacooler, validate, validate_mmetacooler

    staged = session.staged
    if not is_mmetacooler(staged):
        validate(staged, full=True)
        return
    groups = session.touched or {session.group}
    if "" in groups:
        validate_mmetacooler(staged)
        return
    with h5py.File(staged, "r") as f:
        present = sorted(g for g in groups if g in f)
    for name in present:
        validate(f"{staged}::{name}", full=True)


@contextlib.contextmanager
def staged_edit(
    uri: str | Path,
    *,
    output: str | Path | None = None,
    force: bool = False,
    command: str | None = None,
    planned: MetaCooler | Iterable[MetaCooler] | None = None,
) -> Generator[str, None, None]:
    """Edit a map with several writers through one staged copy, written once.

    Under the lock of the map's file (taken unless this thread holds it), the
    maps in *planned* are checked unchanged, and one private copy is made: of
    the whole file, to replace it (*output* ``None``); or, with *output*, of
    the map alone -- an entry of a ``.mmetacool`` becomes a standalone
    ``.metacool`` -- destined for that path.  The block gets the staged map's
    path (a URI into it for an entry edited in place).  Every metacooler writer
    given that path by this thread -- ``annotate(inplace=True)``,
    ``reassign(inplace=True)``, ``store_weights``, ``write_clusters``,
    ``drop_clusters``, ``contact_qc(store=True)``, the ``normalize_*``
    functions, the entry functions -- writes into it directly: no copy, no
    validation, no lock of its own.  When the block ends, the staged map is
    validated in full, given the map's owner, group, mode and extended
    attributes (an output gets a new file's mode), synced and renamed over the
    map's file (a symbolic link to it keeps pointing at it) or onto *output*
    (a link there is replaced, not followed).  On any failure, ``KeyboardInterrupt``
    and ``SystemExit`` included, the staged copy is removed and the map (or
    *output*) is left as it was.  As for any staged write, a map written
    meanwhile by a program that takes no metacooler lock is refused rather
    than lost (:class:`MapChangedError`).

    A session in place joins one this thread already has open on the same map
    (the inner ``with`` gets the same staged map, and the outer one writes it);
    any other nesting on one map is refused.  While a session is open, a write
    of the map itself (rather than of its staged copy) from this thread is
    refused; another thread is refused the staged map.

    Args:
        uri: The map: a ``.metacool`` path, an ``.mmetacool`` entry URI, or a
            whole ``.mmetacool`` (in place only, or copied whole to *output*).
        output: Where the edited map goes instead of over *uri*: a new file
            (refused if it exists, unless *force*), never *uri*'s own file.
        force: Replace *output* if it exists.
        command: The command the map's lock names, if this call takes it.
        planned: The map(s) the edits were planned on: refused if the file
            changed since they were opened.

    Yields:
        The staged map: a ``.metacool`` path, or a URI into a staged
        ``.mmetacool``.

    Raises:
        MapLockedError: If another process or thread holds the map's lock.
        MapChangedError: If a map of *planned* changed since it was opened,
            the map was written meanwhile by another program, or the map
            changed while it was copied to *output*.
        FileExistsError: If *output* exists and *force* is not set.
        PermissionError: In place, if the map is read-only.
        ValueError: If *output* is *uri*'s file, or for a nesting refused.
    """
    from metacooler.uri import parse_uri

    file_text, group = parse_uri(str(uri))
    key = group.strip("/")
    # Nested on one map: a session of this thread -- on its staged map, or on
    # the same map in place -- is joined in place; any other nesting refused.
    joined = _session_of(file_text)
    if joined is None:
        on_map = _session_on(_resolved(file_text))
        if on_map is not None and (output is not None or not on_map.in_place):
            raise ValueError(
                f"{uri} is being copied out, or edited, in an edit session of "
                "this thread: a nested session only joins one that edits the "
                "same map in place (output=None)."
            )
        joined = on_map
    if joined is not None:
        _check_session_thread(joined, file_text)
        if output is not None:
            raise ValueError(
                f"{uri} is the staged map of an edit session of this thread: a "
                "nested session only joins it in place (output=None)."
            )
        check_unchanged(planned)
        joined.touched.add(key)
        yield _session_uri(joined.staged, key)
        return
    source = _resolved(file_text)
    if output is not None:
        with _output_session(
            source,
            key,
            map_output(output),
            force=force,
            command=command,
            planned=planned,
        ) as staged:
            yield staged
        return
    with lock(source, command=command):
        held = _holding(source.with_name(source.name + LOCK_SUFFIX))
        assert held is not None  # taken just above
        check_unchanged(planned)
        refuse_read_only(source)
        if not source.exists():
            raise FileNotFoundError(errno.ENOENT, f"No map to edit: {source}")
        workdir = _new_workdir(source, held)
        staged = workdir / source.name
        session = _Session(
            staged=staged,
            workdir=workdir,
            lock_file=source.with_name(source.name + LOCK_SUFFIX),
            pid=os.getpid(),
            thread=threading.get_ident(),
            thread_name=threading.current_thread().name,
            source=source,
            in_place=True,
            group=key,
        )
        try:
            before = file_state(source)
            copy_fresh(source, staged)
            with _registered(session):
                yield _session_uri(staged, key)
            if os.getpid() != held.pid:
                return
            _validate_session(session)
            _install(
                staged, source, before=before, like=source, what="the edit session"
            )
        except BaseException:
            if os.getpid() == held.pid:
                shutil.rmtree(workdir, ignore_errors=True)
            raise
        shutil.rmtree(workdir, ignore_errors=True)


@contextlib.contextmanager
def _output_session(
    source: Path,
    group: str,
    output: Path,
    *,
    force: bool,
    command: str | None,
    planned: MetaCooler | Iterable[MetaCooler] | None,
) -> Generator[str, None, None]:
    """:func:`staged_edit` with an output: the map copied, standalone, to a
    staging directory beside *output*, then renamed onto it."""
    out = Path(os.path.abspath(output))
    if out == source or (
        out.exists() and not out.is_symlink() and os.path.samefile(out, source)
    ):
        raise ValueError(
            f"output={output} is the file of the map edited: leave output out to "
            "edit the map in place."
        )
    if (out.exists() or out.is_symlink()) and not force:
        raise FileExistsError(
            errno.EEXIST, f"{out} exists: pass force=True to replace it."
        )
    with lock(source, command=command):
        held = _holding(source.with_name(source.name + LOCK_SUFFIX))
        assert held is not None  # taken just above
        check_unchanged(planned)
        workdir = _new_workdir(out, held)
        staged = workdir / out.name
        session = _Session(
            staged=staged,
            workdir=workdir,
            lock_file=source.with_name(source.name + LOCK_SUFFIX),
            pid=os.getpid(),
            thread=threading.get_ident(),
            thread_name=threading.current_thread().name,
            source=source,
        )
        try:
            before = file_state(source)
            copy_fresh(source, staged, group=group or None)
            if file_state(source) != before:
                raise MapChangedError(
                    f"{source} was written by another program while it was "
                    "copied: the copy is not a map. Try again."
                )
            with _registered(session):
                yield str(staged)
            if os.getpid() != held.pid:
                return
            _validate_session(session)
            if not force and (out.exists() or out.is_symlink()):
                raise FileExistsError(
                    errno.EEXIST,
                    f"{out} was created meanwhile: pass force=True to replace it.",
                )
            _install(staged, out, before=None, like=None, what="the edit session")
        except BaseException:
            if os.getpid() == held.pid:
                shutil.rmtree(workdir, ignore_errors=True)
            raise
        shutil.rmtree(workdir, ignore_errors=True)


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


def copy_fresh(
    source: str | Path, target: str | Path, *, group: str | None = None
) -> None:
    """Copy the HDF5 file *source* into a new file *target*, object by object.

    Every object of the root group -- or of *group*, which becomes the root of
    *target*: an entry of a ``.mmetacool`` copied out as a standalone map -- is
    copied as stored (dtype, chunks, filters, attributes, the columns of a
    table in their order), soft links (``/chroms``) as links (one into *group*
    re-pointed at the same object of the new root), then the attributes with
    their dtypes.  The copy holds none of the space that edits left in
    *source*.  *target* is opened with ``"w"``: a file already there, created
    private, keeps its mode.
    """
    key = (group or "").strip("/")
    inside = f"/{key}/" if key else "/"
    with h5py.File(source, "r") as src, h5py.File(target, "w") as dst:
        root = src[key] if key else src
        for name in root:
            link = root.get(name, getlink=True)
            if isinstance(link, h5py.SoftLink):
                path = link.path
                if key and path.startswith(inside):
                    path = "/" + path[len(inside) :]
                dst[name] = h5py.SoftLink(path)
            else:
                src.copy(root[name], dst, name=name)
        for name in root.attrs:
            dst.attrs.create(
                name, root.attrs[name], dtype=root.attrs.get_id(name).dtype
            )


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
    """Give *staged* the owner and group of *path*, or its group alone where
    this process may not give its owner (a member of the map's group editing
    another user's map), naming in a warning what it cannot keep."""
    if not hasattr(os, "chown"):  # pragma: no cover - Windows
        return
    want = os.stat(path)
    with contextlib.suppress(OSError):
        os.chown(staged, want.st_uid, want.st_gid)
        return
    with contextlib.suppress(OSError):
        os.chown(staged, -1, want.st_gid)
    have = os.stat(staged)
    if have.st_gid != want.st_gid:
        logger.warning(
            "%s: its group (gid %d) cannot be kept by this process; the edited "
            "map's group is gid %d.",
            path,
            want.st_gid,
            have.st_gid,
        )
    if have.st_uid != want.st_uid:
        logger.warning(
            "%s: its owner (uid %d) cannot be kept by this process; the edited "
            "map is owned by uid %d.",
            path,
            want.st_uid,
            have.st_uid,
        )


def _copy_xattrs(path: Path, staged: Path) -> None:
    """Give *staged* the extended attributes of *path* (POSIX ACLs included,
    as ``system.posix_acl_*``), where the system has them, naming in a warning
    those it cannot set.  The SELinux label is left to the policy."""
    if not hasattr(os, "listxattr"):  # pragma: no cover - not Linux
        return
    try:
        names = os.listxattr(path)
    except OSError:
        return
    lost = []
    for name in names:
        if name == "security.selinux":
            continue
        try:
            os.setxattr(staged, name, os.getxattr(path, name))
        except OSError:
            lost.append(name)
    if lost:
        logger.warning(
            "%s: its extended attribute(s) %s cannot be kept by this process.",
            path,
            ", ".join(lost),
        )


def _new_workdir(beside: Path, held: _Held) -> Path:
    """Make the private staging directory of an edit of *beside*, under the
    lock *held*: ``.<name>.<uuid>.<n>.d/``, ``0700``, keeping the
    set-group-id bit its directory gives it (so that what it holds takes the
    directory's group)."""
    held.staged += 1
    workdir = beside.with_name(
        f"{_staged_prefix(beside, held.holder['uuid'])}{held.staged}.d"
    )
    os.mkdir(workdir, 0o700)
    inherited = os.stat(beside.parent).st_mode & stat.S_ISGID
    os.chmod(workdir, 0o700 | inherited)
    return workdir


def _install(
    staged: Path,
    target: Path,
    *,
    before: tuple[int, int, int, int] | None,
    like: Path | None,
    what: str,
) -> None:
    """Rename *staged* over *target*: given what *like* has besides its bytes
    (a new file's mode without one), synced, once *target* is checked to be
    as it was (*before*, where *like* is the map it replaces).

    Raises:
        MapChangedError: If *target* was written or replaced meanwhile.
    """
    if like is not None and like.exists():
        _owned_like(like, staged)
        shutil.copymode(like, staged)
        _copy_xattrs(like, staged)
    else:
        os.chmod(staged, 0o666 & ~_umask())
    _fsync(staged)
    if like is not None and file_state(target) != before:
        raise MapChangedError(
            f"{target} was written by another program while {what} was staged "
            "(one that takes no metacooler lock): writing the staged map would "
            "lose that change, so it is not written. Run the edit again."
        )
    os.replace(staged, target)
    _fsync(target.parent)


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
    private directory beside the map (beside the target of a symbolic link),
    ``.<map>.<uuid>.<n>.d/`` (mode ``0700``, named after the lock's uuid): by
    default a fresh copy of the map's file (:func:`copy_fresh`), to edit in
    place; with ``copy=False``, a path where nothing is yet, for a writer that
    writes a whole new map (cooler's temporary files land in that directory
    too).  Once the ``with`` block ends, the staged map is validated in full,
    given the map's owner, group, mode and extended attributes (a new file's
    mode, from the umask, where there is no map yet), synced, checked to
    replace the map as it was when the copy began, and renamed over it.  On any
    failure, ``KeyboardInterrupt`` and ``SystemExit`` included, the staging
    directory is removed and the map is left as it was -- by the process that
    made it only, not by a forked child.

    The staged map of an edit session (:func:`staged_edit`) of this thread is
    written directly: no copy, no validation (the session validates), no lock.

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
        MapChangedError: If a map of *planned* changed since it was opened, or
            the map was written meanwhile by a program that takes no lock.
        PermissionError: If the map has no write permission.
        MapLockedError: If another process or thread holds the lock.
        ValueError: If this thread edits the map in an edit session (write its
            staged map instead).
    """
    from metacooler.schema import is_mmetacooler, validate, validate_mmetacooler
    from metacooler.uri import parse_uri

    file_text, group = parse_uri(str(uri))
    session = _session_of(file_text)
    if session is not None:
        with _in_session(
            session,
            file_text,
            group,
            copy=copy,
            validate_groups=validate_groups,
            planned=planned,
        ) as staged:
            yield staged
        return
    _refuse_the_map_of_a_session(file_text)
    # An in-place edit through a symbolic link edits its target.
    path = _resolved(file_text)
    with lock(path):
        held = _holding(lock_path(path))
        assert held is not None  # taken just above
        check_unchanged(planned)
        refuse_read_only(path)
        workdir = _new_workdir(path, held)
        staged = workdir / path.name
        try:
            before = file_state(path)
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
            _install(staged, path, before=before, like=path, what="this edit")
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
def input_states(uris: Iterable[str | Path]) -> dict[str, tuple | None]:
    """The state of the file of each map in *uris*, to check after reading
    them by path (:func:`check_inputs_unchanged`)."""
    return {str(map_file(uri)): file_state(map_file(uri)) for uri in uris}


def check_inputs_unchanged(states: dict[str, tuple | None]) -> None:
    """Refuse what was made from inputs read by path (``merge``'s
    ``cooler.merge_coolers``, ``collect``'s cooler) when one was written or
    replaced while it was read: it may mix two versions of that input.

    Raises:
        MapChangedError: Naming the inputs that changed.
    """
    changed = [path for path, state in states.items() if file_state(path) != state]
    if changed:
        raise MapChangedError(
            f"{', '.join(changed)} was written or replaced while it was read: "
            "what was read may mix two versions of it, so nothing is written. "
            "Run again."
        )


def file_state(path: str | Path) -> tuple[int, int, int, int] | None:
    """What changes when the file at *path* is replaced (another inode) or
    written (its size, its time of modification): ``None`` when there is no
    file there."""
    try:
        info = os.stat(path)
    except FileNotFoundError:
        return None
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns
