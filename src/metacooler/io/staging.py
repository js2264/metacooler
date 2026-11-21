"""Locked, staged in-place writes of a map.

A writer that edits a map in place -- ``annotate(inplace=True)``,
:func:`~metacooler.io.clusters.write_clusters`,
:func:`~metacooler.io.clusters.drop_clusters`,
:func:`~metacooler.io.weights.store_weights`, ``reassign(inplace=True)`` --
never writes the map itself.  Under the map's lock (:func:`lock`) it writes a
*staged* sibling file (:func:`staged_inplace`): a fresh copy of the map, edited,
validated in full, given the map's mode, and only then renamed over the map.
What that guarantees:

- **The map is whole, always.**  On any failure -- an error, a validation that
  fails, ``KeyboardInterrupt``, ``SystemExit`` -- the staged file and cooler's
  partial files beside it are removed and the map stays as it was.  A process
  killed outright (``kill -9``) leaves its lock and its staged file, both named
  after its lock's ``uuid``: the next writer that takes the stale lock removes
  them.  metacooler installs no signal handler.
- **The map does not grow.**  HDF5 never reclaims the space of what an edit
  replaces, so a map edited in place grew with every edit (129 kB to 272 kB
  over four updates of a small map).  The staged file is a fresh copy of the
  map -- every object copied as stored -- so it carries at most what its own
  edit replaced.
- **One writer at a time.**  The lock, ``<map>.lock``, names its holder (pid,
  host, start, command and a ``uuid``).  A live holder is refused by name; the
  lock of a process that is gone, on this host, is taken over (with a
  warning); a lock of another host is never broken, and the refusal names the
  file to delete.  Only its holder removes it.  It is re-entrant for the
  process that holds it, so a caller that holds a map's lock for a whole run
  (metaTOR does) calls the writers, which take it unless their caller holds it.
- **A symbolic link is replaced, not its target**, and the map keeps its
  mode.

A read is refused when the map is replaced or written while it is read
(:class:`MapChangedError`, raised by :class:`metacooler.MetaCooler`): metacooler
opens the file for each chunk of pixels it reads, so a read that spanned a
replacement would mix the bins of one map with the pixels of the other.
"""

from __future__ import annotations

import contextlib
import dataclasses
import json
import logging
import os
import shutil
import socket
import sys
import time
import uuid
from pathlib import Path
from typing import TYPE_CHECKING

import h5py

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence

logger = logging.getLogger("metacooler")

#: Suffix of a map's lock file, next to the map (next to a symbolic link, not
#: its target: the link is what an in-place write replaces).
LOCK_SUFFIX = ".lock"


class MapLockError(OSError):
    """The lock of a map cannot be taken.

    Attributes:
        lock_file: The lock file.
        holder: What the lock file says of its holder (``pid``, ``host``,
            ``start``, ``command``, ``uuid``): empty when it cannot be read,
            or when the lock could not be created at all.
    """

    def __init__(
        self, message: str, *, lock_file: str | Path, holder: dict | None = None
    ) -> None:
        super().__init__(message)
        self.lock_file = Path(lock_file)
        self.holder = dict(holder or {})

    def __str__(self) -> str:
        return str(self.args[0])


class MapLockedError(MapLockError):
    """Another process holds the lock of a map, or one this one cannot check."""


class MapChangedError(OSError):
    """A map was replaced or written while it was read."""

    def __str__(self) -> str:
        return str(self.args[0])


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------
@dataclasses.dataclass
class _Held:
    """A lock this process holds."""

    holder: dict
    pid: int
    depth: int = 1
    #: Staged files made under the lock so far, which number the next one.
    staged: int = 0


#: The locks this process holds, by the absolute path of their lock file.
_HELD: dict[str, _Held] = {}


def map_file(uri: str | Path) -> Path:
    """The file of a map: the path itself, or the file part of a URI."""
    from metacooler.uri import parse_uri

    file_path, _ = parse_uri(str(uri))
    return Path(file_path)


def lock_path(uri: str | Path) -> Path:
    """The lock file of the map *uri*: ``<file>.lock``, next to its file.

    A URI locks its file: every entry of a ``.mmetacool`` shares the lock of the
    collection, and ``file::/`` is the file itself.
    """
    path = map_file(uri)
    return path.with_name(path.name + LOCK_SUFFIX)


def _key(lock_file: Path) -> str:
    return os.path.abspath(lock_file)


def _holding(lock_file: Path) -> _Held | None:
    """This process's hold of *lock_file*, if it has one."""
    held = _HELD.get(_key(lock_file))
    # A forked child inherits the table, not the lock.
    return held if held is not None and held.pid == os.getpid() else None


def holds_lock(uri: str | Path) -> bool:
    """Whether this process holds the lock of the map *uri*."""
    return _holding(lock_path(uri)) is not None


def _pid_is_running(pid: int) -> bool:
    """Whether the process *pid* runs on this host."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OverflowError:
        return False
    return True


def _create_lock(lock_file: Path, holder: dict) -> bool:
    """Create the lock file, saying who holds it; ``False`` when it exists.

    Any other failure (a directory that is missing, or not writable) is raised
    as it is.
    """
    try:
        handle = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    with os.fdopen(handle, "w") as stream:
        json.dump(holder, stream)
    return True


def _lock_holder(lock_file: Path) -> dict:
    """What a lock file says of its holder: ``{}`` when it cannot be read, or
    holds no JSON object."""
    try:
        with open(lock_file) as stream:
            holder = json.load(stream)
    except (OSError, ValueError):
        return {}
    return holder if isinstance(holder, dict) else {}


def _is_stale(holder: dict) -> bool:
    """Whether a lock's holder is a process of this host that is gone."""
    pid = holder.get("pid")
    return (
        isinstance(pid, int)
        and not isinstance(pid, bool)
        and holder.get("host") == socket.gethostname()
        and not _pid_is_running(pid)
    )


def _take_stale_lock(lock_file: Path, stale: dict, holder: dict) -> bool:
    """Replace a lock whose holder is gone (*stale*, what it says) by *holder*'s,
    unless another process did first: ``True`` when *holder* holds it.

    The stale lock is moved aside first, which only one process can do: one
    that moved another's lock instead (taken in the meantime) puts it back and
    has not taken it.  One case is left open: a third process that creates the
    lock while it is moved aside holds it with the one whose lock was moved,
    which cannot be put back.  Both then write, each through its own staged and
    validated file, the last one replacing the map: it is whole either way.
    """
    aside = lock_file.with_name(f"{lock_file.name}.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        os.rename(lock_file, aside)
    except FileNotFoundError:
        # Another process moved it aside first: the lock is to take.
        return _create_lock(lock_file, holder)
    if _lock_holder(aside) != stale:
        with contextlib.suppress(FileExistsError):
            os.link(aside, lock_file)
        os.remove(aside)
        return False
    os.remove(aside)
    return _create_lock(lock_file, holder)


def _staged_prefix(path: Path, run: str) -> str:
    """The start of the name of every file a holder of uuid *run* stages
    beside the map *path*."""
    return f".{path.name}.{run}."


def _remove_named(directory: Path, prefix: str) -> list[str]:
    """Remove the files of *directory* whose name starts with *prefix*: those
    a writer staged, and the partial files cooler writes beside them."""
    try:
        names = os.listdir(directory)
    except OSError:
        return []
    removed = []
    for name in names:
        if not name.startswith(prefix):
            continue
        try:
            os.remove(directory / name)
        except OSError:
            # What cannot be removed is left: the error being raised, if
            # any, is the one to report.
            continue
        removed.append(name)
    return removed


def _describe(holder: dict) -> str:
    """A lock's holder, as the refusals name it."""
    command = holder.get("command")
    return (
        f"pid {holder.get('pid')} on {holder.get('host')}, since {holder.get('start')}"
        + (f": {command}" if command else "")
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


def _default_command() -> str:
    """The command line of this process, its program by name."""
    argv = list(sys.argv) or ["python"]
    return " ".join([Path(argv[0]).name, *argv[1:]])


@contextlib.contextmanager
def lock(uri: str | Path, *, command: str | None = None) -> Generator[Path, None, None]:
    """Hold the lock of a map for the time of the ``with`` block.

    The lock file, ``<map>.lock`` (see :func:`lock_path`), is created next to
    the map with what says who holds it: ``pid``, ``host``, ``start``,
    ``command`` and a ``uuid`` of its own, which the release compares, so that
    only its holder removes it.  It is re-entrant for the process that holds
    it: a nested ``with lock(...)`` of the same map, in the same process, holds
    it on, and only the outermost one releases it.  Every in-place writer takes
    it (see :func:`staged_inplace`), so a caller that holds it for a whole run
    keeps other writers out from its first check to its last write.

    A lock whose holder is a process of this host that is gone (killed with
    ``SIGKILL``) is taken over, with a warning, and the files that process
    staged beside the map are removed.  A lock of another host is never
    broken: the error names the file to delete.

    Args:
        uri: A map, as a path or a URI: a URI locks its file.
        command: The command the lock names (default: this process's command
            line).

    Yields:
        The lock file.

    Raises:
        MapLockedError: If another process holds the lock, or the lock file
            does not say which (it cannot be read).
        MapLockError: If the lock file cannot be created (the directory of the
            map is missing, or not writable).
    """
    path = map_file(uri)
    lock_file = lock_path(uri)
    held = _holding(lock_file)
    if held is not None:
        held.depth += 1
        try:
            yield lock_file
        finally:
            held.depth -= 1
        return

    mine = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "start": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command": command if command is not None else _default_command(),
        "uuid": uuid.uuid4().hex,
    }
    try:
        taken = _create_lock(lock_file, mine)
        if not taken:
            found = _lock_holder(lock_file)
            if _is_stale(found):
                logger.warning(
                    "Taking over %s: the process that held it (%s) is gone.",
                    lock_file,
                    _describe(found),
                )
                taken = _take_stale_lock(lock_file, found, mine)
                run = found.get("uuid")
                if taken and isinstance(run, str) and run.isalnum():
                    for name in _remove_named(
                        lock_file.parent, _staged_prefix(path, run)
                    ):
                        logger.warning(
                            "Removed %s, which that process left beside %s.",
                            name,
                            path,
                        )
    except OSError as error:
        raise MapLockError(
            f"{lock_file} cannot be created ({error.strerror or error}): an "
            f"in-place write of {path} locks it there, so the directory of the "
            "map must exist and be writable. Write a copy of the map instead.",
            lock_file=lock_file,
        ) from error
    if not taken:
        found = _lock_holder(lock_file)
        raise MapLockedError(
            _refusal(path, lock_file, found), lock_file=lock_file, holder=found
        )
    _HELD[_key(lock_file)] = _Held(holder=mine, pid=os.getpid())
    try:
        yield lock_file
    finally:
        _HELD.pop(_key(lock_file), None)
        # This process's lock only: not one that another process took after
        # it was deleted by hand.
        if _lock_holder(lock_file) == mine:
            with contextlib.suppress(FileNotFoundError):
                os.remove(lock_file)


# ---------------------------------------------------------------------------
# Staged in-place writes
# ---------------------------------------------------------------------------
def _umask() -> int:
    """The umask of the process, which can only be read by setting it."""
    mask = os.umask(0o022)
    os.umask(mask)
    return mask


def copy_fresh(source: str | Path, target: str | Path) -> None:
    """Copy the HDF5 file *source* into a new file *target*, object by object.

    Every object of the root group is copied as stored -- dtype, chunks,
    filters, attributes, the columns of a table in their order -- soft links
    (``/chroms``) as links, then the root's attributes with their dtypes.  The
    copy holds none of the space that edits left in *source*.
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


@contextlib.contextmanager
def staged_inplace(
    uri: str | Path,
    *,
    copy: bool = True,
    validate_groups: Sequence[str] | None = None,
) -> Generator[Path, None, None]:
    """Edit the file of a map through a staged sibling, which then replaces it.

    Under the map's lock (taken unless this process holds it), yields the path
    of the staged file: by default a fresh copy of the map's file
    (:func:`copy_fresh`), to edit in place; with ``copy=False``, a path where
    nothing is yet, for a writer that writes a whole new map.  Once the
    ``with`` block ends, the staged file is validated in full, given the mode
    of the map (of its target, for a symbolic link; a new file's otherwise),
    and renamed over the map's path: a symbolic link there is replaced, not its
    target.  On any failure, ``KeyboardInterrupt`` and ``SystemExit`` included,
    the staged file and every partial file named after it (cooler's
    ``.__writing__``) are removed, and the map is left as it was.

    The staged file is ``.<map name>.<lock uuid>.<n>.staged``, next to the map:
    a process killed while it writes leaves it, and the next writer that takes
    its stale lock removes it (see :func:`lock`).

    Args:
        uri: The map: a ``.metacool`` path, an ``.mmetacool`` path or one of
            its entries (the whole file is staged either way).
        copy: Start from a copy of the map (``False``: from nothing).
        validate_groups: The groups of the staged file to validate in full
            (default: the group *uri* names, or every entry of a bare
            ``.mmetacool``).

    Yields:
        The staged file's path.
    """
    from metacooler.schema import is_mmetacooler, list_metacoolers, validate
    from metacooler.uri import parse_uri

    file_text, group = parse_uri(str(uri))
    path = Path(file_text)
    with lock(path):
        held = _holding(lock_path(path))
        assert held is not None  # taken just above
        held.staged += 1
        prefix = f"{_staged_prefix(path, held.holder['uuid'])}{held.staged}."
        staged = path.with_name(prefix + "staged")
        try:
            if copy:
                copy_fresh(path, staged)
            yield staged
            if validate_groups is not None:
                groups = list(validate_groups)
            elif group.strip("/"):
                groups = [group]
            elif is_mmetacooler(staged):
                groups = list_metacoolers(staged)
            else:
                groups = ["/"]
            for name in groups:
                validate(f"{staged}::{name}", full=True)
            if path.exists():
                shutil.copymode(path, staged)
            else:
                os.chmod(staged, 0o666 & ~_umask())
            os.replace(staged, path)
        except BaseException:
            _remove_named(path.parent, prefix)
            raise


def map_output(output: str | Path) -> Path:
    """The file an output path names: *output*, or the file of ``file::/``.

    An output is a whole file: a URI naming a group inside a file is refused,
    as writing it would make a file named after the URI.

    Raises:
        ValueError: If *output* names a group inside a file.
    """
    text = str(output)
    if "::" not in text:
        return Path(output)
    file_text, _, group = text.partition("::")
    if group.strip().strip("/"):
        raise ValueError(
            f"{text} names the group {group!r} inside {file_text}, but an output "
            "is a whole file: give a path, or edit an entry of a collection in "
            "place."
        )
    return Path(file_text)


# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------
def file_state(path: str | Path) -> tuple[int, int, int, int] | None:
    """What changes when the file at *path* is replaced (another inode) or
    written (its size, its time of modification): ``None`` when there is no
    file there."""
    try:
        stat = os.stat(path)
    except FileNotFoundError:
        return None
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns


class StateCheckedPath(os.PathLike):
    """A file path that checks, each time it is opened, that the file is the
    one it was: what :class:`metacooler.MetaCooler` hands cooler, which opens
    the file for each chunk it reads."""

    def __init__(self, path: str, check) -> None:
        self._path = path
        self._check = check

    def __fspath__(self) -> str:
        self._check()
        return self._path

    def __str__(self) -> str:
        return self._path

    def __repr__(self) -> str:
        return repr(self._path)
