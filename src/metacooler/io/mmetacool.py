"""MMetaCooler I/O operations — create, add, remove, copy, move, list."""

from __future__ import annotations

import logging
from pathlib import Path

import h5py

from metacooler.schema import (
    MMETACOOLER_FORMAT,
    MMETACOOLER_FORMAT_VERSION,
    list_entry_names,
)

logger = logging.getLogger("metacooler")


def init_mmetacool_attrs(f: h5py.File) -> None:
    """Initialise an open HDF5 file as an empty ``.mmetacool`` collection.

    Sets ``format`` to ``HDF5::MMetaCooler`` and ``format-version`` to the
    current :data:`~metacooler.schema.MMETACOOLER_FORMAT_VERSION`, and
    creates the top-level ``/entries`` group.  Both are required by
    :func:`metacooler.schema.is_mmetacooler`, so a file initialised through
    this function is a valid (empty) mmetacool.

    Args:
        f: Open ``h5py.File`` in write or append mode.  Attributes are
            written directly to the file root (``f.attrs``).
    """
    f.attrs["format"] = MMETACOOLER_FORMAT
    f.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
    f.require_group("entries")


#: Deprecated private alias kept so existing importers keep working.
_init_mmetacool_attrs = init_mmetacool_attrs


def _require_mmetacool(path: Path) -> Path:
    """Check that *path* is an existing, valid ``.mmetacool`` file.

    Every mutating helper below opens its target in HDF5 append mode, which
    happily *creates* a missing file — leaving behind a stray file that has
    none of the mmetacool attributes and that
    :func:`metacooler.schema.is_mmetacooler` then rejects.  This guard is
    called first so a bad path fails cleanly instead.

    Args:
        path: Path to the ``.mmetacool`` file to check.

    Returns:
        The path as a :class:`pathlib.Path`.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If *path* exists but is not a valid mmetacool file.
    """
    from metacooler.schema import is_mmetacooler

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"mmetacool file not found: {path}. "
            "Create it first (e.g. with create_mmetacool() or 'metacooler stack')."
        )
    if not is_mmetacooler(path):
        raise ValueError(
            f"Not a valid .mmetacool file: {path} (expected format attribute "
            f"{MMETACOOLER_FORMAT!r} and a top-level 'entries' group)."
        )
    return path


def create_mmetacool(
    output: Path,
    entries: dict[str, Path],
    force: bool = False,
) -> None:
    """Create a new ``.mmetacool`` file from a dict of ``{name: metacool_path}``.

    Each entry is a full metacooler stored as an HDF5 group under
    ``/entries/<name>``.

    Args:
        output: Destination path for the ``.mmetacool`` file.
        entries: Mapping of entry names to ``.metacool`` file paths.
        force: If ``True``, overwrite an existing output file.

    Raises:
        FileExistsError: If *output* already exists and *force* is ``False``.
        FileNotFoundError: If any source metacool file does not exist.
    """
    if output.exists() and not force:
        raise FileExistsError(
            f"Output file already exists: {output}. Pass force=True to overwrite."
        )
    for name, src in entries.items():
        if not Path(src).exists():
            raise FileNotFoundError(f"Source file not found for entry {name!r}: {src}")

    with h5py.File(output, "w") as dst:
        init_mmetacool_attrs(dst)
        for name, src in entries.items():
            logger.info("Adding entry %r from %s", name, src)
            with h5py.File(src, "r") as src_f:
                src_f.copy(src_f["/"], dst["entries"], name=name)

    logger.info("Created mmetacool with %d entries: %s", len(entries), output)


def add_entry(
    mmetacool: Path,
    name: str,
    metacool: Path,
    overwrite: bool = False,
) -> None:
    """Add a metacooler entry to an existing ``.mmetacool`` file.

    Args:
        mmetacool: Path to the ``.mmetacool`` file.
        name: Name for the new entry (stored under ``/entries/<name>``).
        metacool: Path to the source ``.metacool`` file.
        overwrite: If True, replace an existing entry with the same name.

    Raises:
        FileNotFoundError: If *metacool* or *mmetacool* does not exist.
        ValueError: If *mmetacool* is not a valid mmetacool file.
        KeyError: If *name* already exists and *overwrite* is False.
    """
    if not Path(metacool).exists():
        raise FileNotFoundError(f"Source file not found: {metacool}")
    mmetacool = _require_mmetacool(mmetacool)

    with h5py.File(mmetacool, "a") as dst:
        grp = dst.require_group("entries")
        if name in grp:
            if not overwrite:
                raise KeyError(
                    f"Entry {name!r} already exists. Use overwrite=True to replace."
                )
            del grp[name]
        with h5py.File(metacool, "r") as src_f:
            src_f.copy(src_f["/"], grp, name=name)
    logger.info("Added entry %r to %s", name, mmetacool)


def remove_entry(mmetacool: Path, name: str) -> None:
    """Remove an entry from an ``.mmetacool`` file.

    Args:
        mmetacool: Path to the ``.mmetacool`` file.
        name: Name of the entry to remove.

    Raises:
        FileNotFoundError: If *mmetacool* does not exist.
        ValueError: If *mmetacool* is not a valid mmetacool file.
        KeyError: If *name* does not exist.
    """
    mmetacool = _require_mmetacool(mmetacool)
    with h5py.File(mmetacool, "a") as f:
        grp = f.require_group("entries")
        if name not in grp:
            raise KeyError(f"Entry {name!r} not found in {mmetacool}")
        del grp[name]
    logger.info("Removed entry %r from %s", name, mmetacool)


def copy_entry(
    mmetacool: Path,
    src: str,
    dst: str,
    overwrite: bool = False,
) -> None:
    """Copy an entry within the same ``.mmetacool`` file.

    Args:
        mmetacool: Path to the ``.mmetacool`` file.
        src: Source entry name.
        dst: Destination entry name.
        overwrite: If True, replace an existing destination entry.

    Raises:
        FileNotFoundError: If *mmetacool* does not exist.
        KeyError: If *src* does not exist or *dst* exists and *overwrite* is False.
        ValueError: If *src* and *dst* are the same, or *mmetacool* is not a
            valid mmetacool file.
    """
    if src == dst:
        raise ValueError("Source and destination cannot be the same.")
    mmetacool = _require_mmetacool(mmetacool)

    with h5py.File(mmetacool, "a") as f:
        grp = f.require_group("entries")
        if src not in grp:
            raise KeyError(f"Source entry {src!r} not found.")
        if dst in grp:
            if not overwrite:
                raise KeyError(
                    f"Destination entry {dst!r} already exists. Use overwrite=True."
                )
            del grp[dst]
        f.copy(grp[src], grp, name=dst)
    logger.info("Copied entry %r → %r in %s", src, dst, mmetacool)


def move_entry(
    mmetacool: Path,
    src: str,
    dst: str,
    overwrite: bool = False,
) -> None:
    """Move (rename) an entry within the same ``.mmetacool`` file.

    Implemented as an HDF5 link rename (:meth:`h5py.Group.move`), which is
    O(1) and rewrites no data.  Copying then deleting instead would need
    twice the peak disk space and would leave the freed space stranded —
    HDF5 never returns deleted space to the filesystem, so a rename would
    permanently grow the file by the size of the entry.

    Args:
        mmetacool: Path to the ``.mmetacool`` file.
        src: Source entry name.
        dst: Destination entry name.
        overwrite: If True, replace an existing destination entry.

    Raises:
        FileNotFoundError: If *mmetacool* does not exist.
        KeyError: If *src* does not exist or *dst* exists and *overwrite* is False.
        ValueError: If *src* and *dst* are the same, or *mmetacool* is not a
            valid mmetacool file.
    """
    if src == dst:
        raise ValueError("Source and destination cannot be the same.")
    mmetacool = _require_mmetacool(mmetacool)

    with h5py.File(mmetacool, "a") as f:
        grp = f.require_group("entries")
        if src not in grp:
            raise KeyError(f"Source entry {src!r} not found.")
        if dst in grp:
            if not overwrite:
                raise KeyError(
                    f"Destination entry {dst!r} already exists. Use overwrite=True."
                )
            del grp[dst]
        grp.move(src, dst)
    logger.info("Moved entry %r → %r in %s", src, dst, mmetacool)


def list_entries(mmetacool: Path) -> list[str]:
    """List all entry names in an ``.mmetacool`` file.

    Both this and :func:`metacooler.schema.list_metacoolers` are thin
    conventions over the single scan in
    :func:`metacooler.schema.list_entry_names`: this one yields bare entry
    names (``"sample1"``), that one yields URI-ready group paths
    (``"/entries/sample1"``).

    Args:
        mmetacool: Path to the ``.mmetacool`` file.

    Returns:
        A list of entry names (keys under ``/entries``), or an empty list if
        the file has no ``entries`` group.
    """
    return list_entry_names(mmetacool)


def copy_entry_cross_file(
    src_file: Path,
    src_name: str,
    dst_file: Path,
    dst_name: str,
    overwrite: bool = False,
) -> None:
    """Copy an entry from one ``.mmetacool`` file to another.

    Args:
        src_file: Source ``.mmetacool`` file path.
        src_name: Entry name in the source file.
        dst_file: Destination ``.mmetacool`` file path.
        dst_name: Entry name in the destination file.
        overwrite: If True, replace an existing destination entry.

    Raises:
        FileNotFoundError: If *src_file* or *dst_file* does not exist.
        KeyError: If source entry does not exist or destination exists
            and *overwrite* is False.
        ValueError: If source and destination resolve to the same file, or
            either is not a valid mmetacool file.
    """
    if Path(src_file).resolve() == Path(dst_file).resolve():
        raise ValueError(
            f"Source and destination are the same file: {src_file}. "
            "Use rename_entry() to rename entries within a file."
        )
    src_file = _require_mmetacool(src_file)
    dst_file = _require_mmetacool(dst_file)
    with h5py.File(src_file, "r") as sf:
        if "entries" not in sf or src_name not in sf["entries"]:
            raise KeyError(f"Source entry {src_name!r} not found in {src_file}")

        with h5py.File(dst_file, "a") as df:
            grp = df.require_group("entries")
            if dst_name in grp:
                if not overwrite:
                    raise KeyError(
                        f"Destination entry {dst_name!r} already exists in {dst_file}."
                    )
                del grp[dst_name]
            sf.copy(sf[f"entries/{src_name}"], grp, name=dst_name)

    logger.info("Copied %s::%s → %s::%s", src_file, src_name, dst_file, dst_name)
