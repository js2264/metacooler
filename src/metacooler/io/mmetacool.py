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


def _init_mmetacool_attrs(f: h5py.File) -> None:
    """Write the top-level format attributes for a new ``.mmetacool`` file.

    Sets ``format`` to ``HDF5::MMetaCooler`` and ``format-version`` to the
    current :data:`~metacooler.schema.MMETACOOLER_FORMAT_VERSION`.

    Args:
        f: Open ``h5py.File`` in write or append mode.  Attributes are
            written directly to the file root (``f.attrs``).
    """
    f.attrs["format"] = MMETACOOLER_FORMAT
    f.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION


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
        _init_mmetacool_attrs(dst)
        dst.require_group("entries")
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
        FileNotFoundError: If *metacool* does not exist.
        KeyError: If *name* already exists and *overwrite* is False.
    """
    if not metacool.exists():
        raise FileNotFoundError(f"Source file not found: {metacool}")

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
        KeyError: If *name* does not exist.
    """
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
        KeyError: If *src* does not exist or *dst* exists and *overwrite* is False.
        ValueError: If *src* and *dst* are the same.
    """
    if src == dst:
        raise ValueError("Source and destination cannot be the same.")

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

    Args:
        mmetacool: Path to the ``.mmetacool`` file.
        src: Source entry name.
        dst: Destination entry name.
        overwrite: If True, replace an existing destination entry.

    Raises:
        KeyError: If *src* does not exist or *dst* exists and *overwrite* is False.
        ValueError: If *src* and *dst* are the same.
    """
    copy_entry(mmetacool, src, dst, overwrite=overwrite)
    remove_entry(mmetacool, src)
    logger.info("Moved entry %r → %r in %s", src, dst, mmetacool)


def list_entries(mmetacool: Path) -> list[str]:
    """List all entry names in an ``.mmetacool`` file.

    Returns:
        A list of entry names (keys under ``/entries``).
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
        KeyError: If source entry does not exist or destination exists
            and *overwrite* is False.
        ValueError: If source and destination resolve to the same file.
    """
    if Path(src_file).resolve() == Path(dst_file).resolve():
        raise ValueError(
            f"Source and destination are the same file: {src_file}. "
            "Use rename_entry() to rename entries within a file."
        )
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
