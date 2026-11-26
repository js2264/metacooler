"""URI parsing utility for metacooler/mmetacooler files.

Supports the ``file::group_path`` syntax used by cooler for multi-resolution
files.  For plain ``.metacool`` files the group path defaults to ``/``.
"""

from __future__ import annotations

import os
from pathlib import Path

#: The extensions a sample's name leaves out: the maps', then those of the
#: files metacooler writes. The longest that a name ends with goes.
_EXTENSIONS = (
    ".metacool",
    ".mmetacool",
    ".cool",
    ".mcool",
    ".tsv",
    ".tsv.gz",
    ".csv",
    ".csv.gz",
    ".txt",
    ".json",
    ".yaml",
    ".yml",
    ".parquet",
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fna",
    ".fna.gz",
    ".png",
    ".pdf",
    ".svg",
    ".html",
)


def parse_uri(raw: str) -> tuple[Path, str]:
    """Parse a metacooler/mmetacooler URI into ``(path, group_path)``.

    Accepted forms:

    * ``file.metacool``               → ``(Path("file.metacool"), "/")``
    * ``file.mmetacool::entries/e1``   → ``(Path("file.mmetacool"), "/entries/e1")``
    * ``file.mmetacool::``            → ``(Path("file.mmetacool"), "/")``

    Args:
        raw: A file path, optionally followed by ``::group_path``.

    Returns:
        A tuple of ``(file_path, hdf5_group_path)``.  The group path is
        always absolute (starts with ``/``).

    Raises:
        ValueError: If the file part is empty (e.g. ``""`` or ``"::group"``).
    """
    if "::" in raw:
        file_part, group_part = raw.split("::", 1)
        file_part = file_part.strip()
        group_part = group_part.strip()
        if not group_part or group_part == "/":
            group_path = "/"
        else:
            # Ensure leading slash
            group_path = group_part if group_part.startswith("/") else f"/{group_part}"
    else:
        file_part = raw.strip()
        group_path = "/"
    if not file_part:
        raise ValueError(f"URI has an empty file path: {raw!r}")
    return Path(file_part), group_path


def sample_name(uri: str | os.PathLike) -> str:
    """The name of the sample of a map or of an output, which names the log of
    a command and the report of a map, and MultiQC the report's row.

    A map's is its file's name without ``.metacool`` or ``.mmetacool``, or, for
    an entry of a collection (``file.mmetacool::entries/NAME``), the entry's
    name, as metaTOR names a map's sample: one map, one name, in both tools.
    Any other path's is its name without the extension of what metacooler
    writes there (``.cool``, ``.tsv.gz``, ``.png``, ...): a directory's, or a
    prefix's, is its name.

    Args:
        uri: A path, or a map's URI (``file.mmetacool::entries/NAME``).

    Returns:
        The sample's name.

    Raises:
        ValueError: If the URI's file part is empty (see :func:`parse_uri`).
    """
    file_path, group_path = parse_uri(os.fspath(uri))
    if group_path != "/":
        return group_path.rstrip("/").rsplit("/", 1)[-1]
    name = file_path.name
    for extension in sorted(_EXTENSIONS, key=len, reverse=True):
        if name.lower().endswith(extension) and len(name) > len(extension):
            return name[: -len(extension)]
    return name
