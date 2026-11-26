"""URI parsing utility for metacooler/mmetacooler files.

Supports the ``file::group_path`` syntax used by cooler for multi-resolution
files.  For plain ``.metacool`` files the group path defaults to ``/``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path


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
    """The name of a map's sample, which names its report, and MultiQC the
    report's row.

    It is the map's file name without ``.metacool`` or ``.mmetacool``, or, for
    an entry of a collection (``file.mmetacool::entries/NAME``), the entry's
    name: metaTOR's rule, so that one map has one name in both tools.

    Args:
        uri: A map's path, or its URI (``file.mmetacool::entries/NAME``).

    Returns:
        The sample's name.

    Raises:
        ValueError: If the URI's file part is empty (see :func:`parse_uri`).
    """
    file_path, group_path = parse_uri(os.fspath(uri))
    if group_path != "/":
        name = group_path.rstrip("/").rsplit("/", 1)[-1]
    else:
        name = file_path.name
    return re.sub(r"\.m?metacool$", "", name) or name
