"""URI parsing utility for metacooler/mmetacooler files.

Supports the ``file::group_path`` syntax used by cooler for multi-resolution
files.  For plain ``.metacool`` files the group path defaults to ``/``.
"""

from __future__ import annotations

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
