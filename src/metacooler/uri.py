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
    """
    if "::" in raw:
        file_part, group_part = raw.split("::", 1)
        group_part = group_part.strip()
        if not group_part or group_part == "/":
            group_path = "/"
        else:
            # Ensure leading slash
            group_path = group_part if group_part.startswith("/") else f"/{group_part}"
        return Path(file_part), group_path
    return Path(raw), "/"
