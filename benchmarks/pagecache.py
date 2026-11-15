#!/usr/bin/env python
"""Evict files from the page cache, and say how much of each is still cached.

    python pagecache.py evict FILE...      # then fail if any stays > 1 % cached
    python pagecache.py resident FILE...   # the cached share of each, one per line

``bench.sh`` evicts a command's inputs before each of its cold runs, so the
bytes it reads come from the disk rather than from memory. ``evict`` flushes
each file (a page still being written cannot be dropped), asks the kernel to
drop its pages (``posix_fadvise(POSIX_FADV_DONTNEED)``, which any user may do
for a file they can read), then checks with ``mincore`` that they are gone. A
file system that keeps its files in memory (tmpfs) cannot drop them: ``evict``
then exits with status 1 and names the file, so no run is called cold that was
not. Only the named files are evicted; Python and the libraries both tools
import stay cached, as on any machine that ran them recently.
"""

from __future__ import annotations

import ctypes
import mmap
import os
import sys

import numpy as np

#: The largest share of a file that may still be cached after an eviction.
TOLERANCE = 0.01

_libc = ctypes.CDLL(None, use_errno=True)
_libc.mmap.restype = ctypes.c_void_p
_libc.mmap.argtypes = [
    ctypes.c_void_p,
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_long,
]
_libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p]
_libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
_MAP_FAILED = ctypes.c_void_p(-1).value


def resident(path: str) -> float:
    """The share of the file's pages in the page cache (0 for an empty file)."""
    size = os.path.getsize(path)
    if size == 0:
        return 0.0
    fd = os.open(path, os.O_RDONLY)
    try:
        addr = _libc.mmap(None, size, mmap.PROT_READ, mmap.MAP_SHARED, fd, 0)
        if addr in (None, _MAP_FAILED):
            raise OSError(ctypes.get_errno(), f"mmap {path}")
        try:
            pages = -(-size // mmap.PAGESIZE)
            vec = np.zeros(pages, np.uint8)
            if _libc.mincore(addr, size, vec.ctypes.data) != 0:
                raise OSError(ctypes.get_errno(), f"mincore {path}")
            return float((vec & 1).sum()) / pages
        finally:
            _libc.munmap(addr, size)
    finally:
        os.close(fd)


def evict(path: str) -> float:
    """Drop the file's pages from the page cache; the share still cached."""
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)
    return resident(path)


def main(argv: list[str]) -> int:
    if len(argv) < 2 or argv[0] not in ("evict", "resident"):
        print(__doc__, file=sys.stderr)
        return 2
    action, paths = argv[0], argv[1:]
    if action == "resident":
        for path in paths:
            print(f"{resident(path):.4f}\t{path}")
        return 0
    stuck = [(p, share) for p in paths if (share := evict(p)) > TOLERANCE]
    for path, share in stuck:
        print(f"{path}: {share:.0%} still cached after eviction", file=sys.stderr)
    return 1 if stuck else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
