"""Workers of the concurrency tests of ``test_staging.py``, run in spawned
processes: importable on their own, light to import."""

from __future__ import annotations

import logging
import os
from pathlib import Path


def contend(path: str, rounds: int, barrier, results) -> None:
    """Contend for the lock of *path*, round after round.

    Each round: wait for the coordinator to write a stale lock, try to take it,
    and record whether this process held it.  A holder keeps it until every
    contender has tried (the second barrier), so that one exactly may hold it.
    """
    logging.disable(logging.WARNING)
    from metacooler.io.staging import MapLockedError, lock

    for r in range(rounds):
        barrier.wait(60)
        try:
            with lock(path):
                results.put((r, os.getpid(), "held"))
                barrier.wait(60)
        except MapLockedError:
            results.put((r, os.getpid(), "refused"))
            barrier.wait(60)
        barrier.wait(60)


def hold(path: str, ready, release) -> None:
    """Hold the lock of *path* until *release* is set."""
    logging.disable(logging.WARNING)
    from metacooler.io.staging import lock

    with lock(path, command="test holder"):
        ready.set()
        release.wait(60)


def write_in_place(path: str, ready, release, done) -> None:
    """Hold the lock of *path*, then, once released, replace the map through a
    staged copy (a metadata edit)."""
    logging.disable(logging.WARNING)
    from metacooler.io.staging import edit_metadata, lock

    with lock(path, command="test writer"):
        ready.set()
        release.wait(60)
        edit_metadata(path, lambda meta: meta.update({"written": "by another"}))
    done.set()


def exists(path: str) -> bool:
    return Path(path).exists()
