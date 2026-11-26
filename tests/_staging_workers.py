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


def hold(path: str, ready, release, follow_symlinks: bool = True) -> None:
    """Hold the lock of *path* until *release* is set (with *follow_symlinks*
    false, the lock of the path itself: an output's)."""
    logging.disable(logging.WARNING)
    from metacooler.io.staging import lock

    with lock(path, command="test holder", follow_symlinks=follow_symlinks):
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


def try_store(path: str) -> str:
    """Store a weight column in the map *path*: ``"written"``, or what was
    raised, as ``"<type>: <message>"``."""
    logging.disable(logging.WARNING)
    import numpy as np

    from metacooler.api import MetaCooler
    from metacooler.io.weights import store_weights

    try:
        n_bins = int(MetaCooler(path).info["nbins"])
        store_weights(path, "w", np.ones(n_bins), method="test", depends_on="pixels")
    except Exception as error:  # reported to the test
        return f"{type(error).__name__}: {error}"
    return "written"


def try_reassign(source: str, table: str, output: str) -> str:
    """Reassign the map *source* by *table* into *output*, replacing it:
    ``"written"``, or what was raised, as ``"<type>: <message>"``."""
    logging.disable(logging.WARNING)
    from metacooler.io.reassign import reassign

    try:
        reassign(source, table, output=output, force=True)
    except Exception as error:  # reported to the test
        return f"{type(error).__name__}: {error}"
    return "written"
