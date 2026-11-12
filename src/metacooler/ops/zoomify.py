"""Multi-resolution mmetacool generation."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Literal

logger = logging.getLogger("metacooler")


def _geomprog(start: int, mul: int) -> Iterator[int]:
    """Yield an unbounded geometric progression of integers.

    Mirrors ``cooler._reduce.geomprog`` (vendored to avoid importing a cooler
    private).

    Args:
        start: First value of the progression.
        mul: Integer ratio between successive values.

    Yields:
        ``start``, ``start * mul``, ``start * mul ** 2``, …
    """
    start, mul = int(start), int(mul)
    yield start
    while True:
        start *= mul
        yield start


def _niceprog(start: int) -> Iterator[int]:
    """Yield an unbounded "nice" 1-2-5 progression of integers.

    Mirrors ``cooler._reduce.niceprog`` (vendored to avoid importing a cooler
    private): a geometric progression of ratio 10 interspersed with steps of
    ratio 2 and 5.

    Args:
        start: First value of the progression.

    Yields:
        ``start``, ``2 * start``, ``5 * start``, ``10 * start``, …
    """
    start = int(start)
    yield start
    while True:
        for mul in (2, 5, 10):
            yield start * mul
        start *= 10


def _preferred_sequence(
    start: int, stop: int, style: Literal["binary", "nice"] = "nice"
) -> list[int]:
    """Return a sequence of integers with a "preferred" stepping pattern.

    Mirrors ``cooler._reduce.preferred_sequence`` exactly (vendored to avoid
    importing a cooler private), including the inclusive upper bound and the
    empty result when *start* exceeds *stop*.

    Args:
        start: Starting value of the progression.
        stop: Inclusive upper bound; returned values never exceed it.
        style: ``"nice"`` for 1-2-5 x 10^k steps, ``"binary"`` for x2 steps.

    Returns:
        The list of progression values in ascending order.

    Raises:
        ValueError: If *style* is not ``"nice"`` or ``"binary"``.
    """
    if start > stop:
        return []

    if style == "binary":
        gen = _geomprog(start, 2)
    elif style == "nice":
        gen = _niceprog(start)
    else:
        raise ValueError(f"Expected style value of 'binary' or 'nice'; got '{style}'.")

    seq = [next(gen)]
    while True:
        n = next(gen)
        if n > stop:
            break
        seq.append(n)

    return seq


def _parse_resolutions(rstring: str, base_binsize: int, maxres: int) -> list[int]:
    """Expand a cooler-style resolution string into a list of bin sizes.

    Mirrors ``cooler zoomify --resolutions``:

    - a plain integer is used verbatim (``5000``);
    - ``<n>N`` / bare ``N`` gives a "nice" 1-2-5 progression from ``<n>``
      (or from *base_binsize* for bare ``N``) up to *maxres*;
    - ``<n>B`` / bare ``B`` gives a binary (x2) progression;
    - ``4DN`` is an alias for ``1000,2000,5000N``.

    Tokens are comma-separated and case-insensitive (``1000,2000,5000N``).
    """
    out: list[int] = []
    for raw in rstring.split(","):
        token = raw.strip().lower()
        if not token:
            continue
        try:
            if token == "n":
                out.extend(_preferred_sequence(base_binsize, maxres, "nice"))
            elif token == "b":
                out.extend(_preferred_sequence(base_binsize, maxres, "binary"))
            elif token == "4dn":
                out.extend([1000, 2000, *_preferred_sequence(5000, maxres, "nice")])
            elif token.endswith("n"):
                out.extend(_preferred_sequence(int(token[:-1]), maxres, "nice"))
            elif token.endswith("b"):
                out.extend(_preferred_sequence(int(token[:-1]), maxres, "binary"))
            else:
                out.append(int(token))
        except ValueError as exc:
            raise ValueError(
                f"Invalid resolution token {raw!r} in --resolutions "
                f"{rstring!r}. Expected integers or <n>N / <n>B / 4DN."
            ) from exc
    return out


def _check_divisible(targets: list[int], base_binsize: int) -> None:
    """Refuse resolutions that are not exact multiples of the base bin size.

    Every level of the pyramid is produced by re-binning the base, so a
    resolution that is not a multiple of the base bin size cannot nest in it
    (and is rejected by :func:`metacooler.ops.binnify.binnify` anyway).

    Args:
        targets: Requested resolutions, all >= *base_binsize*.
        base_binsize: Bin size of the base resolution.

    Raises:
        ValueError: If any target is not an exact multiple of *base_binsize*.
    """
    bad = sorted({r for r in targets if r % base_binsize})
    if not bad:
        return
    suggestions = ", ".join(
        f"{r} -> {base_binsize * max(1, round(r / base_binsize))}" for r in bad
    )
    raise ValueError(
        f"Requested resolution(s) {bad} are not exact multiples of the base bin "
        f"size ({base_binsize}); such levels would not nest in the base binning. "
        f"Nearest valid resolutions: {suggestions}."
    )


def _bins_per_contig(lengths: np.ndarray, binsize: int) -> np.ndarray:
    """Return the number of bins each contig would get at *binsize*.

    Args:
        lengths: Contig lengths in base pairs.
        binsize: Candidate bin size in base pairs.

    Returns:
        ``ceil(length / binsize)`` for every contig.
    """
    return -(-lengths // binsize)


def zoomify(
    uri: str | Path,
    output: str | Path,
    levels: int = 5,
    resolutions: str | None = None,
    temp_dir: str | Path | None = None,
    force: bool = False,
    nproc: int = 1,
) -> None:
    """Create a multi-resolution ``.mmetacool`` file.

    The input must be a fixed-bin metacooler.  Each resolution is stored as an
    mmetacool entry named ``entries/<binsize>`` (open it as
    ``file.mmetacool::entries/<binsize>``).  The base resolution is streamed
    verbatim into the output; every coarser resolution is produced by
    re-binning, with :func:`metacooler.ops.binnify.binnify`, the coarsest level
    already written whose bin size divides it (10 kb from 5 kb, 2 kb and 5 kb
    from a 1 kb base), as cooler's zoomify does.  Every bin size is a multiple
    of the base, so bins stay aligned to the base and each level holds exactly
    the counts re-binning the base would give, from far fewer pixels.  A level
    is written into *temp_dir*, appended to the output and deleted right away —
    peak temporary disk use is one level, not a whole second pyramid.  Its
    provenance names the base file, as if re-binned from it.

    Resolutions are chosen one of two ways:

    - if *resolutions* is given, it is a cooler-style string
      (``"1000,2000,5000"``, ``"1000N"``, ``"1000B"``, ``"4DN"``, …), expanded
      by :func:`_parse_resolutions`;
    - otherwise *levels* successive doublings of the base bin size are used.

    The base resolution is always included and resolutions finer than the base
    are ignored.  Every remaining resolution must be an exact multiple of the
    base bin size.  Levels that would leave fewer than 2 bins, and levels that
    would be identical to an already-written one (once the bin size exceeds the
    longest contig, every contig collapses to a single bin), are skipped and
    reported in the log.

    Args:
        uri: Path to the source ``.metacool`` file. Must be a standalone file;
            mmetacool entry URIs are not supported.
        output: Destination path for the output ``.mmetacool`` file.
        levels: Number of doublings to use when *resolutions* is not given.
        resolutions: Optional cooler-style resolution string. Takes precedence
            over *levels* when provided.
        temp_dir: Directory for the intermediate re-binned files. Defaults to
            the system temporary directory. Created if missing.
        force: If True, overwrite *output* when it already exists.
        nproc: Levels built at once, each in its own process: with
            ``1000,2000,5000,10000``, 2 kb and 5 kb are built together, then
            10 kb.  Levels are appended in ladder order, so the output is that
            of ``nproc=1``.  A level's file is then kept in *temp_dir* until
            the levels coarsened from it are built.

    Raises:
        FileExistsError: If *output* exists and *force* is False.
        ValueError: If the input is not a fixed-bin metacooler, is an
            mmetacool entry URI, *resolutions* is malformed, or a requested
            resolution is not an exact multiple of the base bin size.
    """
    from metacooler.api import MetaCooler
    from metacooler.io._utils import check_output_path
    from metacooler.io.mmetacool import add_entry, create_mmetacool
    from metacooler.ops.binnify import _binnify

    output = Path(output)
    check_output_path(output, force, inputs=[uri])

    mc = MetaCooler(uri)
    if mc._group_path not in ("/", ""):
        raise ValueError(
            "zoomify does not support mmetacool entry URIs (the base resolution is "
            "copied as a whole file). Extract the entry to a standalone .metacool "
            "first, then zoomify that file."
        )
    base_binsize = mc.binsize
    if base_binsize is None:
        raise ValueError(
            "Input must be a fixed-bin metacooler for zoomify. Got variable-bin input."
        )

    if mc.shape[0] < 2:
        raise ValueError("Input has fewer than 2 bins; cannot zoomify.")

    contig_lengths = mc.contigs["length"].to_numpy(dtype=np.int64)

    # Determine the target resolutions.
    if resolutions is not None:
        maxres = int(contig_lengths.max())
        requested = _parse_resolutions(resolutions, base_binsize, maxres)
    else:
        requested = [base_binsize * (2**i) for i in range(max(levels, 1))]

    below = sorted({r for r in requested if r < base_binsize})
    if below:
        logger.warning(
            "Ignoring resolutions finer than the base bin size %d: %s",
            base_binsize,
            below,
        )
    # Always include the base; keep only base-or-coarser; sort ascending.
    targets = sorted({base_binsize} | {r for r in requested if r >= base_binsize})
    _check_divisible(targets, base_binsize)

    # Predict the bin count of every level from the contig lengths, so that no
    # work is done for a level that is unusable or a duplicate of a coarser one.
    skipped: list[str] = []
    plan: list[int] = []
    collapsed_at: int | None = None
    for i, res in enumerate(targets):
        per_contig = _bins_per_contig(contig_lengths, res)
        if int(per_contig.sum()) < 2:
            skipped += [f"{r} (fewer than 2 bins)" for r in targets[i:]]
            break
        if collapsed_at is not None:
            skipped += [
                f"{r} (identical to resolution {collapsed_at}: one bin per contig)"
                for r in targets[i:]
            ]
            break
        plan.append(res)
        if bool((per_contig == 1).all()):
            collapsed_at = res
    if skipped:
        logger.warning("Skipping %d requested resolution(s): %s", len(skipped), skipped)

    if temp_dir is not None:
        temp_dir = Path(temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Creating mmetacool with %d resolution level(s): %s",
        len(plan),
        [str(r) for r in plan],
    )
    if nproc < 1:
        raise ValueError(f"nproc must be at least 1, got {nproc}")
    source_path = Path(mc._store)
    try:
        with tempfile.TemporaryDirectory(
            dir=None if temp_dir is None else str(temp_dir)
        ) as tmpdir:
            tmpdir_path = Path(tmpdir)
            if nproc > 1 and len(plan) > 2:
                _build_levels_in_parallel(
                    plan, source_path, mc.metadata, output, tmpdir_path, nproc, force
                )
                logger.info("Written multi-resolution mmetacool to %s", output)
                return
            built: list[int] = []
            for res in plan:
                if res == base_binsize:
                    # Base resolution: stream the source straight into the output.
                    create_mmetacool(output, {str(res): source_path}, force=force)
                    built.append(res)
                    continue
                # The coarsest level written so far that nests in this one; the
                # base always does (_check_divisible).
                parent = max(r for r in built if res % r == 0)
                parent_mc = (
                    mc
                    if parent == base_binsize
                    else MetaCooler(f"{output}::entries/{parent}")
                )
                res_path = tmpdir_path / f"res_{res}.metacool"
                logger.info("Coarsening %d -> %d", parent, res)
                _binnify(
                    parent_mc,
                    res_path,
                    res,
                    source_uri=source_path,
                    source_metadata=mc.metadata,
                )
                add_entry(output, str(res), res_path)
                built.append(res)
                # Reclaim the space before computing the next level.
                res_path.unlink(missing_ok=True)
    except Exception:
        # Never leave a half-built pyramid behind.
        if output.exists():
            output.unlink()
        raise

    logger.info("Written multi-resolution mmetacool to %s", output)


def _parents(plan: list[int]) -> dict[int, int]:
    """The level each non-base level of *plan* is coarsened from.

    The coarsest earlier level whose bin size divides it -- the base always
    does (``_check_divisible``).
    """
    return {
        res: max(r for r in plan[:i] if res % r == 0)
        for i, res in enumerate(plan)
        if i > 0
    }


def _build_level(
    parent: str, output: str, binsize: int, source_uri: str, source_metadata: dict
) -> None:
    """Coarsen *parent* into *output* (one zoomify level, in a worker process)."""
    from metacooler.api import MetaCooler
    from metacooler.ops.binnify import _binnify

    _binnify(
        MetaCooler(parent),
        Path(output),
        binsize,
        source_uri=source_uri,
        source_metadata=source_metadata,
    )


def _build_levels_in_parallel(
    plan: list[int],
    source_path: Path,
    source_metadata: dict,
    output: Path,
    tmpdir: Path,
    nproc: int,
    force: bool,
) -> None:
    """Build every level of *plan* on a pool of *nproc* processes.

    A level is submitted as soon as its parent (see :func:`_parents`) exists --
    as a file in *tmpdir*, since the output is only written from this process.
    Finished levels are appended in ladder order, and a level's file is deleted
    once it is appended and every level coarsened from it is built.
    """
    import multiprocessing
    from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

    from metacooler.io.mmetacool import add_entry, create_mmetacool

    base = plan[0]
    parents = _parents(plan)
    create_mmetacool(output, {str(base): source_path}, force=force)
    path_of: dict[int, Path] = {base: source_path}
    done: set[int] = {base}
    appended = 1  # plan[:appended] are in the output
    running: dict = {}

    # spawn: a forked worker would inherit this process's HDF5 state.
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=nproc, mp_context=context) as pool:

        def submit_ready() -> None:
            for res in plan[1:]:
                if res in done or res in running.values() or parents[res] not in done:
                    continue
                path_of[res] = tmpdir / f"res_{res}.metacool"
                logger.info("Coarsening %d -> %d", parents[res], res)
                future = pool.submit(
                    _build_level,
                    str(path_of[parents[res]]),
                    str(path_of[res]),
                    res,
                    str(source_path),
                    source_metadata,
                )
                running[future] = res

        submit_ready()
        while running:
            finished, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for future in finished:
                res = running.pop(future)
                future.result()  # re-raise a worker's error here
                done.add(res)
            while appended < len(plan) and plan[appended] in done:
                add_entry(output, str(plan[appended]), path_of[plan[appended]])
                appended += 1
            for res in plan[1:appended]:
                children = [c for c, p in parents.items() if p == res]
                if all(c in done for c in children):
                    path_of[res].unlink(missing_ok=True)
            submit_ready()
