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
        style: ``"nice"`` for 1-2-5 × 10^k steps, ``"binary"`` for ×2 steps.

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
    - ``<n>B`` / bare ``B`` gives a binary (×2) progression;
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
) -> None:
    """Create a multi-resolution ``.mmetacool`` file.

    The input must be a fixed-bin metacooler.  Each resolution is stored as an
    mmetacool entry named ``entries/<binsize>`` (open it as
    ``file.mmetacool::entries/<binsize>``).  The base resolution is streamed
    verbatim into the output; every coarser resolution is produced by re-binning
    the base with :func:`metacooler.ops.binnify.binnify` (so bins stay aligned
    to the base) into *temp_dir*, appended to the output, and deleted right
    away — peak temporary disk use is one level, not a whole second pyramid.

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

    Raises:
        FileExistsError: If *output* exists and *force* is False.
        ValueError: If the input is not a fixed-bin metacooler, is an
            mmetacool entry URI, *resolutions* is malformed, or a requested
            resolution is not an exact multiple of the base bin size.
    """
    from metacooler.api import MetaCooler
    from metacooler.io._utils import check_output_path
    from metacooler.io.mmetacool import add_entry, create_mmetacool
    from metacooler.ops.binnify import binnify

    output = Path(output)
    check_output_path(output, force)

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
    source_path = Path(mc._store)
    try:
        with tempfile.TemporaryDirectory(
            dir=None if temp_dir is None else str(temp_dir)
        ) as tmpdir:
            tmpdir_path = Path(tmpdir)
            for res in plan:
                if res == base_binsize:
                    # Base resolution: stream the source straight into the output.
                    create_mmetacool(output, {str(res): source_path}, force=force)
                    continue
                res_path = tmpdir_path / f"res_{res}.metacool"
                binnify(uri=source_path, output=res_path, binsize=res)
                add_entry(output, str(res), res_path)
                # Reclaim the space before computing the next level.
                res_path.unlink(missing_ok=True)
    except Exception:
        # Never leave a half-built pyramid behind.
        if output.exists():
            output.unlink()
        raise

    logger.info("Written multi-resolution mmetacool to %s", output)
