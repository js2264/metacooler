"""Multi-resolution mmetacool generation."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

logger = logging.getLogger("metacooler")


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
    from cooler._reduce import preferred_sequence

    out: list[int] = []
    for raw in rstring.split(","):
        token = raw.strip().lower()
        if not token:
            continue
        try:
            if token == "n":
                out.extend(preferred_sequence(base_binsize, maxres, "nice"))
            elif token == "b":
                out.extend(preferred_sequence(base_binsize, maxres, "binary"))
            elif token == "4dn":
                out.extend([1000, 2000, *preferred_sequence(5000, maxres, "nice")])
            elif token.endswith("n"):
                out.extend(preferred_sequence(int(token[:-1]), maxres, "nice"))
            elif token.endswith("b"):
                out.extend(preferred_sequence(int(token[:-1]), maxres, "binary"))
            else:
                out.append(int(token))
        except ValueError as exc:
            raise ValueError(
                f"Invalid resolution token {raw!r} in --resolutions "
                f"{rstring!r}. Expected integers or <n>N / <n>B / 4DN."
            ) from exc
    return out


def zoomify(
    uri: str | Path,
    output: str | Path,
    levels: int = 5,
    resolutions: str | None = None,
) -> None:
    """Create a multi-resolution ``.mmetacool`` file.

    The input must be a fixed-bin metacooler.  Each resolution is stored as an
    mmetacool entry named ``entries/<binsize>`` (open it as
    ``file.mmetacool::entries/<binsize>``).  The base resolution is copied
    verbatim; every coarser resolution is produced by re-binning the base with
    :func:`metacooler.ops.binnify.binnify` (so bins stay aligned to the base).

    Resolutions are chosen one of two ways:

    - if *resolutions* is given, it is a cooler-style string
      (``"1000,2000,5000"``, ``"1000N"``, ``"1000B"``, ``"4DN"``, …), expanded
      by :func:`_parse_resolutions`;
    - otherwise *levels* successive doublings of the base bin size are used.

    The base resolution is always included, resolutions finer than the base are
    ignored, and any resolution that would leave fewer than 2 bins is skipped.

    Args:
        uri: Path to the source ``.metacool`` file. Must be a standalone file;
            mmetacool entry URIs are not supported.
        output: Destination path for the output ``.mmetacool`` file.
        levels: Number of doublings to use when *resolutions* is not given.
        resolutions: Optional cooler-style resolution string. Takes precedence
            over *levels* when provided.

    Raises:
        ValueError: If the input is not a fixed-bin metacooler, is an
            mmetacool entry URI, or *resolutions* is malformed.
    """
    import shutil

    from metacooler.api import MetaCooler
    from metacooler.io.mmetacool import create_mmetacool
    from metacooler.ops.binnify import binnify

    output = Path(output)
    if output.exists():
        raise FileExistsError(
            f"Output file already exists: {output}. Remove it first or choose a "
            "different path."
        )

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

    # Determine the target resolutions.
    if resolutions is not None:
        maxres = int(mc.contigs["length"].max())
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

    entries: dict[str, Path] = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        source_path = Path(mc._store)

        for res in targets:
            res_path = tmpdir_path / f"res_{res}.metacool"
            if res == base_binsize:
                # Base resolution: copy the source file verbatim.
                shutil.copy2(source_path, res_path)
            else:
                try:
                    binnify(uri=source_path, output=res_path, binsize=res)
                except ValueError as exc:
                    logger.warning("Skipping resolution %d: %s", res, exc)
                    continue
                if MetaCooler(res_path).shape[0] < 2:
                    logger.warning(
                        "Skipping resolution %d: fewer than 2 bins after re-binning.",
                        res,
                    )
                    continue
            entries[str(res)] = res_path

        logger.info(
            "Creating mmetacool with %d resolution level(s): %s",
            len(entries),
            sorted(entries.keys(), key=int),
        )
        create_mmetacool(Path(output), entries)

    logger.info("Written multi-resolution mmetacool to %s", output)
