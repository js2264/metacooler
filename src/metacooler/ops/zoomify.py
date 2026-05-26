"""Multi-resolution mmetacool generation."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

logger = logging.getLogger("metacooler")


def zoomify(
    uri: str | Path,
    output: str | Path,
    levels: int = 5,
) -> None:
    """Create a multi-resolution ``.mmetacool`` file.

    The input must be a fixed-bin metacooler.  The base resolution is
    copied as ``resolutions/<base_binsize>``, and each subsequent level
    doubles the bin size until *levels* is reached or fewer than 2 bins
    remain.

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.mmetacool`` file.
        levels: Maximum number of resolution levels to generate.

    Raises:
        ValueError: If the input is not a fixed-bin metacooler.
    """
    import shutil

    from metacooler.api import MetaCooler
    from metacooler.io.mmetacool import create_mmetacool
    from metacooler.ops.binnify import binnify

    output = Path(output)
    if output.exists():
        raise FileExistsError(
            f"Output file already exists: {output}. Remove it first or choose a different path."
        )

    mc = MetaCooler(uri)
    base_binsize = mc.binsize
    if base_binsize is None:
        raise ValueError(
            "Input must be a fixed-bin metacooler for zoomify. Got variable-bin input."
        )

    n_bins = mc.shape[0]
    if n_bins < 2:
        raise ValueError("Input has fewer than 2 bins; cannot zoomify.")

    entries: dict[str, Path] = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Base resolution: copy the original file
        base_path = tmpdir_path / f"res_{base_binsize}.metacool"
        shutil.copy2(uri, base_path)
        entries[str(base_binsize)] = base_path

        current_binsize = base_binsize
        current_uri = Path(uri)
        current_nbins = n_bins

        for i in range(1, levels):
            new_binsize = current_binsize * 2

            # Check if we'd have fewer than 2 bins
            if current_nbins < 4:
                logger.info(
                    "Stopping at level %d: too few bins (%d) for further coarsening.",
                    i,
                    current_nbins,
                )
                break

            new_path = tmpdir_path / f"res_{new_binsize}.metacool"
            try:
                binnify(
                    uri=current_uri,
                    output=new_path,
                    binsize=new_binsize,
                )
            except ValueError:
                logger.info(
                    "Stopping at level %d: binnify failed for binsize %d.",
                    i,
                    new_binsize,
                )
                break

            # Check how many bins the new resolution has
            new_mc = MetaCooler(new_path)
            new_nbins = new_mc.shape[0]

            if new_nbins < 2:
                logger.info(
                    "Stopping at level %d: only %d bins at binsize %d.",
                    i,
                    new_nbins,
                    new_binsize,
                )
                break

            entries[str(new_binsize)] = new_path
            current_binsize = new_binsize
            current_uri = new_path
            current_nbins = new_nbins

        logger.info(
            "Creating mmetacool with %d resolution levels: %s",
            len(entries),
            sorted(entries.keys(), key=int),
        )

        # Create the mmetacool with resolution entries
        # We use "resolutions/<binsize>" as the group structure
        # but create_mmetacool puts them under /entries/<name>
        # So we use the binsize as the entry name
        create_mmetacool(Path(output), entries)

    logger.info("Written multi-resolution mmetacool to %s", output)
