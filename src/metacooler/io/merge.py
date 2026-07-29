"""Merge operations for metacool files."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import h5py
import numpy as np

from metacooler._version import __version__
from metacooler.schema import (
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
    validate,
)

if TYPE_CHECKING:
    from collections.abc import Iterable


def _hash_group(f: h5py.File, group: str) -> str:
    """Return an MD5 hex digest of all datasets in a group.

    Used to compare ``contigs`` and ``genomes`` tables across files for
    merge compatibility checking.

    Args:
        f: Open ``h5py.File`` (or any HDF5 group-like object).
        group: HDF5 group path to hash (e.g. ``"contigs"`` or ``"genomes"``).

    Returns:
        Lowercase hex MD5 digest string of all dataset values, sorted by
        key name for determinism.
    """
    h = hashlib.md5(usedforsecurity=False)
    g = f[group]
    for key in sorted(g.keys()):
        arr = g[key][:]
        # Variable-length string datasets read back as object arrays of Python
        # strings; encode them to bytes for a stable hash.
        if arr.dtype == object or str(arr.dtype).startswith("S"):
            arr = np.array(
                [v.encode() if isinstance(v, str) else bytes(v) for v in arr.ravel()]
            )
        h.update(arr.tobytes())
    return h.hexdigest()


def merge(
    inputs: Iterable[str | Path],
    output: str | Path,
) -> None:
    """Merge multiple ``.metacool`` files with identical axes into one.

    All input files must have identical ``contigs`` and ``genomes`` tables
    and the same ``contig-ordering`` attribute.  Pixels from all inputs are
    concatenated, re-sorted externally by ``cooler.merge_coolers()``, and
    saved to *output*.  Metacooler-specific groups (``genomes``, ``contigs``,
    ``genome_offset``, format attributes) are then copied from the first
    input file.

    Args:
        inputs: Iterable of paths to the input ``.metacool`` files.  At
            least two files are required.
        output: Destination path for the merged ``.metacool`` file.

    Raises:
        ValueError: If fewer than two inputs are provided, if the
            ``contigs`` or ``genomes`` tables differ across files, or if
            the ``contig-ordering`` attribute is inconsistent.
        MetaCoolerFormatError: If any input file fails the metacooler
            format validation check.
    """
    inputs = [Path(p) for p in inputs]
    if len(inputs) < 2:
        raise ValueError("At least two input files required for merge")

    # Validate and compare axes
    first = inputs[0]
    validate(first)

    with h5py.File(first, "r") as f0:
        contigs_hash = _hash_group(f0, "contigs")
        genomes_hash = _hash_group(f0, "genomes")
        ordering = f0.attrs.get("contig-ordering")

    for inp in inputs[1:]:
        validate(inp)
        with h5py.File(inp, "r") as f:
            if _hash_group(f, "contigs") != contigs_hash:
                raise ValueError(f"Contigs table mismatch: {inp}")
            if _hash_group(f, "genomes") != genomes_hash:
                raise ValueError(f"Genomes table mismatch: {inp}")
            if f.attrs.get("contig-ordering") != ordering:
                raise ValueError(f"contig-ordering mismatch: {inp}")

    # Use cooler's merge logic on the pixel tables
    uris = [str(p) for p in inputs]
    # Use cooler's merge logic on the pixel tables.  ``dtypes`` is mandatory:
    # without it the summed counts inherit the int32 of the inputs and wrap.
    cooler.merge_coolers(
        str(output), uris, mergebuf=int(2e7), dtypes={"count": np.int64}
    )

    # Post-process: add back metacooler-specific groups
    with h5py.File(first, "r") as fsrc, h5py.File(output, "a") as fdst:
        # Copy genomes
        if "genomes" in fdst:
            del fdst["genomes"]
        fdst.copy(fsrc["genomes"], "genomes")

        # Copy contigs (overwrites the chroms-compatible one)
        if "contigs" in fdst:
            del fdst["contigs"]
        fdst.copy(fsrc["contigs"], "contigs")

        # Soft link chroms -> contigs
        if "chroms" in fdst:
            del fdst["chroms"]
        fdst["chroms"] = h5py.SoftLink("contigs")

        # Copy genome_offset
        if "genome_offset" in fdst["indexes"]:
            del fdst["indexes"]["genome_offset"]
        fdst["indexes"].copy(fsrc["indexes"]["genome_offset"], "genome_offset")

        # Strip cooler's enum dtype from bins/chrom — cooler.merge_coolers()
        # always writes bins/chrom as an HDF5 enum type, which is inconsistent
        # with the metacooler schema that stores it as plain int32.
        chrom_data = fdst["bins"]["chrom"][:]
        del fdst["bins"]["chrom"]
        fdst["bins"].create_dataset("chrom", data=chrom_data, dtype=np.int32)

        # Set attrs
        fdst.attrs["format"] = COOLER_FORMAT
        fdst.attrs["format-url"] = COOLER_FORMAT_URL
        fdst.attrs["metacooler-format"] = METACOOLER_FORMAT
        fdst.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        fdst.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        fdst.attrs["generated-by"] = f"metacooler-{__version__}"
        fdst.attrs["ngenomes"] = len(fsrc["genomes"]["name"])
        fdst.attrs["contig-ordering"] = ordering

    validate(output)
