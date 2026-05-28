"""HDF5 format constants, dtypes, and validators for metacooler."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

import cooler
import h5py
import numpy as np
import pandas as pd

logger = logging.getLogger("metacooler")

METACOOLER_FORMAT = "HDF5::MetaCooler"
METACOOLER_FORMAT_VERSION = 1
METACOOLER_FORMAT_URL = "https://github.com/js2264/metacooler"
MMETACOOLER_FORMAT = "HDF5::MMetaCooler"
MMETACOOLER_FORMAT_VERSION = 1
COOLER_FORMAT: Literal["HDF5::Cooler"] = cooler.create._constants.MAGIC
COOLER_FORMAT_URL: Literal["https://github.com/open2c/cooler"] = (
    cooler.create._constants.URL
)

# Dense matrix size guard: refuse to allocate matrices larger than this
DENSE_MATRIX_BIN_LIMIT = 50_000

# Dtypes matching cooler conventions
CHROMID_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32

# Column dtypes for metacooler-specific tables
GENOMES_DTYPE = {
    "name": "object",
    "length": "i8",
    "n_contigs": "i4",
    "n_bins": "i8",
}

CONTIGS_DTYPE = {
    "name": "object",
    "length": "i4",
    "genome_id": "i4",
}


class MetaCoolerFormatError(ValueError):
    """Raised when a file does not conform to the metacooler format."""


def is_metacooler(uri: str | Path) -> bool:
    """Return ``True`` if *uri* points to a valid ``.metacool`` file.

    Checks for the ``metacooler-format`` HDF5 attribute equal to
    ``HDF5::MetaCooler`` and the presence of the ``/genomes`` group.
    Note that ``cooler.fileops.is_cooler()`` will *also* return ``True``
    for a metacool file because of the soft-link ``chroms → contigs``
    and the ``format = "HDF5::Cooler"`` attribute.

    Args:
        uri: Path to the file to test.

    Returns:
        ``True`` if the file is a valid metacool file, ``False`` otherwise
        (including if the file does not exist or cannot be opened).
    """
    path = Path(uri)
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as f:
            fmt = f.attrs.get("metacooler-format")
            if fmt != METACOOLER_FORMAT:
                return False
            return "genomes" in f
    except Exception:
        return False


def validate(uri: str | Path) -> None:
    """Validate that *uri* is a well-formed ``.metacool`` file.

    Performs a structural integrity check: required groups, columns,
    indexes, and cross-count consistency.

    Args:
        uri: Path to the ``.metacool`` file to validate.

    Raises:
        MetaCoolerFormatError: If any required group, dataset, or attribute
            is missing or malformed.
    """
    path = Path(uri)
    if not path.exists():
        raise MetaCoolerFormatError(f"File not found: {path}")

    with h5py.File(path, "r") as f:
        # --- format attributes ---
        fmt = f.attrs.get("format")
        if fmt != COOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected format='{COOLER_FORMAT}', got {fmt!r}"
            )

        mfmt = f.attrs.get("metacooler-format")
        if mfmt != METACOOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected metacooler-format='{METACOOLER_FORMAT}', got {mfmt!r}"
            )

        # --- required groups ---
        for grp in ("genomes", "contigs", "bins", "pixels", "indexes"):
            if grp not in f:
                raise MetaCoolerFormatError(f"Missing required group: /{grp}")

        # --- genomes table ---
        g = f["genomes"]
        for col in ("name", "length", "n_contigs", "n_bins"):
            if col not in g:
                raise MetaCoolerFormatError(f"Missing column in /genomes: {col}")

        # --- contigs table ---
        c = f["contigs"]
        for col in ("name", "length", "genome_id"):
            if col not in c:
                raise MetaCoolerFormatError(f"Missing column in /contigs: {col}")

        # --- bins table ---
        b = f["bins"]
        for col in ("chrom", "start", "end"):
            if col not in b:
                raise MetaCoolerFormatError(f"Missing column in /bins: {col}")

        # --- pixels table ---
        p = f["pixels"]
        for col in ("bin1_id", "bin2_id", "count"):
            if col not in p:
                raise MetaCoolerFormatError(f"Missing column in /pixels: {col}")

        # --- indexes ---
        idx = f["indexes"]
        for ds in ("chrom_offset", "bin1_offset"):
            if ds not in idx:
                raise MetaCoolerFormatError(f"Missing index dataset: /indexes/{ds}")

        # --- cross-check counts ---
        n_genomes = len(g["name"])
        n_contigs = len(c["name"])
        n_bins = len(b["chrom"])
        n_pixels = len(p["bin1_id"])
        if n_pixels == 0:
            logger.warning("Empty pixels table")

        if len(idx["chrom_offset"]) != n_contigs + 1:
            raise MetaCoolerFormatError(
                f"chrom_offset length ({len(idx['chrom_offset'])}) != n_contigs+1 ({n_contigs + 1})"
            )
        if len(idx["bin1_offset"]) != n_bins + 1:
            raise MetaCoolerFormatError(
                f"bin1_offset length ({len(idx['bin1_offset'])}) != n_bins+1 ({n_bins + 1})"
            )

        # genome_offset is optional but if present must be consistent
        if "genome_offset" in idx and len(idx["genome_offset"]) != n_genomes + 1:
            raise MetaCoolerFormatError(
                f"genome_offset length ({len(idx['genome_offset'])}) != n_genomes+1 ({n_genomes + 1})"
            )

        # contig genome_ids must be in range
        genome_ids = np.array(c["genome_id"])
        if genome_ids.min() < 0 or genome_ids.max() >= n_genomes:
            raise MetaCoolerFormatError(
                f"genome_id values out of range [0, {n_genomes})"
            )


# ------------------------------------------------------------------
# MMetaCooler helpers
# ------------------------------------------------------------------


def is_mmetacooler(uri: str | Path) -> bool:
    """Return ``True`` if *uri* points to a valid ``.mmetacool`` file.

    Checks for the ``format`` attribute equal to ``HDF5::MMetaCooler`` and
    the presence of the top-level ``entries/`` group.

    Args:
        uri: Path to the file to test.

    Returns:
        ``True`` if the file is a valid mmetacool collection, ``False``
        otherwise (including if the file does not exist or cannot be
        opened).
    """
    path = Path(uri)
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as f:
            fmt = f.attrs.get("format")
            if fmt != MMETACOOLER_FORMAT:
                return False
            return "entries" in f
    except Exception:
        return False


def validate_mmetacooler(uri: str | Path) -> None:
    """Validate that *uri* is a well-formed ``.mmetacool`` file.

    Checks top-level attributes, the ``entries/`` group, and verifies
    that each entry contains the required metacooler tables (``bins``,
    ``pixels``, ``contigs``, ``genomes``, ``indexes``).

    Args:
        uri: Path to the ``.mmetacool`` file to validate.

    Raises:
        MetaCoolerFormatError: If the structure is invalid.
    """
    path = Path(uri)
    if not path.exists():
        raise MetaCoolerFormatError(f"File not found: {path}")

    with h5py.File(path, "r") as f:
        fmt = f.attrs.get("format")
        if fmt != MMETACOOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected format='{MMETACOOLER_FORMAT}', got {fmt!r}"
            )

        fv = f.attrs.get("format-version")
        if fv is None:
            raise MetaCoolerFormatError("Missing 'format-version' attribute.")

        if "entries" not in f:
            raise MetaCoolerFormatError("Missing required key: /entries")

        grp = f["entries"]
        for name in grp:
            entry = grp[name]
            # Each entry must have the core metacooler tables
            for tbl in ("bins", "pixels", "contigs", "genomes", "indexes"):
                if tbl not in entry:
                    raise MetaCoolerFormatError(
                        f"Entry {name!r} missing required group: {tbl}"
                    )


def list_metacoolers(uri: str | Path) -> list[str]:
    """Return HDF5 group paths of all metacooler entries in an mmetacool file.

    Returns paths of the form ``/entries/<name>`` for each entry in the
    ``/entries`` group.

    Args:
        uri: Path to the ``.mmetacool`` file.

    Returns:
        List of HDF5 group paths (e.g. ``["/entries/sample1",
        "/entries/sample2"]``).  Returns an empty list if the file has no
        ``entries`` group.
    """
    path = Path(uri)
    with h5py.File(path, "r") as f:
        if "entries" not in f:
            return []
        return [f"/entries/{name}" for name in f["entries"]]


def read_table(uri: str | Path, group: str) -> pd.DataFrame:
    """Read an HDF5 table group as a pandas DataFrame.

    Matches cooler behaviour: fixed-length byte strings are decoded to
    unicode, and column order follows the HDF5 dataset creation order.

    Args:
        uri: Path to the HDF5 file.
        group: HDF5 group path containing the table datasets (e.g.
            ``"genomes"`` or ``"entries/sample1/contigs"``).

    Returns:
        DataFrame where each HDF5 dataset in *group* becomes a column.
        Byte-string columns are decoded to ``str``; all other dtypes are
        returned as-is.
    """
    path = Path(uri)
    with h5py.File(path, "r", track_order=True) as f:
        g = f[group]
        columns = {}
        for key in g:
            ds = g[key]
            data = ds[:]
            # Decode fixed-length byte strings to str (like cooler)
            if data.dtype.kind == "S":
                data = data.astype("U")
            elif data.dtype == object and len(data) > 0 and isinstance(data[0], bytes):
                data = np.array([x.decode("utf-8") for x in data])
            columns[key] = data
    return pd.DataFrame(columns)
