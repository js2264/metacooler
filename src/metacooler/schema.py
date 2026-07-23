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
            return all(
                g in f for g in ("genomes", "contigs", "bins", "pixels", "indexes")
            )
    except Exception:
        return False


def validate(uri: str | Path) -> None:
    """Validate that *uri* is a well-formed ``.metacool`` file.

    Performs a structural integrity check: required groups, columns,
    indexes, and cross-count consistency.

    Args:
        uri: Path to the ``.metacool`` file, or a URI of the form
            ``file.mmetacool::entries/entry`` to validate a sub-group inside
            an mmetacool file.

    Raises:
        MetaCoolerFormatError: If any required group, dataset, or attribute
            is missing or malformed.
    """
    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(str(uri))
    if not file_path.exists():
        raise MetaCoolerFormatError(f"File not found: {file_path}")

    with h5py.File(file_path, "r") as _froot:
        f = _froot[group_path.lstrip("/")] if group_path != "/" else _froot
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
        if n_contigs == 0:
            raise MetaCoolerFormatError("/contigs table is empty")
        if n_genomes == 0:
            raise MetaCoolerFormatError("/genomes table is empty")

        # --- /chroms must be a soft link to /contigs (cooler compatibility) ---
        link = f.get("chroms", getlink=True)
        if not isinstance(link, h5py.SoftLink):
            raise MetaCoolerFormatError(
                "/chroms must be a soft link to /contigs for cooler compatibility"
            )
        if link.path.lstrip("/").rsplit("/", 1)[-1] != "contigs":
            raise MetaCoolerFormatError(
                f"/chroms soft link must target 'contigs', got {link.path!r}"
            )

        # --- name uniqueness ---
        contig_names = np.array(c["name"])
        genome_names = np.array(g["name"])
        if len(np.unique(contig_names)) != n_contigs:
            raise MetaCoolerFormatError("contig names are not unique")
        if len(np.unique(genome_names)) != n_genomes:
            raise MetaCoolerFormatError("genome names are not unique")

        # --- attribute consistency ---
        nchroms_attr = f.attrs.get("nchroms")
        if nchroms_attr is not None and int(nchroms_attr) != n_contigs:
            raise MetaCoolerFormatError(
                f"attrs['nchroms'] ({int(nchroms_attr)}) != len(/contigs) ({n_contigs})"
            )

        # --- contig lengths must be positive ---
        if (np.array(c["length"]) <= 0).any():
            raise MetaCoolerFormatError("contig lengths must be positive (> 0)")

        # --- bins/chrom range and ordering ---
        bins_chrom = np.array(b["chrom"]).astype(np.int64)
        if n_bins > 0:
            if bins_chrom.min() < 0 or bins_chrom.max() >= n_contigs:
                raise MetaCoolerFormatError(
                    f"bins/chrom values out of range [0, {n_contigs})"
                )
            if np.any(np.diff(bins_chrom) < 0):
                raise MetaCoolerFormatError("bins are not sorted by contig (chrom)")

        # --- index dataset lengths ---
        if len(idx["chrom_offset"]) != n_contigs + 1:
            raise MetaCoolerFormatError(
                f"chrom_offset length ({len(idx['chrom_offset'])}) != n_contigs+1 ({n_contigs + 1})"
            )
        if len(idx["bin1_offset"]) != n_bins + 1:
            raise MetaCoolerFormatError(
                f"bin1_offset length ({len(idx['bin1_offset'])}) != n_bins+1 ({n_bins + 1})"
            )
        if "genome_offset" not in idx:
            raise MetaCoolerFormatError("indexes/genome_offset is missing")
        if len(idx["genome_offset"]) != n_genomes + 1:
            raise MetaCoolerFormatError(
                f"genome_offset length ({len(idx['genome_offset'])}) != n_genomes+1 ({n_genomes + 1})"
            )

        # --- chrom_offset must be the bin-count CSR over the (sorted) bins ---
        chrom_offset = np.array(idx["chrom_offset"]).astype(np.int64)
        expected_chrom = np.zeros(n_contigs + 1, dtype=np.int64)
        if n_bins > 0:
            expected_chrom[1:] = np.cumsum(np.bincount(bins_chrom, minlength=n_contigs))
        if not np.array_equal(chrom_offset, expected_chrom):
            raise MetaCoolerFormatError(
                "chrom_offset is inconsistent with the bins table"
            )

        # --- bin1_offset endpoints must bracket the pixel table ---
        bin1_offset = np.array(idx["bin1_offset"]).astype(np.int64)
        if int(bin1_offset[0]) != 0 or int(bin1_offset[-1]) != n_pixels:
            raise MetaCoolerFormatError(
                f"bin1_offset endpoints ({int(bin1_offset[0])}, {int(bin1_offset[-1])}) "
                f"do not bracket the pixel table (0, {n_pixels})"
            )

        # --- contig genome_ids must be in range ---
        genome_ids = np.array(c["genome_id"]).astype(np.int64)
        if genome_ids.min() < 0 or genome_ids.max() >= n_genomes:
            raise MetaCoolerFormatError(
                f"genome_id values out of range [0, {n_genomes})"
            )

        # --- genome_offset must be the contig-count CSR over genomes ---
        genome_offset = np.array(idx["genome_offset"]).astype(np.int64)
        genome_contig_counts = np.bincount(genome_ids, minlength=n_genomes)
        expected_genome = np.zeros(n_genomes + 1, dtype=np.int64)
        expected_genome[1:] = np.cumsum(genome_contig_counts)
        if not np.array_equal(genome_offset, expected_genome):
            raise MetaCoolerFormatError(
                f"genome_offset {genome_offset.tolist()} is inconsistent with the "
                f"contig genome_id counts; expected {expected_genome.tolist()}"
            )
        # For genome-grouped orderings, contigs must be physically laid out in
        # ascending genome_id order (otherwise genome_offset — though numerically
        # consistent with the counts — would not delimit each genome's contigs).
        ordering = f.attrs.get("contig-ordering")
        if ordering in ("genome-then-length", "genome-then-alpha") and np.any(
            np.diff(genome_ids) < 0
        ):
            raise MetaCoolerFormatError(
                f"contig-ordering={ordering!r} requires contigs grouped by ascending "
                "genome_id, but the genome_id sequence is not non-decreasing, so "
                "genome_offset does not match the physical contig layout"
            )

        # --- genomes/n_contigs and n_bins must match reality ---
        if not np.array_equal(
            np.array(g["n_contigs"]).astype(np.int64), genome_contig_counts
        ):
            raise MetaCoolerFormatError(
                "genomes/n_contigs does not match the contigs table"
            )
        if n_bins > 0:
            genome_bin_counts = np.bincount(genome_ids[bins_chrom], minlength=n_genomes)
        else:
            genome_bin_counts = np.zeros(n_genomes, dtype=np.int64)
        if not np.array_equal(
            np.array(g["n_bins"]).astype(np.int64), genome_bin_counts
        ):
            raise MetaCoolerFormatError("genomes/n_bins does not match the bins table")

        # --- pixel ordering, range, and upper-triangle (sampled first/last rows) ---
        if n_pixels > 0:
            sample_size = min(1000, n_pixels)
            for sl in (
                slice(0, sample_size),
                slice(max(0, n_pixels - sample_size), n_pixels),
            ):
                b1 = np.array(p["bin1_id"][sl])
                b2 = np.array(p["bin2_id"][sl])
                if b1.min() < 0 or b2.max() >= n_bins:
                    raise MetaCoolerFormatError(
                        f"pixel bin ids out of range [0, {n_bins})"
                    )
                if np.any(b1 > b2):
                    raise MetaCoolerFormatError(
                        "pixels are not upper-triangular (bin1_id <= bin2_id)"
                    )
                if len(b1) > 1 and not np.all(
                    (b1[:-1] < b1[1:]) | ((b1[:-1] == b1[1:]) & (b2[:-1] <= b2[1:]))
                ):
                    raise MetaCoolerFormatError(
                        "pixels are not sorted by (bin1_id, bin2_id)"
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
            # Decode fixed-length byte strings to str (like cooler). Use an
            # explicit UTF-8 decode: ``astype("U")`` decodes as ASCII and raises
            # UnicodeDecodeError on any non-ASCII byte (e.g. a "é" in a name),
            # which would break every read of a file that metacooler itself
            # wrote.
            if data.dtype.kind == "S":
                data = np.char.decode(data, "utf-8")
            elif data.dtype == object and len(data) > 0 and isinstance(data[0], bytes):
                data = np.array([x.decode("utf-8") for x in data])
            columns[key] = data
    return pd.DataFrame(columns)
