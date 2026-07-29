"""Shared I/O utilities used across multiple metacooler modules."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd

from metacooler._version import __version__
from metacooler.schema import (
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
)

if TYPE_CHECKING:
    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")


def _str_to_hdf5(arr: np.ndarray) -> np.ndarray:
    """Convert a string/object NumPy array to fixed-length bytes for HDF5 storage.

    h5py cannot write arrays with ``dtype=object`` or pandas ``StringDtype``
    directly.  This function converts them to NumPy fixed-length byte-string
    arrays (``dtype="S<n>"``) suitable for ``h5py.Dataset.create_dataset()``.

    Args:
        arr: Input array.  Arrays with dtype ``"O"``, ``"U"``, ``"S"``, or
            pandas ``"string"`` are converted; all other dtypes are returned
            unchanged via ``numpy.asarray()``.

    Returns:
        A NumPy array with dtype ``S<n>`` where *n* is the length of the
        longest encoded string.  Non-string inputs are returned as-is.
    """
    if arr.dtype.kind in ("O", "U", "S") or str(arr.dtype) == "string":
        str_arr = ["" if pd.isna(v) else str(v) for v in arr]
        encoded = [s.encode("utf-8") for s in str_arr]
        # Size the fixed-length dtype by ENCODED BYTE length, not character
        # count: a multi-byte character (e.g. "é") is longer as UTF-8 than as a
        # Python str, so sizing by len(str) would silently truncate the bytes
        # (and could collide two distinct names).
        max_len = max((len(b) for b in encoded), default=1)
        return np.array(encoded, dtype=f"S{max_len}")
    return np.asarray(arr)


def check_output_path(path: str | Path, force: bool = False) -> None:
    """Refuse to overwrite an existing output file unless *force* is set.

    Single home for the "output exists" guard shared by the CLI commands and
    the ops layer, so every command reports the same message and every command
    that accepts ``--force`` actually honours it.

    Args:
        path: Output path about to be written.
        force: If True, remove an existing file instead of raising.

    Raises:
        FileExistsError: If *path* exists and *force* is False.
        IsADirectoryError: If *path* exists and is a directory.
    """
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        raise IsADirectoryError(f"Output path {str(path)!r} is a directory.")
    if not force:
        raise FileExistsError(
            f"Output file {str(path)!r} already exists. Use --force to overwrite."
        )
    path.unlink()


def resolve_count_dtype(dtype: np.dtype | str, *, context: str = "pixels") -> np.dtype:
    """Pick a lossless on-disk dtype for a ``pixels/count`` column being copied.

    Blindly casting to ``int32`` silently wraps counts written by
    :mod:`metacooler.ops.binnify` (which deliberately writes ``int64``) and
    floors float counts towards zero.  This preserves the source width instead.

    Args:
        dtype: dtype of the source ``count`` column.
        context: Short description of the caller, used in the float warning.

    Returns:
        ``int32`` for narrow integer sources, ``int64`` for wide or unsigned
        ones, and ``float64`` for float sources (with a warning).

    Raises:
        ValueError: If *dtype* is neither integer, boolean, nor float.
    """
    dtype = np.dtype(dtype)
    if dtype.kind == "f":
        logger.warning(
            "%s: source 'count' column is %s; keeping float64 to avoid "
            "truncating fractional counts.",
            context,
            dtype,
        )
        return np.dtype(np.float64)
    if dtype.kind == "b":
        return np.dtype(np.int32)
    if dtype.kind not in ("i", "u"):
        raise ValueError(f"Unsupported 'count' dtype {dtype!r} in {context}.")
    if dtype.itemsize >= 8 or (dtype.kind == "u" and dtype.itemsize >= 4):
        return np.dtype(np.int64)
    return np.dtype(np.int32)


def to_upper_triangle(
    bin1: np.ndarray, bin2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Reorder each ``(bin1, bin2)`` pair so that ``bin1 <= bin2``.

    metacooler stores pixels in ``symmetric-upper`` form, so every writer that
    remaps bin ids has to restore that invariant afterwards.  This used to be
    open-coded five times in two different spellings — a NumPy
    ``a[swap], b[swap] = b[swap], a[swap]`` dance (which relies on fancy
    indexing returning copies) and a pandas ``.loc[...] = ....values`` swap
    (which is only correct because of the trailing ``.values``: without it,
    pandas aligns on column *labels* and the assignment is a silent no-op).
    Both are correct today and both are one careless edit away from not being.
    ``minimum``/``maximum`` says the same thing with nothing to get wrong.

    Args:
        bin1: Row bin ids.
        bin2: Column bin ids, same shape as *bin1*.

    Returns:
        ``(lo, hi)`` — new arrays with ``lo <= hi`` elementwise.  Inputs are
        not modified.
    """
    return np.minimum(bin1, bin2), np.maximum(bin1, bin2)


def enforce_upper_triangle(
    pixels: pd.DataFrame, col1: str = "bin1_id", col2: str = "bin2_id"
) -> pd.DataFrame:
    """Reorder a pixel frame's bin-id columns in place so ``col1 <= col2``.

    DataFrame-shaped wrapper over :func:`to_upper_triangle`.

    Args:
        pixels: Pixel frame to fix up, modified in place.
        col1: Name of the row bin-id column.
        col2: Name of the column bin-id column.

    Returns:
        The same *pixels* object, for chaining.
    """
    lo, hi = to_upper_triangle(pixels[col1].to_numpy(), pixels[col2].to_numpy())
    pixels[col1] = lo
    pixels[col2] = hi
    return pixels


def chrom_codes(bins: pd.DataFrame) -> np.ndarray:
    """Return the integer contig id of every row of a *bins* frame.

    Single spelling for an idiom that was repeated ~12 times across the
    codebase as ``bins["chrom"].cat.codes.values``.  That form assumes cooler
    handed back a *categorical* ``chrom`` column and raises ``AttributeError``
    the moment it did not (a plain integer column read straight from HDF5, or a
    frame that has been through a groupby); this handles both.

    Args:
        bins: A bins DataFrame with a ``chrom`` column — either cooler's
            categorical form or plain integer contig ids.

    Returns:
        A fresh ``int64`` array of contig ids, safe for the caller to mutate.

    Raises:
        KeyError: If *bins* has no ``chrom`` column.
    """
    col = bins["chrom"]
    if isinstance(col.dtype, pd.CategoricalDtype):
        return col.cat.codes.to_numpy(dtype=np.int64)
    return col.to_numpy(dtype=np.int64)


def _finalize_metacool_file(
    path: Path,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    chrom_offset: np.ndarray | None,
    genome_offset: np.ndarray,
    *,
    contig_ordering: str = "genome-then-alpha",
    assembly_name: str | None = None,
) -> None:
    """Add metacooler-specific groups and attributes to a cooler file.

    Opens *path* in HDF5 ``"a"`` (append) mode and performs the following
    in-place modifications:

    - Replaces ``/chroms`` with a soft-link to ``/contigs`` for cooler
      compatibility.
    - Writes the ``/genomes`` and ``/contigs`` HDF5 groups (with
      fixed-length byte-string encoding).
    - Strips the enum dtype from ``/bins/chrom`` and rewrites it as plain
      ``int32``.
    - Writes ``indexes/chrom_offset`` and ``indexes/genome_offset``.
    - Sets all format attributes (``format``, ``metacooler-format``,
      ``ngenomes``, ``contig-ordering``, etc.).

    Args:
        path: Path to the cooler file to finalize (opened in append mode).
        genomes: ``/genomes`` DataFrame with columns ``name``, ``length``,
            ``n_contigs``, and ``n_bins``.
        contigs: ``/contigs`` DataFrame with columns ``name``, ``length``,
            and ``genome_id``.
        chrom_offset: Pre-computed ``chrom_offset`` array of length
            ``n_contigs + 1``.  If ``None``, it is computed on the fly
            from the ``bins/chrom`` dataset.
        genome_offset: Pre-computed ``genome_offset`` array of length
            ``n_genomes + 1``.
        contig_ordering: Value written to the ``contig-ordering`` HDF5
            attribute.  Typically ``"genome-then-alpha"``,
            ``"pairs-header"``, or ``"genomeid-then-alpha"``.
        assembly_name: Optional assembly name written to the
            ``assembly-name`` attribute.  Omitted if ``None``.
    """
    with h5py.File(path, "a") as f:
        # Recreate chroms as a soft link to contigs (cooler compat)
        if "chroms" in f:
            del f["chroms"]
        f["chroms"] = h5py.SoftLink("contigs")

        # Write /genomes
        if "genomes" in f:
            del f["genomes"]
        g = f.create_group("genomes", track_order=True)
        for col in genomes.columns:
            g.create_dataset(col, data=_str_to_hdf5(genomes[col].values))

        # Write /contigs
        if "contigs" in f:
            del f["contigs"]
        c = f.create_group("contigs", track_order=True)
        for col in contigs.columns:
            c.create_dataset(col, data=_str_to_hdf5(contigs[col].values))

        # Strip cooler's enum dtype from bins/chrom so it stores plain int32
        chrom_data = f["bins"]["chrom"][:]
        del f["bins"]["chrom"]
        f["bins"].create_dataset("chrom", data=chrom_data, dtype=np.int32)

        # Compute chrom_offset from bins if not provided
        if chrom_offset is None:
            n_contigs = len(contigs)
            contig_counts = np.bincount(
                chrom_data.astype(np.int64), minlength=n_contigs
            )
            chrom_offset = np.zeros(n_contigs + 1, dtype=np.int64)
            chrom_offset[1:] = np.cumsum(contig_counts)

        # Write index offsets
        idx = f.require_group("indexes")
        for key, data in (
            ("chrom_offset", chrom_offset),
            ("genome_offset", genome_offset),
        ):
            if key in idx:
                del idx[key]
            idx.create_dataset(key, data=data.astype(np.int64))

        # Format attributes
        f.attrs["format"] = COOLER_FORMAT
        f.attrs["format-url"] = COOLER_FORMAT_URL
        f.attrs["metacooler-format"] = METACOOLER_FORMAT
        f.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        f.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        f.attrs["generated-by"] = f"metacooler-{__version__}"
        f.attrs["ngenomes"] = len(genomes)
        f.attrs["contig-ordering"] = contig_ordering
        if assembly_name is not None:
            f.attrs["assembly-name"] = assembly_name


def _filter_by_genomes_and_contigs(
    mc: MetaCooler,
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
) -> set[int]:
    """Return the set of contig IDs matching the genome/contig filters.

    Args:
        mc: Open ``MetaCooler`` instance.
        genomes: Genome names to include (all their contigs are selected).
        contigs: Contig names to include directly.

    Returns:
        Set of integer contig IDs matching the combined filter.

    Raises:
        KeyError: If any genome or contig name is not found in *mc*.
    """
    all_genomes = mc.genomes
    all_contigs = mc.contigs
    selected: set[int] = set()
    known_genomes = set(all_genomes["name"])
    known_contigs = set(all_contigs["name"])

    if genomes:
        for g in genomes:
            if g not in known_genomes:
                available = all_genomes["name"].tolist()
                raise KeyError(
                    f"Genome {g!r} not found. Available ({len(available)} total): {available[:20]}..."
                )
            gid = int(all_genomes.loc[all_genomes["name"] == g].index[0])
            cids = all_contigs.index[all_contigs["genome_id"] == gid].tolist()
            selected.update(cids)

    if contigs:
        for c in contigs:
            if c not in known_contigs:
                available = all_contigs["name"].tolist()
                raise KeyError(
                    f"Contig {c!r} not found. Available ({len(available)} total): {available[:20]}..."
                )
            cid = int(all_contigs.loc[all_contigs["name"] == c].index[0])
            selected.add(cid)

    return selected
