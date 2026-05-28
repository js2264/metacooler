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
        str_arr = [
            "" if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v)
            for v in arr
        ]
        max_len = max((len(s) for s in str_arr), default=1)
        return np.array([s.encode("utf-8") for s in str_arr], dtype=f"S{max_len}")
    return np.asarray(arr)


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
        if assembly_name:
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

    if genomes:
        for g in genomes:
            if g not in all_genomes["name"].values:
                available = all_genomes["name"].tolist()
                raise KeyError(f"Genome {g!r} not found. Available: {available}")
            gid = int(all_genomes.loc[all_genomes["name"] == g].index[0])
            cids = all_contigs.index[all_contigs["genome_id"] == gid].tolist()
            selected.update(cids)

    if contigs:
        for c in contigs:
            if c not in all_contigs["name"].values:
                available = all_contigs["name"].tolist()
                raise KeyError(f"Contig {c!r} not found. Available: {available}")
            cid = int(all_contigs.loc[all_contigs["name"] == c].index[0])
            selected.add(cid)

    return selected
