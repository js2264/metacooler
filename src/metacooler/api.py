"""Primary Python API for metacool files."""

from __future__ import annotations

import contextlib
import functools
from pathlib import Path
from typing import TYPE_CHECKING, cast

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.schema import (
    read_table,
    validate,
)

if TYPE_CHECKING:
    from typing import Literal


class MetaCooler:
    """Primary interface to a ``.metacool`` file.

    Internally holds a ``cooler.Cooler`` instance and delegates standard
    read operations to it. metacooler-specific accessors (``genomes``,
    ``contigs``, ``aggregate_matrix``, etc.) are layered on top.

    Supports URI syntax for mmetacool files:
    ``MetaCooler("file.mmetacool::entries/entry1")`` opens a specific
    entry within an mmetacool collection.
    """

    def __init__(self, store: str | Path):
        """Open a ``.metacool`` (or mmetacool entry) file.

        Args:
            store: File path, optionally with a URI group suffix
                (``file.mmetacool::entries/entry_name``).

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If a group path is given but the file is not an
                mmetacool, or if a plain path points to an mmetacool file
                without specifying an entry.
            KeyError: If the requested HDF5 group does not exist inside
                the mmetacool file.
        """

        from metacooler.uri import parse_uri

        raw = str(store)
        file_path, group_path = parse_uri(raw)
        self._store = file_path
        self._group_path = group_path

        if not self._store.exists():
            raise FileNotFoundError(f"File not found: {self._store}")

        if group_path != "/":
            # Opening a sub-group inside an mmetacool file
            from metacooler.schema import is_mmetacooler

            if not is_mmetacooler(self._store):
                raise ValueError(
                    f"URI group path specified but {self._store} is not an mmetacool file."
                )
            # Validate that the group exists
            with h5py.File(self._store, "r") as f:
                if group_path.lstrip("/") not in f:
                    raise KeyError(f"Group {group_path!r} not found in {self._store}")
            # Use cooler URI to open the specific group
            cooler_uri = f"{self._store}::{group_path}"
            self._clr = cooler.Cooler(cooler_uri)
        else:
            from metacooler.schema import is_mmetacooler, list_metacoolers

            if is_mmetacooler(self._store):
                entries = list_metacoolers(self._store)
                raise KeyError(
                    f"No metacool entry at root of {self._store}. "
                    f"Metacool entries found: {entries}. "
                    f"Use '::' to specify a group path, e.g. "
                    f"{self._store}::{entries[0] if entries else '/entries/<name>'}."
                )
            validate(self._store)
            self._clr = cooler.Cooler(str(self._store))

    @property
    def filename(self) -> str:
        """The file path as a string."""
        return str(self._store)

    @property
    def root(self) -> str:
        """The HDF5 root path (always ``/`` for metacooler)."""
        return "/"

    @property
    def store(self) -> Path:
        """The file path as a ``pathlib.Path``."""
        return self._store

    @property
    def uri(self) -> str:
        """The cooler-style URI (``path::group_path``)."""
        if self._group_path == "/":
            return f"{self._store}::/"
        return f"{self._store}::{self._group_path}"

    # ------------------------------------------------------------------
    # Delegated to cooler.Cooler (pass-through)
    # ------------------------------------------------------------------
    @functools.cached_property
    def _bins_cache(self) -> pd.DataFrame:
        """Cached bin table — bins are immutable once written."""
        return cast(pd.DataFrame, self._clr.bins()[:])

    def bins(self, **kw) -> pd.DataFrame:
        """Return the bins table as a DataFrame.

        Delegates to ``cooler.Cooler.bins()``.  The result is cached on
        the first call (when no keyword arguments are given).

        Args:
            **kw: Keyword arguments forwarded to ``cooler.Cooler.bins()``.

        Returns:
            DataFrame with columns ``chrom`` (categorical contig ID),
            ``start``, and ``end`` (int32), plus any extra columns written
            at creation time.
        """
        if kw:
            return cast(pd.DataFrame, self._clr.bins(**kw)[:])
        return self._bins_cache

    def pixels(self, join: bool = False, **kw) -> pd.DataFrame:
        """Return the pixels table as a DataFrame.

        Delegates to ``cooler.Cooler.pixels()``.

        Args:
            join: If ``True``, annotate each pixel with its bin coordinates
                (adds ``chrom1``, ``start1``, ``end1``, ``chrom2``, ``start2``,
                ``end2`` columns).
            **kw: Additional keyword arguments forwarded to
                ``cooler.Cooler.pixels()``.

        Returns:
            DataFrame with at least ``bin1_id``, ``bin2_id``, and ``count``
            columns.  When *join* is ``True``, bin coordinate columns are
            prepended.
        """
        return cast(pd.DataFrame, self._clr.pixels(join=join, **kw)[:])

    def matrix(self, **kw):
        """Return a cooler matrix selector for slicing the contact matrix.

        Delegates to ``cooler.Cooler.matrix()``.

        Args:
            **kw: Keyword arguments forwarded to ``cooler.Cooler.matrix()``.
                Commonly ``balance=False`` to retrieve raw counts.

        Returns:
            A ``cooler.core.RangeSelector2D`` that can be sliced with bin
            index ranges or region strings.
        """
        return self._clr.matrix(**kw)

    def chroms(self, **kw) -> pd.DataFrame:
        """Return the chroms table as a DataFrame (alias for contigs).

        In a metacool file ``chroms`` is a soft-link to the ``contigs``
        group for cooler compatibility.  This method delegates to
        ``cooler.Cooler.chroms()``.

        Args:
            **kw: Keyword arguments forwarded to ``cooler.Cooler.chroms()``.

        Returns:
            DataFrame with at least ``name`` and ``length`` columns.
        """
        return cast(pd.DataFrame, self._clr.chroms(**kw)[:])

    def extent(self, region: str) -> tuple[int, int]:
        """Return the half-open bin index range for a region string.

        Delegates to ``cooler.Cooler.extent()``.

        Args:
            region: Region string in the format ``contig`` or
                ``contig:start-end`` (1-based, fully closed).

        Returns:
            ``(lo, hi)`` — half-open bin indices such that bins ``[lo, hi)``
            overlap the requested region.
        """
        return self._clr.extent(region)

    def offset(self, region: str) -> int:
        """Return the first bin index for a region string.

        Delegates to ``cooler.Cooler.offset()``.

        Args:
            region: Region string in the format ``contig`` or
                ``contig:start-end`` (1-based, fully closed).

        Returns:
            0-based integer index of the first bin that overlaps *region*.
        """
        return self._clr.offset(region)

    @property
    def binsize(self) -> int | None:
        """Uniform bin size in base pairs, or ``None`` if bins are variable-width."""
        return self._clr.binsize

    @property
    def chromnames(self) -> list[str]:
        """List of contig names in the order they appear in the file."""
        return self._clr.chromnames

    @property
    def chromsizes(self) -> pd.Series:
        """Series mapping each contig name to its length in base pairs."""
        return self._clr.chromsizes

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the square contact matrix ``(n_bins, n_bins)``."""
        return self._clr.shape

    @property
    def info(self) -> dict:
        """Merged file-level attributes in normalised metacooler format.

        Reads the root HDF5 attributes via ``cooler.Cooler.info``, then
        renames metacooler-specific keys (e.g. ``metacooler-format`` →
        ``format``) and casts numeric fields to plain Python ``int``.  The
        ``metadata`` key is removed from the result — use
        :attr:`metadata` instead.

        Returns:
            Dictionary of file-level attributes including ``format``,
            ``format-version``, ``bin-size``, ``nbins``, ``nnz``, ``sum``,
            ``ngenomes``, ``ncontigs``, ``ordering``, and ``assembly``.
        """
        base = dict(self._clr.info)
        # Promote metacooler-specific attrs over their cooler equivalents.
        # Use .pop(key, None) so that files written by older metacooler
        # versions (or partially-written files) don't raise KeyError.
        base["format"] = base.pop("metacooler-format", None)
        base["format-url"] = base.pop("metacooler-format-url", None)
        base["format-version"] = base.pop("metacooler-format-version", None)
        base["ordering"] = base.pop("contig-ordering", None)
        base["ncontigs"] = base.pop("nchroms", None)
        # Remove metadata field if present (accessed via .metadata property)
        base.pop("metadata", None)
        # Cast numeric fields to int (cooler may return numpy scalars or strings)
        for _key in (
            "bin-size",
            "format-version",
            "nbins",
            "nnz",
            "sum",
            "ngenomes",
            "ncontigs",
        ):
            if _key in base:
                with contextlib.suppress(ValueError, TypeError):
                    base[_key] = int(base[_key])
        return base

    @property
    def metadata(self) -> dict:
        """Contents of the optional ``metadata`` group as a nested dict.

        Returns an empty dict if no metadata was stored at creation time.

        Returns:
            Nested dictionary of user-defined metadata values, or ``{}``
            if the file contains no metadata.
        """
        base = dict(self._clr.info)
        meta = base.pop("metadata", None)
        if meta is None:
            return {}
        if isinstance(meta, dict):
            return meta
        return {}

    # ------------------------------------------------------------------
    # metacooler-specific properties
    # ------------------------------------------------------------------
    @functools.cached_property
    def genomes(self) -> pd.DataFrame:
        """The ``/genomes`` table as a DataFrame.

        Loaded once and then cached.  The table has at least the columns
        ``name`` (str), ``length`` (int), and ``n_contigs`` (int).

        Returns:
            DataFrame indexed 0 … n_genomes-1 with one row per genome.
        """
        if self._group_path != "/":
            grp = self._group_path.lstrip("/") + "/genomes"
            return read_table(self._store, grp)
        return read_table(self._store, "genomes")

    @functools.cached_property
    def contigs(self) -> pd.DataFrame:
        """The ``/contigs`` table as a DataFrame.

        Loaded once and then cached.  The table has at least the columns
        ``name`` (str), ``length`` (int), and ``genome_id`` (int, 0-based
        index into the :attr:`genomes` table).

        Returns:
            DataFrame indexed 0 … n_contigs-1 with one row per contig.
        """
        if self._group_path != "/":
            grp = self._group_path.lstrip("/") + "/contigs"
            return read_table(self._store, grp)
        return read_table(self._store, "contigs")

    # ------------------------------------------------------------------
    # Selectors
    # ------------------------------------------------------------------
    def bins_for_genome(self, genome: str) -> pd.DataFrame:
        """Return all bins belonging to *genome*.

        Args:
            genome: Name of the genome as it appears in the ``/genomes``
                table.

        Returns:
            Subset of the bins table (copy) containing only rows whose
            ``chrom`` (contig) belongs to *genome*.  Index values are
            preserved from the original bins table.

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        genomes = self.genomes
        contigs = self.contigs
        if genome not in genomes["name"].values:
            raise KeyError(
                f"Genome {genome!r} not found. Available: {genomes['name'].tolist()}"
            )
        gid = genomes.loc[genomes["name"] == genome].index[0]
        contig_ids = contigs.index[contigs["genome_id"] == gid].tolist()
        bins = self.bins()
        return bins[bins["chrom"].cat.codes.isin(contig_ids)].copy()

    def bins_for_contig(self, contig: str) -> pd.DataFrame:
        """Return all bins belonging to *contig*.

        Args:
            contig: Name of the contig as it appears in the ``/contigs``
                table.

        Returns:
            Subset of the bins table (copy) containing only rows whose
            ``chrom`` equals *contig*.  Index values are preserved from
            the original bins table.

        Raises:
            KeyError: If *contig* is not present in this file.
        """
        contigs = self.contigs
        if contig not in contigs["name"].values:
            raise KeyError(
                f"Contig {contig!r} not found. Available: {contigs['name'].tolist()}"
            )
        cid = contigs.loc[contigs["name"] == contig].index[0]
        bins = self.bins()
        return bins[bins["chrom"].cat.codes == cid].copy()

    def pixels_for_genome_pair(self, g1: str, g2: str) -> pd.DataFrame:
        """Return pixels where one bin is in *g1* and the other is in *g2*.

        Uses ``genome_extent`` index lookups to limit the scan to the
        relevant pixel rows — the full pixel table is never loaded into
        memory.

        Args:
            g1: Name of the first genome.
            g2: Name of the second genome.  May be the same as *g1* to
                retrieve intra-genomic contacts.

        Returns:
            DataFrame with columns ``bin1_id``, ``bin2_id``, and ``count``
            (plus any extra pixel columns).  Returns an empty DataFrame with
            those columns if no matching pixels exist.

        Raises:
            KeyError: If either *g1* or *g2* is not present in this file.
        """
        genomes = self.genomes
        for g in (g1, g2):
            if g not in genomes["name"].values:
                raise KeyError(
                    f"Genome {g!r} not found. Available: {genomes['name'].tolist()}"
                )

        s1, e1 = self.genome_extent(g1)
        s2, e2 = self.genome_extent(g2)

        # Determine the pixel range to scan (contiguous region in the upper triangle)
        lo = min(s1, s2)
        hi = max(e1, e2)

        chunks = []
        chunksize = 1_000_000
        n_pix = self.info["nnz"]
        # Use bin1_offset index to find the pixel range for bins [lo, hi)
        with h5py.File(self._store, "r") as f:
            root = f[self._group_path.lstrip("/")] if self._group_path != "/" else f
            bin1_offset = root["indexes"]["bin1_offset"][:]
        pix_lo = int(bin1_offset[lo])
        pix_hi = int(bin1_offset[hi]) if hi < len(bin1_offset) else n_pix

        bins = self._bins_cache
        bin_to_contig = bins["chrom"].cat.codes.values

        contigs = self.contigs
        gid1 = int(genomes.loc[genomes["name"] == g1].index[0])
        gid2 = int(genomes.loc[genomes["name"] == g2].index[0])
        cids1 = set(contigs.index[contigs["genome_id"] == gid1].tolist())
        cids2 = set(contigs.index[contigs["genome_id"] == gid2].tolist())

        valid1 = np.where(np.isin(bin_to_contig, list(cids1)))[0]
        valid2 = np.where(np.isin(bin_to_contig, list(cids2)))[0]
        set1 = set(valid1.tolist())
        set2 = set(valid2.tolist())

        for i in range(pix_lo, pix_hi, chunksize):
            chunk = self._clr.pixels()[i : i + chunksize]
            b1 = chunk["bin1_id"]
            b2 = chunk["bin2_id"]
            sel = (b1.isin(set1) & b2.isin(set2)) | (b1.isin(set2) & b2.isin(set1))
            if sel.any():
                chunks.append(chunk[sel])

        if not chunks:
            return pd.DataFrame(columns=["bin1_id", "bin2_id", "count"])
        return pd.concat(chunks, ignore_index=True)

    # ------------------------------------------------------------------
    # Extent helpers
    # ------------------------------------------------------------------
    def genome_extent(self, genome: str) -> tuple[int, int]:
        """Return the half-open bin-index range for all bins in *genome*.

        Args:
            genome: Name of the genome as it appears in the ``/genomes``
                table.

        Returns:
            ``(lo, hi)`` — half-open integer interval such that bins
            ``[lo, hi)`` all belong to *genome*.  Returns ``(0, 0)`` if
            *genome* has no bins.

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        genomes = self.genomes
        contigs = self.contigs
        if genome not in genomes["name"].values:
            raise KeyError(
                f"Genome {genome!r} not found. Available: {genomes['name'].tolist()}"
            )
        gid = genomes.loc[genomes["name"] == genome].index[0]
        cids = contigs.index[contigs["genome_id"] == gid].tolist()
        bins = self.bins()
        sel = bins[bins["chrom"].cat.codes.isin(cids)]
        if sel.empty:
            return (0, 0)
        return (int(sel.index[0]), int(sel.index[-1]) + 1)

    def contig_extent(self, contig: str) -> tuple[int, int]:
        """Return the half-open bin-index range for all bins in *contig*.

        Args:
            contig: Name of the contig as it appears in the ``/contigs``
                table.

        Returns:
            ``(lo, hi)`` — half-open integer interval such that bins
            ``[lo, hi)`` all belong to *contig*.  Returns ``(0, 0)`` if
            *contig* has no bins.

        Raises:
            KeyError: If *contig* is not present in this file.
        """
        contigs = self.contigs
        if contig not in contigs["name"].values:
            raise KeyError(
                f"Contig {contig!r} not found. Available: {contigs['name'].tolist()}"
            )
        cid = contigs.loc[contigs["name"] == contig].index[0]
        bins = self.bins()
        sel = bins[bins["chrom"].cat.codes == cid]
        if sel.empty:
            return (0, 0)
        return (int(sel.index[0]), int(sel.index[-1]) + 1)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    def to_cooler(
        self,
        path: str | Path,
        level: Literal["contig", "genome"] = "contig",
    ) -> None:
        """Export an aggregated contact matrix as a standard ``.cool`` file.

        The output is a regular cooler file where the "chromosomes" are
        either contigs (``level="contig"``) or genomes
        (``level="genome"``), making it immediately compatible with the
        cooler / cooltools / HiGlass ecosystem.

        Pixels are streamed in chunks; the full contact matrix is never
        held in memory simultaneously.  ``cooler`` performs an external
        sort to merge cross-chunk duplicates.

        Args:
            path: Destination path for the ``.cool`` file.  The file
                must not already exist.
            level: Aggregation level for the output chromosomes.
                ``"contig"`` keeps one row per contig (default);
                ``"genome"`` sums all contacts between genome pairs,
                producing one chromosome per genome.
        """
        import cooler as clr

        from metacooler.ops.aggregate import _aggregate_pixel_stream

        path = Path(path)

        if level == "contig":
            chroms = self.contigs[["name", "length"]].rename(
                columns={"name": "chrom", "length": "length"}
            )
        elif level == "genome":
            chroms = self.genomes[["name", "length"]].rename(
                columns={"name": "chrom", "length": "length"}
            )
        else:
            raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

        # One bin per chromosome at this level
        bins = pd.DataFrame(
            {
                "chrom": chroms["chrom"].values,
                "start": 0,
                "end": chroms["length"].values,
            }
        )

        clr.create_cooler(
            str(path),
            bins=bins,
            pixels=_aggregate_pixel_stream(self, level=level),
            ordered=False,
            symmetric_upper=True,
        )

    # ------------------------------------------------------------------
    # HDF5 access
    # ------------------------------------------------------------------
    @contextlib.contextmanager
    def open(self, mode: str = "r"):
        """Context manager that opens the underlying HDF5 file.

        Args:
            mode: h5py file-open mode.  ``"r"`` (default) opens read-only;
                ``"r+"`` opens for in-place modification.

        Yields:
            The root ``h5py.File`` object.  When the ``with`` block exits
            the file is closed automatically.
        """
        with h5py.File(self._store, mode) as f:
            yield f
