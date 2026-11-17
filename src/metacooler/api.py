"""Primary Python API for metacool files."""

from __future__ import annotations

import contextlib
import functools
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import TYPE_CHECKING, cast

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.schema import (
    DEFAULT_CHUNKSIZE,
    GENOME_GROUPED_ORDERINGS,
    read_table,
)
from metacooler.schema import (
    validate as _validate,
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

    def __init__(self, store: str | Path, *, validate: bool = True):
        """Open a ``.metacool`` (or mmetacool entry) file.

        Args:
            store: File path, optionally with a URI group suffix
                (``file.mmetacool::entries/entry_name``).
            validate: If ``True`` (default), run the *structural* checks of
                :func:`metacooler.schema.validate` (``full=False``) before
                opening: format attributes, required groups and columns,
                equal column lengths, the ``/chroms`` link, and the length
                and endpoints of every offset index.  They read O(n_columns)
                values, so an open costs milliseconds even on a
                multi-million-contig metagenome.  The O(n_bins) content
                checks (bin geometry, offset monotonicity, genome tables,
                name rules, sampled pixel order) are what ``metacooler
                validate`` and :func:`metacooler.schema.validate` run, and
                every metacooler writer runs them on what it writes.  Pass
                ``False`` to skip even the structural checks; a malformed
                file then surfaces later as an opaque ``IndexError`` or
                ``KeyError`` from cooler.

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
            # Validate that the group exists and points to a well-formed metacooler
            with h5py.File(self._store, "r") as f:
                if group_path.lstrip("/") not in f:
                    raise KeyError(f"Group {group_path!r} not found in {self._store}")
            if validate:
                _validate(f"{self._store}::{group_path}", full=False)
            self._cooler_uri = f"{self._store}::{group_path}"
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
            if validate:
                _validate(self._store, full=False)
            self._cooler_uri = str(self._store)
        self._validate_on_open = validate
        if validate:
            self._check_decodable()

    @functools.cached_property
    def _clr(self) -> cooler.Cooler:
        """The ``cooler.Cooler`` this map delegates to, opened on first use.

        Opening one reads every contig name into a dict (1.5 s for two million
        contigs), which the index-based selectors, ``info``, ``binsize`` and
        ``shape`` never need; so an open, and a genome or contig query, does not
        pay for it.
        """
        return _open_cooler(self._cooler_uri, self._validate_on_open)

    def _check_decodable(self) -> None:
        """Fail now, by name, if cooler could not decode ``/contigs``.

        cooler decodes every byte-string column of ``/chroms`` (= ``/contigs``)
        as ASCII when it opens a map; the open used to fail there.  With cooler
        opened lazily, the same condition -- a byte at or above 0x80 -- is
        checked here, vectorised, and handed to the full validation, whose
        error names the offending values.
        """
        with self.open() as grp:
            columns = grp["contigs"]
            for key in columns:
                column = columns[key]
                if (
                    isinstance(column, h5py.Dataset)
                    and column.dtype.kind == "S"
                    and len(column)
                    and int(np.max(np.ascontiguousarray(column[()]).view(np.uint8)))
                    >= 0x80
                ):
                    _validate(self._cooler_uri)
                    raise UnicodeDecodeError(  # pragma: no cover - validate raised
                        "ascii", b"", 0, 1, f"/contigs/{key} is not ASCII"
                    )

    def _attrs(self) -> dict:
        """The map's attributes as ``cooler.Cooler.info`` reads them, fresh."""
        from cooler.api import info as cooler_info

        with self.open() as grp:
            return cooler_info(grp)

    @functools.cached_property
    def _attrs_at_open(self) -> dict:
        return self._attrs()

    @property
    def filename(self) -> str:
        """The file path as a string."""
        return str(self._store)

    def __repr__(self) -> str:
        return f"<MetaCooler {self.uri!r}>"

    @property
    def root(self) -> str:
        """HDF5 path of the map inside its file, as ``cooler.Cooler.root``.

        ``/`` for a ``.metacool``; ``/entries/<name>`` for an ``.mmetacool``
        entry (it used to claim ``/`` for those too).
        """
        from cooler.util import parse_cooler_uri

        return parse_cooler_uri(self._cooler_uri)[1]

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

        Delegates to ``cooler.Cooler.pixels()``, except for the bin coordinates
        added by *join*. ``/bins/chrom`` is a plain integer column in a metacool
        file, and cooler's join reports such a column as contig *ids*, where
        ``cooler.Cooler.bins()`` converts it to names; so the coordinates are
        joined here, from :meth:`bins`, and name the contigs.

        Args:
            join: If ``True``, annotate each pixel with its bin coordinates
                (adds ``chrom1``, ``start1``, ``end1``, ``chrom2``, ``start2``,
                ``end2`` columns).
            **kw: Additional keyword arguments forwarded to
                ``cooler.Cooler.pixels()``.

        Returns:
            DataFrame with at least ``bin1_id``, ``bin2_id``, and ``count``
            columns.  When *join* is ``True``, bin coordinate columns are
            prepended, with ``chrom1``/``chrom2`` holding contig names (contig
            ids if ``convert_enum=False`` is passed).
        """
        if not join or kw.get("convert_enum", True) is False:
            return cast(pd.DataFrame, self._clr.pixels(join=join, **kw)[:])
        pixels = cast(pd.DataFrame, self._clr.pixels(**kw)[:])
        coords = self._bins_cache[["chrom", "start", "end"]]
        return cooler.api.annotate(pixels, coords, replace=True)

    def iter_pixels(
        self, chunksize: int = DEFAULT_CHUNKSIZE, *, contig_ids: bool = False
    ) -> Iterator[pd.DataFrame]:
        """Yield the pixel table in consecutive chunks, never all of it at once.

        Each chunk is the next *chunksize* rows of ``/pixels``, read with
        cooler's pixel selector and indexed by pixel row as :meth:`pixels` is,
        so the chunks put together are :meth:`pixels` exactly, in file order
        (sorted by ``bin1_id``, then ``bin2_id``).  Memory stays at one chunk
        whatever the size of the map.

        Args:
            chunksize: Most pixel rows per chunk.
            contig_ids: Also give each pixel ``contig1_id`` and ``contig2_id``,
                the contigs (rows of :attr:`contigs`) its two bins belong to,
                looked up through ``indexes/chrom_offset``.  The bin table is
                not read.

        Yields:
            DataFrames with every column of ``/pixels`` (``bin1_id``,
            ``bin2_id``, ``count``, ...), then ``contig1_id`` and ``contig2_id``
            (``int64``) when asked for.  Nothing for a map without pixels.

        Raises:
            ValueError: If *chunksize* is below 1.
        """
        if chunksize < 1:
            raise ValueError(f"chunksize must be at least 1, got {chunksize}")
        pixels = self._clr.pixels()
        contig_of_bin = self._bin_contig_ids if contig_ids else None
        for lo in range(0, len(pixels), chunksize):
            chunk = cast(pd.DataFrame, pixels[lo : lo + chunksize])
            if contig_of_bin is not None:
                chunk["contig1_id"] = contig_of_bin[chunk["bin1_id"].to_numpy()]
                chunk["contig2_id"] = contig_of_bin[chunk["bin2_id"].to_numpy()]
            yield chunk

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
        return self._attrs_at_open["bin-size"]

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
        n_bins = self._attrs_at_open["nbins"]
        return (n_bins, n_bins)

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
        base = dict(self._attrs())
        base["format"] = base.pop("metacooler-format", None)
        base["format-url"] = base.pop("metacooler-format-url", None)
        base["format-version"] = base.pop("metacooler-format-version", None)
        base["ordering"] = base.pop("contig-ordering", None)
        base["ncontigs"] = base.pop("nchroms", None)
        # Normalise the assembly key. The file stores it in cooler's standard
        # "genome-assembly" attribute; "assembly-name" is only read as a
        # fallback so files written by older metacooler versions (which wrote
        # that attribute instead) still report their assembly.
        legacy = base.pop("assembly-name", None)
        assembly = base.pop("genome-assembly", base.pop("assembly", None))
        if assembly in (None, "", "unknown") and legacy is not None:
            assembly = legacy
        base["assembly"] = assembly
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
                    as_int = int(base[_key])
                    # A float-count map has a fractional sum; int() floored it.
                    if as_int == float(base[_key]):
                        base[_key] = as_int
        return base

    @property
    def metadata(self) -> dict:
        """Contents of the optional ``metadata`` group as a nested dict.

        Returns an empty dict if no metadata was stored at creation time.

        Returns:
            Nested dictionary of user-defined metadata values, or ``{}``
            if the file contains no metadata.
        """
        base = dict(self._attrs())
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
    # Index lookups
    #
    # Everything below resolves names to bin ranges through the three offset
    # indexes, never through the bin table: a genome to its contigs via
    # ``indexes/genome_offset`` (or, when the contig ordering does not group
    # genomes, ``/contigs/genome_id``), a contig to its bins via
    # ``indexes/chrom_offset``, and a bin range to its pixel rows via
    # ``indexes/bin1_offset``.  Each of those is O(n_contigs) at worst, so a
    # genome or contig query costs milliseconds whatever the size of the map.
    # ------------------------------------------------------------------
    def _read_index(self, name: str) -> np.ndarray:
        with self.open() as grp:
            return grp["indexes"][name][:].astype(np.int64)

    @functools.cached_property
    def _chrom_offset(self) -> np.ndarray:
        return self._read_index("chrom_offset")

    @functools.cached_property
    def _genome_offset(self) -> np.ndarray:
        return self._read_index("genome_offset")

    @functools.cached_property
    def _bin1_offset(self) -> np.ndarray:
        return self._read_index("bin1_offset")

    @functools.cached_property
    def _bin_contig_ids(self) -> np.ndarray:
        """Contig id of every bin, from ``chrom_offset`` alone (one ``repeat``)."""
        offsets = self._chrom_offset
        return np.repeat(np.arange(len(offsets) - 1, dtype=np.int64), np.diff(offsets))

    @functools.cached_property
    def _contig_genome_ids(self) -> np.ndarray:
        with self.open() as grp:
            return grp["contigs"]["genome_id"][:].astype(np.int64)

    @functools.cached_property
    def _contig_names_raw(self) -> np.ndarray:
        """``/contigs/name`` as stored: fixed-width bytes, never decoded."""
        with self.open() as grp:
            return grp["contigs"]["name"][()]

    def _find_contig(self, contig: str) -> int | None:
        """Row of *contig* in ``/contigs``, or ``None``.

        One vectorised comparison over the stored bytes: a few milliseconds for
        two million names, where a name -> id dict costs half a second to build.
        """
        if "_clr" in self.__dict__:  # cooler already holds that dict
            return self._clr._chromids.get(contig)
        try:
            key = contig.encode("ascii")  # names are ASCII (check_ascii_names)
        except UnicodeEncodeError:
            return None
        hits = np.flatnonzero(self._contig_names_raw == key)
        return int(hits[0]) if hits.size else None

    @functools.cached_property
    def _genome_ids(self) -> dict[str, int]:
        return {str(name): i for i, name in enumerate(self.genomes["name"])}

    @functools.cached_property
    def _genome_grouped(self) -> bool:
        """Whether ``genome_offset`` delimits each genome's contig block."""
        with self.open() as grp:
            ordering = grp.attrs.get("contig-ordering")
        if isinstance(ordering, bytes):
            ordering = ordering.decode()
        return ordering in GENOME_GROUPED_ORDERINGS

    def _genome_id(self, genome: str) -> int:
        gid = self._genome_ids.get(genome)
        if gid is None:
            raise KeyError(
                f"Genome {genome!r} not found. "
                f"Available: {self.genomes['name'].tolist()}"
            )
        return gid

    def _contig_id(self, contig: str) -> int:
        cid = self._find_contig(contig)
        if cid is None:
            raise KeyError(f"Contig {contig!r} not found. Available: {self.chromnames}")
        return cid

    def genome_contig_ids(self, genome: str) -> np.ndarray:
        """Return the sorted contig ids (row numbers of ``/contigs``) of *genome*.

        Under a genome-grouped contig ordering this is one slice of
        ``indexes/genome_offset``; otherwise it is the rows of
        ``/contigs/genome_id`` equal to the genome's id.  The bin table is
        never read.

        Args:
            genome: Name of the genome as it appears in the ``/genomes`` table.

        Returns:
            1-D ``int64`` array of contig ids (possibly empty).

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        gid = self._genome_id(genome)
        if self._genome_grouped:
            lo, hi = self._genome_offset[gid : gid + 2]
            return np.arange(lo, hi, dtype=np.int64)
        return np.flatnonzero(self._contig_genome_ids == gid).astype(np.int64)

    def _contig_bin_ranges(self, contig_ids: np.ndarray) -> np.ndarray:
        """Half-open bin ranges of *contig_ids* (sorted), adjacent ones merged.

        Returns:
            ``(k, 2)`` ``int64`` array of ``[lo, hi)`` rows, sorted, disjoint
            and non-adjacent; empty contigs contribute nothing.
        """
        cids = np.asarray(contig_ids, dtype=np.int64)
        lo = self._chrom_offset[cids]
        hi = self._chrom_offset[cids + 1]
        return _merge_ranges(lo, hi)

    def genome_bin_ranges(self, genome: str) -> np.ndarray:
        """Return *genome*'s bins as sorted, disjoint half-open ranges.

        One range under a genome-grouped ordering; one per run of
        consecutive contigs otherwise.

        Args:
            genome: Name of the genome as it appears in the ``/genomes`` table.

        Returns:
            ``(k, 2)`` ``int64`` array whose rows are ``[lo, hi)`` bin ranges.

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        return self._contig_bin_ranges(self.genome_contig_ids(genome))

    def _bins_in_ranges(self, ranges: np.ndarray) -> pd.DataFrame:
        """The bin table rows in *ranges*, as :meth:`bins` would show them.

        Read with cooler's bin selector, one range at a time (or, for many
        short runs, their span at once), so nothing outside the span is read.
        """
        bins = self._clr.bins()
        if not len(ranges):
            return cast(pd.DataFrame, bins[0:0])
        if len(ranges) <= _MAX_RANGE_READS:
            parts = [cast(pd.DataFrame, bins[int(lo) : int(hi)]) for lo, hi in ranges]
            return pd.concat(parts) if len(parts) > 1 else parts[0]
        # Many short runs: one read of the span beats thousands of small ones,
        # each of which decompresses a whole chunk.
        lo, hi = int(ranges[0, 0]), int(ranges[-1, 1])
        return cast(pd.DataFrame, bins[lo:hi]).iloc[_ranges_to_ids(ranges) - lo]

    def _pixel_offsets(self, bins: np.ndarray) -> np.ndarray:
        """``bin1_offset`` at each of *bins* (int64 array of bin ids)."""
        if len(bins) > _MAX_RANGE_READS or "_bin1_offset" in self.__dict__:
            return self._bin1_offset[bins]
        with self.open() as grp:
            ds = grp["indexes"]["bin1_offset"]
            return np.array([int(ds[int(b)]) for b in bins], dtype=np.int64)

    def _iter_pixel_rows(
        self,
        row_ranges: np.ndarray,
        keep: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None,
        chunksize: int = DEFAULT_CHUNKSIZE,
    ) -> Iterator[pd.DataFrame]:
        """Yield the stored pixels whose ``bin1_id`` lies in *row_ranges*.

        Each range of rows is turned into a block of pixel rows through
        ``indexes/bin1_offset``, and the blocks are read with cooler's pixel
        selector; blocks that nearly touch are read together, and nothing
        outside them is read.  Frames come out in file order, so their
        concatenation is sorted by ``(bin1_id, bin2_id)``.

        Args:
            row_ranges: ``(k, 2)`` sorted, disjoint ``[lo, hi)`` bin ranges.
            keep: Optional ``keep(bin1, bin2) -> bool mask`` applied to every
                block (after the row-range test).
            chunksize: Most pixel rows read at once.

        Yields:
            Non-empty pixel DataFrames, every column of ``/pixels``, indexed
            by pixel row like :meth:`pixels`.
        """
        ranges = np.asarray(row_ranges, dtype=np.int64).reshape(-1, 2)
        ranges = ranges[ranges[:, 1] > ranges[:, 0]]
        if not len(ranges):
            return
        offsets = self._pixel_offsets(ranges.ravel()).reshape(-1, 2)
        spans = _merge_ranges(offsets[:, 0], offsets[:, 1], gap=_PIXEL_READ_GAP)
        exact = len(spans) == len(ranges)  # no gap rows were swept in
        pixels = self._clr.pixels()
        for p0, p1 in spans:
            for lo in range(int(p0), int(p1), chunksize):
                chunk = cast(pd.DataFrame, pixels[lo : min(lo + chunksize, int(p1))])
                if not exact or keep is not None:
                    b1 = chunk["bin1_id"].to_numpy()
                    mask = (
                        np.ones(len(chunk), dtype=bool)
                        if exact
                        else _in_ranges(b1, ranges)
                    )
                    if keep is not None:
                        mask &= keep(b1, chunk["bin2_id"].to_numpy())
                    if not mask.all():
                        chunk = chunk[mask]
                if len(chunk):
                    yield chunk

    def _empty_pixels(self) -> pd.DataFrame:
        """The pixel table's columns, typed, without a row (cooler's selector)."""
        return cast(pd.DataFrame, self._clr.pixels()[:0])

    # ------------------------------------------------------------------
    # Selectors
    # ------------------------------------------------------------------
    def bins_for_genome(self, genome: str) -> pd.DataFrame:
        """Return all bins belonging to *genome*.

        Args:
            genome: Name of the genome as it appears in the ``/genomes``
                table.

        Returns:
            Subset of the bins table containing only rows whose ``chrom``
            (contig) belongs to *genome*.  Index values are preserved from
            the original bins table.

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        return self._bins_in_ranges(self.genome_bin_ranges(genome))

    def bins_for_contig(self, contig: str) -> pd.DataFrame:
        """Return all bins belonging to *contig*.

        Args:
            contig: Name of the contig as it appears in the ``/contigs``
                table.

        Returns:
            Subset of the bins table containing only rows whose ``chrom``
            equals *contig*.  Index values are preserved from the original
            bins table.

        Raises:
            KeyError: If *contig* is not present in this file.
        """
        cid = self._contig_id(contig)
        return self._bins_in_ranges(self._contig_bin_ranges(np.array([cid])))

    def pixels_for_genome_pair(self, g1: str, g2: str) -> pd.DataFrame:
        """Return pixels where one bin is in *g1* and the other is in *g2*.

        Reads only the pixel rows of the two genomes' bins, located through
        ``indexes/bin1_offset``, and keeps those whose other end lies in the
        other genome — milliseconds per genome pair, whatever the size of the
        map.  Correct under any contig ordering, and for ``square`` storage.

        Args:
            g1: Name of the first genome.
            g2: Name of the second genome.  May be the same as *g1* to
                retrieve intra-genomic contacts.

        Returns:
            DataFrame with every column of the pixel table (``bin1_id``,
            ``bin2_id``, ``count``, ...), sorted by ``(bin1_id, bin2_id)``
            and indexed from 0.  Empty (with those columns) if no pixel
            matches.

        Raises:
            KeyError: If either *g1* or *g2* is not present in this file.
        """
        a = self.genome_bin_ranges(g1)
        b = self.genome_bin_ranges(g2)
        if g1 == g2:
            rows = a

            def keep(b1: np.ndarray, b2: np.ndarray) -> np.ndarray:
                return _in_ranges(b2, a)
        else:
            # Two genomes share no contig, so their ranges are disjoint.
            rows = _merge_ranges(
                np.concatenate([a[:, 0], b[:, 0]]),
                np.concatenate([a[:, 1], b[:, 1]]),
            )

            def keep(b1: np.ndarray, b2: np.ndarray) -> np.ndarray:
                return (_in_ranges(b1, a) & _in_ranges(b2, b)) | (
                    _in_ranges(b1, b) & _in_ranges(b2, a)
                )

        chunks = list(self._iter_pixel_rows(rows, keep))
        if not chunks:
            return self._empty_pixels()
        return pd.concat(chunks, ignore_index=True)

    # ------------------------------------------------------------------
    # Extent helpers
    # ------------------------------------------------------------------
    def genome_bin_ids(self, genome: str) -> np.ndarray:
        """Return the sorted array of bin indices belonging to *genome*.

        Selection is by contig membership, so it is correct regardless of the
        file's contig ordering (the bins need not be contiguous).

        Args:
            genome: Name of the genome as it appears in the ``/genomes`` table.

        Returns:
            1-D ``int64`` array of bin indices (possibly empty).

        Raises:
            KeyError: If *genome* is not present in this file.
        """
        return _ranges_to_ids(self.genome_bin_ranges(genome))

    def genome_extent(self, genome: str) -> tuple[int, int]:
        """Return the half-open bin-index range for all bins in *genome*.

        This is only well-defined when the genome's bins form one contiguous
        block — i.e. the file uses a genome-grouped contig ordering
        (``genome-then-length`` / ``genome-then-alpha``, the default).  For
        non-grouped orderings (``fasta``, ``alpha``, ``cool``, ``custom``) a
        genome's contigs may be interleaved with others', so a single bin
        interval would silently include foreign bins; in that case this method
        raises rather than return a misleading range.

        Args:
            genome: Name of the genome as it appears in the ``/genomes``
                table.

        Returns:
            ``(lo, hi)`` — half-open integer interval such that bins
            ``[lo, hi)`` all belong to *genome*.  Returns ``(0, 0)`` if
            *genome* has no bins.

        Raises:
            KeyError: If *genome* is not present in this file.
            ValueError: If the genome's bins are not contiguous under the
                file's contig ordering.
        """
        ranges = self.genome_bin_ranges(genome)
        if len(ranges) == 0:
            return (0, 0)
        if len(ranges) > 1:
            ordering = self.info.get("ordering")
            raise ValueError(
                f"Genome {genome!r} is not contiguous in bin space under "
                f"contig-ordering {ordering!r}; a single bin range cannot "
                f"represent it. Re-create the file with a genome-grouped "
                f"ordering (genome-then-length/genome-then-alpha), or select by "
                f"contig membership (e.g. 'dump --genomes' / 'extract --genomes')."
            )
        return (int(ranges[0, 0]), int(ranges[0, 1]))

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
        cid = self._contig_id(contig)
        lo, hi = (int(x) for x in self._chrom_offset[cid : cid + 2])
        return (0, 0) if lo == hi else (lo, hi)

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
        # Delegate to the shared writer in ops/ so aggregation logic (bin→level
        # mapping, int-overflow handling, create_cooler call) lives in one place.
        from metacooler.ops.aggregate import _write_aggregated_cool

        _write_aggregated_cool(self, Path(path), level=level)

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
            if self._group_path and self._group_path != "/":
                yield f[self._group_path.lstrip("/")]
            else:
                yield f


def _open_cooler(uri: str, validate: bool) -> cooler.Cooler:
    """``cooler.Cooler(uri)``, naming the problem when cooler cannot decode it.

    cooler decodes every ``/chroms`` (= ``/contigs``) column as ASCII, so a
    file holding a non-ASCII name fails here with a bare ``UnicodeDecodeError``.
    The structural checks of an open do not read names (that is O(n_contigs));
    on that failure the full validation runs instead, and its error names the
    offending values.
    """
    try:
        return cooler.Cooler(uri)
    except UnicodeDecodeError:
        if validate:
            _validate(uri)
        raise


#: Up to this many ranges are read one by one; past it, one read of their span.
_MAX_RANGE_READS = 256

#: Pixel row blocks closer than this are read in one go (and the rows between
#: them filtered out): one h5py read per block costs more than a few thousand
#: extra rows of decompression.
_PIXEL_READ_GAP = 1 << 16


def _merge_ranges(lo: np.ndarray, hi: np.ndarray, gap: int = 0) -> np.ndarray:
    """Sort half-open ranges and merge those that overlap or are within *gap*.

    Empty ranges are dropped.  Vectorised: one sort and one cumulative max.

    Returns:
        ``(k, 2)`` ``int64`` array of merged ``[lo, hi)`` rows.
    """
    lo = np.asarray(lo, dtype=np.int64)
    hi = np.asarray(hi, dtype=np.int64)
    keep = hi > lo
    lo, hi = lo[keep], hi[keep]
    if not len(lo):
        return np.empty((0, 2), dtype=np.int64)
    order = np.argsort(lo, kind="stable")
    lo, hi = lo[order], hi[order]
    reach = np.maximum.accumulate(hi)
    starts = np.ones(len(lo), dtype=bool)
    starts[1:] = lo[1:] > reach[:-1] + gap
    first = np.flatnonzero(starts)
    last = np.append(first[1:], len(lo)) - 1
    return np.column_stack([lo[first], reach[last]])


def _in_ranges(values: np.ndarray, ranges: np.ndarray) -> np.ndarray:
    """Boolean mask: which *values* fall in one of the sorted disjoint *ranges*."""
    if not len(ranges):
        return np.zeros(len(values), dtype=bool)
    if len(ranges) == 1:
        return (values >= ranges[0, 0]) & (values < ranges[0, 1])
    i = np.searchsorted(ranges[:, 0], values, side="right") - 1
    return (i >= 0) & (values < ranges[np.maximum(i, 0), 1])


def _ranges_to_ids(ranges: np.ndarray) -> np.ndarray:
    """Concatenate ``arange(lo, hi)`` over *ranges*, without a Python loop."""
    if not len(ranges):
        return np.empty(0, dtype=np.int64)
    lengths = ranges[:, 1] - ranges[:, 0]
    total = int(lengths.sum())
    shift = np.repeat(ranges[:, 0] - (np.cumsum(lengths) - lengths), lengths)
    return shift + np.arange(total, dtype=np.int64)
