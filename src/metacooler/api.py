"""Primary Python API for metacool files."""

from __future__ import annotations

import contextlib
import functools
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, cast

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.schema import (
    CLUSTER_TABLE_COLUMNS,
    CLUSTERS_GROUP,
    DEFAULT_CHUNKSIZE,
    GENOME_GROUPED_ORDERINGS,
    _decode_table_group,
    _decode_text,
    cluster_rounds,
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
    # Clustering rounds (see docs/format.md#clusters)
    # ------------------------------------------------------------------

    @functools.cached_property
    def cluster_rounds(self) -> list[dict]:
        """The clustering rounds of the map, in the order they were made.

        Each is a dict with a ``name`` and what its tool stored with it (by
        convention ``tool``, ``version``, ``parameters`` and ``date``).

        Returns:
            One dict per round; empty when the map has no ``/clusters``.
        """
        with self.open() as grp:
            if CLUSTERS_GROUP not in grp:
                return []
            return cluster_rounds(grp[CLUSTERS_GROUP])

    @functools.cached_property
    def _cluster_table(self) -> pd.DataFrame:
        """``/clusters/table`` with its rounds and parents resolved to names.

        Read once and cached; :attr:`clusters` and the membership lookups name
        their rows from it.
        """
        columns = ["round", "name", "parent_round", "parent", "n_contigs", "length"]
        with self.open() as grp:
            if CLUSTERS_GROUP not in grp:
                return pd.DataFrame(columns=columns)
            table = _decode_table_group(grp[CLUSTERS_GROUP]["table"])
        round_names = np.array([r["name"] for r in self.cluster_rounds], dtype=object)
        rounds = round_names[table["round"].to_numpy(np.int64)]
        parents = table["parent"].to_numpy(np.int64)
        has = parents >= 0
        parent_round = np.full(len(table), None, dtype=object)
        parent = np.full(len(table), None, dtype=object)
        parent_round[has] = rounds[parents[has]]
        parent[has] = table["name"].to_numpy(dtype=object)[parents[has]]
        out = pd.DataFrame(
            {
                "round": rounds,
                "name": table["name"].to_numpy(dtype=object),
                "parent_round": parent_round,
                "parent": parent,
                "n_contigs": table["n_contigs"].to_numpy(),
                "length": table["length"].to_numpy(),
            }
        )
        extra = table.drop(columns=list(CLUSTER_TABLE_COLUMNS))
        return pd.concat([out, extra], axis=1)

    @functools.cached_property
    def clusters(self) -> pd.DataFrame:
        """One row per cluster of every clustering round.

        The columns are ``round`` and ``name``; ``parent_round`` and
        ``parent``, the round and name of the cluster it was split from
        (missing for none); ``n_contigs`` and ``length``, over its members
        (:meth:`cluster_members`); then the columns its tool stored.

        Returns:
            DataFrame indexed by cluster id, the ``cluster_id`` of the stored
            members; empty, with the columns above, when the map has no
            ``/clusters``.
        """
        return self._cluster_table.copy()

    def _check_round(self, round: str) -> None:
        names = [r["name"] for r in self.cluster_rounds]
        if round not in names:
            raise KeyError(
                f"{round!r} is not a clustering round of {self.uri}; its rounds: "
                f"{', '.join(names) or 'none'}."
            )

    def cluster_members(self, round: str | None = None) -> pd.DataFrame:
        """The contigs of each cluster, one row per membership.

        Reads the memberships of *round* only, through
        ``indexes/cluster_offset``: a round's clusters are stored together.

        Args:
            round: The name of one round, to read its clusters only; ``None``
                for every round.

        Returns:
            DataFrame with the columns ``round``, ``cluster`` and ``contig``
            (names), sorted by cluster, then contig row; empty when the map
            has no ``/clusters``.

        Raises:
            KeyError: If *round* is not a round of the map.
        """
        if round is not None:
            self._check_round(round)
        with self.open() as grp:
            if CLUSTERS_GROUP not in grp:
                return self._name_memberships(_NO_ROWS, _NO_ROWS)
            group = grp[CLUSTERS_GROUP]
            lo, hi = 0, len(group["members"]["cluster_id"])
            if round is not None:
                rows = np.flatnonzero(self._cluster_table["round"].to_numpy() == round)
                offsets = group["indexes"]["cluster_offset"]
                if rows.size:
                    lo, hi = int(offsets[rows[0]]), int(offsets[rows[-1] + 1])
                else:  # a round without clusters
                    lo = hi = 0
            cluster_ids = group["members"]["cluster_id"][lo:hi].astype(np.int64)
            contig_ids = group["members"]["contig_id"][lo:hi].astype(np.int64)
        return self._name_memberships(cluster_ids, contig_ids)

    def cluster_contigs(self, round: str, name: str) -> list[str]:
        """The contigs of one cluster, by name.

        One slice of the memberships, through ``indexes/cluster_offset``; only
        the names of those contigs are read.

        Args:
            round: The name of the cluster's round.
            name: The name of the cluster, in its round.

        Returns:
            Its contigs' names, in the order of their rows in ``/contigs``.

        Raises:
            KeyError: If *round* is not a round of the map, or has no cluster
                of that name.
        """
        self._check_round(round)
        table = self._cluster_table
        hit = np.flatnonzero(
            (table["round"].to_numpy() == round) & (table["name"].to_numpy() == name)
        )
        if not hit.size:
            raise KeyError(f"Round {round!r} of {self.uri} has no cluster {name!r}.")
        row = int(hit[0])
        with self.open() as grp:
            group = grp[CLUSTERS_GROUP]
            lo, hi = (int(v) for v in group["indexes"]["cluster_offset"][row : row + 2])
            contig_ids = group["members"]["contig_id"][lo:hi].astype(np.int64)
        return self._contig_names(contig_ids).tolist()

    def contig_clusters(self, contigs: str | Sequence[str]) -> pd.DataFrame:
        """The clusters a contig is in, in every round.

        One slice per contig, through ``indexes/contig_offset``. The clusters
        are numbered round by round, so each contig's come in the order of the
        rounds: its path through a clustering, whatever the method.

        Args:
            contigs: A contig name, or several.

        Returns:
            DataFrame with the columns ``contig``, ``round`` and ``cluster``,
            one row per membership: the contigs in the order given, each one's
            clusters in the order of the rounds. A contig in no cluster has no
            row; a map without ``/clusters`` gives none.

        Raises:
            KeyError: If a name is not a contig of the map.
        """
        names = [contigs] if isinstance(contigs, str) else [str(c) for c in contigs]
        ids = self._contig_ids_of(names)
        cluster_ids, counts = self._memberships_of(ids)
        table = self._cluster_table
        return pd.DataFrame(
            {
                "contig": np.repeat(np.array(names, dtype=object), counts),
                "round": table["round"].to_numpy(dtype=object)[cluster_ids],
                "cluster": table["name"].to_numpy(dtype=object)[cluster_ids],
            },
            columns=["contig", "round", "cluster"],
        )

    def _memberships_of(self, contig_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The clusters of each of *contig_ids*, through the contig index.

        Returns:
            ``(cluster_ids, counts)``: the cluster ids of every membership,
            contig after contig in the order given (each contig's in cluster
            order), and how many each contig has.
        """
        ids = np.asarray(contig_ids, dtype=np.int64)
        with self.open() as grp:
            if CLUSTERS_GROUP not in grp:
                return _NO_ROWS, np.zeros(len(ids), dtype=np.int64)
            indexes = grp[CLUSTERS_GROUP]["indexes"]
            bounds = _read_rows(
                indexes["contig_offset"], np.concatenate([ids, ids + 1])
            )
            lo, counts = bounds[: len(ids)], bounds[len(ids) :] - bounds[: len(ids)]
            # The rows lo[i] .. lo[i] + counts[i] of each contig, end to end.
            rows = np.repeat(lo - np.cumsum(counts) + counts, counts) + np.arange(
                int(counts.sum())
            )
            cluster_ids = _read_rows(indexes["contig_clusters"], rows).astype(np.int64)
        return cluster_ids, counts

    def _name_memberships(
        self, cluster_ids: np.ndarray, contig_ids: np.ndarray
    ) -> pd.DataFrame:
        """Memberships by name: one row per ``(cluster_id, contig_id)``."""
        table = self._cluster_table
        return pd.DataFrame(
            {
                "round": table["round"].to_numpy(dtype=object)[cluster_ids],
                "cluster": table["name"].to_numpy(dtype=object)[cluster_ids],
                "contig": self._contig_names(contig_ids),
            },
            columns=["round", "cluster", "contig"],
        )

    @functools.cached_property
    def _contig_name_strings(self) -> np.ndarray:
        """``/contigs/name`` decoded, as ``str`` objects (cached)."""
        if "contigs" in self.__dict__:  # the table is loaded already
            return self.contigs["name"].to_numpy(dtype=object)
        return _decode_names(self._contig_names_raw)

    def _contig_names(self, contig_ids: np.ndarray) -> np.ndarray:
        """The names of *contig_ids*, as ``str`` objects.

        Never through the whole ``/contigs`` table, whose every column
        :attr:`contigs` decodes: a few names are read from ``/contigs/name``
        alone; more, gathered from the whole column and decoded alone; and
        most of the column, decoded once and kept.
        """
        ids = np.asarray(contig_ids, dtype=np.int64)
        if "_contig_name_strings" in self.__dict__ or "contigs" in self.__dict__:
            return self._contig_name_strings[ids]
        if len(ids) <= _POINT_READS:
            with self.open() as grp:
                return _decode_names(_read_rows(grp["contigs"]["name"], ids))
        if len(ids) < len(self._contig_names_raw) // 2:
            return _decode_names(self._contig_names_raw[ids])
        return self._contig_name_strings[ids]

    def _contig_ids_of(self, names: Sequence[str]) -> np.ndarray:
        """The rows of *names* in ``/contigs``, in one pass over the stored names.

        Raises:
            KeyError: If a name is not a contig of the map.
        """
        if "_clr" in self.__dict__:  # cooler already holds the name -> id dict
            found = {name: self._clr._chromids.get(name) for name in names}
        elif len(names) == 1:
            found = {names[0]: self._find_contig(names[0])}
        else:
            keys = {}
            for name in names:
                with contextlib.suppress(UnicodeEncodeError):  # names are ASCII
                    keys[name] = name.encode("ascii")
            raw = self._contig_names_raw
            hits = np.flatnonzero(np.isin(raw, list(set(keys.values()))))
            row_of = {bytes(raw[i]): int(i) for i in hits}
            found = {name: row_of.get(keys.get(name, b"")) for name in names}
        missing = [name for name in names if found[name] is None]
        if missing:
            shown = ", ".join(map(repr, missing[:5])) + (
                ", ..." if len(missing) > 5 else ""
            )
            raise KeyError(f"Not contigs of {self.uri}: {shown}.")
        return np.array([found[name] for name in names], dtype=np.int64)

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
        A bin table already cached by :meth:`bins` is sliced instead.

        The selector is asked for contig *ids* (``convert_enum=False``) and the
        names are put on once, at the end, exactly as cooler would.  Left to
        convert them itself, cooler reads and decodes the whole ``/chroms``
        table on every read -- once per range -- which made a genome of a few
        hundred contigs scattered over a ``--contig-order alpha`` map cost
        ranges x contigs: a minute at a million contigs.
        """
        if "_bins_cache" in self.__dict__:
            return self._bins_cache.iloc[_ranges_to_ids(ranges)]
        bins = self._clr.bins(convert_enum=False)
        if not len(ranges):
            out = cast(pd.DataFrame, bins[0:0])
        elif len(ranges) <= _MAX_RANGE_READS:
            parts = [cast(pd.DataFrame, bins[int(lo) : int(hi)]) for lo, hi in ranges]
            out = pd.concat(parts) if len(parts) > 1 else parts[0]
        else:
            # Many short runs: one read of the span beats thousands of small
            # ones, each of which decompresses a whole chunk.
            lo, hi = int(ranges[0, 0]), int(ranges[-1, 1])
            out = cast(pd.DataFrame, bins[lo:hi]).iloc[_ranges_to_ids(ranges) - lo]
        out = out.copy()
        out["chrom"] = pd.Categorical.from_codes(
            out["chrom"].to_numpy(), dtype=self._chrom_dtype
        )
        return out

    @functools.cached_property
    def _chrom_dtype(self) -> pd.CategoricalDtype:
        """The dtype cooler gives ``bins.chrom``: the contig names, ordered.

        Read once with cooler's own chroms selector, so the categories match
        what ``cooler.Cooler.bins()`` builds on each read, dtype included.
        """
        names = cast(pd.Series, self._clr.chroms()["name"][:])
        return pd.CategoricalDtype(pd.Index(names), ordered=True)

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


#: An empty row selection.
_NO_ROWS = np.zeros(0, dtype=np.int64)

#: Up to this many rows, a lookup reads them one by one, each in the chunk that
#: holds it; past it, reading the whole dataset once is faster.
_POINT_READS = 1024


def _decode_names(raw: np.ndarray) -> np.ndarray:
    """Contig names as stored, decoded to ``str`` objects.

    Names are ASCII (``check_ascii_names`` on write, ``_check_decodable`` on
    open), which numpy decodes in C, several times faster than the UTF-8
    decode of :func:`metacooler.schema._decode_text`, kept for what is not.
    """
    try:
        return raw.astype(str).astype(object)
    except UnicodeDecodeError:
        return _decode_text(raw).astype(object)


def _read_rows(dataset: h5py.Dataset, rows: np.ndarray) -> np.ndarray:
    """The values of *dataset* at *rows*, in that order (repeats allowed)."""
    rows = np.asarray(rows, dtype=np.int64)
    if len(rows) > _POINT_READS:
        return np.asarray(dataset[()])[rows]
    unique, inverse = np.unique(rows, return_inverse=True)
    if not len(unique):
        return np.asarray(dataset[0:0])
    return np.asarray(dataset[unique])[inverse]
