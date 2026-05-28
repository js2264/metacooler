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
        if kw:
            return cast(pd.DataFrame, self._clr.bins(**kw)[:])
        return self._bins_cache

    def pixels(self, join: bool = False, **kw) -> pd.DataFrame:
        return cast(pd.DataFrame, self._clr.pixels(join=join, **kw)[:])

    def matrix(self, **kw):
        return self._clr.matrix(**kw)

    def chroms(self, **kw) -> pd.DataFrame:
        return cast(pd.DataFrame, self._clr.chroms(**kw)[:])

    def extent(self, region: str) -> tuple[int, int]:
        return self._clr.extent(region)

    def offset(self, region: str) -> int:
        return self._clr.offset(region)

    @property
    def binsize(self) -> int | None:
        return self._clr.binsize

    @property
    def chromnames(self) -> list[str]:
        return self._clr.chromnames

    @property
    def chromsizes(self) -> pd.Series:
        return self._clr.chromsizes

    @property
    def shape(self) -> tuple[int, int]:
        return self._clr.shape

    @property
    def info(self) -> dict:
        """Merged metadata in metacool format."""
        base = dict(self._clr.info)
        # Promote metacooler-specific attrs over their cooler equivalents
        base["format"] = base.pop("metacooler-format")
        base["format-url"] = base.pop("metacooler-format-url")
        base["format-version"] = base.pop("metacooler-format-version")
        base["ordering"] = base.pop("contig-ordering")
        base["ncontigs"] = base.pop("nchroms")
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
        """The contents of the optional ``metadata`` group as a nested dict."""
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
        """The ``/genomes`` table as a DataFrame."""
        if self._group_path != "/":
            grp = self._group_path.lstrip("/") + "/genomes"
            return read_table(self._store, grp)
        return read_table(self._store, "genomes")

    @functools.cached_property
    def contigs(self) -> pd.DataFrame:
        """The ``/contigs`` table as a DataFrame (includes ``genome_id``)."""
        if self._group_path != "/":
            grp = self._group_path.lstrip("/") + "/contigs"
            return read_table(self._store, grp)
        return read_table(self._store, "contigs")

    # ------------------------------------------------------------------
    # Selectors
    # ------------------------------------------------------------------
    def bins_for_genome(self, genome: str) -> pd.DataFrame:
        """Return all bins belonging to *genome*."""
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
        """Return all bins belonging to *contig*."""
        contigs = self.contigs
        if contig not in contigs["name"].values:
            raise KeyError(
                f"Contig {contig!r} not found. Available: {contigs['name'].tolist()}"
            )
        cid = contigs.loc[contigs["name"] == contig].index[0]
        bins = self.bins()
        return bins[bins["chrom"].cat.codes == cid].copy()

    def pixels_for_genome_pair(self, g1: str, g2: str) -> pd.DataFrame:
        """Return pixels where one bin is in *g1* and the other in *g2*.

        Uses ``genome_extent`` index lookups instead of loading the full pixel
        table — safe for large datasets.
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
        """Return the inclusive bin-id range for *genome*."""
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
        """Return the inclusive bin-id range for *contig*."""
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
        """Export an aggregated matrix as a standard ``.cool`` file.

        The output is a regular cooler file where the "chromosomes" are
        either contigs or genomes, making it immediately compatible with
        the cooler/cooltools/HiGlass ecosystem.

        Pixels are streamed chunk-by-chunk; the full contact matrix is never
        materialised in memory.  ``cooler`` performs an external sort to merge
        cross-chunk duplicates.
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
        """Context manager for the root HDF5 file."""
        with h5py.File(self._store, mode) as f:
            yield f
