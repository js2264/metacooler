"""Primary Python API for metacool files."""

from __future__ import annotations

import contextlib
from pathlib import Path
from typing import TYPE_CHECKING, cast

import cooler
import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp

from metacooler.schema import (
    METACOOLER_FORMAT_URL,
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
        """The cooler-style URI (``path::/``)."""
        return f"{self._store}::/"

    # ------------------------------------------------------------------
    # Delegated to cooler.Cooler (pass-through)
    # ------------------------------------------------------------------
    def bins(self, **kw) -> pd.DataFrame:
        return cast(pd.DataFrame, self._clr.bins(**kw)[:])

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
        with h5py.File(self._store, "r") as f:
            root = f[self._group_path.lstrip("/")] if self._group_path != "/" else f
            base["format"] = "HDF5::MetaCooler"
            base["format-url"] = root.attrs.get(
                "metacooler-format-url", METACOOLER_FORMAT_URL
            )
            base["ordering"] = root.attrs.get("contig-ordering", "unknown")
            base["ncontigs"] = len(root["contigs"]["name"])
            base["ngenomes"] = len(root["genomes"]["name"])
            # Remove cooler-specific keys that are redundant
            base.pop("metacooler-format", None)
            base.pop("contig-ordering", None)
            base.pop("n_genomes", None)
            base.pop("n_contigs", None)
            base.pop("nchroms", None)
            base.pop("metadata", None)
        # Cast numeric fields to int (cooler may return them as strings or numpy scalars)
        for _key in ("nbins", "nnz", "sum"):
            if _key in base:
                with contextlib.suppress(ValueError, TypeError):
                    base[_key] = int(base[_key])
        return base

    # ------------------------------------------------------------------
    # metacooler-specific properties
    # ------------------------------------------------------------------
    @property
    def genomes(self) -> pd.DataFrame:
        """The ``/genomes`` table as a DataFrame."""
        if self._group_path != "/":
            grp = self._group_path.lstrip("/") + "/genomes"
            return read_table(self._store, grp)
        return read_table(self._store, "genomes")

    @property
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
        """Return pixels where one bin is in *g1* and the other in *g2*."""
        genomes = self.genomes
        contigs = self.contigs
        for g in (g1, g2):
            if g not in genomes["name"].values:
                raise KeyError(
                    f"Genome {g!r} not found. Available: {genomes['name'].tolist()}"
                )

        gid1 = genomes.loc[genomes["name"] == g1].index[0]
        gid2 = genomes.loc[genomes["name"] == g2].index[0]
        cids1 = set(contigs.index[contigs["genome_id"] == gid1].tolist())
        cids2 = set(contigs.index[contigs["genome_id"] == gid2].tolist())

        bins = self.bins()
        bin_to_contig = bins["chrom"].cat.codes.values
        mask1 = np.isin(bin_to_contig, list(cids1))
        mask2 = np.isin(bin_to_contig, list(cids2))

        # Build a set of valid bin ids for each genome
        bins1 = set(bins.index[mask1].tolist())
        bins2 = set(bins.index[mask2].tolist())

        pix = self.pixels()
        sel = (pix["bin1_id"].isin(bins1) & pix["bin2_id"].isin(bins2)) | (
            pix["bin1_id"].isin(bins2) & pix["bin2_id"].isin(bins1)
        )
        return pix[sel].copy()

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
    # Aggregation
    # ------------------------------------------------------------------
    def aggregate_matrix(
        self,
        level: Literal["contig", "genome"] = "contig",
        sparse: bool = True,
    ) -> sp.csr_matrix | np.ndarray:
        """Aggregate the contact matrix to *contig* or *genome* level.

        Returns a scipy sparse CSR matrix by default, or a dense numpy
        array when *sparse* is ``False``.
        """
        n_contigs = len(self.contigs)
        n_genomes = len(self.genomes)

        if level == "contig":
            n = n_contigs
            bins = self.bins()[["chrom"]]
            mapping = bins["chrom"].cat.codes.values
        elif level == "genome":
            n = n_genomes
            contigs = self.contigs
            bins = self.bins()[["chrom"]]
            # Map bin -> contig -> genome
            bin_to_contig = bins["chrom"].cat.codes.values
            mapping = contigs["genome_id"].values[bin_to_contig]
        else:
            raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

        # Chunked accumulation — stream directly from the cooler selector
        rows, cols, data = [], [], []
        chunk_size = 1_000_000
        n_pix = len(self._clr.pixels())
        for i in range(0, n_pix, chunk_size):
            chunk = self._clr.pixels()[i : i + chunk_size]
            c1 = mapping[chunk["bin1_id"].values]
            c2 = mapping[chunk["bin2_id"].values]
            vals = chunk["count"].values

            # Enforce upper-triangular at the aggregated level
            swap = c1 > c2
            c1[swap], c2[swap] = c2[swap], c1[swap]

            rows.extend(c1.tolist())
            cols.extend(c2.tolist())
            data.extend(vals.tolist())

        agg = sp.coo_matrix((data, (rows, cols)), shape=(n, n))
        agg = agg.tocsr()
        # Symmetrise
        agg = (
            agg
            + agg.T
            - sp.diags(
                agg.diagonal().astype(np.int64), offsets=0, format="csr", dtype=np.int64
            )
        )
        if sparse:
            return agg
        return agg.toarray()

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
        """
        import cooler as clr

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

        # Aggregate pixels
        agg = cast(sp.csr_matrix, self.aggregate_matrix(level=level, sparse=True))
        coo = agg.tocoo()
        pixels = pd.DataFrame(
            {
                "bin1_id": coo.row.astype(np.int64),
                "bin2_id": coo.col.astype(np.int64),
                "count": coo.data.astype(np.int32),
            }
        )
        # Keep only upper triangle
        pixels = pixels[pixels["bin1_id"] <= pixels["bin2_id"]]

        clr.create_cooler(
            str(path),
            bins=bins,
            pixels=pixels,
            ordered=True,
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
