"""Aggregation operations for metacool files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import numpy as np
import pandas as pd
from scipy import sparse

from metacooler.api import MetaCooler

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal


def _aggregate_sparse(
    mc: MetaCooler,
    level: Literal["contig", "genome"] = "contig",
    chunksize: int = 1_000_000,
) -> sparse.csr_matrix:
    """Return the aggregated contact matrix as a sparse CSR matrix.

    Args:
        mc: Open MetaCooler instance.
        level: Aggregation level — ``"contig"`` or ``"genome"``.
        chunksize: Number of pixels to process per chunk.

    Returns:
        Symmetrised sparse CSR matrix of shape ``(n, n)`` where ``n`` is
        the number of contigs or genomes at the requested level.
    """
    n_contigs = len(mc.contigs)
    n_genomes = len(mc.genomes)

    if level == "contig":
        n = n_contigs
        bins = mc.bins()[["chrom"]]
        mapping = bins["chrom"].cat.codes.values
    elif level == "genome":
        n = n_genomes
        contigs = mc.contigs
        bins = mc.bins()[["chrom"]]
        bin_to_contig = bins["chrom"].cat.codes.values
        mapping = contigs["genome_id"].values[bin_to_contig]
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    rows, cols, data = [], [], []
    n_pix = len(mc._clr.pixels())
    for i in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[i : i + chunksize]
        c1 = mapping[chunk["bin1_id"].values]
        c2 = mapping[chunk["bin2_id"].values]
        vals = chunk["count"].values

        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap], c1[swap]

        rows.extend(c1.tolist())
        cols.extend(c2.tolist())
        data.extend(vals.tolist())

    agg = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    agg = agg.tocsr()
    agg = (
        agg
        + agg.T
        - sparse.diags(
            agg.diagonal().astype(np.int64), offsets=0, format="csr", dtype=np.int64
        )
    )
    return agg


def _aggregate_pixel_stream(
    mc: MetaCooler,
    level: Literal["contig", "genome"] = "contig",
    chunksize: int = 1_000_000,
) -> Generator[pd.DataFrame, None, None]:
    """Yield aggregated pixel chunks without building a sparse matrix.

    Maps bin IDs to contig or genome IDs, enforces the upper triangle, and
    aggregates within each chunk.  Cross-chunk duplicates are resolved by the
    caller (pass to ``cooler.create_cooler`` with ``ordered=False``).

    Args:
        mc: Open MetaCooler instance.
        level: Aggregation level — ``"contig"`` or ``"genome"``.
        chunksize: Number of pixels to process per chunk.

    Yields:
        DataFrames with columns ``bin1_id``, ``bin2_id``, ``count``
        (upper-triangle only, aggregated within each chunk).
    """
    n_contigs = len(mc.contigs)
    n_genomes = len(mc.genomes)

    if level == "contig":
        bins = mc.bins()[["chrom"]]
        mapping = bins["chrom"].cat.codes.values.copy()
    elif level == "genome":
        contigs = mc.contigs
        bins = mc.bins()[["chrom"]]
        bin_to_contig = bins["chrom"].cat.codes.values
        mapping = contigs["genome_id"].values[bin_to_contig].copy()
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    del n_contigs, n_genomes  # not needed here

    n_pix = len(mc._clr.pixels())
    for i in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[i : i + chunksize]
        c1 = mapping[chunk["bin1_id"].values].copy()
        c2 = mapping[chunk["bin2_id"].values].copy()
        vals = chunk["count"].values

        # Enforce upper triangle
        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap], c1[swap]

        df = pd.DataFrame(
            {
                "bin1_id": c1.astype(np.int64),
                "bin2_id": c2.astype(np.int64),
                "count": vals.astype(np.int32),
            }
        )
        # Aggregate within-chunk duplicates to reduce cooler's external sort load
        df = df.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()
        yield df[["bin1_id", "bin2_id", "count"]]


def aggregate(
    uri: str | Path,
    output: str | Path,
    level: Literal["contig", "genome"] = "contig",
    count_col: str = "count",
    chunksize: int = 1_000_000,
) -> None:
    """Aggregate a metacool file to contig or genome level.

    The output is a standard ``.cool`` file where the "chromosomes" are
    either contigs or genomes.
    """
    mc = MetaCooler(uri)

    if level == "contig":
        chroms = mc.contigs[["name", "length"]].rename(
            columns={"name": "chrom", "length": "length"}
        )
    elif level == "genome":
        chroms = mc.genomes[["name", "length"]].rename(
            columns={"name": "chrom", "length": "length"}
        )
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    agg = _aggregate_sparse(mc, level=level, chunksize=chunksize)
    coo = agg.tocoo()

    pixels = pd.DataFrame(
        {
            "bin1_id": coo.row.astype(np.int64),
            "bin2_id": coo.col.astype(np.int64),
            "count": coo.data.astype(np.int32),
        }
    )
    pixels = pixels[pixels["bin1_id"] <= pixels["bin2_id"]]

    bins_out = pd.DataFrame(
        {
            "chrom": chroms["chrom"].values,
            "start": 0,
            "end": chroms["length"].values,
        }
    )

    cooler.create_cooler(
        str(output),
        bins=bins_out,
        pixels=pixels,
        ordered=True,
        symmetric_upper=True,
    )
