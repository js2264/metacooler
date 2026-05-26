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
    from typing import Literal


def aggregate(
    uri: str | Path,
    output: str | Path,
    level: Literal["contig", "genome"] = "contig",
    count_col: str = "count",
    chunksize: int = 10_000_000,
) -> None:
    """Aggregate a metacool file to contig or genome level.

    The output is a standard ``.cool`` file where the "chromosomes" are
    either contigs or genomes.
    """
    mc = MetaCooler(uri)
    n_contigs = len(mc.contigs)
    n_genomes = len(mc.genomes)

    if level == "contig":
        n = n_contigs
        bins = mc.bins()[["chrom"]]
        mapping = bins["chrom"].cat.codes.values
        chroms = mc.contigs[["name", "length"]].rename(
            columns={"name": "chrom", "length": "length"}
        )
    elif level == "genome":
        n = n_genomes
        contigs = mc.contigs
        bins = mc.bins()[["chrom"]]
        bin_to_contig = bins["chrom"].cat.codes.values
        mapping = contigs["genome_id"].values[bin_to_contig]
        chroms = mc.genomes[["name", "length"]].rename(
            columns={"name": "chrom", "length": "length"}
        )
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    # Chunked accumulation — stream directly from the cooler selector
    rows, cols, data = [], [], []
    n_pix = len(mc._clr.pixels())
    for i in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[i : i + chunksize]
        c1 = mapping[chunk["bin1_id"].values]
        c2 = mapping[chunk["bin2_id"].values]
        vals = chunk[count_col].values

        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap], c1[swap]

        rows.extend(c1.tolist())
        cols.extend(c2.tolist())
        data.extend(vals.tolist())

    agg = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    agg = agg.tocsr()
    # Symmetrise: subtract diagonal to avoid double-counting (matches api.aggregate_matrix)
    agg = (
        agg
        + agg.T
        - sparse.diags(
            agg.diagonal().astype(np.int64), offsets=0, format="csr", dtype=np.int64
        )
    )
    coo = agg.tocoo()

    pixels = pd.DataFrame(
        {
            "bin1_id": coo.row.astype(np.int64),
            "bin2_id": coo.col.astype(np.int64),
            "count": coo.data.astype(np.int32),
        }
    )
    pixels = pixels[pixels["bin1_id"] <= pixels["bin2_id"]]

    # One bin per chromosome at this level
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
