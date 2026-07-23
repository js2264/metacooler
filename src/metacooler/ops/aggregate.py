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

    rows_list, cols_list, data_list = [], [], []
    n_pix = len(mc._clr.pixels())
    for i in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[i : i + chunksize]
        c1 = mapping[chunk["bin1_id"].values]
        c2 = mapping[chunk["bin2_id"].values]
        # int64 so scipy's duplicate-summing in .tocsr() cannot overflow int32
        # (a genome pair can exceed 2**31 contacts on real metagenomic data).
        vals = chunk["count"].values.astype(np.int64)

        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap], c1[swap]

        rows_list.append(c1)
        cols_list.append(c2)
        data_list.append(vals)

    rows = np.concatenate(rows_list) if rows_list else np.array([], dtype=np.int64)
    cols = np.concatenate(cols_list) if cols_list else np.array([], dtype=np.int64)
    data = np.concatenate(data_list) if data_list else np.array([], dtype=np.int64)
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
    count_col: str = "count",
) -> Generator[pd.DataFrame, None, None]:
    """Yield aggregated pixel chunks without building a sparse matrix.

    Maps bin IDs to contig or genome IDs, enforces the upper triangle, and
    aggregates within each chunk.  Cross-chunk duplicates are resolved by the
    caller (pass to ``cooler.create_cooler`` with ``ordered=False``).

    Args:
        mc: Open MetaCooler instance.
        level: Aggregation level — ``"contig"`` or ``"genome"``.
        chunksize: Number of pixels to process per chunk.
        count_col: Pixel column to use as contact count.

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
        vals = chunk[count_col].values
        # Promote integer counts to int64 (a genome pair can exceed 2**31
        # contacts on real data); preserve float weights (e.g. balanced
        # values) instead of truncating them to zero.
        if np.issubdtype(vals.dtype, np.integer):
            vals = vals.astype(np.int64)

        # Enforce upper triangle
        swap = c1 > c2
        c1[swap], c2[swap] = c2[swap], c1[swap]

        df = pd.DataFrame(
            {
                "bin1_id": c1.astype(np.int64),
                "bin2_id": c2.astype(np.int64),
                "count": vals,
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
    either contigs or genomes.  Pixels are streamed chunk-by-chunk; the
    full contact matrix is never materialised in memory.

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.cool`` file.
        level: Aggregation level — ``"contig"`` or ``"genome"``.
        count_col: Pixel column to use as contact count (default: ``"count"``).
        chunksize: Number of pixels to process per chunk.
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

    bins_out = pd.DataFrame(
        {
            "chrom": chroms["chrom"].values,
            "start": 0,
            "end": chroms["length"].values,
        }
    )

    # Store an int64 count (integer contacts) or float64 (weighted columns) so
    # cooler's default int32 cannot silently saturate large genome-pair totals.
    src_dtype = mc._clr.pixels()[:0][count_col].dtype
    count_dtype = np.int64 if np.issubdtype(src_dtype, np.integer) else np.float64

    cooler.create_cooler(
        str(output),
        bins=bins_out,
        pixels=_aggregate_pixel_stream(
            mc, level=level, chunksize=chunksize, count_col=count_col
        ),
        dtypes={"count": count_dtype},
        ordered=False,
        symmetric_upper=True,
        assembly=mc.info.get("assembly"),
        metadata=mc.metadata if mc.metadata else None,
    )
