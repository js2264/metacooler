"""Aggregation operations for metacool files.

All aggregation shares two primitives:

- :func:`_bin_to_level_mapping` maps each bin id to a contig or genome id;
- :func:`_write_aggregated_cool` streams the aggregated pixels into a ``.cool``.

``_aggregate_sparse`` (dense/sparse, used by ``show``) and the streaming
``aggregate`` / :meth:`MetaCooler.to_cooler` all build on these so the
int-overflow handling and the bin→level mapping live in exactly one place.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import sparse

from metacooler.api import MetaCooler
from metacooler.io._utils import chrom_codes, to_upper_triangle, write_cooler
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal


def _bin_to_level_mapping(
    mc: MetaCooler, level: Literal["contig", "genome"]
) -> np.ndarray:
    """Map each bin id to a contig id (``contig``) or genome id (``genome``).

    Returns a fresh ``int`` array of length ``n_bins`` (safe for callers to
    mutate).  Raises ``ValueError`` for an unknown *level*.
    """
    codes = chrom_codes(mc.bins()[["chrom"]])
    if level == "contig":
        return codes.copy()
    if level == "genome":
        return mc.contigs["genome_id"].values[codes].copy()
    raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")


def _level_chroms(mc: MetaCooler, level: Literal["contig", "genome"]) -> pd.DataFrame:
    """Return the output "chromosomes" (``name``, ``length``) for *level*."""
    if level == "contig":
        return mc.contigs[["name", "length"]]
    if level == "genome":
        return mc.genomes[["name", "length"]]
    raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")


def _aggregate_sparse(
    mc: MetaCooler,
    level: Literal["contig", "genome"] = "contig",
    chunksize: int = DEFAULT_CHUNKSIZE,
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
    mapping = _bin_to_level_mapping(mc, level)  # also validates *level*
    n = len(mc.contigs) if level == "contig" else len(mc.genomes)

    # Reduce each chunk into a running aggregate instead of buffering every
    # chunk: peak memory stays bounded by the output matrix plus one chunk,
    # rather than growing with the whole pixel table (~24 bytes/pixel).
    agg = sparse.csr_matrix((n, n), dtype=np.int64)
    pixels = mc._clr.pixels()
    n_pix = len(pixels)
    for i in range(0, n_pix, chunksize):
        chunk = pixels[i : i + chunksize]
        c1 = mapping[chunk["bin1_id"].values]
        c2 = mapping[chunk["bin2_id"].values]
        # int64 so scipy's duplicate-summing in .tocsr() cannot overflow int32
        # (a genome pair can exceed 2**31 contacts on real metagenomic data).
        vals = chunk["count"].values.astype(np.int64)

        c1, c2 = to_upper_triangle(c1, c2)

        # .tocsr() sums the duplicates within this chunk; the addition sums
        # across chunks, so the total is identical to concatenating first.
        agg = agg + sparse.coo_matrix((vals, (c1, c2)), shape=(n, n)).tocsr()

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
    chunksize: int = DEFAULT_CHUNKSIZE,
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
    mapping = _bin_to_level_mapping(mc, level)  # also validates *level*

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
        c1, c2 = to_upper_triangle(c1, c2)

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


def _write_aggregated_cool(
    mc: MetaCooler,
    output: str | Path,
    level: Literal["contig", "genome"] = "contig",
    count_col: str = "count",
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> None:
    """Stream an aggregated ``.cool`` for an already-open MetaCooler.

    Shared by :func:`aggregate` (opens a URI) and
    :meth:`metacooler.api.MetaCooler.to_cooler` (already holds ``mc``) so the
    bins construction, int-dtype handling and ``create_cooler`` call live once.
    """
    chroms = _level_chroms(mc, level)  # also validates *level*
    bins_out = pd.DataFrame(
        {
            "chrom": chroms["name"].values,
            "start": 0,
            "end": chroms["length"].values,
        }
    )

    # Store an int64 count (integer contacts) or float64 (weighted columns) so
    # cooler's default int32 cannot silently saturate large genome-pair totals.
    src_dtype = mc._clr.pixels()[:0][count_col].dtype
    count_dtype = np.int64 if np.issubdtype(src_dtype, np.integer) else np.float64

    write_cooler(
        output,
        bins_out,
        _aggregate_pixel_stream(
            mc, level=level, chunksize=chunksize, count_col=count_col
        ),
        dtypes={"count": count_dtype},
        ordered=False,
        symmetric_upper=True,
        assembly=mc.info.get("assembly"),
        metadata=mc.metadata if mc.metadata else None,
    )


def aggregate(
    uri: str | Path,
    output: str | Path,
    level: Literal["contig", "genome"] = "contig",
    count_col: str = "count",
    chunksize: int = DEFAULT_CHUNKSIZE,
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
    _write_aggregated_cool(
        MetaCooler(uri),
        output,
        level=level,
        count_col=count_col,
        chunksize=chunksize,
    )
