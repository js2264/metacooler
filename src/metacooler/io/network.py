"""Network edge list export from metacool files."""

from __future__ import annotations

import contextlib
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, TextIO, cast

import numpy as np
import pandas as pd

from metacooler.io._utils import (
    aggregate_pixels,
    check_output_path,
    to_upper_triangle,
)
from metacooler.io.weights import read_weights, spurious_cutoff
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")

#: The columns of an edge list, at every level.
_EDGE_COLUMNS = ("node1", "node2", "value")

#: Number of per-chunk groupby frames buffered before they are folded into the
#: running aggregate.  Keeps peak memory at O(FOLD_EVERY x chunk edges) instead
#: of growing one frame per chunk for the whole pixel table.
_FOLD_EVERY = 32


def _fold(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge per-chunk edge frames, re-summing node pairs seen more than once.

    The frames are pixel-shaped (``bin1_id``/``bin2_id``/``count`` standing for
    ``node1``/``node2``/``value``), so this is :func:`aggregate_pixels`: a
    key-based merge whose output is sorted by node pair.

    Args:
        frames: Non-empty list of sorted, duplicate-free edge frames.

    Returns:
        One frame with a single row per node pair, sorted.
    """
    return aggregate_pixels(frames)


def _weighted(
    chunk: pd.DataFrame, weight: str, weights: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(bin1, bin2, value)`` of a pixel chunk, normalized when *weights* is set.

    A normalized value is ``value * weights[bin1] * weights[bin2]``; a pixel of
    a masked bin (``NaN`` weight) is left out.
    """
    b1 = chunk["bin1_id"].to_numpy()
    b2 = chunk["bin2_id"].to_numpy()
    values = chunk[weight].to_numpy()
    if weights is None:
        return b1, b2, values
    values = values * weights[b1] * weights[b2]
    ok = np.isfinite(values)
    return b1[ok], b2[ok], values[ok]


def _chunk_edges(
    chunk: pd.DataFrame,
    bin_to_node: np.ndarray,
    weight: str,
    weights: np.ndarray | None,
) -> pd.DataFrame:
    """The edges of one pixel chunk, summed by node pair (pixel-shaped, as
    :func:`_fold` takes them)."""
    b1, b2, values = _weighted(chunk, weight, weights)
    n1, n2 = to_upper_triangle(bin_to_node[b1], bin_to_node[b2])
    # A node pair sums every pixel between two contigs or genomes, and
    # aggregate_pixels keeps the dtype it is given: an int32 count column
    # (what create writes) wrapped past 2**31 to a negative edge weight.
    # Widen integers as ops/aggregate does; float weights pass unchanged.
    if np.issubdtype(values.dtype, np.integer):
        values = values.astype(np.int64, copy=False)
    return aggregate_pixels(
        [pd.DataFrame({"bin1_id": n1, "bin2_id": n2, "count": values}, copy=False)]
    )


def _as_edges(pixels: pd.DataFrame) -> pd.DataFrame:
    """A pixel-shaped frame of node pairs, its columns named as edges."""
    return pixels.rename(
        columns={"bin1_id": "node1", "bin2_id": "node2", "count": "value"}
    )


def _bin_edges(
    pixels, n_pix: int, chunksize: int, weight: str, weights: np.ndarray | None
) -> pd.DataFrame:
    """Every pixel as an edge, in one frame: ``network(level="bin")``'s return.

    The three columns are allocated once, the pixel count being known, and
    filled chunk by chunk; a frame per chunk, concatenated at the end, held
    every edge twice.  They are allocated once the first chunk is read, so
    that the read's own peak does not come on top of every edge's.  Their
    dtypes are a chunk's -- ``int64`` bin ids, and the weight column's own
    dtype, or ``float64`` when normalized -- and a map without a pixel gets
    those of a chunk of none.  A pixel of a masked bin is left out, so the
    columns may end short of the pixel count.

    Args:
        pixels: The map's pixel selector (``mc._clr.pixels()``).
        n_pix: Its number of pixels.
        chunksize: Pixels read per chunk.
        weight: Pixels column used as the edge weight.
        weights: Normalization weights of the bins (see :func:`_weighted`).

    Returns:
        ``node1``, ``node2``, ``value``, one row per pixel kept, in the
        pixels' order.
    """
    columns: list[np.ndarray] = []
    n = 0
    for lo in range(0, n_pix, chunksize):
        chunk = _weighted(
            cast(pd.DataFrame, pixels[lo : lo + chunksize]), weight, weights
        )
        if not columns:
            columns = [np.empty(n_pix, dtype=values.dtype) for values in chunk]
        k = len(chunk[0])
        for column, values in zip(columns, chunk, strict=True):
            column[n : n + k] = values
        n += k
    if not columns:
        columns = list(_weighted(cast(pd.DataFrame, pixels[0:0]), weight, weights))
    return pd.DataFrame(
        dict(zip(_EDGE_COLUMNS, (column[:n] for column in columns), strict=True)),
        copy=False,
    )


def _aggregate_pairs(
    mc,
    bin_to_node: np.ndarray,
    weight: str,
    n_pix: int,
    chunksize: int,
    weights: np.ndarray | None = None,
) -> pd.DataFrame | None:
    """Sum pixel weights per node pair, folding partial results as we go.

    Args:
        mc: Open :class:`~metacooler.api.MetaCooler`.
        bin_to_node: Array mapping each bin ID to its node ID (contig or genome).
        weight: Pixels column used as the edge weight.
        n_pix: Total number of pixels in the source table.
        chunksize: Number of pixels read per chunk.
        weights: Normalization weights of the bins (see :func:`_weighted`).

    Returns:
        A frame with columns ``node1``/``node2``/``value`` sorted by node pair,
        where ``node1 <= node2``, or ``None`` if the pixel table is empty.
    """
    pending: list[pd.DataFrame] = []
    running: pd.DataFrame | None = None

    for lo in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[lo : lo + chunksize]
        pending.append(_chunk_edges(chunk, bin_to_node, weight, weights))
        if len(pending) >= _FOLD_EVERY:
            running = _fold(pending if running is None else [running, *pending])
            pending = []

    if running is None and not pending:
        return None
    return _as_edges(_fold(pending if running is None else [running, *pending]))


def _spurious_rule(mc, balance: str | None, level: str) -> tuple[float, str]:
    """The spurious-contact cutoff ``network --drop-spurious`` applies.

    Raises:
        ValueError: Without *balance*, at bin level, or when the normalization
            records no usable cutoff.
    """
    if balance is None:
        raise ValueError(
            "--drop-spurious applies the spurious-contact cutoff of a "
            "normalization: pass --balance NAME as well."
        )
    if level == "bin":
        raise ValueError(
            "--drop-spurious drops contig pairs, which bins are not: use "
            "--level contig or --level genome."
        )
    found = spurious_cutoff(mc, balance)
    if found is None or not np.isfinite(found[0]):
        raise ValueError(
            f"/bins/{balance} records no spurious-contact cutoff; "
            "'metacooler normalize normcc' and 'normalize hiczin' store one."
        )
    return found


def _without_spurious(
    edges: pd.DataFrame | None, cutoff: tuple[float, str]
) -> pd.DataFrame | None:
    """Drop the contig pairs below the spurious-contact cutoff.

    The cutoff applies between two *different* contigs (NormCC and HiCzin
    leave a contig's contacts with itself alone).  Under the ``<`` rule a pair
    below the cutoff is spurious; under ``<=``, one at it too.
    """
    if edges is None:
        return None
    threshold, rule = cutoff
    value = edges["value"].to_numpy()
    spurious = value <= threshold if rule == "<=" else value < threshold
    spurious &= edges["node1"].to_numpy() != edges["node2"].to_numpy()
    logger.info(
        "Dropped %d spurious contig pair(s) (normalized contact %s %.6g)",
        int(spurious.sum()),
        rule,
        threshold,
    )
    return edges[~spurious].reset_index(drop=True)


def _regroup(edges: pd.DataFrame, contig_to_genome: np.ndarray) -> pd.DataFrame | None:
    """Sum contig-pair edges into genome-pair edges."""
    if edges.empty:
        return None
    g1, g2 = to_upper_triangle(
        contig_to_genome[edges["node1"].to_numpy()],
        contig_to_genome[edges["node2"].to_numpy()],
    )
    summed = aggregate_pixels(
        [
            pd.DataFrame(
                {"bin1_id": g1, "bin2_id": g2, "count": edges["value"].to_numpy()},
                copy=False,
            )
        ]
    )
    return _as_edges(summed)


def network(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "contig",
    weight: str = "count",
    out: str | Path | TextIO | None = None,
    header: bool = False,
    chunksize: int = DEFAULT_CHUNKSIZE,
    force: bool = False,
    balance: str | None = None,
    drop_spurious: bool = False,
    ids: bool = False,
) -> pd.DataFrame | None:
    """Export a network edge list from a metacool file, or return it.

    The edges are ``node1``, ``node2``, ``value``: contig or genome names at
    ``contig`` and ``genome`` level (their ids with *ids*), one row per pair
    of nodes with ``node1 <= node2`` in the map's order (a node's contacts
    with itself included, once), and bin ids at ``bin`` level, one row per
    pixel.

    Args:
        uri: Path to the source ``.metacool`` file.
        level: Aggregation level — ``bin``, ``contig``, or ``genome``.
        weight: Column name in the pixels table to use as edge weight.
        out: Where to write the edges as TSV: a file path, or an open text
            stream (the ``network`` command passes stdout), which is written
            to and left open.  If ``None``, nothing is written: the edges are
            returned (before metacooler 0.4.0 they were printed, as
            ``out=sys.stdout`` still does).  At ``bin`` level that holds every
            pixel in memory, where a path or a stream gets them chunk by chunk.
        header: If ``True``, write a ``node1\\tnode2\\tvalue`` header line.
        chunksize: Number of pixels to process per chunk.
        force: Overwrite *out* if it already exists.
        balance: Normalize the edges with this weight column of ``/bins``
            (``metacooler normalize``, or any tool's weights stored with
            :func:`metacooler.io.weights.store_weights`): every pixel is
            multiplied by the weights of its two bins before it is summed, and
            pixels of masked bins are left out.
        drop_spurious: Drop the contig pairs whose normalized contact the
            *balance* normalization deems spurious (its stored
            ``spurious_cutoff``; NormCC and HiCzin record one) before the
            edges are summed to genomes.
        ids: Give each contig or genome as its id instead of its name: its
            row in ``/contigs`` or ``/genomes``, the index of
            ``MetaCooler.contigs`` or ``MetaCooler.genomes`` (a row of this
            map, which a ``reassign`` may renumber), as ``int64``, the type of
            the bin ids at ``bin`` level.  The rows and their order are those
            of the names, which ``contigs["name"].to_numpy()[edges["node1"]]``
            gives back.  The edges are summed by id either way: *ids* only
            leaves the names off, and with them a string per node of every
            edge.  At ``bin`` level, whose nodes are bin ids already, it
            changes nothing.

    Raises:
        ValueError: If *weight* column does not exist in the pixels table, if
            *level* is not one of ``bin``/``contig``/``genome``, or if
            *drop_spurious* is set without a *balance* that records a cutoff
            or at ``level="bin"``.
        KeyError: If *balance* names no weight column.
        FileExistsError: If *out* exists and *force* is False.

    Returns:
        The edges, when *out* is ``None``; otherwise ``None``.  With *ids*,
        ``node1`` and ``node2`` are ``int64`` even without an edge.  A map
        without a pixel gives no row, in the dtypes a map with some gives:
        the names as ``/contigs`` or ``/genomes`` hold them (``int64`` ids
        with *ids*), or ``int64`` bin ids; and a ``value`` of the weight
        column's dtype at ``bin`` level, ``int64`` for an integer weight
        summed at the other two, ``float64`` with *balance*.
    """
    from metacooler.api import MetaCooler

    mc = MetaCooler(uri)

    # Validate weight column
    with mc.open() as f:
        pixel_cols = list(f["pixels"].keys())
    if weight not in pixel_cols:
        raise ValueError(
            f"Weight column {weight!r} not found in pixels table. "
            f"Available columns: {pixel_cols}"
        )
    if level not in ("bin", "contig", "genome"):
        raise ValueError(f"level must be 'bin', 'contig', or 'genome', got {level!r}")

    # A file to write, or a stream the caller opened (and closes).
    path: str | Path | None = None
    stream: TextIO | None = None
    if isinstance(out, (str, os.PathLike)):
        path = Path(out)
        check_output_path(path, force, inputs=[uri])
    elif out is not None:
        stream = out

    weights = read_weights(mc, balance) if balance is not None else None
    cutoff = _spurious_rule(mc, balance, level) if drop_spurious else None
    n_pix = len(mc._clr.pixels())

    if level == "bin":
        if out is None:
            edges = _bin_edges(mc._clr.pixels(), n_pix, chunksize, weight, weights)
            logger.info("Generated %d edges at %s level", len(edges), level)
            return edges

        # Raw pixel edge list, chunk by chunk.
        def chunk_edges(lo: int) -> pd.DataFrame:
            chunk = cast(pd.DataFrame, mc._clr.pixels()[lo : lo + chunksize])
            b1, b2, values = _weighted(chunk, weight, weights)
            return pd.DataFrame({"node1": b1, "node2": b2, "value": values})

        frames = (chunk_edges(lo) for lo in range(0, n_pix, chunksize))
        # Streamed to the file.  A stream the caller gave is left open:
        # closing sys.stdout would break any later output in the process.
        cm = open(path, "w") if path is not None else contextlib.nullcontext(stream)  # noqa: SIM115
        with cm as fout:
            assert fout is not None
            if header:
                fout.write("node1\tnode2\tvalue\n")
            for df in frames:
                fout.write(df.to_csv(sep="\t", index=False, header=False))
        return None

    # Contig and genome levels: map every bin to its node, then reduce, in ids
    # throughout -- a bin's contig from indexes/chrom_offset, a contig's genome
    # from /contigs/genome_id, neither the bin table nor a name -- so that the
    # names are read and put on at the end, and only without *ids*.
    bin_to_contig = mc._bin_contig_ids
    contig_to_genome = mc._contig_genome_ids
    bin_to_node = (
        bin_to_contig if level == "contig" else contig_to_genome[bin_to_contig]
    )

    if cutoff is None:
        agg = _aggregate_pairs(mc, bin_to_node, weight, n_pix, chunksize, weights)
    else:
        agg = _without_spurious(
            _aggregate_pairs(mc, bin_to_contig, weight, n_pix, chunksize, weights),
            cutoff,
        )
        if agg is not None and level == "genome":
            agg = _regroup(agg, contig_to_genome)
    if agg is None:
        # No edge: an empty pixel chunk's, so that the columns keep the types
        # edges would have.
        empty = cast(pd.DataFrame, mc._clr.pixels()[:0])
        agg = _as_edges(_chunk_edges(empty, bin_to_node, weight, weights))
    if ids:
        edges = agg
    else:
        names = (mc.contigs if level == "contig" else mc.genomes)["name"].values
        edges = pd.DataFrame(
            {
                "node1": names[agg["node1"].values],
                "node2": names[agg["node2"].values],
                "value": agg["value"].values,
            }
        )

    logger.info("Generated %d edges at %s level", len(edges), level)

    if out is None:
        return edges
    # Write TSV (header optional for compatibility with igraph, networkx, etc.)
    edges.to_csv(
        path if path is not None else stream, sep="\t", index=False, header=header
    )
    return None
