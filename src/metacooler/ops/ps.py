"""P(s): contact probability as a function of genomic distance.

The shape of P(s) is the single most informative quality plot for a Hi-C
library — a proper experiment decays smoothly over orders of magnitude, while a
library dominated by religation or by undigested fragments shows a flat or
bumpy short-range head — and it is measurable on a metagenome without knowing
anything about the genomes involved, because it only uses contacts whose two
ends sit on the *same contig*, where a genomic distance is defined at all.
Contacts between two contigs (of the same MAG or not) have no distance and are
counted by :mod:`metacooler.qc` instead.

Two numbers per distance bin make the curve comparable across libraries:

``n_contacts``
    Contacts observed at that separation.
``n_bin_pairs``
    Pairs of bins that *could* have contributed one, given the contig lengths.
    Without it, a curve reports the shrinking number of long-range bin pairs
    (a contig of *n* bins offers ``n - d`` pairs at separation *d*, so the tail
    is starved by construction) as if it were a decay of contact probability.

``avg_contacts`` is their ratio, and ``prob`` normalises it so the curve
integrates to 1 over the sampled range — which is what puts two libraries of
different depth on the same axes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.api import MetaCooler
from metacooler.io._utils import (
    check_output_path,
    chrom_codes,
    resolve_contig_selection,
)
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

logger = logging.getLogger("metacooler")

# Distance bins per decade when the caller does not say.  Ten is the usual
# resolution for a log-log P(s): fine enough to see the shoulder where
# short-range polymer contacts give way to genuine 3D signal, coarse enough
# that a shallow library still has counts in every bin.
DEFAULT_BINS_PER_DECADE = 10

# Curves drawn with a legend.  Past this many genomes the legend covers the
# plot, so the lines are drawn without one.
_MAX_LEGEND_ENTRIES = 20

# Column layout of the returned table, in order.
PS_COLUMNS = (
    "group",
    "dist_min",
    "dist_max",
    "dist_mid",
    "n_contacts",
    "n_bin_pairs",
    "avg_contacts",
    "prob",
)


def _distance_edges(
    max_distance: int,
    min_distance: int,
    bins_per_decade: int,
) -> np.ndarray:
    """Build log-spaced distance-bin edges, with a leading ``[0, min)`` bin.

    The leading bin is what holds the diagonal (both ends in the same bin,
    separation 0), which no log-spaced grid can express and which is the single
    largest value on most curves.

    Args:
        max_distance: Largest separation to cover, in bp (exclusive upper edge
            is one past it).
        min_distance: Upper edge of the leading bin, i.e. where the log grid
            starts.
        bins_per_decade: Log-grid resolution.

    Returns:
        Strictly increasing ``int64`` edges, ``edges[0] == 0``.

    Raises:
        ValueError: If the arguments do not describe at least one bin.
    """
    if max_distance < 0:
        raise ValueError(f"max_distance must be >= 0, got {max_distance}")
    if min_distance < 1:
        raise ValueError(f"min_distance must be >= 1, got {min_distance}")
    if bins_per_decade < 1:
        raise ValueError(f"bins_per_decade must be >= 1, got {bins_per_decade}")
    top = max_distance + 1
    if top <= min_distance:
        return np.array([0, top], dtype=np.int64)
    n_decades = np.log10(top / min_distance)
    n_edges = max(int(np.ceil(n_decades * bins_per_decade)), 1) + 1
    grid = np.logspace(np.log10(min_distance), np.log10(top), n_edges)
    edges = np.unique(np.concatenate([[0.0], np.round(grid)])).astype(np.int64)
    # np.round can collapse two neighbouring log edges onto the same integer at
    # short distances; np.unique above already dropped the duplicates.
    return edges


def _bin_pair_counts(
    key: np.ndarray,
    groups: np.ndarray,
    n_groups: int,
    edges: np.ndarray,
) -> np.ndarray:
    """Count, per group and distance bin, the bin pairs that *could* contact.

    Uses the same folded key trick as the rest of the codebase: with the contig
    id multiplied above every coordinate, ``(contig, midpoint)`` ordering is
    plain integer ordering, so one global ``searchsorted`` per distance edge
    answers "how many later bins of *my own* contig are within this distance"
    for every bin at once — no Python loop over contigs, and no per-contig
    O(n²) enumeration of pairs.

    Only the *selected* bins are passed in.  A selection is always whole
    contigs, so no partner of a selected bin is ever missing from the array.

    Args:
        key: ``contig_id * K + midpoint`` of every selected bin, ascending.
        groups: Output-row group of every selected bin.
        n_groups: Number of output groups.
        edges: Distance-bin edges from :func:`_distance_edges`.

    Returns:
        ``(n_groups, len(edges) - 1)`` float64 array of possible bin pairs.
    """
    index = np.arange(len(key), dtype=np.int64)

    # Cumulative counts at each edge: pairs strictly closer than the edge.
    cumulative = np.zeros((n_groups, len(edges)), dtype=np.float64)
    for i, edge in enumerate(edges):
        if edge == 0:
            continue  # no pair is closer than 0 bp
        # First bin whose key is >= key_j + edge; the fold guarantees it never
        # reaches into the next contig.
        upper = np.searchsorted(key, key + np.int64(edge), side="left")
        partners = (upper - index - 1).astype(np.float64)
        cumulative[:, i] = np.bincount(
            groups,
            weights=partners,
            minlength=n_groups,
        )[:n_groups]

    counts = np.diff(cumulative, axis=1)
    # The diagonal: one "pair" per bin, at separation 0, which the strictly-
    # increasing-key count above cannot see (a bin is not its own partner).
    counts[:, 0] += np.bincount(groups, minlength=n_groups)[:n_groups]
    return counts


def contact_probability(
    uri: str | Path,
    *,
    by: Literal["all", "genome"] = "all",
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
    filters: Sequence[str] | None = None,
    bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
    min_distance: int | None = None,
    max_distance: int | None = None,
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> pd.DataFrame:
    """Compute P(s) from the intra-contig contacts of a metacool file.

    Args:
        uri: ``file.metacool`` path or ``file.mmetacool::entries/name`` URI.
        by: ``"all"`` pools every selected contig into one curve;
            ``"genome"`` emits one curve per genome (MAG).
        genomes: Restrict to these genomes' contigs.
        contigs: Restrict to these contigs.
        filters: ``--filter`` expressions narrowing the selection further
            (see :mod:`metacooler.filters`).
        bins_per_decade: Log-grid resolution.
        min_distance: Upper edge of the leading ``[0, min)`` bin, in bp.
            Defaults to the file's bin size (or the median bin width on a
            variable-width file), since nothing finer is resolvable.
        max_distance: Largest separation to report, in bp.  Defaults to the
            longest intra-contig separation available in the selection.
        chunksize: Pixels read per iteration.

    Returns:
        A long-format DataFrame with the :data:`PS_COLUMNS` columns, one row
        per (group, distance bin), sorted by group then distance.  Distance
        bins with no possible bin pair are omitted.

    Raises:
        ValueError: If *by* is not a known grouping, or if the selection holds
            no bins.
    """
    if by not in ("all", "genome"):
        raise ValueError(f"by must be 'all' or 'genome', got {by!r}")

    mc = MetaCooler(uri)
    bins = mc.bins()
    contig_of_bin = chrom_codes(bins)
    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)
    mids = (starts + ends) // 2

    selected_contigs = resolve_contig_selection(
        mc, genomes=genomes, contigs=contigs, filters=filters
    )
    if selected_contigs is None:
        selected = np.ones(len(bins), dtype=bool)
    else:
        selected = np.isin(contig_of_bin, np.fromiter(selected_contigs, dtype=np.int64))
    if not selected.any():
        raise ValueError("The selection holds no bins; nothing to compute P(s) from.")

    genome_of_bin = mc.contigs["genome_id"].to_numpy(dtype=np.int64)[contig_of_bin]
    if by == "genome":
        group_of_bin = genome_of_bin
        group_names = mc.genomes["name"].astype(str).to_numpy()
    else:
        group_of_bin = np.zeros(len(bins), dtype=np.int64)
        group_names = np.array(["all"])
    n_groups = len(group_names)

    # Distance grid. The widest separation any contig can offer is the span of
    # its own midpoints, so anything beyond it would be an empty bin.
    span = _max_intra_contig_span(contig_of_bin, mids, selected)
    if min_distance is None:
        binsize = mc.binsize
        widths = ends - starts
        min_distance = int(binsize) if binsize else max(int(np.median(widths)), 1)
    if max_distance is None:
        max_distance = int(span)
    edges = _distance_edges(int(max_distance), int(min_distance), bins_per_decade)

    # Fold the contig id in above every coordinate so one global searchsorted
    # answers a per-contig question (see _bin_pair_counts).
    fold = int(mids.max()) + int(edges[-1]) + 1 if len(mids) else 1
    key = (contig_of_bin.astype(np.int64) * fold + mids)[selected]

    logger.info(
        "P(s) over %s selected bin(s) in %d group(s), on %d distance bin(s) "
        "up to %s bp.",
        f"{int(selected.sum()):,}",
        n_groups,
        len(edges) - 1,
        f"{int(edges[-1]):,}",
    )
    possible = _bin_pair_counts(key, group_of_bin[selected], n_groups, edges)
    observed = _observed_counts(
        mc, contig_of_bin, mids, group_of_bin, selected, edges, n_groups, chunksize
    )

    return _assemble(group_names, edges, observed, possible)


def _max_intra_contig_span(
    contig_of_bin: np.ndarray, mids: np.ndarray, selected: np.ndarray
) -> int:
    """Largest midpoint span within a single selected contig, in bp."""
    if not selected.any():
        return 0
    cids = contig_of_bin[selected]
    values = mids[selected]
    order = np.lexsort((values, cids))
    cids, values = cids[order], values[order]
    unique = np.unique(cids)
    first = np.searchsorted(cids, unique, side="left")
    last = np.searchsorted(cids, unique, side="right") - 1
    return int((values[last] - values[first]).max(initial=0))


def _observed_counts(
    mc: MetaCooler,
    contig_of_bin: np.ndarray,
    mids: np.ndarray,
    group_of_bin: np.ndarray,
    selected: np.ndarray,
    edges: np.ndarray,
    n_groups: int,
    chunksize: int,
) -> np.ndarray:
    """Sum pixel counts per (group, distance bin) over the intra-contig pixels."""
    n_dbins = len(edges) - 1
    totals = np.zeros(n_groups * n_dbins, dtype=np.float64)
    pixels = mc._clr.pixels()
    n_pixels = len(pixels)
    for lo in range(0, n_pixels, chunksize):
        chunk = pixels[lo : lo + chunksize]
        b1 = chunk["bin1_id"].to_numpy(dtype=np.int64)
        b2 = chunk["bin2_id"].to_numpy(dtype=np.int64)
        keep = (contig_of_bin[b1] == contig_of_bin[b2]) & selected[b1]
        if not keep.any():
            continue
        b1, b2 = b1[keep], b2[keep]
        distance = np.abs(mids[b2] - mids[b1])
        dbin = np.searchsorted(edges, distance, side="right") - 1
        # A separation at or beyond the top edge only happens when the caller
        # capped --max-distance; those contacts are outside the reported range.
        inside = (dbin >= 0) & (dbin < n_dbins)
        if not inside.any():
            continue
        flat = group_of_bin[b1[inside]] * n_dbins + dbin[inside]
        totals += np.bincount(
            flat,
            weights=chunk["count"].to_numpy()[keep][inside].astype(np.float64),
            minlength=n_groups * n_dbins,
        )[: n_groups * n_dbins]
    return totals.reshape(n_groups, n_dbins)


def _assemble(
    group_names: np.ndarray,
    edges: np.ndarray,
    observed: np.ndarray,
    possible: np.ndarray,
) -> pd.DataFrame:
    """Turn the two count matrices into the long-format output table."""
    n_groups, n_dbins = observed.shape
    lo = edges[:-1]
    hi = edges[1:]
    # Geometric centre for a log axis; the leading [0, min) bin has no geometric
    # centre, so it gets the arithmetic one.
    mid = np.where(lo > 0, np.sqrt(np.maximum(lo, 1) * hi), hi / 2.0)
    width = (hi - lo).astype(np.float64)

    frame = pd.DataFrame(
        {
            "group": np.repeat(group_names, n_dbins),
            "dist_min": np.tile(lo, n_groups),
            "dist_max": np.tile(hi, n_groups),
            "dist_mid": np.tile(mid, n_groups),
            "n_contacts": observed.ravel(),
            "n_bin_pairs": possible.ravel(),
        }
    )
    frame["_width"] = np.tile(width, n_groups)
    frame = frame[frame["n_bin_pairs"] > 0].copy()
    frame["avg_contacts"] = frame["n_contacts"] / frame["n_bin_pairs"]

    # Normalise each curve to unit area so libraries of different depth are
    # comparable. A group with no contact at all keeps NaN rather than a
    # fabricated flat line.
    area = (frame["avg_contacts"] * frame["_width"]).groupby(frame["group"]).sum()
    frame["prob"] = frame["avg_contacts"] / frame["group"].map(
        area.replace(0.0, np.nan)
    )
    frame = frame.drop(columns="_width")

    frame["n_contacts"] = frame["n_contacts"].round().astype(np.int64)
    frame["n_bin_pairs"] = frame["n_bin_pairs"].round().astype(np.int64)
    return frame[list(PS_COLUMNS)].reset_index(drop=True)


def plot_ps(
    table: pd.DataFrame,
    out: str | Path,
    *,
    value: str = "prob",
    title: str | None = None,
    dpi: int = 300,
    force: bool = False,
) -> None:
    """Render a P(s) table as a log-log curve, one line per group.

    Args:
        table: Output of :func:`contact_probability`.
        out: Image path; the extension picks the format (``.png``, ``.pdf``,
            ``.svg``).
        value: Column to plot on the y axis — ``"prob"`` (unit-area normalised)
            or ``"avg_contacts"`` (raw average per bin pair).
        title: Figure title.
        dpi: Raster resolution.
        force: Overwrite an existing file.

    Raises:
        FileExistsError: If *out* exists and *force* is False.
        ValueError: If *value* is not a column of *table*.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if value not in table.columns:
        raise ValueError(f"{value!r} is not a column of the P(s) table.")
    check_output_path(out, force)

    groups = list(dict.fromkeys(table["group"]))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    handles: list = []
    labels: list[str] = []
    for name in groups:
        curve = table[table["group"] == name]
        # A log axis cannot show a zero, and an empty distance bin is a gap in
        # the curve rather than a probability of zero.
        drawable = curve[(curve[value] > 0) & (curve["dist_mid"] > 0)]
        if drawable.empty:
            continue
        (line,) = ax.plot(
            drawable["dist_mid"], drawable[value], label=str(name), linewidth=1.2
        )
        handles.append(line)
        labels.append(str(name))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genomic distance (bp)")
    ax.set_ylabel("P(s)" if value == "prob" else "Mean contacts per bin pair")
    ax.set_title(title or "Contact probability vs genomic distance")
    if 1 < len(handles) <= _MAX_LEGEND_ENTRIES:
        # Handles are passed explicitly: matplotlib hides any label starting
        # with an underscore from an automatic legend, and ``__unbinned__`` —
        # the genome every unassigned contig lands in — is exactly that.
        ax.legend(handles, labels, fontsize="small", frameon=False)
    elif len(handles) > _MAX_LEGEND_ENTRIES:
        logger.info(
            "%d curves drawn without a legend (more than %d).",
            len(groups),
            _MAX_LEGEND_ENTRIES,
        )
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    logger.info("Wrote P(s) plot to %s", out)


def ps(
    uri: str | Path,
    *,
    out: str | Path | None = None,
    plot: str | Path | None = None,
    by: Literal["all", "genome"] = "all",
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
    filters: Sequence[str] | None = None,
    bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
    min_distance: int | None = None,
    max_distance: int | None = None,
    value: str = "prob",
    header: bool = True,
    dpi: int = 300,
    chunksize: int = DEFAULT_CHUNKSIZE,
    force: bool = False,
) -> pd.DataFrame:
    """Compute P(s), optionally writing the table and/or the plot.

    Args:
        uri: Source metacool URI.
        out: TSV destination.  ``None`` writes nothing (the table is still
            returned); the CLI passes stdout explicitly.
        plot: Image destination, or ``None`` for no figure.
        by: ``"all"`` or ``"genome"``.
        genomes: Restrict to these genomes.
        contigs: Restrict to these contigs.
        filters: ``--filter`` expressions.
        bins_per_decade: Log-grid resolution.
        min_distance: Upper edge of the leading bin, in bp.
        max_distance: Largest separation to report, in bp.
        value: Column plotted on the y axis.
        header: Write a header row into the TSV.
        dpi: Plot resolution.
        chunksize: Pixels read per iteration.
        force: Overwrite existing outputs.

    Returns:
        The P(s) table.
    """
    table = contact_probability(
        uri,
        by=by,
        genomes=genomes,
        contigs=contigs,
        filters=filters,
        bins_per_decade=bins_per_decade,
        min_distance=min_distance,
        max_distance=max_distance,
        chunksize=chunksize,
    )
    if out is not None:
        check_output_path(out, force)
        table.to_csv(out, sep="\t", index=False, header=header)
        logger.info("Wrote P(s) table (%d rows) to %s", len(table), out)
    if plot is not None:
        plot_ps(
            table,
            Path(plot),
            value=value,
            title=f"P(s) — {Path(str(uri).split('::')[0]).name}",
            dpi=dpi,
            force=force,
        )
    return table
