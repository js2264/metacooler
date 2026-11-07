"""P(s) for a metacool file that already exists.

The grid, the accumulator, the denominators and the output table live in
:mod:`metacooler.ps`; this module is the file-facing half — it opens a map,
walks its pixel table (or re-reads a pairs file), and writes the result.

Two sources, one table:

- **pixels** (default): what the map itself holds.  Intra-contig distances are
  gaps between bin midpoints, so nothing below the bin size is resolvable.
- **pairs** (``--pairs``): the original contacts, at base-pair resolution.  The
  map is still needed — it supplies the contig→genome assignment, the contig
  lengths, and any ``--genomes``/``--contigs``/``--filter`` selection — but the
  curve is built from read positions, so it runs from 1 bp up.
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
from metacooler.ps import (
    DEFAULT_BINS_PER_DECADE,
    PS_COLUMNS,
    PSCounts,
    assemble,
    bin_pair_denominators,
    distance_edges,
    pairs_distance_grid,
    position_pair_denominators,
    tally_intra_contig,
    write_ps_table,
)
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

logger = logging.getLogger("metacooler")

# Curves drawn with a legend.  Past this many genomes the legend covers the
# plot, so the lines are drawn without one.
_MAX_LEGEND_ENTRIES = 20

__all__ = [
    "DEFAULT_BINS_PER_DECADE",
    "PS_COLUMNS",
    "contact_probability",
    "plot_ps",
    "ps",
]


class _Selection:
    """What a P(s) run is computed over: which contigs, grouped how.

    Resolving this once keeps the pixel and pairs paths from drifting apart on
    the two questions that decide what the curves mean — which contigs take
    part, and which row each one contributes to.
    """

    def __init__(
        self,
        mc: MetaCooler,
        by: Literal["all", "genome"],
        genomes: Sequence[str] | None,
        contigs: Sequence[str] | None,
        filters: Sequence[str] | None,
    ) -> None:
        if by not in ("all", "genome"):
            raise ValueError(f"by must be 'all' or 'genome', got {by!r}")
        self.mc = mc
        contig_table = mc.contigs
        self.n_contigs = len(contig_table)
        self.contig_genome = contig_table["genome_id"].to_numpy(dtype=np.int64)
        self.contig_length = contig_table["length"].to_numpy(dtype=np.int64)

        selected_ids = resolve_contig_selection(
            mc, genomes=genomes, contigs=contigs, filters=filters
        )
        if selected_ids is None:
            self.contig_selected = np.ones(self.n_contigs, dtype=bool)
        else:
            self.contig_selected = np.zeros(self.n_contigs, dtype=bool)
            self.contig_selected[np.fromiter(selected_ids, dtype=np.int64)] = True
        if not self.contig_selected.any():
            raise ValueError(
                "The selection holds no contigs; nothing to compute P(s) from."
            )

        if by == "genome":
            self.contig_group = self.contig_genome
            self.group_names = mc.genomes["name"].astype(str).to_numpy()
        else:
            self.contig_group = np.zeros(self.n_contigs, dtype=np.int64)
            self.group_names = np.array(["all"])
        self.n_groups = len(self.group_names)


def contact_probability(
    uri: str | Path,
    *,
    pairs: str | Path | None = None,
    by: Literal["all", "genome"] = "all",
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
    filters: Sequence[str] | None = None,
    bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
    min_distance: int | None = None,
    max_distance: int | None = None,
    chunksize: int = DEFAULT_CHUNKSIZE,
    pairs_kwargs: dict | None = None,
) -> pd.DataFrame:
    """Compute P(s) from the intra-contig contacts of a metacool file.

    Args:
        uri: ``file.metacool`` path or ``file.mmetacool::entries/name`` URI.
        pairs: Optional pairs file to measure from instead of the pixel table
            (plain text, ``.gz``, ``.bgz``) — read through the very same
            streamer ``create pairs`` uses, so compression and column handling
            are identical.  Gives base-pair distance resolution instead of bin
            resolution.
        by: ``"all"`` pools every selected contig into one curve;
            ``"genome"`` emits one curve per genome (MAG).
        genomes: Restrict to these genomes' contigs.
        contigs: Restrict to these contigs.
        filters: ``--filter`` expressions narrowing the selection further
            (see :mod:`metacooler.filters`).
        bins_per_decade: Log-grid resolution.
        min_distance: Upper edge of the leading ``[0, min)`` bin, in bp.
            Defaults to 1 bp with *pairs*, and to the file's bin size (or the
            median bin width on a variable-width file) without, since nothing
            finer is resolvable there.
        max_distance: Largest separation to report, in bp.  Defaults to the
            longest separation the selected contigs can offer.
        chunksize: Pixels (or pairs) read per iteration.
        pairs_kwargs: Extra keyword arguments for the pairs reader (column
            indices), as accepted by ``create pairs``.

    Returns:
        A long-format DataFrame with the :data:`PS_COLUMNS` columns, one row
        per (genome, distance bin), sorted by genome then distance.

    Raises:
        ValueError: If *by* is not a known grouping, or the selection is empty.
    """
    mc = MetaCooler(uri)
    selection = _Selection(mc, by, genomes, contigs, filters)
    if pairs is not None:
        return _from_pairs_file(
            selection,
            Path(pairs),
            bins_per_decade=bins_per_decade,
            min_distance=min_distance,
            max_distance=max_distance,
            chunksize=chunksize,
            pairs_kwargs=pairs_kwargs or {},
        )
    return _from_pixels(
        selection,
        bins_per_decade=bins_per_decade,
        min_distance=min_distance,
        max_distance=max_distance,
        chunksize=chunksize,
    )


def _from_pixels(
    selection: _Selection,
    *,
    bins_per_decade: int,
    min_distance: int | None,
    max_distance: int | None,
    chunksize: int,
) -> pd.DataFrame:
    """P(s) from the binned pixel table: bin-midpoint distances, bin-pair denominators."""
    mc = selection.mc
    bins = mc.bins()
    contig_of_bin = chrom_codes(bins)
    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)
    mids = (starts + ends) // 2
    selected = selection.contig_selected[contig_of_bin]

    if min_distance is None:
        binsize = mc.binsize
        min_distance = (
            int(binsize) if binsize else max(int(np.median(ends - starts)), 1)
        )
    if max_distance is None:
        max_distance = _max_intra_contig_span(contig_of_bin, mids, selected)
    edges = distance_edges(int(max_distance), int(min_distance), bins_per_decade)

    # Fold the contig id in above every coordinate so one global searchsorted
    # answers a per-contig question (see ps.bin_pair_denominators).
    fold = int(mids.max()) + int(edges[-1]) + 1 if len(mids) else 1
    key = (contig_of_bin.astype(np.int64) * fold + mids)[selected]
    group_of_bin = selection.contig_group[contig_of_bin]

    logger.info(
        "P(s) from pixels: %s selected bin(s) in %d group(s), %d distance bin(s) "
        "up to %s bp.",
        f"{int(selected.sum()):,}",
        selection.n_groups,
        len(edges) - 1,
        f"{int(edges[-1]):,}",
    )
    possible = bin_pair_denominators(
        key, group_of_bin[selected], selection.n_groups, edges
    )

    tally = PSCounts(selection.n_groups, edges)
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
        tally.tally(
            group_of_bin[b1],
            np.abs(mids[b2] - mids[b1]),
            chunk["count"].to_numpy()[keep],
        )
    return assemble(selection.group_names, edges, tally.matrix, possible)


def _from_pairs_file(
    selection: _Selection,
    pairs: Path,
    *,
    bins_per_decade: int,
    min_distance: int | None,
    max_distance: int | None,
    chunksize: int,
    pairs_kwargs: dict,
) -> pd.DataFrame:
    """P(s) from a pairs file: exact distances, position-pair denominators."""
    from metacooler.io.create import (
        _decode_pairs_chunk,
        _PairsDropStats,
        _stream_pairs,
    )

    mc = selection.mc
    lengths = selection.contig_length[selection.contig_selected]
    groups = selection.contig_group[selection.contig_selected]
    edges = pairs_distance_grid(lengths, bins_per_decade, min_distance, max_distance)

    logger.info(
        "P(s) from %s: %s selected contig(s) in %d group(s), %d distance bin(s) "
        "up to %s bp.",
        pairs,
        f"{int(selection.contig_selected.sum()):,}",
        selection.n_groups,
        len(edges) - 1,
        f"{int(edges[-1]):,}",
    )
    possible = position_pair_denominators(lengths, groups, selection.n_groups, edges)

    # Contig id per name, then the group of each contig — with the unselected
    # ones marked -1 so tally_intra_contig drops them in one comparison.
    name_to_id = {str(n): i for i, n in enumerate(mc.contigs["name"])}
    contig_group = np.where(selection.contig_selected, selection.contig_group, -1)

    tally = PSCounts(selection.n_groups, edges)
    stats = _PairsDropStats()
    contig_length = mc.contigs["length"].to_numpy(dtype=np.int64)
    for chunk in _stream_pairs(pairs, chunksize=chunksize, **pairs_kwargs):
        decoded = _decode_pairs_chunk(chunk, name_to_id, stats, contig_length)
        if decoded is None:
            continue
        c1, c2, p1, p2 = decoded
        tally_intra_contig(tally, contig_group, c1, c2, np.abs(p2 - p1))
    stats.report()
    return assemble(selection.group_names, edges, tally.matrix, possible)


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


def plot_ps(
    table: pd.DataFrame,
    out: str | Path,
    *,
    value: str = "prob",
    title: str | None = None,
    dpi: int = 300,
    force: bool = False,
) -> None:
    """Render a P(s) table as a log-log curve, one line per genome.

    Args:
        table: Output of :func:`contact_probability`.
        out: Image path; the extension picks the format (``.png``, ``.pdf``,
            ``.svg``).
        value: Column to plot on the y axis — ``"prob"`` (unit-area normalised)
            or ``"avg_contacts"`` (raw rate per opportunity).
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

    groups = list(dict.fromkeys(table["genome"]))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    handles: list = []
    labels: list[str] = []
    for name in groups:
        curve = table[table["genome"] == name]
        # A log axis cannot show a zero, and an empty distance bin is a gap in
        # the curve rather than a probability of zero.
        drawable = curve[(curve[value] > 0) & (curve["dist"] > 0)]
        if drawable.empty:
            continue
        (line,) = ax.plot(
            drawable["dist"], drawable[value], label=str(name), linewidth=1.2
        )
        handles.append(line)
        labels.append(str(name))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genomic distance (bp)")
    ax.set_ylabel("P(s)" if value == "prob" else "Mean contacts per opportunity")
    ax.set_title(title or "Contact probability vs genomic distance")
    if 1 < len(handles) <= _MAX_LEGEND_ENTRIES:
        # Handles are passed explicitly: matplotlib hides any label starting
        # with an underscore from an automatic legend, and ``__unbinned__`` —
        # the genome every unassigned contig lands in — is exactly that.
        ax.legend(handles, labels, fontsize="small", frameon=False)
    elif len(handles) > _MAX_LEGEND_ENTRIES:
        logger.info(
            "%d curves drawn without a legend (more than %d).",
            len(handles),
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
    pairs: str | Path | None = None,
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
            returned); the CLI writes to stdout in that case.
        plot: Image destination, or ``None`` for no figure.
        pairs: Measure from this pairs file instead of the pixel table.
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
        chunksize: Rows read per iteration.
        force: Overwrite existing outputs.

    Returns:
        The P(s) table.
    """
    table = contact_probability(
        uri,
        pairs=pairs,
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
        # Through the shared writer, so a table written here is byte-identical
        # to the one `create pairs` drops beside the file.
        write_ps_table(out, table, header=header, force=force)
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
