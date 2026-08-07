"""P(s) — the shared engine: distance grid, accumulator, denominators, table.

Pure numpy/pandas apart from the one writer at the end, which owns how the
table is named and formatted so that two producers cannot emit two dialects of
the same file.  Both producers go through this module, so a curve means the
same thing however it was measured:

- :func:`metacooler.io.create.from_pairs` accumulates while it is already
  streaming the pairs, at **base-pair** resolution, and writes the table beside
  the ``.metacool`` as ``<name>.ps.tsv``;
- :mod:`metacooler.ops.ps` recomputes from a file's pixel table, at **bin**
  resolution, or re-reads a pairs file to get the fine one back.

Only contacts whose two ends sit on the same contig take part — those are the
ones a genomic distance is defined for.

The two resolutions differ in what the denominator counts, and that is the
whole reason to prefer pairs:

``pixels``
    Pairs of *bins* at that separation.  A contig of *n* bins offers ``n - d``
    of them at separation *d*.  Nothing below the bin size is resolvable: every
    sub-bin contact lands on the diagonal.
``pairs``
    Pairs of *positions* at that separation.  A contig of *L* bp offers
    ``L - s`` of them at separation *s*, so the curve runs from 1 bp up, which
    is where the short-range artifacts a metaHi-C library is judged on live
    (undigested fragments, self-religation).

Both are "how many opportunities existed", so ``avg_contacts`` is a contact
rate per opportunity in both cases, and ``prob`` normalises it to unit area so
that libraries of different depth — and curves measured either way — are
comparable.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("metacooler")

# Distance bins per decade when the caller does not say.  Ten is the usual
# resolution for a log-log P(s): fine enough to see the shoulder where
# short-range polymer contacts give way to genuine 3D signal, coarse enough
# that a shallow library still has counts in every bin.
DEFAULT_BINS_PER_DECADE = 10

# Column layout of the P(s) table, in order.  ``genome`` holds ``"all"`` when
# the curves were pooled rather than split per MAG.
PS_COLUMNS = (
    "genome",
    "dist",
    "n_contacts",
    "n_bin_pairs",
    "avg_contacts",
    "prob",
)


def distance_edges(
    max_distance: int,
    min_distance: int,
    bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
) -> np.ndarray:
    """Build log-spaced distance-bin edges, with a leading ``[0, min)`` bin.

    The leading bin is what holds separation 0 — the diagonal of a binned map,
    or two reads at the same position — which no log grid can express and which
    is the largest single value on most curves.

    Args:
        max_distance: Largest separation to cover, in bp.
        min_distance: Upper edge of the leading bin, i.e. where the log grid
            starts.  1 bp on the pairs path; the bin size on the pixel path.
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
    # np.round can collapse two neighbouring log edges onto the same integer at
    # short distances (below ~1/bins_per_decade of a bp); np.unique drops the
    # duplicates rather than leaving a zero-width bin.
    return np.unique(np.concatenate([[0.0], np.round(grid)])).astype(np.int64)


class PSCounts:
    """Running ``(group, distance bin)`` tally of intra-contig contacts.

    Accumulated in ``float64``, which is exact for integer counts below 2**53 —
    six orders of magnitude above the contact count of any real library — and
    also carries the fractional weights a balanced pixel column would bring.

    Args:
        n_groups: Number of output rows per distance bin (1 when pooled, else
            one per genome).
        edges: Distance-bin edges from :func:`distance_edges`.
    """

    def __init__(self, n_groups: int, edges: np.ndarray) -> None:
        self.n_groups = int(n_groups)
        self.edges = np.asarray(edges, dtype=np.int64)
        self.n_dbins = len(self.edges) - 1
        self.counts = np.zeros(self.n_groups * self.n_dbins, dtype=np.float64)

    def tally(
        self,
        groups: np.ndarray,
        distances: np.ndarray,
        weights: np.ndarray | None = None,
    ) -> None:
        """Add a chunk of intra-contig contacts to the tally.

        Args:
            groups: Group id of each contact — a genome id, or 0 when pooled.
            distances: Separation of each contact, in bp.
            weights: Contact count per record (the pixel ``count`` column), or
                ``None`` for one contact per record, which is what the pairs
                path hands in.
        """
        if len(distances) == 0:
            return
        dbin = np.searchsorted(self.edges, distances, side="right") - 1
        # A separation at or past the top edge only happens when the caller
        # capped max_distance; those contacts are outside the reported range.
        inside = (dbin >= 0) & (dbin < self.n_dbins)
        if not inside.any():
            return
        flat = groups[inside] * self.n_dbins + dbin[inside]
        size = self.n_groups * self.n_dbins
        self.counts += np.bincount(
            flat,
            weights=None if weights is None else weights[inside].astype(np.float64),
            minlength=size,
        )[:size]

    @property
    def matrix(self) -> np.ndarray:
        """The tally as a ``(n_groups, n_dbins)`` array."""
        return self.counts.reshape(self.n_groups, self.n_dbins)


def pairs_distance_grid(
    lengths: np.ndarray,
    bins_per_decade: int = DEFAULT_BINS_PER_DECADE,
    min_distance: int | None = None,
    max_distance: int | None = None,
) -> np.ndarray:
    """Distance-bin edges for a curve measured from raw pairs.

    Defaults differ from the pixel path in both directions: the grid starts at
    1 bp, because read positions resolve that far and the short-range head is
    the part a library is judged on, and it stops at the longest contig, which
    is the widest separation any intra-contig pair can have.

    Args:
        lengths: Lengths of the contigs taking part, in bp.
        bins_per_decade: Log-grid resolution.
        min_distance: Override for the leading bin's upper edge.
        max_distance: Override for the largest separation reported.

    Returns:
        Edges as returned by :func:`distance_edges`.
    """
    if max_distance is None:
        max_distance = int(np.max(lengths)) if len(lengths) else 0
    return distance_edges(
        int(max_distance),
        1 if min_distance is None else int(min_distance),
        bins_per_decade,
    )


def tally_intra_contig(
    tally: PSCounts,
    contig_group: np.ndarray,
    contig1: np.ndarray,
    contig2: np.ndarray,
    distance: np.ndarray,
) -> None:
    """Add the intra-contig contacts of one decoded pairs chunk to *tally*.

    The single classification step for every pairs-derived curve, shared by
    ``create pairs`` (which is already holding these arrays for its own
    binning) and ``ps --pairs``.

    Args:
        tally: Accumulator to update.
        contig_group: Output group of every contig, with ``-1`` marking a
            contig outside the selection.  ``create`` passes genome ids, where
            nothing is excluded and the check costs one comparison.
        contig1: Contig id of the first end of every pair.
        contig2: Contig id of the second end.
        distance: ``|pos2 - pos1|`` of every pair, read only where the two ends
            share a contig.
    """
    same_contig = contig1 == contig2
    if not same_contig.any():
        return
    groups = contig_group[contig1[same_contig]]
    keep = groups >= 0
    if keep.all():
        tally.tally(groups, distance[same_contig])
        return
    tally.tally(groups[keep], distance[same_contig][keep])


def bin_pair_denominators(
    key: np.ndarray,
    groups: np.ndarray,
    n_groups: int,
    edges: np.ndarray,
) -> np.ndarray:
    """Count the pairs of *bins* that could contact, per group and distance bin.

    Uses the same folded key trick as the rest of the codebase: with the contig
    id multiplied above every coordinate, ``(contig, midpoint)`` ordering is
    plain integer ordering, so one global ``searchsorted`` per distance edge
    answers "how many later bins of *my own* contig are within this distance"
    for every bin at once — no Python loop over contigs, and no per-contig
    O(n²) enumeration of pairs.

    Only the selected bins are passed in.  A selection is always whole contigs,
    so no partner of a selected bin is ever missing from the array.

    Args:
        key: ``contig_id * K + midpoint`` of every selected bin, ascending.
        groups: Output group of every selected bin.
        n_groups: Number of output groups.
        edges: Distance-bin edges.

    Returns:
        ``(n_groups, len(edges) - 1)`` float64 array.
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
        cumulative[:, i] = np.bincount(groups, weights=partners, minlength=n_groups)[
            :n_groups
        ]

    counts = np.diff(cumulative, axis=1)
    # The diagonal: one "pair" per bin, at separation 0, which the strictly-
    # increasing-key count above cannot see (a bin is not its own partner).
    counts[:, 0] += np.bincount(groups, minlength=n_groups)[:n_groups]
    return counts


def position_pair_denominators(
    lengths: np.ndarray,
    groups: np.ndarray,
    n_groups: int,
    edges: np.ndarray,
) -> np.ndarray:
    """Count the pairs of *positions* that could contact, per group and bin.

    The base-pair analogue of :func:`bin_pair_denominators`, and the reason a
    pairs-derived curve can start at 1 bp: a contig of *L* bp offers ``L - s``
    ordered position pairs at separation *s* (and ``L`` at ``s = 0``, one per
    position), so a distance bin ``[lo, hi)`` offers

    .. code-block:: text

        sum_{s=lo}^{min(hi, L) - 1} (L - s)

    which is an arithmetic series and therefore closed-form per contig — one
    vectorised expression per distance edge over the whole contig-length array,
    rather than any iteration over separations.

    Args:
        lengths: Length in bp of every selected contig.
        groups: Output group of every selected contig.
        n_groups: Number of output groups.
        edges: Distance-bin edges.

    Returns:
        ``(n_groups, len(edges) - 1)`` float64 array.
    """
    lengths = np.asarray(lengths, dtype=np.float64)
    counts = np.zeros((n_groups, len(edges) - 1), dtype=np.float64)
    for i in range(len(edges) - 1):
        lo = float(edges[i])
        hi = float(edges[i + 1])
        upper = np.minimum(hi, lengths)  # exclusive: separations >= L do not exist
        n_terms = np.maximum(upper - lo, 0.0)
        # Sum of the run (L - lo), (L - lo - 1), … , (L - upper + 1).
        per_contig = n_terms * (2.0 * lengths - lo - upper + 1.0) / 2.0
        counts[:, i] = np.bincount(groups, weights=per_contig, minlength=n_groups)[
            :n_groups
        ]
    return counts


def assemble(
    group_names: np.ndarray | list[str],
    edges: np.ndarray,
    observed: np.ndarray,
    possible: np.ndarray,
) -> pd.DataFrame:
    """Turn the two count matrices into the :data:`PS_COLUMNS` table.

    Args:
        group_names: One name per row of *observed* / *possible*.
        edges: Distance-bin edges.
        observed: ``(n_groups, n_dbins)`` contacts.
        possible: ``(n_groups, n_dbins)`` opportunities.

    Returns:
        Long-format DataFrame, one row per (group, distance bin), sorted by
        group then distance.  Distance bins no pair could populate are omitted
        rather than reported as a contact rate of zero.
    """
    group_names = np.asarray(group_names)
    n_groups, n_dbins = observed.shape
    lo = edges[:-1]
    hi = edges[1:]
    # Geometric centre, the natural x on a log axis; the leading [0, min) bin
    # has none, so it gets the arithmetic one.
    dist = np.where(lo > 0, np.sqrt(np.maximum(lo, 1) * hi), hi / 2.0)
    width = (hi - lo).astype(np.float64)

    frame = pd.DataFrame(
        {
            "genome": np.repeat(group_names, n_dbins),
            "dist": np.tile(dist, n_groups),
            "n_contacts": observed.ravel(),
            "n_bin_pairs": possible.ravel(),
            "_width": np.tile(width, n_groups),
        }
    )
    frame = frame[frame["n_bin_pairs"] > 0].copy()
    frame["avg_contacts"] = frame["n_contacts"] / frame["n_bin_pairs"]

    # Normalise each curve to unit area so libraries of different depth are
    # comparable. A group with no contact at all keeps NaN rather than a
    # fabricated flat line.
    area = (frame["avg_contacts"] * frame["_width"]).groupby(frame["genome"]).sum()
    frame["prob"] = frame["avg_contacts"] / frame["genome"].map(
        area.replace(0.0, np.nan)
    )

    frame["n_contacts"] = frame["n_contacts"].round().astype(np.int64)
    frame["n_bin_pairs"] = frame["n_bin_pairs"].round().astype(np.int64)
    return frame[list(PS_COLUMNS)].reset_index(drop=True)


def ps_table_path(uri: str | Path) -> Path:
    """Where ``create pairs`` writes the P(s) table for an output file.

    Beside the map and named after it — ``gut_T0.metacool`` yields
    ``gut_T0.ps.tsv`` — so the curve travels with the file it describes without
    being locked inside the HDF5, where a plotting script would need metacooler
    to read it back.

    Args:
        uri: Output ``.metacool`` path.

    Returns:
        The companion table's path.
    """
    return Path(uri).with_suffix(".ps.tsv")


def write_ps_table(
    path: str | Path,
    table: pd.DataFrame,
    *,
    header: bool = True,
    force: bool = False,
) -> None:
    """Write a P(s) table as TSV, in the one format every producer emits.

    Args:
        path: Destination.
        table: Frame with the :data:`PS_COLUMNS` columns.
        header: Write the column-name row.
        force: Overwrite an existing file.

    Raises:
        FileExistsError: If *path* exists and *force* is False.
    """
    from metacooler.io._utils import check_output_path

    check_output_path(path, force)
    table.to_csv(path, sep="\t", index=False, header=header)
    logger.info("Wrote P(s) table (%d rows) to %s", len(table), path)
