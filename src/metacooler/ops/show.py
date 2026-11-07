"""Visualization operations for metacool files."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from metacooler.api import MetaCooler
from metacooler.io._utils import check_output_path, chrom_codes
from metacooler.schema import DENSE_MATRIX_BIN_LIMIT

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

logger = logging.getLogger("metacooler")

# One tick per row/column is only legible up to this many rows/columns; beyond
# it the labels overlap into a smear, so the axis falls back to one tick per
# selected genome, or to none at all.
_MAX_CATEGORICAL_TICKS = 50

# Base-pair position ticks placed on a bin-level axis. Bin-level axes routinely
# hold thousands of bins, so positions are sampled at evenly spaced bins rather
# than labelled one per bin.
_MAX_POSITION_TICKS = 6

# Names spelled out in a figure title before it switches to "(+N more)".
_MAX_TITLE_NAMES = 4

# Names listed in an "unknown name" error before it is truncated.
_MAX_ERROR_NAMES = 20

# Boundary lines: contig (or genome-row) boundaries are hairlines; the blocks of
# a multi-name selection get a heavier line, since that is the structure the
# reader is there to see.
_SEP_COLOR = "#666666"
_SEP_WIDTH = 0.5
_GROUP_SEP_COLOR = "#222222"
_GROUP_SEP_WIDTH = 1.0

_LEVEL_NOUN = {"genome": "genomes", "contig": "contigs", "bin": "bins"}


class _Selection(NamedTuple):
    """What the user asked to put on both axes.

    Attributes:
        kind: ``"all"`` (the whole file), ``"genomes"``, or ``"contigs"``.
        names: The selected names, in the order given — which is the order the
            blocks appear in, so the caller controls the layout.
    """

    kind: Literal["all", "genomes", "contigs"]
    names: list[str]


class _Axis(NamedTuple):
    """One axis (row or column dimension) of the rendered heatmap.

    The tick placement is decided by whichever builder created the axis — it is
    the only place that knows whether the rows are genomes, contigs, or bins of
    a single contig — so rendering just applies what is here.

    Attributes:
        scope: What the axis is restricted to — a single name, ``"N genomes"``,
            a range string, or ``None`` for "everything in the file".  Used in
            the figure title.
        n: Number of rows (or columns) along this axis.
        label: The X/Y axis label text.
        tick_pos: Tick positions, in row/column coordinates.  Empty for an
            unlabelled axis (too many rows/columns to label legibly).
        tick_labels: One label per entry of *tick_pos*.
        seps: Positions of the contig (or genome-row) boundary lines.
        group_seps: Positions of the selection-block boundary lines, drawn
            heavier.
    """

    scope: str | None
    n: int
    label: str
    tick_pos: list[int]
    tick_labels: list[str]
    seps: list[float]
    group_seps: list[float]


class _Panel(NamedTuple):
    """The heatmap to draw: a matrix, the two axes describing it, and a title."""

    mat: np.ndarray
    row: _Axis
    col: _Axis
    title: str


def _check_dense_size(n_rows: int, n_cols: int, hint: str) -> None:
    """Refuse to materialise a dense matrix larger than the shared limit.

    The limit is expressed as an element count (``DENSE_MATRIX_BIN_LIMIT``
    squared) so that rectangular selections (e.g. a ``-r``/``-r2`` block) are
    checked on what is actually about to be allocated.

    Args:
        n_rows: Number of rows of the dense matrix to be allocated.
        n_cols: Number of columns of the dense matrix to be allocated.
        hint: Suggestion appended to the error message.

    Raises:
        ValueError: If ``n_rows * n_cols`` exceeds the limit.
    """
    if n_rows * n_cols > DENSE_MATRIX_BIN_LIMIT**2:
        gb = n_rows * n_cols * 8 / 1e9
        raise ValueError(
            f"Refusing to allocate a {n_rows}×{n_cols} dense matrix "
            f"({gb:.1f} GB). {hint}"
        )


def _genome_id(mc: MetaCooler, name: str) -> int:
    """Return the integer genome id of *name*."""
    return int(mc.genomes.loc[mc.genomes["name"] == name].index[0])


def _contig_id(mc: MetaCooler, name: str) -> int:
    """Return the integer contig id of *name*."""
    return int(mc.contigs.loc[mc.contigs["name"] == name].index[0])


def _genome_contig_ids(mc: MetaCooler, name: str) -> list[int]:
    """Return the contig ids belonging to genome *name*."""
    gid = _genome_id(mc, name)
    return [int(i) for i in mc.contigs.index[mc.contigs["genome_id"] == gid]]


def _plural(n: int, noun: str) -> str:
    """Return ``"<n> <noun>"``, de-pluralising *noun* when ``n == 1``."""
    if n == 1 and noun.endswith("s"):
        noun = noun[:-1]
    return f"{n} {noun}"


def _format_bp(pos: int) -> str:
    """Format a base-pair coordinate compactly (``0 bp``, ``120 kb``, ``1.5 Mb``)."""
    for scale, unit in ((1_000_000, "Mb"), (1_000, "kb")):
        if abs(pos) >= scale:
            return f"{pos / scale:.2f}".rstrip("0").rstrip(".") + f" {unit}"
    return f"{pos} bp"


def _join_names(names: Sequence[str], limit: int = _MAX_TITLE_NAMES) -> str:
    """Join *names* for a title, truncating past *limit* with a ``(+N more)`` tail."""
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + f" (+{len(names) - limit} more)"


def _scoped(scope: str | None, detail: str) -> str:
    """Prefix an axis-label *detail* with the selection it describes."""
    if scope is None:
        return f"{detail[0].upper()}{detail[1:]}"
    return f"{scope} — {detail}"


def _dedupe_seps(seps: list[float], group_seps: Sequence[float]) -> list[float]:
    """Drop hairlines that a heavier selection-block line already draws.

    At contig level every block boundary *is* a genome boundary, so the two
    lists coincide entirely; at bin level they overlap only where a block
    happens to end at a contig boundary — which is always.
    """
    heavy = set(group_seps)
    return [pos for pos in seps if pos not in heavy]


def _group_bounds(sizes: Sequence[int]) -> tuple[list[int], list[float]]:
    """Return the start position and the boundary line of each block.

    Args:
        sizes: Number of rows/columns contributed by each block, in order.

    Returns:
        ``(starts, seps)`` — the first row/column of every block, and the
        positions of the lines *between* blocks (so one fewer than *starts*).
    """
    starts: list[int] = []
    seps: list[float] = []
    offset = 0
    for size in sizes:
        if starts:
            seps.append(offset - 0.5)
        starts.append(offset)
        offset += size
    return starts, seps


def _apply_ticks(ax, axis: _Axis, which: Literal["x", "y"]) -> None:
    """Apply the axis label and the ticks *axis* chose to one axis of *ax*."""
    if which == "x":
        ax.set_xlabel(axis.label, fontsize=8)
        set_ticks, set_labels, rotation = ax.set_xticks, ax.set_xticklabels, 90
    else:
        ax.set_ylabel(axis.label, fontsize=8)
        set_ticks, set_labels, rotation = ax.set_yticks, ax.set_yticklabels, 0
    if not axis.tick_pos:
        return
    set_ticks(axis.tick_pos)
    set_labels(axis.tick_labels, rotation=rotation, fontsize=6)


# ------------------------------------------------------------------
# Selection
# ------------------------------------------------------------------


def _truncate_names(names: Sequence[str], limit: int = _MAX_ERROR_NAMES) -> str:
    """Render a name list for an error message, truncating a long one."""
    if len(names) <= limit:
        return str(list(names))
    return f"{list(names[:limit])} ... ({len(names) - limit} more)"


def _validate_names(names: Sequence[str], known: Sequence[str], noun: str) -> list[str]:
    """Validate a ``--genomes``/``--contigs`` list against the file.

    Args:
        names: The names the user gave, in order.
        known: Every name of that kind present in the file.
        noun: ``"genomes"`` or ``"contigs"``, used in the messages.

    Returns:
        *names* as a list, unchanged.

    Raises:
        ValueError: If the list is empty, repeats a name, or names something
            that is not in the file.
    """
    singular = noun[:-1]
    if not names:
        raise ValueError(
            f"--{noun} was given but is empty; omit it to show every genome."
        )
    dupes = sorted(n for n, k in Counter(names).items() if k > 1)
    if dupes:
        raise ValueError(
            f"--{noun} lists {', '.join(repr(d) for d in dupes)} more than once; "
            f"each {singular} may appear only once."
        )
    known_set = set(known)
    missing = [n for n in names if n not in known_set]
    if missing:
        named = (
            f"{singular}: {missing[0]!r}" if len(missing) == 1 else f"{noun}: {missing}"
        )
        raise ValueError(
            f"Unknown {named}. Available {noun}: {_truncate_names(sorted(known_set))}"
        )
    return list(names)


def _resolve_selection(
    mc: MetaCooler,
    genomes: Sequence[str] | None,
    contigs: Sequence[str] | None,
    filters: Sequence[str] | None = None,
) -> _Selection:
    """Validate the axis selection and return it in normalised form.

    Args:
        mc: Open MetaCooler instance.
        genomes: The ``--genomes`` values, if given.
        contigs: The ``--contigs`` values, if given.
        filters: The ``--filter`` expressions, if given.  They narrow whatever
            the names selected — or the whole file, when no name is given.

    Raises:
        ValueError: If both name options are given, if either fails
            validation, or if nothing survives the filters.
    """
    if genomes is not None and contigs is not None:
        raise ValueError(
            "--genomes and --contigs are mutually exclusive; both put their "
            "names on *both* axes, so combining them has no single meaning. "
            "Use --range/--range2 for a different selection per axis."
        )
    if genomes is not None:
        base = _Selection(
            "genomes",
            _validate_names(
                [str(g) for g in genomes], mc.genomes["name"].tolist(), "genomes"
            ),
        )
    elif contigs is not None:
        base = _Selection(
            "contigs",
            _validate_names(
                [str(c) for c in contigs], mc.contigs["name"].tolist(), "contigs"
            ),
        )
    else:
        base = _Selection("all", [])
    return base if not filters else _apply_filters(mc, base, filters)


def _apply_filters(
    mc: MetaCooler, base: _Selection, filters: Sequence[str]
) -> _Selection:
    """Narrow *base* by the ``--filter`` predicates, keeping the given order.

    A contig-level filter has no meaning in genome id space (a genome is one
    row, and half of it cannot be dropped), so any contig filter resolves the
    selection to ``"contigs"`` — the same upgrade ``--contigs`` already
    triggers.  A genome-only filter leaves the kind alone.

    Args:
        mc: Open MetaCooler instance.
        base: The selection the names produced.
        filters: Raw ``--filter`` expressions.

    Returns:
        The narrowed selection, with names in their original order.

    Raises:
        metacooler.filters.FilterError: If a filter is malformed, names an
            unknown column, or leaves nothing selected.
    """
    from metacooler.filters import (
        FilterError,
        explain_empty,
        filter_mask,
        parse_filters,
    )
    from metacooler.io._utils import contigs_with_genome_names

    parsed = parse_filters(filters)
    genomes_df = mc.genomes
    contigs_df = contigs_with_genome_names(mc)
    kept_genomes = set(
        genomes_df.loc[filter_mask(genomes_df, parsed, table="genomes"), "name"]
    )
    contig_mask = filter_mask(contigs_df, parsed, table="contigs")
    kept_contig_ids = set(contigs_df.index[contig_mask])

    if base.kind == "contigs" or any(f.table == "contigs" for f in parsed):
        # Candidate contigs, in the order the caller asked for.
        if base.kind == "contigs":
            candidates = [_contig_id(mc, c) for c in base.names]
        elif base.kind == "genomes":
            candidates = [cid for g in base.names for cid in _genome_contig_ids(mc, g)]
        else:
            candidates = list(contigs_df.index)
        genome_of = contigs_df["genome_id"].map(genomes_df["name"])
        names = [
            str(contigs_df.at[cid, "name"])
            for cid in candidates
            if cid in kept_contig_ids and genome_of[cid] in kept_genomes
        ]
        kind: Literal["all", "genomes", "contigs"] = "contigs"
        n_candidates = len(candidates)
    else:
        candidates_g = (
            base.names if base.kind == "genomes" else genomes_df["name"].tolist()
        )
        names = [str(g) for g in candidates_g if g in kept_genomes]
        kind = "genomes"
        n_candidates = len(candidates_g)

    logger.info("--filter kept %d of %d %s.", len(names), n_candidates, kind)
    if not names:
        raise FilterError(
            explain_empty(
                {"genomes": genomes_df, "contigs": contigs_df},
                parsed,
                context=(
                    f"...and nothing named by --{base.kind} survived them."
                    if base.kind != "all"
                    else ""
                ),
            )
        )
    return _Selection(kind, names)


def _selection_scope(sel: _Selection) -> str | None:
    """The axis scope of *sel*: the single name, ``"N genomes"``, or ``None``."""
    if sel.kind == "all":
        return None
    if len(sel.names) == 1:
        return sel.names[0]
    return _plural(len(sel.names), sel.kind)


def _selection_title(sel: _Selection, level: str) -> str:
    """The figure title for a genome/contig selection (or the whole file)."""
    if sel.kind == "all":
        head = "All genomes"
    elif len(sel.names) == 1:
        head = sel.names[0]
    else:
        head = f"{_plural(len(sel.names), sel.kind)}: {_join_names(sel.names)}"
    return f"{head} — {level} level"


def _range_title(row: _Axis, col: _Axis, level: str) -> str:
    """The figure title for a ``--range``/``--range2`` selection."""
    head = str(row.scope) if row.scope == col.scope else f"{row.scope} vs {col.scope}"
    return f"{head} — {level} level"


def _selection_blocks(
    mc: MetaCooler, level: Literal["contig", "genome"], sel: _Selection
) -> list[tuple[str, list[int]]]:
    """Return ``(block name, ids)`` for every name in *sel*, in the given order."""
    if sel.kind == "genomes":
        if level == "genome":
            return [(g, [_genome_id(mc, g)]) for g in sel.names]
        return [(g, _genome_contig_ids(mc, g)) for g in sel.names]
    # A contig selection is meaningless in genome id space; ``show`` coerces
    # ``level="genome"`` to ``"contig"`` before we get here.
    return [(c, [_contig_id(mc, c)]) for c in sel.names]


def _selection_extents(mc: MetaCooler, sel: _Selection) -> list[tuple[int, int]]:
    """Return the half-open bin range of every name in *sel*, in the given order."""
    if sel.kind == "genomes":
        return [mc.genome_extent(g) for g in sel.names]
    return [mc.contig_extent(c) for c in sel.names]


# ------------------------------------------------------------------
# Axes
# ------------------------------------------------------------------


def _level_axis(
    mc: MetaCooler,
    agg_level: Literal["contig", "genome"],
    ids: list[int],
    scope: str | None,
    group_spans: Sequence[tuple[str, int]] | None = None,
) -> _Axis:
    """Describe an axis of an aggregated (genome- or contig-level) heatmap.

    Each row/column is one genome or contig, so each gets its own tick — until
    there are too many for that to be legible, at which point a multi-name
    selection is ticked once per block instead of going unlabelled.

    Args:
        mc: Open MetaCooler instance.
        agg_level: ``"genome"`` or ``"contig"`` — the id space of *ids*.
        ids: The genome or contig ids on this axis, in display order.
        scope: What the axis is restricted to, for the label and title.
        group_spans: ``(block name, number of ids)`` per selection block, when
            the axis is a concatenation of several.
    """
    names = mc.genomes["name"] if agg_level == "genome" else mc.contigs["name"]
    if agg_level == "genome":
        # Every genome is its own row/column, so every boundary is a separator.
        seps = [i + 0.5 for i in range(len(ids) - 1)]
    else:
        gid = mc.contigs.loc[ids, "genome_id"].to_numpy()
        seps = (np.where(np.diff(gid) != 0)[0] + 0.5).tolist()

    group_seps: list[float] = []
    if 0 < len(ids) <= _MAX_CATEGORICAL_TICKS:
        tick_pos = list(range(len(ids)))
        tick_labels = [str(x) for x in names.loc[ids]]
    else:
        tick_pos, tick_labels = [], []
    if group_spans is not None:
        starts, group_seps = _group_bounds([size for _, size in group_spans])
        if not tick_pos and len(group_spans) <= _MAX_CATEGORICAL_TICKS:
            # Too many contigs to name one by one; name the blocks instead.
            tick_pos = starts
            tick_labels = [name for name, _ in group_spans]

    count = _plural(len(ids), _LEVEL_NOUN[agg_level])
    return _Axis(
        scope=scope,
        n=len(ids),
        # ``scope == count`` for a many-genome selection at genome level, where
        # the genomes *are* the rows: "4 genomes", not "4 genomes (4 genomes)".
        label=count if scope is None or scope == count else f"{scope} ({count})",
        tick_pos=tick_pos,
        tick_labels=tick_labels,
        seps=_dedupe_seps(seps, group_seps),
        group_seps=group_seps,
    )


def _bin_axis(
    mc: MetaCooler,
    all_bins,
    extents: Sequence[tuple[int, int]],
    scope: str | None,
    group_names: Sequence[str] | None = None,
) -> _Axis:
    """Describe an axis of a bin-level heatmap covering the given bin ranges.

    Ticks are base-pair coordinates while the selection stays inside one contig
    — a bin index is not a position, and it starts at 0 even when the selection
    does not, so ``--range contig:20000-50000`` used to read as though it began
    at the start of the contig.  A selection spanning several contigs has no
    single coordinate system, so there the ticks name the contigs instead, one
    at the first bin of each (matching the separator lines) — or the selection
    blocks, when the axis is a concatenation of several.

    Args:
        mc: Open MetaCooler instance.
        all_bins: The full bins table.
        extents: Half-open ``(lo, hi)`` bin ranges, in display order.
        scope: What the axis is restricted to, for the label and title.
        group_names: One name per entry of *extents*, when the axis is a
            concatenation of several.
    """
    positions = (
        np.concatenate([np.arange(lo, hi) for lo, hi in extents])
        if extents
        else np.empty(0, dtype=int)
    )
    block = all_bins.iloc[positions]
    codes = chrom_codes(block)
    starts = block["start"].to_numpy()
    n = int(positions.size)
    contig_starts = [0, *(int(i) for i in np.flatnonzero(np.diff(codes) != 0) + 1)]
    n_contigs = len(contig_starts) if n else 0
    group_starts, group_seps = _group_bounds([hi - lo for lo, hi in extents])
    group_seps = group_seps if len(extents) > 1 else []

    if n_contigs > 1:
        contig_count = _plural(n_contigs, "contigs")
        detail = f"{contig_count} ({_plural(n, 'bins')})"
        # "3 contigs — 3 contigs (30 bins)" says it twice.
        label = detail if scope == contig_count else _scoped(scope, detail)
        if group_names is not None and len(extents) > 1:
            labelled = len(extents) <= _MAX_CATEGORICAL_TICKS
            tick_pos = group_starts if labelled else []
            tick_labels = list(group_names) if labelled else []
        else:
            labelled = n_contigs <= _MAX_CATEGORICAL_TICKS
            tick_pos = contig_starts if labelled else []
            contig_names = mc.contigs["name"]
            tick_labels = [str(contig_names.iloc[codes[i]]) for i in tick_pos]
    else:
        label = _scoped(scope, "position (bp)")
        n_ticks = min(_MAX_POSITION_TICKS, n)
        tick_pos = (
            np.unique(np.linspace(0, n - 1, n_ticks).round().astype(int)).tolist()
            if n
            else []
        )
        tick_labels = [_format_bp(int(starts[i])) for i in tick_pos]

    return _Axis(
        scope=scope,
        n=n,
        label=label,
        tick_pos=tick_pos,
        tick_labels=tick_labels,
        seps=_dedupe_seps([pos - 0.5 for pos in contig_starts[1:]], group_seps),
        group_seps=group_seps,
    )


# ------------------------------------------------------------------
# Matrix
# ------------------------------------------------------------------


def _dense_block(
    selector,
    row_extents: Sequence[tuple[int, int]],
    col_extents: Sequence[tuple[int, int]],
) -> np.ndarray:
    """Materialise the dense matrix spanned by the given bin ranges.

    Each genome's (or contig's) bins form one contiguous range, but the ranges
    of several are neither adjacent nor necessarily in the requested order, so a
    multi-name selection is not one slice.  Fetching the bounding block instead
    would pull in everything that happens to lie between the two extremes — up
    to the whole matrix — so the blocks are fetched pair by pair and placed into
    the output.

    Args:
        selector: A cooler matrix selector (``mc.matrix(balance=False)``).
        row_extents: Half-open ``(lo, hi)`` bin ranges of the row axis.
        col_extents: Half-open ``(lo, hi)`` bin ranges of the column axis.

    Returns:
        The dense matrix, rows and columns ordered as the extents were given.
    """
    if len(row_extents) == 1 and len(col_extents) == 1:
        (i0, i1), (j0, j1) = row_extents[0], col_extents[0]
        return selector[i0:i1, j0:j1]

    n_rows = sum(hi - lo for lo, hi in row_extents)
    n_cols = sum(hi - lo for lo, hi in col_extents)
    mat = np.zeros((n_rows, n_cols), dtype=np.float64)
    row_at = 0
    for i0, i1 in row_extents:
        col_at = 0
        for j0, j1 in col_extents:
            mat[row_at : row_at + (i1 - i0), col_at : col_at + (j1 - j0)] = selector[
                i0:i1, j0:j1
            ]
            col_at += j1 - j0
        row_at += i1 - i0
    return mat


def _transform_matrix(
    mat: np.ndarray, transform: Literal["log1p", "log10", "none"]
) -> np.ndarray:
    """Apply the requested value transform to a (dense) matrix."""
    arr = np.asarray(mat, dtype=float)
    if transform == "log1p":
        return np.log1p(arr)
    if transform == "log10":
        return np.where(arr > 0, np.log10(arr), 0.0)
    return arr


# ------------------------------------------------------------------
# Panels
# ------------------------------------------------------------------


def _aggregated_panel(
    mc: MetaCooler,
    level: Literal["contig", "genome"],
    sel: _Selection,
) -> _Panel:
    """Build the heatmap of a genome- or contig-level figure."""
    from metacooler.ops.aggregate import _aggregate_sparse

    hint = (
        "Use a smaller file."
        if level == "genome"
        else "Use --level genome, --genomes/--contigs, or a smaller file."
    )

    spans: list[tuple[str, int]] | None = None
    if sel.kind == "all":
        ids = list(range(len(mc.genomes) if level == "genome" else len(mc.contigs)))
    else:
        blocks = _selection_blocks(mc, level, sel)
        ids = [i for _, block in blocks for i in block]
        spans = [(name, len(block)) for name, block in blocks]
        # At contig level a --contigs selection is one row per block, so blocking
        # every row would just draw a heavy line between all of them.
        if sel.kind == "contigs" and level == "contig":
            spans = None

    _check_dense_size(len(ids), len(ids), hint)

    agg = _aggregate_sparse(mc, level=level)
    # Only "all" may take the whole-matrix fast path: a selection naming every
    # genome in a *different* order has the same length but needs the permutation.
    mat = agg.toarray() if sel.kind == "all" else agg[ids, :][:, ids].toarray()

    axis = _level_axis(mc, level, ids, _selection_scope(sel), spans)
    if level == "genome" and mat.shape == (1, 1):
        logger.warning(
            "Selection at --level genome is a single genome (1×1 image); "
            "use --level contig or --level bin for a more informative view."
        )
    return _Panel(mat=mat, row=axis, col=axis, title=_selection_title(sel, level))


def _bin_panel(
    mc: MetaCooler,
    sel: _Selection,
    range1: str | None,
    range2: str | None,
) -> _Panel:
    """Build the heatmap of a bin-level figure."""
    all_bins = mc.bins()

    if range1 is not None:
        from metacooler.ops.range import parse_metacooler_range

        # The one asymmetric selection: a different range per axis.
        row_extents = [parse_metacooler_range(mc, range1)]
        col_extents = (
            [parse_metacooler_range(mc, range2)] if range2 else list(row_extents)
        )
        _check_dense_size(
            row_extents[0][1] - row_extents[0][0],
            col_extents[0][1] - col_extents[0][0],
            "Use a narrower --range/--range2.",
        )
        mat = _dense_block(mc.matrix(balance=False), row_extents, col_extents)
        row = _bin_axis(mc, all_bins, row_extents, range1)
        col = _bin_axis(mc, all_bins, col_extents, range2 or range1)
        return _Panel(mat=mat, row=row, col=col, title=_range_title(row, col, "bin"))

    if sel.kind == "all":
        raise ValueError(
            "At bin level, --genomes, --contigs, or --range is required to avoid "
            "rendering an impossibly large matrix."
        )

    extents = _selection_extents(mc, sel)
    size = sum(hi - lo for lo, hi in extents)
    _check_dense_size(
        size, size, "Use --range to select a smaller region, or --level contig."
    )
    mat = _dense_block(mc.matrix(balance=False), extents, extents)
    axis = _bin_axis(
        mc,
        all_bins,
        extents,
        _selection_scope(sel),
        sel.names if len(extents) > 1 else None,
    )
    return _Panel(mat=mat, row=axis, col=axis, title=_selection_title(sel, "bin"))


def show(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "genome",
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
    filters: Sequence[str] | None = None,
    range1: str | None = None,
    range2: str | None = None,
    transform: Literal["log1p", "log10", "none"] = "log1p",
    cmap: str = "afmhot_r",
    vmax: float | None = None,
    separators: bool = True,
    out: str | Path | None = None,
    dpi: int = 500,
    force: bool = False,
) -> None:
    """Render a contact map heatmap.

    There are two ways to choose what the axes show, and they cannot be
    combined:

    * *genomes* / *contigs* — a **symmetric** selection: every name goes on
      **both** axes, so the figure holds all the intra- *and* inter-block
      contacts between them, laid out in the order the names were given.  This
      is uniform in the number of names: one name is that genome (or contig)
      against itself, and N names are all N against all N.
    * *range1* / *range2* — the **asymmetric** selection: a different region per
      axis, at bin level.

    Both are honoured at every *level*: at ``level="genome"`` a genome is a
    single row/column, at ``"contig"`` it is its contigs, at ``"bin"`` its bins.
    A *contigs* selection has no meaning in genome id space, so it upgrades
    ``level="genome"`` to ``"contig"``; likewise *range1* implies
    ``level="bin"``.

    Args:
        uri: Path to a ``.metacool`` file.
        level: Resolution level to display.
        genomes: Genome names to put on both axes.
        contigs: Contig names to put on both axes.
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`) that
            narrow *genomes*/*contigs* — or the whole file when neither is
            given.  A ``contigs.*`` filter resolves the axes to contigs, just
            as *contigs* does.  Cannot be combined with *range1*/*range2*.
        range1: Row-dimension range: genome name, contig name, or
            ``contig:start-end``.  Implies ``level="bin"`` when set.
        range2: Column-dimension range (defaults to *range1* when omitted).
            Implies ``level="bin"`` when set.
        transform: Value transform (``log1p``, ``log10``, or ``none``).
        cmap: Matplotlib colormap name.
        vmax: Upper bound for the colormap. Defaults to the 99th percentile of
            non-zero values in the (transformed) matrix.
        separators: Draw lines at genome (or contig) boundaries.
        out: Output file path (display interactively if ``None``).
        dpi: DPI for the output image.
        force: Overwrite *out* if it already exists.

    Raises:
        ValueError: If *level*, *transform* or a name is unknown, if
            incompatible selection arguments are combined, or if the requested
            matrix is too large to materialise densely.
        FileExistsError: If *out* exists and *force* is False.
    """
    if level not in _LEVEL_NOUN:
        raise ValueError(f"Unknown level: {level!r}")
    if transform not in ("log1p", "log10", "none"):
        raise ValueError(f"Unknown transform: {transform!r}")

    mc = MetaCooler(uri)
    if (range1 is not None or range2 is not None) and (genomes or contigs or filters):
        raise ValueError(
            "--range/--range2 cannot be combined with --genomes, --contigs or "
            "--filter; they are two ways of selecting the same axes."
        )
    sel = _resolve_selection(mc, genomes, contigs, filters)

    # A selector that only makes sense at one level upgrades the level rather
    # than being silently dropped. Done here, not in the CLI, so the Python API
    # behaves identically.
    if range1 is not None and level != "bin":
        logger.info("--range implies --level bin; ignoring --level %s.", level)
        level = "bin"
    elif sel.kind == "contigs" and level == "genome":
        logger.info(
            "%s has no meaning at --level genome (a genome is one row/column, "
            "and part of it cannot be dropped); showing --level contig instead.",
            "--contigs" if contigs else "A contigs.* --filter",
        )
        level = "contig"

    if out is not None:
        check_output_path(out, force, inputs=[uri])
        matplotlib.use("Agg")

    panel = (
        _bin_panel(mc, sel, range1, range2)
        if level == "bin"
        else _aggregated_panel(mc, level, sel)
    )

    mat = _transform_matrix(panel.mat, transform)
    if vmax is None:
        nonzero = mat[mat > 0]
        vmax_val = float(np.percentile(nonzero, 99)) if nonzero.size > 0 else 1.0
    else:
        vmax_val = float(vmax)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=dpi, layout="constrained")
    im = ax.imshow(
        mat,
        cmap=cmap,
        aspect="equal",
        interpolation="none",
        vmin=0,
        vmax=vmax_val,
    )
    ax.set_title(panel.title, fontsize=10)
    _apply_ticks(ax, panel.col, "x")
    _apply_ticks(ax, panel.row, "y")
    if separators:
        for pos in panel.row.seps:
            ax.axhline(pos, color=_SEP_COLOR, linewidth=_SEP_WIDTH)
        for pos in panel.col.seps:
            ax.axvline(pos, color=_SEP_COLOR, linewidth=_SEP_WIDTH)
        for pos in panel.row.group_seps:
            ax.axhline(pos, color=_GROUP_SEP_COLOR, linewidth=_GROUP_SEP_WIDTH)
        for pos in panel.col.group_seps:
            ax.axvline(pos, color=_GROUP_SEP_COLOR, linewidth=_GROUP_SEP_WIDTH)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    if out:
        plt.savefig(out, dpi=dpi)
    else:
        plt.show()
    plt.close(fig)
