"""Visualization operations for metacool files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from metacooler.api import MetaCooler
from metacooler.io._utils import check_output_path, chrom_codes
from metacooler.schema import DENSE_MATRIX_BIN_LIMIT

if TYPE_CHECKING:
    from typing import Literal

logger = logging.getLogger("metacooler")


def _check_dense_size(n_rows: int, n_cols: int, hint: str) -> None:
    """Refuse to materialise a dense matrix larger than the shared limit.

    The limit is expressed as an element count (``DENSE_MATRIX_BIN_LIMIT``
    squared) so that rectangular selections (e.g. an inter-genome block) are
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


def _member_ids(
    mc: MetaCooler, level: Literal["contig", "genome"], name: str
) -> list[int]:
    """Return the ids of *name* in the id space of *level*.

    At ``level="genome"`` a genome is a single row/column, so the result is a
    one-element list; at ``level="contig"`` it is the list of the genome's
    contig ids.
    """
    gid = _genome_id(mc, name)
    if level == "genome":
        return [gid]
    return [int(i) for i in mc.contigs.index[mc.contigs["genome_id"] == gid]]


def show(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "genome",
    genome: str | None = None,
    genomes: tuple[str, str] | None = None,
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

    Selection is uniform across the three levels: ``genome`` restricts both
    axes to that genome (its single row/column at ``level="genome"``, its
    contigs at ``level="contig"``, its bins at ``level="bin"``), and
    ``genomes=(g1, g2)`` puts *g1* on the row axis and *g2* on the column
    axis (the inter-genome block).  ``range1``/``range2`` are an alternative,
    bin-level way of selecting the two axes and cannot be combined with
    ``genome``/``genomes``.

    Args:
        uri: Path to a ``.metacool`` file.
        level: Resolution level to display.
        genome: Show only intra-genome contacts for this genome (all levels).
        genomes: Show inter-genome contacts between a pair of genomes, the
            first on the row axis and the second on the column axis (all
            levels).
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
        ValueError: If *level*, *transform* or a genome name is unknown, if
            incompatible selection arguments are combined, or if the requested
            matrix is too large to materialise densely.
        FileExistsError: If *out* exists and *force* is False.
    """
    mc = MetaCooler(uri)
    known_genomes = set(mc.genomes["name"].tolist())

    if genome is not None and genome not in known_genomes:
        raise ValueError(
            f"Unknown genome: {genome!r}. Available genomes: {sorted(known_genomes)}"
        )
    if genomes is not None:
        missing = [g for g in genomes if g not in known_genomes]
        if missing:
            raise ValueError(
                f"Unknown genome(s): {missing}. "
                f"Available genomes: {sorted(known_genomes)}"
            )
    if genome is not None and genomes is not None:
        raise ValueError("--genome and --genomes are mutually exclusive.")
    if (range1 is not None or range2 is not None) and (
        genome is not None or genomes is not None
    ):
        raise ValueError(
            "--range/--range2 cannot be combined with --genome or --genomes; "
            "they are two ways of selecting the same axes."
        )

    if out is not None:
        check_output_path(out, force)

    if out:
        matplotlib.use("Agg")

    labels: list | None = None
    col_labels: list | None = None
    row_sep: list[float] = []
    col_sep: list[float] = []

    if level in ("genome", "contig"):
        from metacooler.ops.aggregate import _aggregate_sparse

        agg_level: Literal["contig", "genome"] = (
            "genome" if level == "genome" else "contig"
        )
        n_level = len(mc.genomes) if agg_level == "genome" else len(mc.contigs)
        # The selection is resolved (and size-checked) before aggregating, and
        # applied in sparse space: only the requested block is ever densified,
        # so a 50-contig request is not refused because the file holds 100k
        # contigs.
        if genome is not None:
            row_ids = col_ids = _member_ids(mc, agg_level, genome)
        elif genomes is not None:
            row_ids = _member_ids(mc, agg_level, genomes[0])
            col_ids = _member_ids(mc, agg_level, genomes[1])
        else:
            row_ids = col_ids = list(range(n_level))

        hint = (
            "Use a smaller file."
            if level == "genome"
            else "Use --level genome, --genome/--genomes, or a smaller file."
        )
        _check_dense_size(len(row_ids), len(col_ids), hint)

        agg = _aggregate_sparse(mc, level=agg_level)
        if len(row_ids) == agg.shape[0] and len(col_ids) == agg.shape[1]:
            mat = agg.toarray()
        else:
            mat = agg[row_ids, :][:, col_ids].toarray()

        if agg_level == "genome":
            names = mc.genomes["name"]
            row_sep = [i + 0.5 for i in range(len(row_ids) - 1)]
            col_sep = [i + 0.5 for i in range(len(col_ids) - 1)]
            if mat.shape == (1, 1) and (genome is not None or genomes is not None):
                logger.warning(
                    "Selection at --level genome is a single genome pair "
                    "(1×1 image); use --level contig or --level bin for "
                    "a more informative view."
                )
        else:
            names = mc.contigs["name"]
            row_gid = mc.contigs.loc[row_ids, "genome_id"].values
            col_gid = mc.contigs.loc[col_ids, "genome_id"].values
            row_sep = (np.where(np.diff(row_gid) != 0)[0] + 0.5).tolist()
            col_sep = (np.where(np.diff(col_gid) != 0)[0] + 0.5).tolist()
        labels = names.loc[row_ids].tolist()
        col_labels = names.loc[col_ids].tolist()
    elif level == "bin":
        _all_bins = mc.bins()
        if range1 is not None:
            from metacooler.ops.range import parse_metacooler_range

            i0, i1 = parse_metacooler_range(mc, range1)
            j0, j1 = parse_metacooler_range(mc, range2) if range2 else (i0, i1)
            hint = "Use a narrower --range/--range2."
        elif genome is None and genomes is None:
            raise ValueError(
                "At bin level, --genome, --genomes, or --range is required to avoid "
                "rendering an impossibly large matrix."
            )
        elif genome is not None:
            i0, i1 = mc.genome_extent(genome)
            j0, j1 = i0, i1
            hint = "Use --range to select a smaller region, or --level contig."
        else:
            assert genomes is not None
            i0, i1 = mc.genome_extent(genomes[0])
            j0, j1 = mc.genome_extent(genomes[1])
            hint = "Use --range to select a smaller region, or --level contig."

        # mc.matrix()[...] materialises a dense ndarray, so guard it exactly
        # like the aggregated levels above.
        _check_dense_size(i1 - i0, j1 - j0, hint)
        mat = mc.matrix(balance=False)[i0:i1, j0:j1]

        row_chrom = chrom_codes(_all_bins.iloc[i0:i1])
        col_chrom = chrom_codes(_all_bins.iloc[j0:j1])
        row_sep = (np.where(np.diff(row_chrom) != 0)[0] + 0.5).tolist()
        col_sep = (np.where(np.diff(col_chrom) != 0)[0] + 0.5).tolist()
        if genomes is None:
            row_bins = _all_bins.iloc[i0:i1]["chrom"].astype(str).tolist()
            col_bins = _all_bins.iloc[j0:j1]["chrom"].astype(str).tolist()
            labels = row_bins if len(row_bins) <= 50 else None
            col_labels = col_bins if len(col_bins) <= 50 else None
    else:
        raise ValueError(f"Unknown level: {level!r}")

    mat_array: np.ndarray = np.asarray(mat, dtype=float)
    if transform == "log1p":
        mat_array = np.log1p(mat_array)
    elif transform == "log10":
        mat_array = np.where(mat_array > 0, np.log10(mat_array), 0.0)
    elif transform != "none":
        raise ValueError(f"Unknown transform: {transform!r}")

    if vmax is None:
        nonzero = mat_array[mat_array > 0]
        vmax_val = float(np.percentile(nonzero, 99)) if nonzero.size > 0 else 1.0
    else:
        vmax_val = float(vmax)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=dpi)
    im = ax.imshow(
        mat_array,
        cmap=cmap,
        aspect="equal",
        interpolation="none",
        vmin=0,
        vmax=vmax_val,
    )
    ax.set_xlabel(
        "Genome" if level == "genome" else "Contig" if level == "contig" else "Bin"
    )
    ax.set_ylabel(
        "Genome" if level == "genome" else "Contig" if level == "contig" else "Bin"
    )
    if labels and len(labels) <= 50:
        _col_labels = col_labels if col_labels is not None else labels
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=6)
        if _col_labels and len(_col_labels) <= 50:
            ax.set_xticks(range(len(_col_labels)))
            ax.set_xticklabels(_col_labels, rotation=90, fontsize=6)
    if separators:
        lw = 0.5
        lc = "#666666"
        for pos in row_sep:
            ax.axhline(pos, color=lc, linewidth=lw)
        for pos in col_sep:
            ax.axvline(pos, color=lc, linewidth=lw)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()

    if out:
        plt.savefig(out, dpi=dpi)
    else:
        plt.show()
    plt.close(fig)
