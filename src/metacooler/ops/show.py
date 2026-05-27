"""Visualization operations for metacool files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from metacooler.api import MetaCooler
from metacooler.schema import DENSE_MATRIX_BIN_LIMIT

if TYPE_CHECKING:
    from typing import Literal


def show(
    uri: str | Path,
    level: Literal["bin", "contig", "genome"] = "genome",
    genome: str | None = None,
    genomes: tuple[str, str] | None = None,
    range1: str | None = None,
    range2: str | None = None,
    transform: Literal["log1p", "none"] = "log1p",
    cmap: str = "YlOrRd",
    out: str | Path | None = None,
    dpi: int = 500,
) -> None:
    """Render a contact map heatmap.

    Args:
        uri: Path to a ``.metacool`` file.
        level: Resolution level to display.
        genome: Show only intra-genome contacts for this genome (bin level).
        genomes: Show inter-genome contacts between a pair of genomes (bin level).
        range1: Row-dimension range: genome name, contig name, or
            ``contig:start-end``.  Implies ``level="bin"`` when set.
        range2: Column-dimension range (defaults to *range1* when omitted).
            Implies ``level="bin"`` when set.
        transform: Value transform.
        cmap: Matplotlib colormap name.
        out: Output file path (display interactively if ``None``).
        dpi: DPI for the output image.
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

    if out:
        matplotlib.use("Agg")

    col_labels: list | None = None

    if level == "genome":
        from metacooler.ops.aggregate import _aggregate_sparse

        n = len(mc.genomes)
        if n > DENSE_MATRIX_BIN_LIMIT:
            raise ValueError(
                f"Refusing to allocate a {n}\u00d7{n} dense matrix "
                f"({n**2 * 8 / 1e9:.1f} GB). Use a smaller file."
            )
        mat = _aggregate_sparse(mc, level="genome").toarray()
        labels = mc.genomes["name"].tolist()
    elif level == "contig":
        from metacooler.ops.aggregate import _aggregate_sparse

        n = len(mc.contigs)
        if n > DENSE_MATRIX_BIN_LIMIT:
            raise ValueError(
                f"Refusing to allocate a {n}\u00d7{n} dense matrix "
                f"({n**2 * 8 / 1e9:.1f} GB). Use --level genome or a smaller file."
            )
        mat = _aggregate_sparse(mc, level="contig").toarray()
        labels = mc.contigs["name"].tolist()
    elif level == "bin":
        if range1 is not None:
            from metacooler.ops.range import parse_metacooler_range

            i0, i1 = parse_metacooler_range(mc, range1)
            j0, j1 = parse_metacooler_range(mc, range2) if range2 else (i0, i1)
            mat = mc.matrix(balance=False)[i0:i1, j0:j1]
            row_bins = mc.bins().iloc[i0:i1]["chrom"].astype(str).tolist()
            col_bins = mc.bins().iloc[j0:j1]["chrom"].astype(str).tolist()
            labels = row_bins if len(row_bins) <= 50 else None
            col_labels = col_bins if len(col_bins) <= 50 else None
        elif genome is None and genomes is None:
            raise ValueError(
                "At bin level, --genome, --genomes, or --range is required to avoid "
                "rendering an impossibly large matrix."
            )
        elif genome is not None:
            start, end = mc.genome_extent(genome)
            mat = mc.matrix(balance=False)[start:end, start:end]
            labels = mc.bins().iloc[start:end]["chrom"].astype(str).tolist()
        else:
            assert genomes is not None
            g1, g2 = genomes
            s1, e1 = mc.genome_extent(g1)
            s2, e2 = mc.genome_extent(g2)
            mat = mc.matrix(balance=False)[s1:e1, s2:e2]
            labels = None
    else:
        raise ValueError(f"Unknown level: {level!r}")

    mat_array: np.ndarray = np.asarray(mat)
    if transform == "log1p":
        mat_array = np.log1p(mat_array)
    elif transform != "none":
        raise ValueError(f"Unknown transform: {transform!r}")

    fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
    im = ax.imshow(
        mat_array,
        cmap=cmap,
        aspect="auto",
        interpolation="none",
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
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()

    if out:
        plt.savefig(out, dpi=dpi)
    else:
        plt.show()
    plt.close(fig)
