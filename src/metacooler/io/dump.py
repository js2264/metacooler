"""Dump helpers for metacooler tables."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.api import MetaCooler
from metacooler.io._utils import _filter_by_genomes_and_contigs
from metacooler.schema import DENSE_MATRIX_BIN_LIMIT

if TYPE_CHECKING:
    from typing import Literal


def _decode_chrom(series: pd.Series) -> pd.Series:
    """Ensure the chrom column contains str, not bytes (cooler version quirk)."""
    if hasattr(series, "cat"):
        cats = series.cat.categories
        if len(cats) > 0 and isinstance(cats[0], bytes):
            series = series.cat.rename_categories([c.decode("utf-8") for c in cats])
    elif (
        series.dtype == object and len(series) > 0 and isinstance(series.iloc[0], bytes)
    ):
        series = series.str.decode("utf-8")
    return series


def dump_table(
    uri: str | Path,
    table: Literal["bins", "chroms", "contigs", "genomes", "pixels", "indexes", "all"],
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
    range1: str | None = None,
    range2: str | None = None,
    header: bool = False,
    out: str | Path | None = None,
    prefix: str | None = None,
    pixels_format: Literal["long", "wide"] = "long",
    weight: str = "count",
) -> None:
    """Dump a metacooler table to TSV.

    Args:
        uri: Path to the metacool file.
        table: Which table to dump. Use ``"all"`` to dump all tables.
        genomes: Multiple genome filter.
        contigs: Multiple contig filter.
        range1: Row-dimension range selector (genome, contig, or
            ``contig:start-end``).  Applied only when ``table="pixels"``.
        range2: Column-dimension range selector.  Defaults to *range1* when
            omitted.  Applied only when ``table="pixels"``.
        header: Include header row.
        out: Output file path (default: stdout).
        prefix: Required when ``table="all"``. File prefix for output files.
        pixels_format: ``"long"`` (default) or ``"wide"`` (2D matrix).
        weight: Column name for wide format values (default: ``"count"``).
    """
    # Merge single/multi genome/contig filters
    genome_list = list(genomes) if genomes else []
    contig_list = list(contigs) if contigs else []
    mc = MetaCooler(uri)

    if table == "all":
        if not prefix:
            raise ValueError("--prefix is required when --table all is specified.")
        for tbl in ("genomes", "contigs", "bins", "pixels"):
            out_path = Path(f"{prefix}.{tbl}.tsv")
            dump_table(
                uri=uri,
                table=tbl,  # type: ignore[arg-type]
                genomes=genome_list or None,
                contigs=contig_list or None,
                range1=range1,
                range2=range2,
                header=header,
                out=out_path,
                pixels_format=pixels_format,
                weight=weight,
            )
        return

    # Resolve contig filter
    selected_contig_ids: set[int] | None = None
    if genome_list or contig_list:
        selected_contig_ids = _filter_by_genomes_and_contigs(
            mc,
            genomes=genome_list or None,
            contigs=contig_list or None,
        )

    if table == "bins":
        df = mc.bins()
        if selected_contig_ids is not None:
            bin_chrom_codes = df["chrom"].cat.codes.values
            mask = np.isin(bin_chrom_codes, list(selected_contig_ids))
            df = df[mask]
        df = df.copy()
        df["chrom"] = _decode_chrom(df["chrom"])
    elif table == "chroms":
        df = mc.contigs[["name", "length"]]
    elif table == "contigs":
        df = mc.contigs
        if selected_contig_ids is not None:
            df = df.loc[df.index.isin(selected_contig_ids)]
        core = ["name", "length", "genome_id"]
        extras = [c for c in df.columns if c not in core]
        df = df[core + extras]
    elif table == "genomes":
        df = mc.genomes
        core = ["name", "length", "n_contigs"]
        extras = [c for c in df.columns if c not in core]
        df = df[core + extras]
    elif table == "pixels":
        df = mc.pixels()
        if selected_contig_ids is not None:
            bins = mc.bins()
            bin_chrom_codes = bins["chrom"].cat.codes.values
            valid_bin_ids = np.where(
                np.isin(bin_chrom_codes, list(selected_contig_ids))
            )[0]
            mask = np.isin(df["bin1_id"].values, valid_bin_ids) & np.isin(
                df["bin2_id"].values, valid_bin_ids
            )
            df = df[mask]

        if range1 is not None:
            from metacooler.ops.range import parse_metacooler_range

            i0, i1 = parse_metacooler_range(mc, range1)
            j0, j1 = parse_metacooler_range(mc, range2) if range2 else (i0, i1)
            # Keep pixels intersecting the query box (symmetric-upper aware)
            in_row = df["bin1_id"].between(i0, i1 - 1)
            in_col = df["bin2_id"].between(j0, j1 - 1)
            in_row2 = df["bin1_id"].between(j0, j1 - 1)
            in_col2 = df["bin2_id"].between(i0, i1 - 1)
            df = df[(in_row & in_col) | (in_row2 & in_col2)]

        if pixels_format == "wide":
            # Validate weight column
            if weight not in df.columns:
                raise ValueError(
                    f"Weight column {weight!r} not found. Available: {list(df.columns)}"
                )
            n_bins = mc.shape[0]
            if n_bins > DENSE_MATRIX_BIN_LIMIT:
                raise ValueError(
                    f"Wide format refuses to allocate a {n_bins}\u00d7{n_bins} matrix "
                    f"({n_bins**2 * 8 / 1e9:.1f} GB). Use --pixels-format long."
                )
            # Build dense matrix (vectorised scatter)
            mat = np.zeros((n_bins, n_bins), dtype=np.float64)
            b1 = df["bin1_id"].values.astype(int)
            b2 = df["bin2_id"].values.astype(int)
            vals = df[weight].values.astype(np.float64)
            mat[b1, b2] = vals
            mat[b2, b1] = vals  # symmetrise
            df = pd.DataFrame(mat)
    elif table == "indexes":
        import h5py

        rows = []
        with h5py.File(mc._store, "r") as f:
            idx = f["indexes"]
            for key in idx:
                arr = idx[key][:]
                for i, val in enumerate(arr):
                    rows.append({"index": key, "position": i, "value": val})
        df = pd.DataFrame(rows)
    else:
        raise ValueError(f"Unknown table: {table!r}")

    if out:
        with open(out, "w") as out_stream:
            df.to_csv(out_stream, sep="\t", index=False, header=header)
    else:
        df.to_csv(sys.stdout, sep="\t", index=False, header=header)
