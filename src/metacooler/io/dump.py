"""Dump helpers for metacooler tables."""

from __future__ import annotations

import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.api import MetaCooler, _merge_ranges, _ranges_to_ids
from metacooler.io._utils import (
    check_output_path,
    chrom_codes,
    resolve_contig_selection,
)
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


def _selected_pixels(
    mc: MetaCooler,
    valid_bin_ids: np.ndarray | None,
    range_bounds: tuple[int, int, int, int] | None,
) -> Iterator[pd.DataFrame]:
    """Yield the stored pixels a ``dump -t pixels`` selection keeps, in order.

    Only the pixel rows that can match are read, located through
    ``indexes/bin1_offset``: the row blocks of the two ``--range`` dimensions,
    or of the selected contigs' bins.  A pixel ``(bin1, bin2)`` is kept when it
    links the row range to the column range in either orientation, and both its
    bins are selected.

    Args:
        mc: Open map.
        valid_bin_ids: Sorted bins of the ``--genomes``/``--contigs``/
            ``--filter`` selection, or ``None`` for every bin.
        range_bounds: ``(i0, i1, j0, j1)`` half-open bin ranges of ``--range``
            and ``--range2``, or ``None``.

    Yields:
        Non-empty pixel frames, in file order.
    """
    n_bins = mc.shape[0]
    is_valid: np.ndarray | None = None
    if valid_bin_ids is not None:
        is_valid = np.zeros(n_bins, dtype=bool)
        is_valid[valid_bin_ids] = True

    if range_bounds is not None:
        i0, i1, j0, j1 = range_bounds
        rows = _merge_ranges(np.array([i0, j0]), np.array([i1, j1]))
    elif valid_bin_ids is not None:
        rows = _merge_ranges(valid_bin_ids, valid_bin_ids + 1)
    else:
        rows = np.array([[0, n_bins]], dtype=np.int64)

    def keep(b1: np.ndarray, b2: np.ndarray) -> np.ndarray:
        mask = np.ones(len(b1), dtype=bool)
        if range_bounds is not None:
            i0, i1, j0, j1 = range_bounds
            mask = ((b1 >= i0) & (b1 < i1) & (b2 >= j0) & (b2 < j1)) | (
                (b1 >= j0) & (b1 < j1) & (b2 >= i0) & (b2 < i1)
            )
        if is_valid is not None:
            mask &= is_valid[b1] & is_valid[b2]
        return mask

    selective = range_bounds is not None or is_valid is not None
    yield from mc._iter_pixel_rows(rows, keep if selective else None)


def dump_table(
    uri: str | Path,
    table: Literal[
        "bins",
        "chroms",
        "contigs",
        "genomes",
        "pixels",
        "indexes",
        "clusters",
        "members",
        "all",
    ],
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
    filters: list[str] | None = None,
    range1: str | None = None,
    range2: str | None = None,
    header: bool = False,
    out: str | Path | None = None,
    prefix: str | None = None,
    pixels_format: Literal["long", "wide"] = "long",
    weight: str = "count",
    join: bool = False,
    force: bool = False,
) -> None:
    """Dump a metacooler table to TSV.

    Args:
        uri: Path to the metacool file.
        table: Which table to dump. Use ``"all"`` to dump all tables.
        genomes: Multiple genome filter.
        contigs: Multiple contig filter.
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`), e.g.
            ``["contigs.length >= 10kb"]``.  ANDed with each other and applied
            *within* the *genomes*/*contigs* selection rather than beside it.
        range1: Row-dimension range selector (genome, contig, or
            ``contig:start-end``).  Applied only when ``table="pixels"``.
        range2: Column-dimension range selector.  Defaults to *range1* when
            omitted.  Applied only when ``table="pixels"``.
        header: Include header row.
        out: Output file path (default: stdout).
        prefix: Required when ``table="all"``. File prefix for output files.
        pixels_format: ``"long"`` (default) or ``"wide"`` (2D matrix).
        weight: Column name for wide format values (default: ``"count"``).
        join: Annotate pixel output with contig names and coordinates.  For
            ``pixels_format="long"`` this replaces ``bin1_id``/``bin2_id`` with
            ``contig1/start1/end1``/``contig2/start2/end2`` (cooler ``--join``
            semantics); for ``"wide"`` it labels the matrix rows/columns with
            ``contig:start-end``.  Applied only when ``table="pixels"``.
        force: Overwrite existing output file(s).  With ``table="all"`` every
            ``<prefix>.<table>.tsv`` is guarded.

    Raises:
        FileExistsError: If an output file exists and *force* is False.
        ValueError: If *table* is unknown, if ``table="all"`` is used without
            *prefix*, or if the wide-format selection is too large to
            materialise densely.
    """
    # Merge single/multi genome/contig filters
    genome_list = list(genomes) if genomes else []
    contig_list = list(contigs) if contigs else []
    filter_list = list(filters) if filters else []
    mc = MetaCooler(uri)

    if table == "all":
        if not prefix:
            raise ValueError("--prefix is required when --table all is specified.")
        # ``chroms`` is intentionally omitted: it is the cooler soft-link alias
        # of ``contigs[["name", "length"]]`` and would duplicate the contigs file.
        tables = ("genomes", "contigs", "bins", "pixels", "indexes")
        if mc.cluster_rounds:
            tables += ("clusters", "members")
        # Guard every destination up front so a clash on the last table does
        # not leave the first four already overwritten.
        for tbl in tables:
            check_output_path(Path(f"{prefix}.{tbl}.tsv"), force, inputs=[uri])
        for tbl in tables:
            out_path = Path(f"{prefix}.{tbl}.tsv")
            dump_table(
                uri=uri,
                table=tbl,  # type: ignore[arg-type]
                genomes=genome_list or None,
                contigs=contig_list or None,
                filters=filter_list or None,
                range1=range1,
                range2=range2,
                header=header,
                out=out_path,
                pixels_format=pixels_format,
                weight=weight,
                join=join,
                force=True,
            )
        return

    if out is not None:
        check_output_path(out, force, inputs=[uri])

    # Resolve contig filter
    selected_contig_ids = resolve_contig_selection(
        mc,
        genomes=genome_list or None,
        contigs=contig_list or None,
        filters=filter_list or None,
    )

    write_index = False
    if table == "bins":
        df = mc.bins()
        if selected_contig_ids is not None:
            bin_chrom_codes = chrom_codes(df)
            mask = np.isin(bin_chrom_codes, list(selected_contig_ids))
            df = df[mask]
        df = df.copy()
        df["chrom"] = _decode_chrom(df["chrom"])
    elif table == "chroms":
        df = mc.contigs[["name", "length"]]
        if selected_contig_ids is not None:
            df = df.loc[df.index.isin(selected_contig_ids)]
    elif table == "contigs":
        df = mc.contigs.copy()
        if selected_contig_ids is not None:
            df = df.loc[df.index.isin(selected_contig_ids)]
        # Surface the genome/MAG name alongside the internal integer genome_id.
        df["genome"] = df["genome_id"].map(mc.genomes["name"])
        core = ["name", "length", "genome", "genome_id"]
        extras = [c for c in df.columns if c not in core]
        df = df[core + extras]
    elif table == "genomes":
        df = mc.genomes
        if selected_contig_ids is not None:
            sel_genome_ids = set(mc.contigs.loc[list(selected_contig_ids), "genome_id"])
            df = df.loc[df.index.isin(sel_genome_ids)]
        core = ["name", "length", "n_contigs"]
        extras = [c for c in df.columns if c not in core]
        df = df[core + extras]
    elif table in ("clusters", "members"):
        # The clustering rounds: every cluster, or every membership, that
        # holds a selected contig -- found through the contig index, so a
        # selection reads its own memberships only.
        if selected_contig_ids is None:
            df = mc.cluster_members() if table == "members" else mc.clusters
        else:
            contig_ids = np.array(sorted(selected_contig_ids), dtype=np.int64)
            cluster_ids, counts = mc._memberships_of(contig_ids)
            if table == "members":
                # In the order of the stored members: cluster, then contig.
                contig_ids = np.repeat(contig_ids, counts)
                order = np.lexsort((contig_ids, cluster_ids))
                df = mc._name_memberships(cluster_ids[order], contig_ids[order])
            else:
                df = mc.clusters.loc[np.unique(cluster_ids)]
    elif table == "pixels":
        # Pre-compute filter parameters that are independent of chunk
        valid_bin_ids: np.ndarray | None = None
        if selected_contig_ids is not None:
            # Through chrom_offset, not the bin table.
            valid_bin_ids = _ranges_to_ids(
                mc._contig_bin_ranges(np.array(sorted(selected_contig_ids)))
            )

        range_bounds: tuple[int, int, int, int] | None = None
        if range1 is not None:
            from metacooler.ops.range import parse_metacooler_range

            i0, i1 = parse_metacooler_range(mc, range1)
            j0, j1 = parse_metacooler_range(mc, range2) if range2 else (i0, i1)
            range_bounds = (i0, i1, j0, j1)

        if pixels_format == "wide":
            pixel_columns = list(mc._empty_pixels().columns)
            if weight not in pixel_columns:
                raise ValueError(
                    f"Weight column {weight!r} not found. Available: {pixel_columns}"
                )

            # Resolve the *selected* bin set before touching the pixel table:
            # the dense matrix is allocated over that set only, and the size
            # guard must fire before anything is materialised.
            n_bins = mc.shape[0]
            selected_bins = np.arange(n_bins, dtype=np.int64)
            if valid_bin_ids is not None:
                selected_bins = np.asarray(valid_bin_ids, dtype=np.int64)
            if range_bounds is not None:
                i0, i1, j0, j1 = range_bounds
                # A range selection keeps pixels in the row block, the column
                # block, or their transpose, so the axis is the union of both.
                in_range = np.union1d(
                    np.arange(i0, i1, dtype=np.int64),
                    np.arange(j0, j1, dtype=np.int64),
                )
                selected_bins = np.intersect1d(
                    selected_bins, in_range, assume_unique=True
                )
            n_sel = int(selected_bins.size)
            if n_sel > DENSE_MATRIX_BIN_LIMIT:
                raise ValueError(
                    f"Wide format refuses to allocate a {n_sel}\u00d7{n_sel} matrix "
                    f"({n_sel**2 * 8 / 1e9:.1f} GB). Use --pixels-format long."
                )

            # Global bin id -> position along the selected axis (-1 = dropped).
            local_of_bin = np.full(n_bins, -1, dtype=np.int64)
            local_of_bin[selected_bins] = np.arange(n_sel, dtype=np.int64)

            mat = np.zeros((n_sel, n_sel), dtype=np.float64)
            for chunk in _selected_pixels(mc, valid_bin_ids, range_bounds):
                b1 = local_of_bin[chunk["bin1_id"].to_numpy(dtype=np.int64)]
                b2 = local_of_bin[chunk["bin2_id"].to_numpy(dtype=np.int64)]
                keep = (b1 >= 0) & (b2 >= 0)
                if not keep.any():
                    continue
                b1, b2 = b1[keep], b2[keep]
                vals = chunk[weight].to_numpy()[keep].astype(np.float64)
                mat[b1, b2] = vals
                mat[b2, b1] = vals
            if join:
                bdf = mc.bins()[["chrom", "start", "end"]].iloc[selected_bins]
                chrom = _decode_chrom(bdf["chrom"]).astype(str).to_numpy()
                labels = [
                    f"{c}:{s}-{e}"
                    for c, s, e in zip(
                        chrom,
                        bdf["start"].to_numpy(),
                        bdf["end"].to_numpy(),
                        strict=False,
                    )
                ]
                df = pd.DataFrame(mat, index=labels, columns=labels)
                write_index = True
            else:
                df = pd.DataFrame(mat)
        else:
            # Long format: stream pixel chunks to output without full materialisation
            ann: pd.DataFrame | None = None
            joined_empty_cols: list[str] | None = None
            if join:
                import cooler

                ann = mc.bins()[["chrom", "start", "end"]].copy()
                ann["chrom"] = _decode_chrom(ann["chrom"])
                pixel_cols = [
                    c
                    for c in mc._empty_pixels().columns
                    if c not in ("bin1_id", "bin2_id")
                ]
                joined_empty_cols = [
                    "contig1",
                    "start1",
                    "end1",
                    "contig2",
                    "start2",
                    "end2",
                    *pixel_cols,
                ]
            out_stream = open(out, "w") if out else sys.stdout  # noqa: SIM115
            try:
                first_chunk = True
                for chunk in _selected_pixels(mc, valid_bin_ids, range_bounds):
                    if join:
                        assert ann is not None
                        chunk = cooler.api.annotate(chunk, ann, replace=True).rename(
                            columns={"chrom1": "contig1", "chrom2": "contig2"}
                        )
                    chunk.to_csv(
                        out_stream,
                        sep="\t",
                        index=False,
                        header=header and first_chunk,
                    )
                    first_chunk = False
                if first_chunk and header:
                    # Nothing matched: still honour --header.
                    empty = (
                        pd.DataFrame(columns=joined_empty_cols)
                        if join
                        else mc._empty_pixels()
                    )
                    empty.to_csv(out_stream, sep="\t", index=False, header=True)
            finally:
                if out:
                    out_stream.close()
            return
    elif table == "indexes":
        rows = []
        with mc.open() as f:
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
            df.to_csv(
                out_stream, sep="\t", index=write_index, header=header or write_index
            )
    else:
        df.to_csv(sys.stdout, sep="\t", index=write_index, header=header or write_index)
