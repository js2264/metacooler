"""Extract a genome-level cool file from a metacool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import _filter_by_genomes_and_contigs, resolve_count_dtype

logger = logging.getLogger("metacooler")


def extract(
    uri: str | Path,
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
    output: str | Path = "extracted.cool",
    split: bool = False,
    force: bool = False,
) -> None:
    """Extract selected genomes and/or contigs from a metacool file.

    In single-file mode, ``genomes`` and ``contigs`` can be combined: the
    union of all contigs from the specified genomes plus the directly
    specified contigs is extracted.  In ``split`` mode only ``genomes`` is
    allowed (one file per genome); combining it with ``contigs`` is an error.

    Args:
        uri: Path to the source ``.metacool`` file.
        genomes: List of genome names to extract.
        contigs: List of contig names to extract (single-file mode only).
        output: Destination path. If *split* is True, treated as a
            filename prefix (one ``<prefix>_<genome>.cool`` per genome).
        split: If True, write one ``.cool`` per genome instead of a
            single merged file. Requires ``genomes``; rejects ``contigs``.
        force: Overwrite existing per-genome output files (split mode).

    Raises:
        KeyError: If any genome or contig is not found.
        ValueError: If neither genomes nor contigs are specified, or if
            ``split`` is combined with ``contigs`` / used without ``genomes``.
    """
    from metacooler.api import MetaCooler

    if not genomes and not contigs:
        raise ValueError("At least one of --genomes or --contigs must be specified.")
    if split and contigs:
        raise ValueError(
            "--contigs cannot be combined with --split. Split writes one file per "
            "genome; pass --genomes only, or drop --split to extract the contigs "
            "into a single file."
        )

    mc = MetaCooler(uri)

    if split:
        # One cool per genome (contigs are excluded above; genomes is non-empty
        # here because split+contigs and the empty-selection case both raised).
        all_genomes_df = mc.genomes
        all_contigs = mc.contigs
        output_path = Path(output)
        for g in genomes or []:
            gid = int(all_genomes_df.loc[all_genomes_df["name"] == g].index[0])
            cids = all_contigs.index[all_contigs["genome_id"] == gid].tolist()
            out_file = output_path.parent / f"{output_path.stem}_{g}.cool"
            if out_file.exists() and not force:
                raise FileExistsError(
                    f"Output file already exists: {out_file}. Use --force to overwrite."
                )
            _extract_contig_ids(mc, sorted(cids), out_file, assembly=g)
    else:
        # Single output file with all selected contigs
        selected_contig_ids = _filter_by_genomes_and_contigs(
            mc, genomes=genomes, contigs=contigs
        )
        _extract_contig_ids(
            mc, sorted(selected_contig_ids), Path(output), assembly="multi"
        )


def _extract_contig_ids(
    mc,
    contig_ids: list[int],
    output: Path,
    assembly: str = "unknown",
) -> None:
    """Extract bins/pixels for a set of contig IDs and write a cool file."""
    all_bins = mc.bins()
    bin_chrom_codes = all_bins["chrom"].cat.codes.values

    # Boolean mask: bins belonging to selected contigs
    genome_bin_mask = np.isin(bin_chrom_codes, contig_ids)
    genome_bin_indices = np.where(genome_bin_mask)[0]

    global_to_local_bin = {
        int(gbin): local for local, gbin in enumerate(genome_bin_indices)
    }

    genome_bins = all_bins.iloc[genome_bin_indices].copy()
    genome_bins_chrom_names = genome_bins["chrom"].values

    bins_df = pd.DataFrame(
        {
            "chrom": genome_bins_chrom_names,
            "start": genome_bins["start"].values.astype(np.int32),
            "end": genome_bins["end"].values.astype(np.int32),
        }
    )

    # Preserve the width of the source ``count`` column: a binnified source
    # writes int64 and would silently wrap if narrowed back to int32 here.
    count_dtype = resolve_count_dtype(
        mc._clr.pixels()[:0]["count"].dtype, context="extract"
    )

    genome_bin_set = set(genome_bin_indices.tolist())
    genome_bin_array = np.array(sorted(genome_bin_set), dtype=np.int64)
    valid_min, valid_max = genome_bin_array[0], genome_bin_array[-1]
    chunksize = 1_000_000

    # Use bin1_offset index to find pixel range for this genome's bins
    with h5py.File(str(mc.store), "r") as _f:
        _root = _f[mc._group_path.lstrip("/")] if mc._group_path != "/" else _f
        bin1_offset = _root["indexes"]["bin1_offset"][:]
    pix_lo = int(bin1_offset[valid_min])
    n_pix = len(mc._clr.pixels())
    pix_hi = (
        int(bin1_offset[valid_max + 1]) if valid_max + 1 < len(bin1_offset) else n_pix
    )

    def _pixel_stream():
        for _i in range(pix_lo, pix_hi, chunksize):
            chunk = mc._clr.pixels()[_i : _i + chunksize].copy()
            mask = chunk["bin1_id"].isin(genome_bin_set) & chunk["bin2_id"].isin(
                genome_bin_set
            )
            if not mask.any():
                continue
            chunk = chunk[mask].copy()
            chunk["bin1_id"] = chunk["bin1_id"].map(global_to_local_bin)
            chunk["bin2_id"] = chunk["bin2_id"].map(global_to_local_bin)
            chunk["bin1_id"] = chunk["bin1_id"].astype(np.int64)
            chunk["bin2_id"] = chunk["bin2_id"].astype(np.int64)
            chunk["count"] = chunk["count"].astype(count_dtype)
            swap = chunk["bin1_id"] > chunk["bin2_id"]
            chunk.loc[swap, ["bin1_id", "bin2_id"]] = chunk.loc[
                swap, ["bin2_id", "bin1_id"]
            ].values
            yield chunk[["bin1_id", "bin2_id", "count"]]

    pixel_stream = _pixel_stream()
    try:
        first_chunk = next(pixel_stream)
    except StopIteration:
        # No pixels for the selected contigs: write an empty cool file
        empty_pixels = pd.DataFrame(
            {
                "bin1_id": pd.Series(dtype=np.int64),
                "bin2_id": pd.Series(dtype=np.int64),
                "count": pd.Series(dtype=count_dtype),
            }
        )
        cooler.create_cooler(
            str(output),
            bins=bins_df,
            pixels=empty_pixels,
            dtypes={"count": count_dtype},
            ordered=True,
            symmetric_upper=True,
            assembly=assembly,
            metadata={"source-file": str(mc.store)},
        )
    else:

        def _all_chunks():
            yield first_chunk
            yield from pixel_stream

        cooler.create_cooler(
            str(output),
            bins=bins_df,
            pixels=_all_chunks(),
            dtypes={"count": count_dtype},
            ordered=False,
            symmetric_upper=True,
            assembly=assembly,
            metadata={"source-file": str(mc.store)},
        )

    c = cooler.Cooler(str(output))
    n_bins = len(genome_bins)
    n_contigs = len(set(genome_bins["chrom"]))
    nnz = int(c.info["nnz"])
    total = int(c.info["sum"])
    logger.info(
        f"Created {output}: bins={n_bins}, contigs={n_contigs}. "
        f"NNZ={nnz}, sum={total}.",
    )
