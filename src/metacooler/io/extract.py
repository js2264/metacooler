"""Extract a genome-level cool file from a metacool file."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from metacooler.io._utils import (
    enforce_upper_triangle,
    resolve_contig_selection,
    resolve_count_dtype,
    write_cooler,
)

logger = logging.getLogger("metacooler")


def extract(
    uri: str | Path,
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
    filters: list[str] | None = None,
    output: str | Path = "extracted.cool",
    split: bool = False,
    force: bool = False,
) -> None:
    """Extract selected genomes and/or contigs from a metacool file.

    In single-file mode, ``genomes`` and ``contigs`` can be combined: the
    union of all contigs from the specified genomes plus the directly
    specified contigs is extracted.  In ``split`` mode only ``genomes`` is
    allowed (one file per genome); combining it with ``contigs`` is an error.

    ``filters`` composes with both modes: it narrows whatever the names
    selected (or the whole file, when no name is given), so ``--split
    --filter 'genomes.length > 1Mb'`` writes one cool per genome above 1 Mb.

    Args:
        uri: Path to the source ``.metacool`` file.
        genomes: List of genome names to extract.
        contigs: List of contig names to extract (single-file mode only).
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`).
        output: Destination path. If *split* is True, treated as a
            filename prefix (one ``<prefix>_<genome>.cool`` per genome).
        split: If True, write one ``.cool`` per genome instead of a
            single merged file. Rejects ``contigs``.
        force: Overwrite existing per-genome output files (split mode).

    Raises:
        KeyError: If any genome or contig is not found.
        ValueError: If no selector at all is specified, or if ``split`` is
            combined with ``contigs``.
    """
    from metacooler.api import MetaCooler
    from metacooler.io.staging import writing_output

    if not genomes and not contigs and not filters:
        raise ValueError(
            "At least one of --genomes, --contigs or --filter must be specified."
        )
    if split and contigs:
        raise ValueError(
            "--contigs cannot be combined with --split. Split writes one file per "
            "genome; pass --genomes only, or drop --split to extract the contigs "
            "into a single file."
        )

    mc = MetaCooler(uri)
    selected_contig_ids = resolve_contig_selection(
        mc, genomes=genomes, contigs=contigs, filters=filters
    )
    assert selected_contig_ids is not None  # a selector was required above

    if split:
        # One cool per genome, over the selected contigs only: with --genomes
        # that is exactly the named genomes, and with --filter it is whichever
        # genomes still have a contig left.
        all_contigs = mc.contigs
        genome_names = mc.genomes["name"]
        output_path = Path(output)
        # Strip only a trailing ".cool": Path.stem cut at the *last* dot, so
        # "-o sample.rep1" wrote "sample_<genome>.cool", and a second run with
        # "-o sample.rep2" collided with (or under --force, overwrote) it.
        name = output_path.name
        prefix = name[: -len(".cool")] if name.endswith(".cool") else name
        by_genome: dict[int, list[int]] = {}
        for cid in sorted(selected_contig_ids):
            by_genome.setdefault(int(all_contigs.at[cid, "genome_id"]), []).append(cid)
        for gid, cids in sorted(by_genome.items()):
            g = str(genome_names[gid])
            out_file = output_path.parent / f"{prefix}_{g}.cool"
            with writing_output(out_file, force=force):
                if out_file.exists() and not force:
                    raise FileExistsError(
                        f"Output file already exists: {out_file}. Use --force to "
                        "overwrite."
                    )
                _extract_contig_ids(mc, cids, out_file, assembly=g)
    else:
        # Single output file with all selected contigs
        with writing_output(output):
            _extract_contig_ids(
                mc, sorted(selected_contig_ids), Path(output), assembly="multi"
            )


def _extract_contig_ids(
    mc,
    contig_ids: list[int],
    output: Path,
    assembly: str = "unknown",
) -> None:
    """Extract bins/pixels for a set of contig IDs and write a cool file.

    The selection is resolved to bin ranges through ``indexes/chrom_offset``
    and its pixels are read by row block through ``indexes/bin1_offset``, so
    neither the whole bin table nor the pixels outside the selected rows are
    read.  Renumbering the kept bins in order is monotonic, so a
    ``symmetric-upper`` source streams out sorted and is written without a
    merge.
    """
    from metacooler.api import _in_ranges

    ranges = mc._contig_bin_ranges(np.asarray(sorted(contig_ids), dtype=np.int64))
    genome_bins = mc._bins_in_ranges(ranges)
    genome_bin_indices = genome_bins.index.to_numpy(dtype=np.int64)

    bins_df = pd.DataFrame(
        {
            "chrom": genome_bins["chrom"].values,
            "start": genome_bins["start"].to_numpy().astype(np.int32),
            "end": genome_bins["end"].to_numpy().astype(np.int32),
        }
    )

    # Preserve the width of the source ``count`` column: a binnified source
    # writes int64 and would silently wrap if narrowed back to int32 here.
    count_dtype = resolve_count_dtype(
        mc._empty_pixels()["count"].dtype, context="extract"
    )
    square = mc.info.get("storage-mode") == "square"

    def _pixel_stream():
        rows = mc._iter_pixel_rows(ranges, lambda b1, b2: _in_ranges(b2, ranges))
        for chunk in rows:
            local = pd.DataFrame(
                {
                    "bin1_id": np.searchsorted(
                        genome_bin_indices, chunk["bin1_id"].to_numpy()
                    ).astype(np.int64),
                    "bin2_id": np.searchsorted(
                        genome_bin_indices, chunk["bin2_id"].to_numpy()
                    ).astype(np.int64),
                    "count": chunk["count"].to_numpy().astype(count_dtype),
                }
            )
            if square:
                enforce_upper_triangle(local)
            yield local

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
        write_cooler(
            output,
            bins_df,
            empty_pixels,
            ordered=True,
            symmetric_upper=True,
            assembly=assembly,
            dtypes={"count": count_dtype},
            metadata={"source-file": str(mc.store)},
        )
    else:

        def _all_chunks():
            yield first_chunk
            yield from pixel_stream

        # Kept bins are renumbered in order, so an upper-triangle source stays
        # sorted; a square one is folded per chunk and needs the merge.
        write_cooler(
            output,
            bins_df,
            _all_chunks(),
            ordered=not square,
            symmetric_upper=True,
            assembly=assembly,
            dtypes={"count": count_dtype},
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
