"""Re-binning operation for metacool files."""

from __future__ import annotations

import logging
from itertools import pairwise
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from metacooler.io._utils import (
    check_output_path,
    chrom_codes,
    pixel_keys,
    pixels_from_keys,
    reduce_sorted_keys,
    resolve_count_dtype,
    to_upper_triangle,
)
from metacooler.qc import (
    is_qc_column,
    strip_contig_qc_columns,
    strip_qc_columns,
    strip_qc_metadata,
)
from metacooler.schema import CLUSTERS_GROUP, DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Literal

    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

# Default on-disk dtype for the aggregated ``pixels/count`` column.
_DEFAULT_COUNT_DTYPE = np.dtype(np.int64)


def binnify(
    uri: str | Path,
    output: str | Path,
    binsize: int,
    level: Literal["contig", "genome"] = "contig",
    chunksize: int = DEFAULT_CHUNKSIZE,
    force: bool = False,
) -> None:
    """Re-bin a metacool file to a new, coarser fixed bin size.

    The target *binsize* must be **strictly coarser** than the source binning:

    - for a fixed-bin source, *binsize* must be a strict, **exact multiple** of
      the source bin size.  A non-multiple target is rejected: the midpoint
      rule would pack an oscillating number of source bins into each output bin
      (e.g. 1 kb → 3.5 kb alternates 3 and 4 source bins, a built-in ±14 %
      checkerboard) while the file would advertise a uniform bin size;
    - for a variable-bin source (fragment-level files, contig-level files) there
      is no exact multiple, so the *median* source bin width is used as the
      effective source resolution and *binsize* must be strictly larger than it.
      A target that is larger than the median but smaller than the widest source
      bin is accepted with a warning (those widest bins each collapse into a
      single output bin, so the local resolution is uneven).

    Args:
        uri: Path to the source ``.metacool`` file.
        output: Destination path for the output ``.metacool`` file.
        binsize: New bin size in base pairs.
        level: Tiling level. ``contig`` (default) tiles each contig
            independently, preserving the contig hierarchy. ``genome``
            artificially ligates each genome's contigs (ordered by decreasing
            length) into one sequence and tiles that uniformly; the output's
            contigs become the ligated genomes (one per genome), so bins may
            span the original contig boundaries.
        chunksize: Number of source pixels loaded, remapped and aggregated per
            chunk.  Pixels are never fully materialised in memory.
        force: Overwrite *output* if it already exists.

    Raises:
        ValueError: If *binsize* is not strictly coarser than the source
            resolution, if it is not an exact multiple of a fixed source bin
            size, or if *level* is unknown.
        FileExistsError: If *output* exists and *force* is False.
    """
    from metacooler.api import MetaCooler
    from metacooler.io.staging import writing_output

    output = Path(output)
    check_output_path(output, force=force, inputs=[uri])
    if binsize <= 0:
        raise ValueError(f"binsize must be a positive integer, got {binsize}.")
    # Checked before the write begins: a refusal fails no edit session whose
    # staged map is the output.
    with writing_output(output, force=force):
        mc = MetaCooler(uri)
        _binnify(
            mc,
            output,
            binsize,
            level=level,
            chunksize=chunksize,
            source_uri=uri,
            source_metadata=mc.metadata,
        )


def _binnify(
    mc: MetaCooler,
    output: Path,
    binsize: int,
    *,
    level: Literal["contig", "genome"] = "contig",
    chunksize: int = DEFAULT_CHUNKSIZE,
    source_uri: str | Path,
    source_metadata: dict,
) -> None:
    """:func:`binnify` of an open map, with its provenance given separately.

    ``zoomify`` coarsens a level from a finer level it has already written,
    but records the *base* map as the source, as it did when every level was
    re-binned from the base; hence *source_uri* and *source_metadata*.
    """
    uri = source_uri
    contigs = mc.contigs
    genomes = mc.genomes
    old_bins = mc.bins()

    _check_target_binsize(mc.binsize, old_bins, binsize)

    # Get contig codes for each bin
    bin_chrom_codes = chrom_codes(old_bins)

    if level == "contig":
        # Tile each contig independently. Built vectorised and already in
        # canonical (contig-id, start) order, so bin rows match /contigs.
        new_bins = _tile_uniform(contigs["length"].values, binsize)
        if new_bins.empty:
            raise ValueError("No bins generated. Check input data.")
        old_to_new = _map_old_bins_per_contig(
            old_bins, bin_chrom_codes, len(contigs), new_bins
        )
    elif level == "genome":
        # Ligate each genome's contigs (ordered by decreasing length) into a
        # single artificial sequence, then tile that sequence uniformly. The
        # output's "contigs" ARE the ligated genomes (one contig per genome), so
        # bins genuinely span contig boundaries within a genome.
        new_bins, old_to_new, contigs, genomes = _tile_per_genome_ligated(
            contigs, genomes, old_bins, bin_chrom_codes, binsize
        )
        if new_bins.empty:
            raise ValueError("No bins generated. Check input data.")
    else:
        raise ValueError(f"level must be 'contig' or 'genome', got {level!r}")

    # Filter out unmapped bins
    valid_mask = old_to_new >= 0

    # Build output bins with contig names
    # Contig *ids*, as create writes them: cooler.create only needs a chrom
    # value per bin, and the /chroms table it derives from them is replaced by
    # the link to /contigs (_finalize_metacool_file).  Names would make cooler
    # build and iterate a column of two million Arrow strings, several times,
    # on every level of a zoomify.
    new_bins_out = pd.DataFrame(
        {
            "chrom": new_bins["chrom"].to_numpy(dtype=np.int32),
            "start": new_bins["start"].values.astype(np.int32),
            "end": new_bins["end"].values.astype(np.int32),
        }
    )

    # Recompute n_bins per genome from the new bins
    new_bin_genome_ids = contigs["genome_id"].values[new_bins_out["chrom"].to_numpy()]
    new_genome_bin_counts = np.bincount(
        new_bin_genome_ids.astype(np.int64), minlength=len(genomes)
    ).astype(np.int64)
    genomes = genomes.copy()
    genomes["n_bins"] = new_genome_bin_counts

    # A coarse bin pair sums many source pixels, so integer counts are always
    # widened to int64; float counts (e.g. balanced values) stay float64.
    src_count_dtype = resolve_count_dtype(
        mc._clr.pixels()[:0]["count"].dtype, context="binnify"
    )
    count_dtype = (
        _DEFAULT_COUNT_DTYPE if src_count_dtype.kind == "i" else src_count_dtype
    )

    # The row-block path needs every source pixel's new row to be the new row
    # of its bin1: true when the map is monotonic (contig order is kept and
    # bins only coarsen within a contig, i.e. --level contig) and the source
    # stores the upper triangle (a square source's lower half folds over).
    ordered = mc._clr.info.get("storage-mode") != "square" and _is_monotonic(
        old_to_new[valid_mask]
    )
    if ordered:
        pixels = _coarsen_row_blocks(
            mc, old_to_new, chunksize=chunksize, count_dtype=count_dtype
        )
    else:
        pixels = _remap_pixel_chunks(
            mc, old_to_new, valid_mask, chunksize=chunksize, count_dtype=count_dtype
        )

    # Write as a new metacool file using cooler + h5py so the metacooler
    # structure is preserved.
    _write_rebinned_metacool(
        output=output,
        bins=new_bins_out,
        pixels=pixels,
        ordered=ordered,
        genomes=genomes,
        contigs=contigs,
        binsize=binsize,
        level=level,
        source_uri=uri,
        source_metadata=source_metadata,
        count_dtype=count_dtype,
        assembly_name=mc.info.get("assembly"),
        # api.info always inserts "ordering" (possibly None), so .get()'s
        # default never fires; `or` is what actually supplies the fallback.
        contig_ordering=mc.info.get("ordering") or "genome-then-alpha",
    )
    # The rounds were copied from the source by path, after its pixels were
    # read: the map read is the map copied (MapChangedError otherwise).
    mc._check_unchanged()

    logger.info(
        "Re-binned %s → %s: binsize=%d, %d bins",
        uri,
        output,
        binsize,
        len(new_bins_out),
    )


def _check_target_binsize(
    source_binsize: int | None, old_bins: pd.DataFrame, binsize: int
) -> None:
    """Refuse a target bin size that is not strictly coarser than the source.

    Args:
        source_binsize: Fixed source bin size, or ``None`` for a variable-bin
            source.
        old_bins: Source bins table (used to derive an effective resolution for
            variable-bin sources).
        binsize: Requested target bin size in base pairs.

    Raises:
        ValueError: If *binsize* is finer than (or equal to) the source
            resolution, or is not an exact multiple of a fixed *source_binsize*.
    """
    if source_binsize is not None:
        if binsize <= source_binsize:
            raise ValueError(
                f"New binsize ({binsize}) must be strictly larger than current "
                f"binsize ({source_binsize})."
            )
        remainder = binsize % source_binsize
        if remainder:
            lower = (binsize // source_binsize) * source_binsize
            upper = lower + source_binsize
            raise ValueError(
                f"New binsize ({binsize}) must be an exact multiple of the source "
                f"binsize ({source_binsize}) — {binsize} % {source_binsize} = "
                f"{remainder}. A non-multiple target puts an oscillating number of "
                "source bins into each output bin (a systematic checkerboard) while "
                f"the file would advertise a uniform {binsize} bp bin size. Nearest "
                f"valid multiples: {lower} and {upper}."
            )
        return

    # Variable-bin source: there is no exact multiple to check, so use the
    # median bin width as the effective source resolution.
    widths = old_bins["end"].values.astype(np.int64) - old_bins["start"].values.astype(
        np.int64
    )
    if widths.size == 0:
        raise ValueError("Source file contains no bins.")
    median_width = float(np.median(widths))
    max_width = int(widths.max())
    if binsize <= median_width:
        raise ValueError(
            f"New binsize ({binsize}) must be strictly larger than the effective "
            f"resolution of this variable-bin source (median bin width "
            f"{median_width:g} bp; widest bin {max_width} bp). Re-binning to a "
            "finer size would leave most output bins empty."
        )
    if binsize < max_width:
        logger.warning(
            "Target binsize %d is coarser than the median source bin width (%g bp) "
            "but finer than the widest source bin (%d bp): those widest bins each "
            "collapse into a single output bin, so the effective resolution is "
            "uneven across the file.",
            binsize,
            median_width,
            max_width,
        )


def _tile_uniform(lengths: np.ndarray, binsize: int) -> pd.DataFrame:
    """Tile each sequence of *lengths* with fixed-width bins, vectorised.

    Args:
        lengths: Sequence lengths in base pairs (one entry per contig/genome).
        binsize: Bin width in base pairs.

    Returns:
        DataFrame with columns ``chrom`` (0-based sequence index), ``start`` and
        ``end``, already in canonical ``(chrom, start)`` order with a
        ``RangeIndex`` giving the bin id.  The last bin of each sequence is
        truncated to its length.
    """
    lengths = np.asarray(lengths).astype(np.int64)
    counts = -(-lengths // binsize)  # ceil division, 0 for zero-length entries
    total = int(counts.sum())
    if total == 0:
        return pd.DataFrame(
            {
                "chrom": np.array([], dtype=np.int64),
                "start": np.array([], dtype=np.int64),
                "end": np.array([], dtype=np.int64),
            }
        )
    chrom = np.repeat(np.arange(len(lengths), dtype=np.int64), counts)
    first_bin = np.cumsum(counts) - counts
    within = np.arange(total, dtype=np.int64) - np.repeat(first_bin, counts)
    start = within * binsize
    end = np.minimum(start + binsize, np.repeat(lengths, counts))
    return pd.DataFrame({"chrom": chrom, "start": start, "end": end})


def _remap_pixel_chunks(
    mc: MetaCooler,
    old_to_new: np.ndarray,
    valid_mask: np.ndarray,
    chunksize: int = DEFAULT_CHUNKSIZE,
    count_dtype: np.dtype = _DEFAULT_COUNT_DTYPE,
) -> Generator[pd.DataFrame, None, None]:
    """Stream source pixels, remap their bin ids and aggregate per chunk.

    Each chunk is remapped through *old_to_new*, folded into the upper triangle
    and aggregated on ``(bin1_id, bin2_id)`` before being yielded, so no more
    than one chunk of pixels is ever resident.  Cross-chunk duplicates are
    resolved by ``cooler.create_cooler(..., ordered=False)``, whose merge pass
    (``cooler._reduce.CoolerMerger``) sums duplicate bin pairs.

    Args:
        mc: Open source MetaCooler.
        old_to_new: Map from source bin id to output bin id (``-1`` = drop).
        valid_mask: Boolean mask of source bins that map to an output bin.
        chunksize: Number of source pixels loaded per chunk.
        count_dtype: dtype for the accumulated ``count`` column.

    Yields:
        One ``pd.DataFrame`` per chunk with columns ``bin1_id``, ``bin2_id``,
        ``count``, sorted by ``(bin1_id, bin2_id)`` and free of duplicate pairs.
    """
    n_pix = len(mc._clr.pixels())
    emitted = False
    for lo in range(0, n_pix, chunksize):
        chunk = mc._clr.pixels()[lo : lo + chunksize]
        b1 = chunk["bin1_id"].values
        b2 = chunk["bin2_id"].values
        # Widened up front so the per-chunk and cross-chunk sums cannot overflow
        # the source's (typically int32) count dtype.
        counts = chunk["count"].values.astype(count_dtype)

        pixel_mask = valid_mask[b1] & valid_mask[b2]
        if not pixel_mask.any():
            continue
        new_b1 = old_to_new[b1[pixel_mask]]
        new_b2 = old_to_new[b2[pixel_mask]]
        new_counts = counts[pixel_mask]

        # Enforce upper triangle
        new_b1, new_b2 = to_upper_triangle(new_b1, new_b2)

        # Sorted and duplicate-free, as cooler's unordered writer requires.
        emitted = True
        yield _sum_by_key(new_b1, new_b2, new_counts)

    if not emitted:
        # cooler's unordered path needs at least one chunk to merge.
        yield pd.DataFrame(
            {
                "bin1_id": np.array([], dtype=np.int64),
                "bin2_id": np.array([], dtype=np.int64),
                "count": np.array([], dtype=count_dtype),
            }
        )


def _is_monotonic(values: np.ndarray) -> bool:
    return bool(len(values) < 2 or (values[1:] >= values[:-1]).all())


def _sum_by_key(b1: np.ndarray, b2: np.ndarray, counts: np.ndarray) -> pd.DataFrame:
    """Sum *counts* per ``(b1, b2)``; the result is sorted (see aggregate_pixels).

    Re-mapped source pixels arrive as one sorted run per source row (a bin map
    that keeps contig order within each contig keeps each row's ``bin2`` order),
    which a stable sort merges about twice as fast as quicksort sorts them.
    """
    keys = pixel_keys(b1, b2)
    order = np.argsort(keys, kind="stable")
    return pixels_from_keys(*reduce_sorted_keys(keys[order], counts[order]))


def _row_blocks(
    old_to_new: np.ndarray, bin1_offset: np.ndarray, chunksize: int
) -> list[tuple[int, int]]:
    """Cut the source rows into blocks of whole *output* rows, ~*chunksize* pixels.

    With a monotonic *old_to_new*, the source bins of one output row are
    contiguous, so the source can be cut wherever the output row changes.  The
    cuts are chosen so each block holds about *chunksize* source pixels (a
    single output row larger than that is one block on its own).

    Returns:
        ``[(lo, hi), ...]`` half-open source-bin ranges, in order, covering
        every source bin.
    """
    n_old = len(old_to_new)
    if n_old == 0:
        return []
    # A dropped bin (-1) inside an output row must not start a new row:
    # carry the previous row over it (a running max is a forward fill here,
    # since the kept ids never decrease).
    rows = np.maximum.accumulate(old_to_new)
    starts = np.flatnonzero(np.concatenate(([True], rows[1:] != rows[:-1])))
    pixels_before = bin1_offset[starts]
    targets = np.arange(chunksize, int(bin1_offset[-1]), chunksize, dtype=np.int64)
    picks = np.searchsorted(pixels_before, targets, side="right") - 1
    cuts = np.unique(np.concatenate(([0], starts[np.maximum(picks, 0)], [n_old])))
    return [(int(lo), int(hi)) for lo, hi in pairwise(cuts)]


def _coarsen_row_blocks(
    mc: MetaCooler,
    old_to_new: np.ndarray,
    chunksize: int = DEFAULT_CHUNKSIZE,
    count_dtype: np.dtype = _DEFAULT_COUNT_DTYPE,
) -> Generator[pd.DataFrame, None, None]:
    """Stream the re-binned pixels in final order, one block of output rows at a time.

    Requires a monotonic *old_to_new* and an upper-triangular source (see
    :func:`binnify`): then a source pixel ``(b1, b2)`` lands on ``(old_to_new[b1],
    old_to_new[b2])`` with the new ``bin1 <= bin2``, and every pixel of an
    output row comes from one contiguous run of source rows.  Each block's rows
    are located through ``indexes/bin1_offset`` and read with cooler's pixel
    selector, then mapped, sorted on their pixel key and reduced;
    blocks follow each other in output-row order, so their concatenation is the
    finished, sorted pixel table, written without cooler's external merge, in
    memory proportional to one block.

    Args:
        mc: Open source map.
        old_to_new: Map from source bin id to output bin id (``-1`` = drop).
        chunksize: Source pixels per block (whole output rows are kept together).
        count_dtype: dtype of the accumulated ``count`` column.

    Yields:
        Sorted, duplicate-free pixel frames, in order (at least one, possibly
        empty, so cooler's ordered writer always has a frame to write).
    """
    bin1_offset = mc._bin1_offset
    pixels = mc._clr.pixels()[["bin1_id", "bin2_id", "count"]]
    emitted = False
    for lo, hi in _row_blocks(old_to_new, bin1_offset, chunksize):
        p0, p1 = int(bin1_offset[lo]), int(bin1_offset[hi])
        if p0 == p1:
            continue
        chunk = pixels[p0:p1]
        b1 = old_to_new[chunk["bin1_id"].to_numpy()]
        b2 = old_to_new[chunk["bin2_id"].to_numpy()]
        counts = chunk["count"].to_numpy().astype(count_dtype)
        keep = (b1 >= 0) & (b2 >= 0)
        if not keep.all():
            b1, b2, counts = b1[keep], b2[keep], counts[keep]
        if not len(b1):
            continue
        emitted = True
        yield _sum_by_key(b1, b2, counts)
    if not emitted:
        yield pd.DataFrame(
            {
                "bin1_id": np.array([], dtype=np.int64),
                "bin2_id": np.array([], dtype=np.int64),
                "count": np.array([], dtype=count_dtype),
            }
        )


def _map_old_bins_per_contig(
    old_bins: pd.DataFrame,
    bin_chrom_codes: np.ndarray,
    n_contigs: int,
    new_bins: pd.DataFrame,
) -> np.ndarray:
    """Map each old bin to a new bin within the SAME contig (midpoint rule).

    The target bin size is guaranteed (by :func:`_check_target_binsize`) to be
    an exact multiple of a fixed source bin size, so every output bin receives
    the same number of source bins.

    One global ``searchsorted``, not one per contig.  The per-contig form also
    rebuilt an ``old_chrom == cid`` mask on every iteration, which is a full
    pass over the source bins per contig -- O(n_contigs x n_bins), so a
    metagenome with millions of contigs never finished.  Both output bins and
    source bins are already in ``(contig, start)`` order, so folding the contig
    id into the sort key above every coordinate turns that lexicographic order
    into plain integer order and the whole mapping becomes one sorted lookup.
    """
    old_chrom = bin_chrom_codes.astype(np.int64)
    # Integer midpoints: every bin boundary is an integer, so flooring cannot
    # move a midpoint across one, and it keeps the key arithmetic exact.
    old_mid = (
        old_bins["start"].values.astype(np.int64)
        + old_bins["end"].values.astype(np.int64)
    ) // 2
    new_chrom = new_bins["chrom"].values.astype(np.int64)
    new_start = new_bins["start"].values.astype(np.int64)
    new_ids = new_bins.index.values.astype(np.int64)

    old_to_new = np.full(len(old_bins), -1, dtype=np.int64)
    if new_start.size == 0 or old_mid.size == 0:
        return old_to_new

    # Where each contig's run of output bins starts, so a hit that fell short of
    # its own contig can be told from one inside it.
    new_count = np.bincount(new_chrom, minlength=n_contigs)
    new_first = np.concatenate([[0], np.cumsum(new_count)])[:n_contigs]

    scale = int(max(int(new_start.max()), int(old_mid.max()))) + 1
    idx = (
        np.searchsorted(
            new_chrom * scale + new_start,
            old_chrom * scale + old_mid,
            side="right",
        )
        - 1
    )
    # The two exclusions the per-contig loop made: a contig with no output bin
    # at all, and a midpoint preceding its contig's first bin (which lands the
    # search in an earlier contig's run, hence below new_first).  A midpoint
    # past a stale /contigs/length still lands on that contig's last bin,
    # because the next key belongs to the next contig.
    ok = (new_count[old_chrom] > 0) & (idx >= new_first[old_chrom])
    old_to_new[ok] = new_ids[idx[ok]]
    return old_to_new


def _ligated_contig_offsets(
    contig_len: np.ndarray, contig_gid: np.ndarray, n_genomes: int
) -> tuple[np.ndarray, np.ndarray]:
    """Lay out each genome's contigs end-to-end, longest first.

    Args:
        contig_len: Contig lengths.
        contig_gid: Genome id of each contig.
        n_genomes: Total number of genomes.

    Returns:
        ``(contig_offset, genome_length)``: the start position of each contig
        inside its ligated genome, and the total ligated length per genome.
    """
    genome_length = np.bincount(
        contig_gid, weights=contig_len, minlength=n_genomes
    ).astype(np.int64)
    # Sort by (genome, decreasing length, contig id) — deterministic ties.
    order = np.lexsort(
        (np.arange(len(contig_len), dtype=np.int64), -contig_len, contig_gid)
    )
    sorted_len = contig_len[order]
    sorted_gid = contig_gid[order]
    # Exclusive cumulative sum, made genome-local by subtracting the cumulative
    # length of all preceding genomes (the sort keeps each genome contiguous).
    genome_start = np.concatenate([[0], np.cumsum(genome_length)])
    within = (np.cumsum(sorted_len) - sorted_len) - genome_start[sorted_gid]
    contig_offset = np.empty(len(contig_len), dtype=np.int64)
    contig_offset[order] = within
    return contig_offset, genome_length


def _tile_per_genome_ligated(
    contigs: pd.DataFrame,
    genomes: pd.DataFrame,
    old_bins: pd.DataFrame,
    bin_chrom_codes: np.ndarray,
    binsize: int,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Ligate contigs per genome (decreasing length) and tile at genome level.

    Genomes with no contigs are dropped (a zero-length contig is rejected by
    :func:`metacooler.schema.validate`) and the remaining genome ids are
    renumbered.  Extra ``/genomes`` columns are carried onto the output contigs,
    which are the ligated genomes themselves.

    Returns ``(new_bins, old_to_new, out_contigs, out_genomes)`` where the
    output contigs are the ligated genomes (one per genome, named after the
    genome) so that a re-binned bin can span the original contig boundaries.
    """
    from metacooler.schema import CONTIGS_DTYPE, GENOMES_DTYPE

    n_genomes = len(genomes)
    contig_len = contigs["length"].values.astype(np.int64)
    contig_gid = contigs["genome_id"].values.astype(np.int64)

    contig_offset, genome_length = _ligated_contig_offsets(
        contig_len, contig_gid, n_genomes
    )

    # Drop genomes with no contigs: they would become a zero-length contig,
    # which schema.validate rejects ("contig lengths must be positive").
    keep = np.flatnonzero(genome_length > 0)
    if len(keep) < n_genomes:
        dropped = genomes["name"].values[
            np.setdiff1d(np.arange(n_genomes, dtype=np.int64), keep)
        ]
        logger.warning(
            "Dropping %d genome(s) with no contigs at --level genome: %s",
            n_genomes - len(keep),
            ", ".join(str(name) for name in dropped),
        )
    new_gid = np.full(n_genomes, -1, dtype=np.int64)
    new_gid[keep] = np.arange(len(keep), dtype=np.int64)
    kept_length = genome_length[keep]

    # Uniform tiling of each ligated genome (already in (gid, start) order).
    new_bins = _tile_uniform(kept_length, binsize)
    if new_bins.empty:
        return new_bins, np.array([], dtype=np.int64), contigs, genomes
    n_new_per_genome = -(-kept_length // binsize)
    genome_bin_start = np.concatenate([[0], np.cumsum(n_new_per_genome)]).astype(
        np.int64
    )

    # Map each old bin to the ligated-genome bin containing its midpoint.
    old_chrom = bin_chrom_codes.astype(np.int64)
    old_mid = (
        old_bins["start"].values.astype(np.int64)
        + old_bins["end"].values.astype(np.int64)
    ) // 2
    ligated_pos = contig_offset[old_chrom] + old_mid
    gid_of_old = new_gid[contig_gid[old_chrom]]
    raw = genome_bin_start[gid_of_old] + (ligated_pos // binsize)
    # Clamp inside the genome's own bin range: a stale /contigs/length (never
    # checked by validate) would otherwise attribute pixels to the NEXT
    # genome's first bin.
    upper = genome_bin_start[gid_of_old + 1] - 1
    old_to_new = np.clip(raw, genome_bin_start[gid_of_old], upper)
    n_clamped = int((raw > upper).sum())
    if n_clamped:
        logger.warning(
            "%d source bin(s) fell past the end of their genome (stale contig "
            "lengths?); clamped to the genome's last bin.",
            n_clamped,
        )

    # Output "contigs" are the ligated genomes (one contig per genome).
    reserved = {"name", "length", "n_contigs", "n_bins"}
    out_contigs = pd.DataFrame(
        {
            "name": genomes["name"].values[keep],
            "length": _to_contig_length(kept_length),
            "genome_id": np.arange(len(keep), dtype=CONTIGS_DTYPE["genome_id"]),
        }
    )
    # Extra /genomes columns describe the ligated sequences too; keep them so
    # this branch does not silently drop columns the contig branch preserves --
    # but not the QC ones, the source's per-genome counts: in /contigs they
    # would pass for user metadata that no writer strips.
    for col in genomes.columns:
        if col not in reserved and not is_qc_column("genomes", col):
            out_contigs[col] = genomes[col].values[keep]

    out_genomes = genomes.iloc[keep].copy().reset_index(drop=True)
    out_genomes["length"] = kept_length.astype(GENOMES_DTYPE["length"])
    out_genomes["n_contigs"] = np.ones(len(keep), dtype=GENOMES_DTYPE["n_contigs"])
    return new_bins, old_to_new, out_contigs, out_genomes


def _to_contig_length(lengths: np.ndarray) -> np.ndarray:
    """Cast ligated genome lengths to the ``/contigs/length`` dtype.

    Raises:
        ValueError: If a ligated length overflows the declared dtype (casting
            silently would write a negative contig length).
    """
    from metacooler.schema import CONTIGS_DTYPE

    dtype = np.dtype(CONTIGS_DTYPE["length"])
    limit = int(np.iinfo(str(CONTIGS_DTYPE["length"])).max)
    if lengths.size and int(lengths.max()) > limit:
        raise ValueError(
            f"A ligated genome is {int(lengths.max())} bp, which exceeds the "
            f"maximum contig length representable in the metacooler format "
            f"({limit} bp). Use --level contig for this file."
        )
    return lengths.astype(dtype)


def _write_rebinned_metacool(
    output: Path,
    bins: pd.DataFrame,
    pixels,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    binsize: int,
    ordered: bool = False,
    level: str = "contig",
    source_uri: str | Path = "",
    source_metadata: dict | None = None,
    count_dtype: np.dtype = _DEFAULT_COUNT_DTYPE,
    assembly_name: str | None = None,
    contig_ordering: str = "genome-then-alpha",
) -> None:
    """Write a re-binned metacool file preserving the full metacooler schema.

    Args:
        output: Destination path for the output ``.metacool`` file.
        bins: Bins DataFrame whose ``chrom`` column holds contig ids (row
            numbers of *contigs*), sorted.
        pixels: Iterable of pixel DataFrames, consumed one chunk at a time.
        genomes: Genomes DataFrame.
        contigs: Contigs DataFrame.
        ordered: Whether *pixels* arrive in final order (row blocks of a
            monotonic re-binning), so cooler appends them without its merge.
        binsize: New bin size in base pairs.  Written explicitly to the output's
            ``bin-size``/``bin-type`` attributes: cooler infers ``None`` (hence
            ``bin-type="variable"``) as soon as every contig collapses to a
            single bin, which would then let a later ``binnify`` to a *finer*
            size slip past the resolution guard.
        level: Tiling level used, recorded in the provenance metadata.
        source_uri: Source file URI, recorded in the provenance metadata.
        source_metadata: Metadata of the source file, propagated to the output.
        count_dtype: On-disk dtype of the ``pixels/count`` column.
        assembly_name: Optional assembly name to propagate from source file.
        contig_ordering: Contig ordering string to propagate from source file.
    """
    import h5py

    from metacooler.io._utils import _finalize_metacool_file, write_cooler
    from metacooler.io.staging import recorded_uri
    from metacooler.schema import validate

    # Pre-compute index offsets from the bins DataFrame
    contig_names = contigs["name"].values
    bin_chrom_codes = bins["chrom"].to_numpy(dtype=np.int64)
    chrom_counts = np.bincount(bin_chrom_codes, minlength=len(contig_names))
    chrom_offset = np.concatenate([[0], np.cumsum(chrom_counts)]).astype(np.int64)

    # genome_offset indexes the /contigs table (CSR over contigs), so it must
    # count CONTIGS per genome — not bins. Counting bins would produce offsets
    # running up to n_bins, i.e. out-of-range indices into the contigs table.
    contig_genome_ids = contigs["genome_id"].values
    genome_counts = np.bincount(
        contig_genome_ids.astype(np.int64), minlength=len(genomes)
    )
    genome_offset = np.concatenate([[0], np.cumsum(genome_counts)]).astype(np.int64)

    # Provenance: keep the source file's metadata and append this re-binning.
    # The QC block is the exception — its short/long intra-contig split is a
    # statement about the source resolution, so carrying it onto a coarser grid
    # would attach numbers to a file they were never measured on. Same for the
    # qc_* columns of the tables below, per genome and per contig.
    metadata = strip_qc_metadata(source_metadata)
    history = list(metadata.get("binnify", []))
    history.append(
        {"source": recorded_uri(source_uri), "binsize": int(binsize), "level": level}
    )
    metadata["binnify"] = history
    genomes = strip_qc_columns(genomes)
    contigs = strip_contig_qc_columns(contigs)

    # Pixel chunks arrive sorted and de-duplicated.  Row blocks of a monotonic
    # re-binning are also globally ordered; otherwise the chunk sequence is not,
    # and cooler's merge pass sums duplicate bin pairs across chunks. int64 /
    # float64 counts so a coarse bin pair summing >2**31 source contacts is not
    # saturated by cooler's default int32.
    write_cooler(
        output,
        bins,
        pixels,
        binsize=binsize,
        dtypes={"count": count_dtype},
        ordered=ordered,
        symmetric_upper=True,
        metadata=metadata,
    )

    _finalize_metacool_file(
        output,
        genomes,
        contigs,
        chrom_offset,
        genome_offset,
        contig_ordering=contig_ordering,
        assembly_name=assembly_name,
    )

    # Record the REQUESTED resolution: cooler's get_binsize() returns None when
    # a contig holds a single (truncated) bin, which would write bin-size=null.
    with h5py.File(output, "r+") as f:
        f.attrs["bin-type"] = "fixed"
        f.attrs["bin-size"] = int(binsize)

    # Binned by contig, the contigs are the source's, row for row: so are its
    # clustering rounds. Binned by genome, they are not, and the rounds, which
    # group the source's contigs, are dropped.
    if source_uri:
        from metacooler.io.clusters import copy_clusters
        from metacooler.uri import parse_uri

        source_path, source_group = parse_uri(str(source_uri))
        key = source_group.strip("/")
        with h5py.File(source_path, "r") as fs, h5py.File(output, "r+") as f:
            source = fs[key] if key else fs
            if level == "contig":
                copy_clusters(source, f)
            elif CLUSTERS_GROUP in source:
                logger.warning(
                    "%s holds clustering rounds of its contigs, which %s, "
                    "binned by genome, has not: they are not carried.",
                    source_uri,
                    output,
                )

    validate(output)
    logger.info("Written rebinned metacool to %s", output)
