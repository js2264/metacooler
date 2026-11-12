"""Contact-quality statistics: the shared tally engine.

A metagenomic Hi-C library is judged before anything else on *where* its
contacts land, and the four places they can land are the four categories
below.  This module owns that classification and nothing else: it is pure
numpy, it never opens a file, and both producers of the statistics go through
it so they cannot drift apart.

- :func:`metacooler.io.create.from_pairs` tallies while it is already
  streaming the pairs file, from the raw base-pair positions
  (``source="pairs"``).
- :func:`metacooler.ops.qc.qc` tallies from the pixel table of a file that
  already exists (``source="pixels"``), which is the only option for a map
  imported with ``create cool`` — at the cost of measuring intra-contig
  distances between *bin midpoints* rather than between read positions.

The categories partition every contact exactly once:

``intra_contig_short``
    Both ends on the same contig, less than ``distance_threshold`` apart.
``intra_contig_long``
    Both ends on the same contig, at least ``distance_threshold`` apart.
``inter_contig_intra_genome``
    Two different contigs of the same genome (MAG).
``inter_genome``
    Two contigs assigned to different genomes.

``__unbinned__`` is not a genome: it collects every contig that no binner
assigned, so two of its contigs share nothing but the absence of an
assignment.  A contact between two *different* unbinned contigs is therefore
``inter_genome`` -- each unbinned contig counts as a genome of its own.  Treating
the bucket as one genome used to inflate ``intra_genome_fraction`` by every
contact among the unassigned contigs, which on a real assembly are often most
of them.

The headline number, ``trans_contig_ratio`` (the "3D ratio"), is the share of
contacts joining two *different* contigs — the signal that assembles a
metagenome — and is by construction
``(inter_contig_intra_genome + inter_genome) / n_contacts``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

logger = logging.getLogger("metacooler")

# Version of the ``metadata["qc"]`` block written by this module.  Bump it if
# the meaning of a key changes, so a reader can tell an old block apart from a
# new one rather than silently comparing two different quantities.
#   2: contacts between two different __unbinned__ contigs moved from
#      inter_contig_intra_genome to inter_genome.
QC_VERSION = 2

# Name of the genome that collects unassigned contigs (see module docstring).
UNBINNED_GENOME = "__unbinned__"

# Key the block is stored under inside a file's ``metadata`` attribute.
QC_METADATA_KEY = "qc"

# Default split between "short-range" and "long-range" intra-contig signal.
# 10 kb is the convention in the metaHi-C literature: below it the signal is
# dominated by the polymer's local contacts (and by undigested/religated
# fragments), above it by genuine long-range 3D proximity.
DEFAULT_DISTANCE_THRESHOLD = 10_000

# The four mutually exclusive categories, in report order.
QC_CATEGORIES = (
    "intra_contig_short",
    "intra_contig_long",
    "inter_contig_intra_genome",
    "inter_genome",
)

# Prefix of every per-genome column this module writes into ``/genomes``.
# Anything a derived file inherits under this prefix describes the *source*
# file's pixels and is stripped rather than carried over (see
# :func:`strip_qc_columns`).
QC_COLUMN_PREFIX = "qc_"

# Per-genome columns written into ``/genomes``, in order.
QC_GENOME_COLUMNS = tuple(f"{QC_COLUMN_PREFIX}{c}" for c in QC_CATEGORIES) + (
    f"{QC_COLUMN_PREFIX}intra_genome_ratio",
)


def _jsonable(value: Any) -> Any:
    """Return *value* as a plain Python scalar where it is a NumPy one.

    The block is written straight into a JSON attribute, and ``json.dumps``
    refuses a ``numpy.int64`` — which is what ``cooler``'s ``binsize`` and
    every count derived from an array actually are.  Left unconverted, they
    reach the file as the *string* ``"500"`` (via a ``default=str`` fallback)
    or not at all.
    """
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


class BinMaps(NamedTuple):
    """Per-bin lookup arrays needed to classify a pixel.

    Attributes:
        contig: Contig id of every bin.
        genome: Genome id of every bin.
        mid: Midpoint coordinate of every bin, in bp.  Non-decreasing within
            a contig, which is what makes an intra-contig distance a plain
            subtraction.
    """

    contig: np.ndarray
    genome: np.ndarray
    mid: np.ndarray


def bin_maps(bins: pd.DataFrame, contigs: pd.DataFrame) -> BinMaps:
    """Build the per-bin contig / genome / midpoint arrays for a pixel tally.

    Args:
        bins: Bins table with ``chrom`` (contig id, cooler categorical or
            plain integer), ``start`` and ``end``.
        contigs: ``/contigs`` table, used for its ``genome_id`` column.

    Returns:
        A :class:`BinMaps` triple, each array of length ``len(bins)``.
    """
    from metacooler.io._utils import chrom_codes

    contig_of_bin = chrom_codes(bins).astype(np.int64, copy=False)
    genome_of_bin = (
        contigs["genome_id"].to_numpy(dtype=np.int64)[contig_of_bin]
        if len(contig_of_bin)
        else np.zeros(0, dtype=np.int64)
    )
    mid = (
        bins["start"].to_numpy(dtype=np.int64) + bins["end"].to_numpy(dtype=np.int64)
    ) // 2
    return BinMaps(contig_of_bin, genome_of_bin, mid)


def unbinned_genome_id(names: Sequence[str] | np.ndarray) -> int | None:
    """Return the row of ``__unbinned__`` in a ``/genomes`` name list, if any."""
    hits = np.flatnonzero(np.asarray(names).astype(str) == UNBINNED_GENOME)
    return int(hits[0]) if hits.size else None


class ContactQC:
    """Running tally of contacts per category, in total and per genome.

    Counts are accumulated as ``int64`` and promoted to ``float64`` only if a
    float weight column (e.g. balanced values) is ever handed in, so an
    ordinary integer ``count`` column is summed exactly.

    Args:
        n_genomes: Number of rows of the ``/genomes`` table.
        distance_threshold: Intra-contig short/long split, in bp.
        source: ``"pairs"`` when the tally comes from read positions,
            ``"pixels"`` when it comes from binned pixels.
        unbinned_genome: Row of ``__unbinned__`` in ``/genomes`` (see
            :func:`unbinned_genome_id`), whose contigs are each classified as
            a genome of their own; ``None`` if the file has none.
    """

    def __init__(
        self,
        n_genomes: int,
        distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
        source: Literal["pairs", "pixels"] = "pixels",
        unbinned_genome: int | None = None,
    ) -> None:
        if distance_threshold < 0:
            raise ValueError(
                f"distance_threshold must be >= 0, got {distance_threshold}"
            )
        self.n_genomes = int(n_genomes)
        self.distance_threshold = int(distance_threshold)
        self.source = source
        self.unbinned_genome = unbinned_genome
        self.counts = np.zeros(len(QC_CATEGORIES), dtype=np.int64)
        # int64 until `_promote_for` widens it to float64 for fractional weights.
        self.per_genome: np.ndarray = np.zeros(
            (self.n_genomes, len(QC_CATEGORIES)), dtype=np.int64
        )

    # ------------------------------------------------------------------
    # Accumulation
    # ------------------------------------------------------------------
    def _promote_for(self, weights: np.ndarray | None) -> None:
        """Widen the accumulators to float64 the first time a float weight arrives."""
        if weights is None or weights.dtype.kind != "f":
            return
        if self.counts.dtype.kind != "f":
            self.counts = self.counts.astype(np.float64)
            self.per_genome = self.per_genome.astype(np.float64)

    def _add(
        self,
        category: int,
        mask: np.ndarray,
        genome_ids: np.ndarray,
        weights: np.ndarray | None,
        *,
        add_total: bool = True,
    ) -> None:
        """Add the masked records to one category, in total and per genome.

        ``np.bincount`` only accepts float64 weights, which is exact for
        integer sums below 2**53 — six orders of magnitude above the contact
        count of any real library — so the per-genome totals are cast back to
        the accumulator dtype without loss.

        Args:
            category: Index into :data:`QC_CATEGORIES`.
            mask: Which records of the current chunk fall in this category.
            genome_ids: Genome id of the *side* of each record being credited.
            weights: Contact count per record, or ``None`` for one each.
            add_total: Whether to also add to the file-level total.  ``False``
                for the second side of an inter-genome contact, which is
                credited to both genomes but counted once overall.
        """
        if not mask.any():
            return
        if weights is None:
            if add_total:
                self.counts[category] += np.int64(mask.sum())
            hits = np.bincount(genome_ids[mask], minlength=self.n_genomes)
        else:
            w = weights[mask]
            if add_total:
                self.counts[category] += w.sum()
            hits = np.bincount(
                genome_ids[mask],
                weights=w.astype(np.float64, copy=False),
                minlength=self.n_genomes,
            )
        # bincount sizes its output by the largest id seen; a genome id can
        # never exceed n_genomes-1, so this only ever trims the zero tail.
        self.per_genome[:, category] += hits[: self.n_genomes].astype(
            self.per_genome.dtype
        )

    def tally(
        self,
        contig1: np.ndarray,
        contig2: np.ndarray,
        genome1: np.ndarray,
        genome2: np.ndarray,
        distance: np.ndarray,
        weights: np.ndarray | None = None,
    ) -> None:
        """Classify a chunk of contacts and add them to the tally.

        Args:
            contig1: Contig id of the first end of every contact.
            contig2: Contig id of the second end.
            genome1: Genome id of the first end.
            genome2: Genome id of the second end.
            distance: Separation in bp.  Only read where the two ends share a
                contig, so callers may pass anything (typically the raw
                difference) elsewhere.
            weights: Contact count per record — the pixel ``count`` column.
                ``None`` means one contact per record, which is what the pairs
                path hands in.
        """
        if len(contig1) == 0:
            return
        self._promote_for(weights)

        same_contig = contig1 == contig2
        short = same_contig & (distance < self.distance_threshold)
        long_range = same_contig & ~short
        same_genome = genome1 == genome2
        if self.unbinned_genome is not None:
            # Two different unassigned contigs share no genome.
            same_genome &= same_contig | (genome1 != self.unbinned_genome)
        inter_contig_cis = ~same_contig & same_genome
        inter_genome = ~same_genome

        # The first three categories are intra-genome, so both ends carry the
        # same genome id and one credit is the whole story.
        self._add(0, short, genome1, weights)
        self._add(1, long_range, genome1, weights)
        self._add(2, inter_contig_cis, genome1, weights)
        # An inter-genome contact belongs to two genomes: credit both, but
        # count it once in the file-level total.
        self._add(3, inter_genome, genome1, weights)
        self._add(3, inter_genome, genome2, weights, add_total=False)

    def tally_pixels(
        self,
        bin1: np.ndarray,
        bin2: np.ndarray,
        counts: np.ndarray,
        maps: BinMaps,
    ) -> None:
        """Classify a chunk of *pixels* using the per-bin lookup arrays.

        Args:
            bin1: ``bin1_id`` column of the chunk.
            bin2: ``bin2_id`` column of the chunk.
            counts: ``count`` column of the chunk.
            maps: Output of :func:`bin_maps` for the same file.
        """
        b1 = np.asarray(bin1, dtype=np.int64)
        b2 = np.asarray(bin2, dtype=np.int64)
        self.tally(
            maps.contig[b1],
            maps.contig[b2],
            maps.genome[b1],
            maps.genome[b2],
            np.abs(maps.mid[b2] - maps.mid[b1]),
            np.asarray(counts),
        )

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    @property
    def n_contacts(self) -> int | float:
        """Total number of classified contacts."""
        return self.counts.sum().item()

    def summary(self, **extra: Any) -> dict:
        """Assemble the ``metadata["qc"]`` block.

        Args:
            **extra: Additional keys to record verbatim — ``binsize`` for a
                pixel tally, ``n_pairs_read`` / ``n_pairs_discarded`` for a
                pairs tally.  ``None`` values are dropped.

        Returns:
            A flat, JSON-serialisable dict.  Every ``*_fraction`` and the
            ``trans_contig_ratio`` are ``None`` on an empty file rather than
            ``0.0``, which would claim a measurement that was never made.
        """
        short, long_range, inter_contig_cis, inter_genome = (
            v.item() for v in self.counts
        )
        intra_contig = short + long_range
        inter_contig = inter_contig_cis + inter_genome
        intra_genome = intra_contig + inter_contig_cis
        total = intra_contig + inter_contig

        def frac(value: float) -> float | None:
            return round(value / total, 6) if total else None

        block: dict[str, Any] = {
            "version": QC_VERSION,
            "source": self.source,
            "distance_threshold": self.distance_threshold,
            "n_contacts": total,
            "intra_contig_short": short,
            "intra_contig_long": long_range,
            "inter_contig_intra_genome": inter_contig_cis,
            "inter_genome": inter_genome,
            "intra_contig": intra_contig,
            "inter_contig": inter_contig,
            "intra_genome": intra_genome,
            "trans_contig_ratio": frac(inter_contig),
            "intra_contig_short_fraction": frac(short),
            "intra_contig_long_fraction": frac(long_range),
            "inter_contig_intra_genome_fraction": frac(inter_contig_cis),
            "intra_genome_fraction": frac(intra_genome),
            "inter_genome_fraction": frac(inter_genome),
        }
        block.update(
            {k: _jsonable(v) for k, v in extra.items() if v is not None},
        )
        return block

    def genome_table(self, names: Sequence[str] | np.ndarray) -> pd.DataFrame:
        """Per-genome counts plus the intra-genome ratio, one row per genome.

        The ratio answers "how self-contained is this MAG?": the share of a
        genome's contacts that stay inside it.  A bin whose contigs really
        belong together sits near 1; a chimeric or badly split one leaks
        contacts to its neighbours and drops.

        Args:
            names: Genome names, in ``/genomes`` order.

        Returns:
            DataFrame with a ``name`` column and the
            :data:`QC_GENOME_COLUMNS` columns.

        Raises:
            ValueError: If *names* does not have one entry per tallied genome.
        """
        if len(names) != self.n_genomes:
            raise ValueError(
                f"expected {self.n_genomes} genome name(s), got {len(names)}"
            )
        per_genome = self.per_genome
        totals = per_genome.sum(axis=1)
        intra = per_genome[:, 0] + per_genome[:, 1] + per_genome[:, 2]
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(totals > 0, intra / np.maximum(totals, 1), np.nan)
        frame = pd.DataFrame({"name": np.asarray(names)})
        for i, column in enumerate(QC_GENOME_COLUMNS[: len(QC_CATEGORIES)]):
            frame[column] = per_genome[:, i]
        frame[QC_GENOME_COLUMNS[-1]] = ratio
        return frame


def strip_qc_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop every ``qc_*`` column from a genomes/contigs table.

    Called by the writers that derive a new file from an existing one
    (``binnify``, ``merge``).  Their output has a different bin grid or a
    different pixel table, so counts copied across would describe contacts
    the output does not contain — a stale number is worse than no number,
    because nothing about it looks stale.

    Args:
        frame: Table about to be written to a derived file.

    Returns:
        *frame* unchanged when it carries no QC column, else a narrower copy.
    """
    drop = [c for c in frame.columns if str(c).startswith(QC_COLUMN_PREFIX)]
    return frame.drop(columns=drop) if drop else frame


def strip_qc_metadata(metadata: dict | None) -> dict:
    """Drop the ``qc`` block from metadata being propagated to a derived file.

    Companion to :func:`strip_qc_columns`; see there for why.

    Args:
        metadata: Source file's metadata, or ``None``.

    Returns:
        A copy without the ``qc`` key.
    """
    meta = dict(metadata or {})
    meta.pop(QC_METADATA_KEY, None)
    return meta
