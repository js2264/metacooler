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

The first two categories do not depend on the contig→genome assignment: a
contact inside one contig stays inside it whichever genome the contig is given
to.  They are therefore also counted per *contig*
(:data:`QC_CONTIG_COLUMNS`, stored in ``/contigs``), which is what lets
:func:`metacooler.io.reassign.reassign` keep a ``source="pairs"`` split after a
re-binning: it sums them per new genome (:meth:`ContactQC.set_intra_contig`)
and re-measures only the two inter-contig categories, which the pixels
classify exactly at any bin size.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any

logger = logging.getLogger("metacooler")

# Version of the ``metadata["qc"]`` block written by this module.  Bump it if
# the meaning of a key changes, so a reader can tell an old block apart from a
# new one rather than silently comparing two different quantities.  Bump it
# too for a new key that adds up over libraries, as the four categories and
# the pairs counts do: ``merge`` sums only the counts it knows of, and carries
# any other key once, as every input records it, so a release that did not
# know the key would keep one library's value where the sum belongs.  It sums
# blocks of its own version only, so the bump has it drop them instead.
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

# Per-contig columns written into ``/contigs``: the two categories that do not
# depend on the genome assignment, under the names their per-genome sums carry
# in ``/genomes``.
QC_CONTIG_COLUMNS = QC_GENOME_COLUMNS[:2]

# The extra keys a ``create pairs`` block records that count something: the
# pairs read and the pairs discarded while binning.  Like the four categories,
# they add up over libraries; every other extra key (``binsize``) says how the
# contacts were counted.
QC_PAIRS_COUNT_KEYS = ("n_pairs_read", "n_pairs_discarded")

# Records classified per pass.  Every temporary of one block -- a handful of
# int64 and bool arrays this long -- stays in a core's L2 cache, where a whole
# 50 M-row chunk at once streams gigabytes of them through main memory: the
# blocking alone halves the tally's time, and caps its scratch memory at a few
# MB instead of ~3 GB per chunk.
_TALLY_BLOCK = 1 << 18


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


def _is_integer(value: object) -> bool:
    """Whether a value read back from a block is an integer (a bool is not)."""
    return isinstance(value, int | np.integer) and not isinstance(value, bool)


def _is_count(value: object) -> bool:
    """Whether a value read back from a block is a finite number, not a bool."""
    if _is_integer(value):
        return True
    return isinstance(value, float | np.floating) and math.isfinite(value)


def _counts_problem(
    block: dict, keys: Sequence[str], *, whole: bool = False
) -> str | None:
    """Why the values of *keys* in *block* are not counts, or ``None``.

    A count is a finite number, not negative; with *whole*, an integer too,
    as a count of pairs (:data:`QC_PAIRS_COUNT_KEYS`) is.
    """
    not_numbers = [key for key in keys if not _is_count(block.get(key))]
    if not_numbers:
        return f"not a finite number: {', '.join(not_numbers)}"
    negative = [key for key in keys if block[key] < 0]
    if negative:
        return f"negative: {', '.join(negative)}"
    fractional = [key for key in keys if whole and not _is_integer(block[key])]
    if fractional:
        return f"not an integer, as a count of pairs is: {', '.join(fractional)}"
    return None


def _sums_match(columns: np.ndarray, totals: Sequence[object]) -> bool:
    """Whether the columns of *columns* sum to *totals*, one total per column.

    Exactly for integers; to a relative ``1e-9`` when a column or a total is a
    float, since floats summed in another order differ in their last bits.
    ``False`` when *columns* holds anything but numbers -- columns read
    straight from a file can hold text -- or a total is not a finite number.
    """
    array = np.asarray(columns)
    if (
        array.ndim != 2
        or array.dtype.kind not in "iuf"
        or array.shape[1] != len(totals)
        or not all(_is_count(total) for total in totals)
    ):
        return False
    sums = array.sum(axis=0)
    if sums.dtype.kind == "f" or not all(_is_integer(total) for total in totals):
        expected = np.asarray(totals, dtype=np.float64)
        return bool(np.allclose(sums, expected, rtol=1e-9, atol=0))
    return sums.tolist() == list(totals)


def _refuse_unlike_counts(
    name: str, array: np.ndarray, totals: list, *, integers: bool, note: str = ""
) -> None:
    """Refuse *array*, stored counts, unless they are a tally's: what it sums to.

    Args:
        name: What *array* is, for the message.
        array: The counts, a 2-D array of numbers.
        totals: What its columns must sum to, from the block.
        integers: Whether the block counts in integers, which its columns
            must then do too.
        note: Said after the totals when they do not match.
    """
    if integers and array.dtype.kind == "f":
        raise ValueError(f"{name} holds floats, where the block counts in integers")
    if array.dtype.kind == "f" and not np.isfinite(array).all():
        raise ValueError(f"{name} holds a count that is not finite")
    if (array < 0).any():
        raise ValueError(f"{name} holds a negative count")
    if not _sums_match(array, totals):
        sums = [_jsonable(value) for value in array.sum(axis=0)]
        raise ValueError(
            f"{name} does not add up to the block: its columns sum to {sums}, "
            f"not {totals}{note}"
        )


def _count_matrix(values: np.ndarray, n_columns: int, name: str) -> np.ndarray:
    """*values* as a 2-D array of numbers with *n_columns* columns, or refuse it."""
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] != n_columns or array.dtype.kind not in "iuf":
        raise ValueError(
            f"{name} must be a 2-D array of numbers with {n_columns} columns, "
            f"got {array.dtype} of shape {array.shape}"
        )
    return array


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
        n_contigs: Number of rows of the ``/contigs`` table, to also count the
            two intra-contig categories per contig (see :meth:`contig_table`);
            ``None`` (the default) counts per genome only.
    """

    def __init__(
        self,
        n_genomes: int,
        distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
        source: Literal["pairs", "pixels"] = "pixels",
        unbinned_genome: int | None = None,
        n_contigs: int | None = None,
    ) -> None:
        if distance_threshold < 0:
            raise ValueError(
                f"distance_threshold must be >= 0, got {distance_threshold}"
            )
        self.n_genomes = int(n_genomes)
        self.distance_threshold = int(distance_threshold)
        self.source = source
        self.unbinned_genome = unbinned_genome
        # Both int64 until `_promote_for` widens them to float64 for fractional
        # weights.
        self.counts: np.ndarray = np.zeros(len(QC_CATEGORIES), dtype=np.int64)
        self.per_genome: np.ndarray = np.zeros(
            (self.n_genomes, len(QC_CATEGORIES)), dtype=np.int64
        )
        # (short, long) per contig, when asked for.  Intra-contig records are
        # queued and summed in batches (see `_queue_intra_contig`).
        self._per_contig: np.ndarray | None = (
            None
            if n_contigs is None
            else np.zeros((int(n_contigs), len(QC_CONTIG_COLUMNS)), dtype=np.int64)
        )
        self._contig_queue: list[tuple[np.ndarray, np.ndarray | None]] = []
        self._contig_queued = 0

    # ------------------------------------------------------------------
    # Accumulation
    # ------------------------------------------------------------------
    def _promote_for(self, weights: np.ndarray | None) -> None:
        """Widen the accumulators to float64 the first time a float weight arrives."""
        if weights is None or weights.dtype.kind != "f":
            return
        if self.counts.dtype.kind != "f":
            # What is queued was weighted by integers: summed into the integer
            # accumulator first, it stays exact.
            self._flush_contigs()
            self.counts = self.counts.astype(np.float64)
            self.per_genome = self.per_genome.astype(np.float64)
            if self._per_contig is not None:
                self._per_contig = self._per_contig.astype(np.float64)

    def _queue_intra_contig(
        self,
        contig: np.ndarray,
        same_contig: np.ndarray,
        long_range: np.ndarray,
        weights: np.ndarray | None,
    ) -> None:
        """Queue one block's intra-contig records for the per-contig tally.

        The per-contig table has two cells per contig -- millions on a
        metagenome -- so a ``bincount`` sized to it for every cache-sized block
        would cost more than the block itself.  The records are queued instead
        and summed once at least as many have gathered as the table has cells
        (:meth:`_flush_contigs`), which bounds that fixed cost by the work it
        does while the queue stays a few times the size of the table.
        """
        rows = np.flatnonzero(same_contig)
        key = np.multiply(contig[rows], len(QC_CONTIG_COLUMNS), dtype=np.int64)
        key += long_range[rows]
        self._contig_queue.append((key, None if weights is None else weights[rows]))
        self._contig_queued += len(key)
        assert self._per_contig is not None
        if self._contig_queued >= max(_TALLY_BLOCK, self._per_contig.size):
            self._flush_contigs()

    def _flush_contigs(self) -> None:
        """Add the queued intra-contig records to the per-contig table."""
        if self._per_contig is None or not self._contig_queue:
            return
        queue, self._contig_queue, self._contig_queued = self._contig_queue, [], 0
        keys = np.concatenate([key for key, _ in queue])
        w = None
        if any(weights is not None for _, weights in queue):
            w = np.concatenate(
                [
                    np.ones(len(key)) if weights is None else weights
                    for key, weights in queue
                ]
            ).astype(np.float64, copy=False)
        # As in `_tally_block`: bincount sums in float64, exact for integers
        # below 2**53, and the cast back to the accumulator dtype is lossless;
        # a contig id never exceeds n_contigs-1, so the slice only trims the
        # zero tail.
        cells = np.bincount(keys, weights=w, minlength=self._per_contig.size)
        self._per_contig += (
            cells[: self._per_contig.size]
            .reshape(self._per_contig.shape)
            .astype(self._per_contig.dtype, copy=False)
        )

    @property
    def per_contig(self) -> np.ndarray | None:
        """``(n_contigs, 2)`` short/long intra-contig counts, or ``None``.

        ``None`` unless the tally was created with *n_contigs*.
        """
        self._flush_contigs()
        return self._per_contig

    def _blocks(self, n: int) -> Iterator[slice]:
        """Cut a chunk of *n* records into the slices :meth:`_tally_block` takes.

        An integer tally is cut into cache-sized blocks (see
        :data:`_TALLY_BLOCK`): integer sums are exact in any order.  Once a
        float weight has turned the accumulators to float64 they are not --
        float addition does not associate -- so each chunk is then classified
        in one pass, which adds every value up in the order an unblocked tally
        always has.  A block is never smaller than the per-genome table it is
        summed into, which bounds the fixed cost of a block by the work it
        does.
        """
        if self.counts.dtype.kind == "f":
            step = max(n, 1)
        else:
            step = max(_TALLY_BLOCK, self.n_genomes * len(QC_CATEGORIES))
        for lo in range(0, n, step):
            yield slice(lo, lo + step)

    def _tally_block(
        self,
        contig1: np.ndarray,
        contig2: np.ndarray,
        genome1: np.ndarray,
        genome2: np.ndarray,
        distance: np.ndarray,
        weights: np.ndarray | None,
    ) -> None:
        """Classify one block of records and add it to the tally.

        Every record gets the index of its category in :data:`QC_CATEGORIES`,
        and one ``bincount`` keyed on ``(genome of the first end, category)``
        sums the whole block -- every category of every genome at once, in a
        single pass, where five masked passes used to select each category's
        records first.  Within each ``(genome, category)`` cell it adds the
        same values in the same order as those did, so even float weights sum
        to the same bits.

        ``np.bincount`` only accepts float64 weights, which is exact for
        integer sums below 2**53 -- six orders of magnitude above the contact
        count of any real library -- so the sums are cast back to the
        accumulator dtype without loss.
        """
        n_categories = len(QC_CATEGORIES)
        same_contig = contig1 == contig2
        # Two genomes, or two different contigs that no binner assigned.
        inter_genome = genome1 != genome2
        if self.unbinned_genome is not None:
            inter_genome |= genome1 == self.unbinned_genome
        inter_genome &= ~same_contig
        # 0/1: intra-contig short/long; 2: inter-contig, same genome;
        # 3: inter-genome.  Bit 1 says "two contigs", bit 0 "two genomes" or,
        # for one contig, "long-range".
        category = (~same_contig).view(np.uint8) << 1
        category |= inter_genome.view(np.uint8)
        long_range = distance >= self.distance_threshold
        long_range &= same_contig
        category |= long_range.view(np.uint8)
        if self._per_contig is not None:
            self._queue_intra_contig(contig1, same_contig, long_range, weights)

        key = np.multiply(genome1, n_categories, dtype=np.int64)
        key += category
        n_cells = self.n_genomes * n_categories
        w = None if weights is None else weights.astype(np.float64, copy=False)
        cells = np.bincount(key, weights=w, minlength=n_cells)
        # An inter-genome contact belongs to two genomes: the cells above
        # credit the first, this credits the second (and the file-level total
        # counts it once).
        second = np.bincount(
            genome2[inter_genome],
            weights=None if w is None else w[inter_genome],
            minlength=self.n_genomes,
        )
        # bincount sizes its output by the largest id seen; a genome id can
        # never exceed n_genomes-1, so this only ever trims the zero tail.
        per_genome = cells[:n_cells].reshape(self.n_genomes, n_categories)

        if weights is not None and weights.dtype.kind == "f":
            # Summed category by category in the weights' own dtype, as the
            # masked passes did, so that a float total keeps its exact bits.
            for k in range(n_categories):
                selected = category == k
                if selected.any():
                    self.counts[k] += weights[selected].sum()
        else:
            self.counts += per_genome.sum(axis=0).astype(np.int64, copy=False)
        dtype = self.per_genome.dtype
        self.per_genome += per_genome.astype(dtype, copy=False)
        self.per_genome[:, n_categories - 1] += second[: self.n_genomes].astype(
            dtype, copy=False
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
        for block in self._blocks(len(contig1)):
            self._tally_block(
                contig1[block],
                contig2[block],
                genome1[block],
                genome2[block],
                distance[block],
                None if weights is None else weights[block],
            )

    def tally_pixels(
        self,
        bin1: np.ndarray,
        bin2: np.ndarray,
        counts: np.ndarray,
        maps: BinMaps,
    ) -> None:
        """Classify a chunk of *pixels* using the per-bin lookup arrays.

        The per-bin lookups are made block by block too (see :meth:`_blocks`),
        so that none of the five arrays they build outgrows the cache.

        Args:
            bin1: ``bin1_id`` column of the chunk.
            bin2: ``bin2_id`` column of the chunk.
            counts: ``count`` column of the chunk.
            maps: Output of :func:`bin_maps` for the same file.
        """
        bin1, bin2, counts = np.asarray(bin1), np.asarray(bin2), np.asarray(counts)
        if len(counts) == 0:
            return
        self._promote_for(counts)
        for block in self._blocks(len(counts)):
            b1 = bin1[block].astype(np.int64, copy=False)
            b2 = bin2[block].astype(np.int64, copy=False)
            self.tally(
                maps.contig[b1],
                maps.contig[b2],
                maps.genome[b1],
                maps.genome[b2],
                np.abs(maps.mid[b2] - maps.mid[b1]),
                counts[block],
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

    def contig_table(self, names: Sequence[str] | np.ndarray) -> pd.DataFrame:
        """Per-contig short/long intra-contig counts, one row per contig.

        Args:
            names: Contig names, in ``/contigs`` order.

        Returns:
            DataFrame with a ``name`` column and the :data:`QC_CONTIG_COLUMNS`
            columns.

        Raises:
            ValueError: If the tally was created without *n_contigs*, or
                *names* does not have one entry per tallied contig.
        """
        per_contig = self.per_contig
        if per_contig is None:
            raise ValueError(
                "this tally does not count per contig; create it with n_contigs"
            )
        if len(names) != len(per_contig):
            raise ValueError(
                f"expected {len(per_contig)} contig name(s), got {len(names)}"
            )
        frame = pd.DataFrame({"name": np.asarray(names)})
        for i, column in enumerate(QC_CONTIG_COLUMNS):
            frame[column] = per_contig[:, i]
        return frame

    def set_intra_contig(
        self, per_contig: np.ndarray, contig_genome: np.ndarray
    ) -> None:
        """Take both intra-contig categories from per-contig counts.

        Neither category depends on the genome assignment, so a tally of the
        pixels under a *new* assignment can take them from the per-contig
        counts of an earlier, finer tally -- the read positions ``create
        pairs`` measured -- and keep only its own two inter-contig categories,
        which a pixel classifies exactly at any bin size.  The two categories'
        totals and per-genome cells are rebuilt from *per_contig*, summed over
        each genome's contigs, and *per_contig* becomes the per-contig table.

        Args:
            per_contig: ``(n_contigs, 2)`` short/long counts, one row per
                contig of the tallied file, in its ``/contigs`` order (the
                :data:`QC_CONTIG_COLUMNS` of ``/contigs``).
            contig_genome: Genome id of every contig, same order.

        Raises:
            ValueError: If the arrays do not describe the same contigs, or a
                genome id is not a row of the tally.
        """
        per_contig = np.asarray(per_contig)
        contig_genome = np.asarray(contig_genome, dtype=np.int64)
        if (
            per_contig.ndim != 2
            or per_contig.shape[1] != len(QC_CONTIG_COLUMNS)
            or len(per_contig) != len(contig_genome)
        ):
            raise ValueError(
                f"expected one (short, long) row per contig for "
                f"{len(contig_genome)} contig(s), got shape {per_contig.shape}"
            )
        if len(contig_genome) and (
            contig_genome.min() < 0 or contig_genome.max() >= self.n_genomes
        ):
            raise ValueError(
                f"genome ids must lie in [0, {self.n_genomes}), got "
                f"[{contig_genome.min()}, {contig_genome.max()}]"
            )
        self._promote_for(per_contig)
        dtype = self.counts.dtype
        per_contig = per_contig.astype(dtype)
        for k in range(len(QC_CONTIG_COLUMNS)):
            # Exact for integers below 2**53, as in `_tally_block`.
            sums = np.bincount(
                contig_genome,
                weights=per_contig[:, k].astype(np.float64),
                minlength=self.n_genomes,
            )
            self.per_genome[:, k] = sums[: self.n_genomes].astype(dtype)
            self.counts[k] = per_contig[:, k].sum()
        self._contig_queue, self._contig_queued = [], 0
        self._per_contig = per_contig

    # ------------------------------------------------------------------
    # Stored tallies
    # ------------------------------------------------------------------
    @classmethod
    def from_stored(
        cls,
        block: dict,
        per_genome: np.ndarray,
        per_contig: np.ndarray | None = None,
        *,
        unbinned_genome: int | None = None,
    ) -> ContactQC:
        """Rebuild the tally a stored QC block and its count columns describe.

        The inverse of storing one: :meth:`summary` of the result gives back
        the block's counts and every key it derives from them,
        :meth:`genome_table` the ``/genomes`` ``qc_*`` columns, and
        :meth:`contig_table` the ``/contigs`` ones.  Keys the block records
        besides the tally (``binsize``, ``n_pairs_read``, ...) are not part of
        it: pass them to :meth:`summary` again.  The counts stay ``int64``
        unless one of them is a float, as a float weight would have made them.

        What is stored is checked first, since nothing else would catch a
        column edited after its block was stored, or a block from another
        map: the counts must be finite and not negative; an integer block's
        columns must hold integers too; and the columns must add up to the
        block.  Per genome, the first three sum to its first three categories
        and the fourth to twice ``inter_genome``, since an inter-genome
        contact counts for both its genomes; per contig, the two sum to the
        two intra-contig categories.  A block measured on pairs may count in
        floats: ``reassign`` keeps one and re-measures its two inter-contig
        categories from the pixels, which a rescaled map holds as floats.

        Args:
            block: A ``metadata["qc"]`` block of version :data:`QC_VERSION`.
            per_genome: ``(n_genomes, 4)`` counts, one row per genome in
                ``/genomes`` order: its ``qc_<category>`` columns, in
                :data:`QC_CATEGORIES` order.
            per_contig: ``(n_contigs, 2)`` short/long counts, one row per
                contig in ``/contigs`` order (its :data:`QC_CONTIG_COLUMNS`),
                or ``None`` for a tally that does not count per contig.
            unbinned_genome: Row of ``__unbinned__`` in ``/genomes`` (see
                :func:`unbinned_genome_id`), or ``None`` if there is none.

        Returns:
            A new tally; *per_genome* and *per_contig* are copied, not shared.

        Raises:
            ValueError: If *block* is not a JSON object of version
                :data:`QC_VERSION` with a ``source`` of ``"pairs"`` or
                ``"pixels"``, a non-negative integer ``distance_threshold``
                and a count for each category, if an array is not shaped and
                typed as above, or if the counts are not those of one tally:
                negative, not finite, or not adding up.
        """
        if not isinstance(block, dict):
            raise ValueError(
                f"a QC block is a JSON object, not a {type(block).__name__}"
            )
        version = block.get("version")
        if version != QC_VERSION:
            raise ValueError(f"a version {version!r} QC block, not {QC_VERSION}")
        source = block.get("source")
        if source not in ("pairs", "pixels"):
            raise ValueError(f"source must be 'pairs' or 'pixels', not {source!r}")
        threshold = block.get("distance_threshold")
        if not isinstance(threshold, int) or isinstance(threshold, bool):
            raise ValueError(
                f"distance_threshold must be an integer, not {threshold!r}"
            )
        problem = _counts_problem(block, QC_CATEGORIES)
        if problem:
            raise ValueError(problem)
        per_genome = _count_matrix(per_genome, len(QC_CATEGORIES), "per_genome")
        if per_contig is not None:
            per_contig = _count_matrix(per_contig, len(QC_CONTIG_COLUMNS), "per_contig")
        short, long_range, inter_contig_cis, inter_genome = (
            block[key] for key in QC_CATEGORIES
        )
        integers = all(_is_integer(block[key]) for key in QC_CATEGORIES)
        _refuse_unlike_counts(
            "per_genome",
            per_genome,
            [short, long_range, inter_contig_cis, 2 * inter_genome],
            integers=integers,
            note=" (the last is twice inter_genome: an inter-genome contact "
            "counts for both its genomes)",
        )
        if per_contig is not None:
            _refuse_unlike_counts(
                "per_contig", per_contig, [short, long_range], integers=integers
            )

        tally = cls(
            len(per_genome),
            threshold,
            source=source,
            unbinned_genome=unbinned_genome,
            n_contigs=None if per_contig is None else len(per_contig),
        )
        counts = np.asarray([block[key] for key in QC_CATEGORIES])
        for array in (counts, per_genome, per_contig):
            tally._promote_for(array)
        dtype = tally.counts.dtype
        tally.counts = counts.astype(dtype)
        tally.per_genome = per_genome.astype(dtype)
        if per_contig is not None:
            tally._per_contig = per_contig.astype(dtype)
        return tally

    def _measurement(self) -> str:
        """How the contacts were counted, as :meth:`add` compares tallies."""
        return f"source={self.source!r}, distance_threshold={self.distance_threshold}"

    def _layout(self) -> str:
        """The genomes and contigs counted, as :meth:`add` compares tallies."""
        unbinned = (
            "no __unbinned__"
            if self.unbinned_genome is None
            else f"__unbinned__ at row {self.unbinned_genome}"
        )
        per_contig = self.per_contig
        contigs = (
            "no per-contig counts"
            if per_contig is None
            else f"{len(per_contig)} contig(s)"
        )
        return f"{self.n_genomes} genome(s) ({unbinned}), {contigs}"

    def add(self, other: ContactQC) -> None:
        """Add another tally of the same genomes and contigs to this one.

        For another library of the same map: every cell adds up -- the four
        totals, every genome's counts and, when both count per contig, every
        contig's -- so the sum of the libraries' tallies is the tally of all
        their contacts at once, exactly for integer counts.  Float counts add
        up to the same value, though not always to the same bits: float
        addition does not associate.  An integer tally becomes a float one
        when *other* is float, as it would on a float weight.

        Args:
            other: The tally to add, left as it was.

        Raises:
            ValueError: If *other* was counted differently (another ``source``
                or ``distance_threshold``), or counts other genomes or contigs
                (another number of either, another ``__unbinned__`` row, or
                per contig where this tally does not, or the reverse).
        """
        for describe in (ContactQC._measurement, ContactQC._layout):
            if describe(other) != describe(self):
                raise ValueError(
                    f"cannot add a tally of {describe(other)} to one of "
                    f"{describe(self)}"
                )
        self._promote_for(other.counts)
        dtype = self.counts.dtype
        self.counts += other.counts.astype(dtype, copy=False)
        self.per_genome += other.per_genome.astype(dtype, copy=False)
        theirs = other.per_contig
        if self._per_contig is not None and theirs is not None:
            self._per_contig += theirs.astype(dtype, copy=False)


# Every key ContactQC.summary() writes itself: the tally's parameters (version,
# source, distance_threshold), its four counts and what it derives from them.
# Any other key of a stored block is an ``**extra`` it was given, recording how
# the counts were measured (binsize, n_pairs_read, n_pairs_discarded, ...).
QC_SUMMARY_KEYS = frozenset(ContactQC(0).summary())


def _contig_counts_add_up(block: dict, per_contig: np.ndarray) -> bool:
    """Whether per-contig QC counts add up to a block's intra-contig totals.

    They do when the block and the ``/contigs`` columns were stored together.
    A ``qc --store`` from a release without per-contig columns (0.2.0 and
    earlier) rewrites the block and leaves older columns behind, which no
    longer describe it.

    Args:
        block: A ``metadata["qc"]`` block.
        per_contig: ``(n_contigs, 2)`` short/long counts: the
            :data:`QC_CONTIG_COLUMNS` of ``/contigs``.

    Returns:
        Whether the two column totals are the block's ``intra_contig_short``
        and ``intra_contig_long``: exactly for integer counts, to a relative
        ``1e-9`` for float ones.  ``False`` when a total is missing or not a
        number, or when the columns are not numbers.
    """
    totals = [block.get("intra_contig_short"), block.get("intra_contig_long")]
    return _sums_match(per_contig, totals)


def is_qc_column(table: str, column: object) -> bool:
    """Whether a column of ``/genomes`` or ``/contigs`` belongs to the QC block.

    In ``/genomes``, every ``qc_*`` column does: metacooler writes those names
    there, and a user's own are stripped with them (:func:`strip_qc_columns`).
    In ``/contigs``, only the two per-contig columns, :data:`QC_CONTIG_COLUMNS`,
    do: any other ``qc_*`` column there (``qc_pass``, from a ``--contigs``
    table) is a user's contig metadata, which every writer keeps.

    Args:
        table: ``"genomes"`` or ``"contigs"``.
        column: The column name.

    Raises:
        ValueError: If *table* is neither.
    """
    if table == "contigs":
        return column in QC_CONTIG_COLUMNS
    if table == "genomes":
        return str(column).startswith(QC_COLUMN_PREFIX)
    raise ValueError(f"table must be 'genomes' or 'contigs', got {table!r}")


def strip_qc_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop every ``qc_*`` column from a genomes table.

    Called by the writers that derive a new file from an existing one
    (``binnify``, ``reassign``).  Their output has a different bin grid or
    different genomes, so counts copied across would describe contacts the
    output does not contain — a stale number is worse than no number, because
    nothing about it looks stale.  (``reassign`` then stores a block for its
    own output, carrying the per-contig :data:`QC_CONTIG_COLUMNS` over when it
    can.  ``merge`` drops each input's columns too, and stores their sum: see
    :meth:`ContactQC.add`.)

    For ``/contigs``, see :func:`strip_contig_qc_columns`.

    Args:
        frame: Genomes table about to be written to a derived file.

    Returns:
        *frame* unchanged when it carries no QC column, else a narrower copy.
    """
    drop = [c for c in frame.columns if is_qc_column("genomes", c)]
    return frame.drop(columns=drop) if drop else frame


def strip_contig_qc_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop the per-contig QC columns, :data:`QC_CONTIG_COLUMNS`, from a contigs table.

    The ``/contigs`` counterpart of :func:`strip_qc_columns`, and narrower on
    purpose: these two are the only ``qc_*`` columns metacooler writes to
    ``/contigs``, while a user's ``--contigs`` table may bring a ``qc_*`` column
    of its own (``qc_pass``) -- contig metadata, which a derived file keeps.

    Args:
        frame: Contigs table about to be written to a derived file.

    Returns:
        *frame* unchanged when it has neither column, else a narrower copy.
    """
    drop = [c for c in frame.columns if is_qc_column("contigs", c)]
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
