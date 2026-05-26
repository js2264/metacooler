"""Range resolution for metacooler region selectors.

Provides :func:`parse_metacooler_range`, which maps a range string to a
half-open ``(bin_i0, bin_i1)`` interval.  The range string may refer to:

* a **genome name** — expands to all bins across every contig in that genome;
* a **contig name** — expands to all bins in that contig;
* a **contig range** (``contig:start-end``) — the bins that overlap the
  given base-pair interval (delegated to ``cooler.Cooler.extent``).
"""

from __future__ import annotations

from metacooler.api import MetaCooler


def parse_metacooler_range(mc: MetaCooler, range_str: str) -> tuple[int, int]:
    """Resolve a range string to a ``(bin_i0, bin_i1)`` half-open interval.

    Resolution order:
    1. If the name (part before ``:``) is a **genome name**, return
       :py:meth:`MetaCooler.genome_extent`.  Any ``:start-end`` suffix is
       ignored — genome ranges are always whole-genome.
    2. Otherwise treat the string as a **contig** (optionally with
       ``contig:start-end`` coordinates) and delegate to
       ``mc._clr.extent()``, which reuses cooler's ``parse_region`` logic.

    Args:
        mc: An open :class:`~metacooler.api.MetaCooler` instance.
        range_str: Range specifier, e.g. ``"MAG_1"``, ``"contig_A"``, or
            ``"contig_A:100-500"``.

    Returns:
        ``(bin_i0, bin_i1)`` half-open bin index interval.

    Raises:
        ValueError: If the label (part before ``:``) is neither a known
            genome nor a known contig.
    """
    label = range_str.split(":")[0]

    if label in mc.genomes["name"].values:
        return mc.genome_extent(label)

    if label in mc.contigs["name"].values:
        # Delegate to cooler — handles both "contig" and "contig:start-end"
        return mc._clr.extent(range_str)

    raise ValueError(
        f"Range {range_str!r}: {label!r} is not a known genome or contig.\n"
        f"  Genomes: {mc.genomes['name'].tolist()}\n"
        f"  Contigs: {mc.contigs['name'].tolist()}"
    )
