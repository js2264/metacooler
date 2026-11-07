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

    genome_names = mc.genomes["name"].tolist()
    contig_names = mc.contigs["name"].tolist()

    # A coordinate suffix only makes sense on a contig, so when the label is
    # both (``binnify --level genome`` names each ligated contig after its
    # genome) the suffix decides: "G_a:1-6000" is a contig range, which used to
    # be refused because the genome branch was tried first.
    if ":" in range_str and label in contig_names:
        return mc._clr.extent(range_str)

    if label in genome_names:
        if ":" in range_str:
            raise ValueError(
                f"Range {range_str!r}: {label!r} is a genome name — coordinate "
                f"suffixes are not allowed for genome ranges. Use just {label!r}."
            )
        return mc.genome_extent(label)

    if label in contig_names:
        # Delegate to cooler — handles both "contig" and "contig:start-end"
        return mc._clr.extent(range_str)

    def _truncate(names: list, limit: int = 20) -> str:
        if len(names) <= limit:
            return str(names)
        return str(names[:limit]) + f" ... ({len(names) - limit} more)"

    raise ValueError(
        f"Range {range_str!r}: {label!r} is not a known genome or contig.\n"
        f"  Genomes: {_truncate(genome_names)}\n"
        f"  Contigs: {_truncate(contig_names)}"
    )
