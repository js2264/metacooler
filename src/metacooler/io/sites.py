"""Restriction-site counts per contig, read from the assembly FASTA.

NormCC and HiCzin correct for how many ligation opportunities a contig offers:
its number of restriction sites.  A contig table can carry those counts
already; otherwise they are counted here from the FASTA and the enzymes of the
library.

The sites of every enzyme, ``N`` read as any base, are joined into one
regular expression and its non-overlapping matches counted along each
sequence, so a position cut by two enzymes counts once.  (MetaCC and HiCBin
count with Biopython, enzyme by enzyme, and so count such a position once per
enzyme; the two agree wherever no two of the enzymes' sites overlap.)  A site
that is not its own reverse complement is searched on both strands; every
common Hi-C enzyme's site is palindromic, where this changes nothing.
Sequences are upper-cased first, so soft-masked bases still count.
"""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np

logger = logging.getLogger("metacooler")

#: Recognition sites of the restriction enzymes used for Hi-C, by the names NEB
#: and Biopython give them.  ``Arima`` is the DpnII + HinfI cocktail of the
#: Arima kits, as hicstuff names it.
ENZYME_SITES: dict[str, tuple[str, ...]] = {
    "AluI": ("AGCT",),
    "Arima": ("GATC", "GANTC"),
    "BamHI": ("GGATCC",),
    "BglII": ("AGATCT",),
    "Csp6I": ("GTAC",),
    "CviQI": ("GTAC",),
    "DdeI": ("CTNAG",),
    "DpnI": ("GATC",),
    "DpnII": ("GATC",),
    "EcoRI": ("GAATTC",),
    "HaeIII": ("GGCC",),
    "HindIII": ("AAGCTT",),
    "HinfI": ("GANTC",),
    "HpaII": ("CCGG",),
    "MboI": ("GATC",),
    "MluCI": ("AATT",),
    "MseI": ("TTAA",),
    "MspI": ("CCGG",),
    "NcoI": ("CCATGG",),
    "NlaIII": ("CATG",),
    "RsaI": ("GTAC",),
    "Sau3AI": ("GATC",),
}

#: IUPAC nucleotide codes, as regular-expression classes.
_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[CG]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": ".",
}

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _reverse_complement(site: str) -> str:
    return site.translate(_COMPLEMENT)[::-1]


def enzyme_sites(enzymes: Sequence[str]) -> list[str]:
    """Return the recognition sites of *enzymes*.

    Each entry is an enzyme name (case-insensitive, see :data:`ENZYME_SITES`)
    or a recognition site written in IUPAC codes (``GATC``, ``GANTC``).

    Raises:
        ValueError: If an entry is neither, or no entry is given.
    """
    by_name = {name.lower(): sites for name, sites in ENZYME_SITES.items()}
    sites: list[str] = []
    for raw in enzymes:
        entry = raw.strip()
        if not entry:
            continue
        if entry.lower() in by_name:
            sites.extend(by_name[entry.lower()])
        elif set(entry.upper()) <= set(_IUPAC):
            sites.append(entry.upper())
        else:
            known = ", ".join(sorted(ENZYME_SITES))
            raise ValueError(
                f"Unknown enzyme {entry!r}: give an enzyme name ({known}) or its "
                "recognition site in IUPAC letters (e.g. GATC)."
            )
    if not sites:
        raise ValueError("No enzyme given to count restriction sites with.")
    return sites


def site_pattern(sites: Sequence[str]) -> re.Pattern[str]:
    """One regular expression matching any of *sites*, on either strand."""
    variants = {s for site in sites for s in (site, _reverse_complement(site))}
    parts = ["".join(_IUPAC[base] for base in site) for site in variants]
    return re.compile("|".join(sorted(set(parts))))


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(name, sequence)`` for each record of a FASTA file.

    The name is the header up to its first whitespace.  A ``.gz`` file is read
    through gzip.  One record is held in memory at a time.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    name: str | None = None
    chunks: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                fields = line[1:].split()
                name = fields[0] if fields else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def count_sites(
    fasta: str | Path,
    contig_names: Sequence[str],
    contig_lengths: Sequence[int] | np.ndarray,
    enzymes: Sequence[str],
) -> np.ndarray:
    """Count the restriction sites of every contig of a map.

    Args:
        fasta: The assembly the map was built from (plain or ``.gz``).
        contig_names: The map's contigs, in ``/contigs`` order.
        contig_lengths: Their lengths, checked against the FASTA so a map is
            never counted against the wrong assembly.
        enzymes: Enzyme names or IUPAC sites (see :func:`enzyme_sites`).

    Returns:
        One count per contig, in *contig_names* order.

    Raises:
        ValueError: If a contig is missing from the FASTA or has another
            length there.
    """
    pattern = site_pattern(enzyme_sites(enzymes))
    row = {str(name): i for i, name in enumerate(contig_names)}
    lengths = np.asarray(contig_lengths, dtype=np.int64)
    counts = np.full(len(row), -1, dtype=np.int64)
    mismatched: list[str] = []
    for name, sequence in iter_fasta(fasta):
        i = row.get(name)
        if i is None:
            continue
        if len(sequence) != lengths[i]:
            mismatched.append(
                f"{name} ({len(sequence)} bp in the FASTA, {lengths[i]} in the map)"
            )
            continue
        counts[i] = len(pattern.findall(sequence.upper()))
    if mismatched:
        raise ValueError(
            f"{len(mismatched)} contig(s) have another length in {fasta}, so it is "
            f"not the assembly of this map: {', '.join(mismatched[:5])}"
            + (" ..." if len(mismatched) > 5 else "")
        )
    missing = [name for name, i in row.items() if counts[i] < 0]
    if missing:
        raise ValueError(
            f"{len(missing)} contig(s) of the map are not in {fasta}: "
            f"{', '.join(missing[:5])}" + (" ..." if len(missing) > 5 else "")
        )
    logger.info(
        "Counted %s restriction site(s) of %s on %d contig(s)",
        f"{int(counts.sum()):,}",
        ", ".join(enzymes),
        len(counts),
    )
    return counts
