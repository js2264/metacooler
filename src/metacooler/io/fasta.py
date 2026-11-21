"""Write the sequences of a map's genomes, or of a round's clusters, as FASTA.

``metacooler fasta`` (:func:`write_fasta`) reads the assembly a map was made
from and writes, for each genome of ``/genomes`` -- or each cluster of one
clustering round of ``/clusters``, or each group of contigs a caller names --
the records of its contigs:

- **as they are**, by default: each contig's record copied byte for byte, its
  header line and its line breaks included, exactly as ``pyfastx`` gives it
  (``raw``), so a tool that wrote its FASTA files record by record gets the
  same bytes from here. A genome's or a round's contigs come in map order (the
  rows of ``/contigs``), a caller's group's in the order given;
- **stitched** (``stitch=True``): one record per genome or cluster, named
  after it, its contigs' bases joined by ``spacer`` N and wrapped at
  :data:`LINE_WIDTH` columns. A genome's or a round's contigs come longest
  first, ties in map order -- a cluster may span genomes, whose order in the
  map follows the current assignment -- and a caller's group's in the order
  given.

One file per genome or cluster (``<name>.fa`` in a directory), or every record
in one file (``single_file=True``), where a contig record's header line is
prefixed with its genome or cluster (``<genome>|<header>``). The contigs in no
genome, or in no cluster of the round, go to ``__unbinned__`` with
``unbinned=True``: as they are, in map order, even when stitching.

A record is named, and measured, as the map's creator reads it: its name is
the header up to the first whitespace, past any space after the ``>``
(``> c1 desc`` is ``c1``, as for ``create`` and metaTOR), and a contig is the
map's when the map's length counts either its sequence's bytes but line
breaks (as pyfaidx does for ``create --fasta``) or its bases alone.

The assembly is read once, to find its records, and each record is then copied
from it: a plain file is mapped in memory; a gzip-compressed one (gzip or
bgzip), or one that is no regular file (a pipe, ``<(zcat ...)``), is first
copied, plain, into a private temporary directory. Nothing is written next to
it.
"""

from __future__ import annotations

import contextlib
import gzip
import logging
import mmap
import os
import re
import shutil
import stat
import tempfile
import zlib
from collections import Counter
from collections.abc import Callable, Generator, Mapping, Sequence
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING, BinaryIO, NoReturn

import numpy as np
import pandas as pd

from metacooler.io._utils import (
    check_output_path,
    clusters_filter_frame,
    contigs_filter_frame,
    genomes_filter_frame,
)
from metacooler.io.staging import _umask
from metacooler.qc import UNBINNED_GENOME, unbinned_genome_id

if TYPE_CHECKING:
    from metacooler.api import MetaCooler
    from metacooler.filters import Filter

logger = logging.getLogger("metacooler")

#: Width of the sequence lines of a stitched record.
LINE_WIDTH = 80
#: Number of N between two contigs of a stitched record, by default.
DEFAULT_SPACER = 100
#: Between the genome (or cluster) and a contig's header, with ``single_file``.
SINGLE_FILE_SEPARATOR = "|"
#: The extension of the file of each genome or cluster.
FASTA_SUFFIX = ".fa"

_GZIP_MAGIC = b"\x1f\x8b"
# What a sequence line may hold besides its bases.
_WHITESPACE = b" \t\r\n\v\f"
# A record starts with a line that starts with ">"; its name runs to the first
# whitespace, past any space or tab after the ">", as metacooler's readers and
# metaTOR read it. The first record is the file's first text; every other
# follows a line break, which a search for the literal "\n>" finds several
# times faster than a multi-line "^>".
_FIRST_TEXT = re.compile(rb"\S")
_FIRST_HEADER = re.compile(rb">[ \t]*(\S*)")
_HEADER = re.compile(rb"\n>[ \t]*(\S*)")
# How many names an error shows, before counting the rest.
_SHOWN = 5
_COPY_BUFFER = 1 << 24


@dataclass(frozen=True)
class _Group:
    """What one genome, cluster or group writes: its contigs (rows of
    ``/contigs``), in map order or in the caller's (``given``)."""

    name: str
    contig_ids: np.ndarray
    #: The contigs in no genome or cluster: never stitched.
    unbinned: bool = False
    #: The order is the caller's: kept when stitching too.
    given: bool = False


@dataclass(frozen=True)
class _Record:
    """Where a record lies in the assembly: from its ``>`` (``start``) to the
    next record (``end``), its sequence from ``sequence``."""

    start: int
    sequence: int
    end: int


def write_fasta(
    uri: str | Path,
    fasta: str | Path,
    out: str | Path,
    *,
    genomes: str | Sequence[str] | None = None,
    filters: str | Sequence[str] | None = None,
    round: str | None = None,
    groups: Mapping[str, Sequence[str]] | None = None,
    unbinned: bool = False,
    stitch: bool = False,
    spacer: int | None = None,
    single_file: bool = False,
    force: bool = False,
    temp_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Write the contigs of each genome, or of each cluster of a round, as FASTA.

    By default every genome of ``/genomes`` but ``__unbinned__`` gets a file,
    ``<out>/<genome>.fa``, holding its contigs' records exactly as they are in
    *fasta* (header and line breaks included), in map order.

    Each call reads the assembly once (and copies a gzip-compressed one
    first): give every group of a run to one call.

    Args:
        uri: A ``.metacool`` path, or an entry of an ``.mmetacool``
            (``file.mmetacool::entries/<name>``).
        fasta: The assembly the map was made from: FASTA, plain or compressed
            with gzip (or bgzip); a pipe is read too. Every contig written must
            be in it, once, with the map's length. Nothing is written next to
            it.
        out: The directory to write ``<name>.fa`` into, created if missing
            (files of other names in it are left alone); with *single_file*,
            the FASTA file to write (its directory created if missing).
        genomes: The genome(s) to write, by name; with *round*, the clusters
            of that round. Default: every one (but ``__unbinned__``).
        filters: ``--filter`` expression(s) (see :mod:`metacooler.filters`),
            narrowing the selection: ``genomes.<column>`` selects genomes,
            ``contigs.<column>`` the contigs written, and with *round*,
            ``clusters.<column>`` selects clusters, over the round's table as
            ``dump -t clusters`` prints it.
        round: A clustering round of ``/clusters`` (``mc.cluster_rounds``):
            write its clusters instead of the genomes.
        groups: Name -> contig names: write these groups of contigs instead,
            each in the order given (the contigs must be contigs of the map).
            For clusters the map does not hold yet. Not combined with
            *genomes*, *filters* or *round*.
        unbinned: Also write the contigs in no genome (in no cluster of
            *round*, in no group of *groups*) as ``__unbinned__``: as they
            are, in map order, never stitched. Names or genome (cluster)
            filters that keep no genome (cluster) are refused all the same.
        stitch: Write each genome or cluster as one record named after it:
            its contigs' bases, longest first (ties in map order; a group's
            in the order given), joined by *spacer* N, wrapped at
            :data:`LINE_WIDTH` columns.
        spacer: Number of N between two stitched contigs, with *stitch* only;
            :data:`DEFAULT_SPACER` by default.
        single_file: Write every record to the one file *out*: a contig
            record's header line becomes ``<genome>|<header line>``; a
            stitched record stays named after its genome.
        force: Overwrite output files that exist. An output that is the map
            or the assembly is refused all the same.
        temp_dir: Where to copy a gzip-compressed or piped assembly, plain
            (created if missing). Default: the system temporary directory
            (``TMPDIR``).

    Returns:
        Each genome, cluster or group written (``__unbinned__`` included),
        in the order written, mapped to its file.

    Raises:
        KeyError: If a genome, a round, a cluster or a group's contig is not
            one of the map.
        ValueError: If the arguments do not go together, a name cannot be
            that of a file or a record (or two names that of one file, case
            aside), the assembly is not FASTA, lacks a contig, names one twice
            or gives one another length than the map, or a compressed one
            cannot be decompressed (a :class:`metacooler.filters.FilterError`
            for a filter that does not apply, or selects nothing).
        FileExistsError: If an output exists and *force* is False.
        IsADirectoryError: If an output file is a directory.
        NotADirectoryError: If *out* is a file, without *single_file*.
        OSError: If the plain copy of the assembly cannot be written (the
            message names it: see *temp_dir*).
    """
    from metacooler.api import MetaCooler
    from metacooler.filters import TABLES_WITH_CLUSTERS, parse_filters

    if spacer is None:
        spacer = DEFAULT_SPACER
    elif not stitch:
        raise ValueError("spacer separates stitched contigs: give stitch=True too.")
    if isinstance(spacer, bool) or not isinstance(spacer, (int, np.integer)):
        raise ValueError(f"spacer is a number of N: got {spacer!r}.")
    if spacer < 0:
        raise ValueError(f"spacer is a number of N, 0 or more: got {spacer}.")
    if groups is not None and (genomes or filters or round is not None):
        raise ValueError(
            "groups names the contigs to write: it cannot be combined with "
            "genomes, filters or round."
        )
    if isinstance(genomes, str):
        genomes = [genomes]
    if isinstance(filters, str):
        filters = [filters]
    parsed = parse_filters(filters, tables=TABLES_WITH_CLUSTERS)
    fasta, out = Path(fasta), Path(out)
    mc = MetaCooler(uri)
    contigs = mc.contigs
    if groups is not None:
        selection = _given_groups(mc, contigs, groups, unbinned=unbinned)
    elif round is not None:
        selection = _round_groups(mc, contigs, round, genomes, parsed, unbinned)
    else:
        selection = _genome_groups(mc, contigs, genomes, parsed, unbinned)
    _check_names(selection, files=not single_file, headers=stitch or single_file)
    destinations = _destinations(
        selection, out, single_file=single_file, force=force, inputs=[uri, fasta]
    )

    names = contigs["name"].to_numpy(dtype=object)
    lengths = contigs["length"].to_numpy(dtype=np.int64)
    wanted = np.unique(np.concatenate([group.contig_ids for group in selection]))
    temp = None if temp_dir is None else Path(temp_dir)
    with _assembly(fasta, temp) as data:
        records = _find_records(data, {str(names[i]) for i in wanted}, fasta)
        _check_lengths(data, records, names[wanted], lengths[wanted], fasta, mc.uri)
        writer = _Writer(data, records, names, lengths, stitch=stitch, spacer=spacer)
        if single_file:
            out.parent.mkdir(parents=True, exist_ok=True)
            plan = [(out, lambda handle: writer.write_all(handle, selection))]
        else:
            out.mkdir(parents=True, exist_ok=True)
            plan = [
                (destinations[group.name], _bind(writer, group)) for group in selection
            ]
        _write_files(plan)
    written = [group for group in selection if not group.unbinned]
    kind = "genome(s)" if round is None else "cluster(s)"
    if groups is not None:
        kind = "group(s)"
    logger.info(
        "Wrote %d FASTA file(s) to %s: %d %s%s, %d record(s), %s bp of %d contig(s).",
        len(set(destinations.values())),
        out,
        len(written),
        kind,
        f" and {UNBINNED_GENOME}" if len(written) < len(selection) else "",
        writer.records,
        f"{writer.bases:,}",
        writer.contigs,
    )
    return destinations


def _bind(writer: _Writer, group: _Group) -> Callable[[BinaryIO], None]:
    """What writes *group*'s file."""
    return lambda handle: writer.write(handle, group)


# ----------------------------------------------------------------------
# Selection
# ----------------------------------------------------------------------
def _shown(items: Sequence[object]) -> str:
    """The first items, then how many more."""
    shown = ", ".join(str(item) for item in items[:_SHOWN])
    return shown + (f" and {len(items) - _SHOWN} more" if len(items) > _SHOWN else "")


def _refuse_tables(parsed: Sequence[Filter], table: str, reason: str) -> None:
    """Refuse a filter on *table*, which this selection does not read."""
    from metacooler.filters import FilterError

    for filt in parsed:
        if filt.table == table:
            raise FilterError(f"Filter {filt.raw!r}: {reason}")


def _check_known(asked: Sequence[str], known: np.ndarray, what: str) -> None:
    """Refuse names that are none of *known*."""
    present = set(known.tolist())
    missing = [name for name in asked if name not in present]
    if missing:
        raise KeyError(
            f"Not {what}: {_shown(missing)}. There are {len(known)}: "
            f"{_shown(known.tolist())}."
        )


def _runs(keys: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    """The positions of each run of equal *keys*: its key, and its positions."""
    if not len(keys):
        return keys, []
    cuts = np.flatnonzero(keys[1:] != keys[:-1]) + 1
    return keys[np.concatenate([[0], cuts])], np.split(np.arange(len(keys)), cuts)


def _unbinned_group(rows: np.ndarray, taken: Sequence[str]) -> list[_Group]:
    """The ``__unbinned__`` group of *rows*, if any; never alongside a genome,
    cluster or group of that name."""
    if UNBINNED_GENOME in taken:
        raise ValueError(
            f"A cluster or group is named {UNBINNED_GENOME!r}, which names the "
            "contigs in none: rename it to write the unbinned contigs too."
        )
    if not rows.size:
        logger.info("No contig is in no genome or cluster: no %s.", UNBINNED_GENOME)
        return []
    return [_Group(UNBINNED_GENOME, rows, unbinned=True)]


def _nothing_selected(
    parsed: Sequence[Filter],
    frames: Mapping[str, pd.DataFrame],
    named: Sequence[str] | None,
    what: str,
    otherwise: str,
    hint: str = "",
) -> NoReturn:
    """Refuse a selection that keeps nothing: which filter emptied it, or
    *otherwise* without a filter."""
    from metacooler.filters import FilterError, explain_empty

    if not parsed:
        raise ValueError(otherwise)
    context = [hint] if hint else []
    if named:
        context.insert(
            0,
            f"The {len(named)} {what} named by --genomes were then intersected "
            "with the above, leaving none.",
        )
    raise FilterError(explain_empty(frames, parsed, context="\n".join(context)))


def _genome_groups(
    mc: MetaCooler,
    contigs: pd.DataFrame,
    genomes: Sequence[str] | None,
    parsed: Sequence[Filter],
    unbinned: bool,
) -> list[_Group]:
    """Every genome (but ``__unbinned__``), or those named, narrowed by the
    filters; a genome whose every contig a filter drops writes nothing."""
    from metacooler.filters import filter_mask

    _refuse_tables(
        parsed,
        "clusters",
        "it selects the clusters of a round: give the round (--round NAME), or "
        "filter genomes.",
    )
    table = mc.genomes
    genome_names = table["name"].to_numpy(dtype=object)
    by_genome = genomes_filter_frame(table)
    by_contig = contigs_filter_frame(contigs, table["name"])
    keep_genome = filter_mask(by_genome, parsed, table="genomes")
    if genomes:
        asked = [str(name) for name in genomes]
        if UNBINNED_GENOME in asked:
            raise ValueError(
                f"{UNBINNED_GENOME} is no genome: it holds the contigs in none. "
                "Write them with --unbinned (unbinned=True)."
            )
        _check_known(asked, genome_names, f"genomes of {mc.uri}")
        keep_genome &= np.isin(genome_names, asked)
    unbinned_id = unbinned_genome_id(genome_names)
    # The filters kept __unbinned__, which is no genome: --unbinned writes it.
    kept_unbinned = unbinned_id is not None and bool(keep_genome[unbinned_id])
    if unbinned_id is not None:
        keep_genome[unbinned_id] = False
    keep_contig = filter_mask(by_contig, parsed, table="contigs")
    genome_ids = contigs["genome_id"].to_numpy(dtype=np.int64)

    rows = np.flatnonzero(keep_contig & keep_genome[genome_ids])
    rows = rows[np.argsort(genome_ids[rows], kind="stable")]
    ids, positions = _runs(genome_ids[rows])
    selection = [
        _Group(str(genome_names[gid]), rows[at])
        for gid, at in zip(ids, positions, strict=True)
    ]
    # Every genome holds a contig: only a contig filter empties one.
    emptied = int(keep_genome.sum()) - len(selection)
    if emptied:
        logger.info(
            "%d genome(s) have no contig left after --filter: no file for them.",
            emptied,
        )
    restricted = bool(genomes) or any(f.table == "genomes" for f in parsed)
    frames = {"genomes": by_genome, "contigs": by_contig}
    if not selection and (not unbinned or (restricted and not kept_unbinned)):
        _nothing_selected(
            parsed,
            frames,
            genomes,
            "genome(s)",
            f"{mc.uri} has no genome to write: every contig is in "
            f"{UNBINNED_GENOME}. Write them with --unbinned (unbinned=True).",
            hint=(
                f"{UNBINNED_GENOME} is no genome: write its contigs with --unbinned."
                if kept_unbinned and not unbinned
                else ""
            ),
        )
    if unbinned:
        # A map without __unbinned__ has every contig in a genome.
        in_none = genome_ids == (-1 if unbinned_id is None else unbinned_id)
        selection += _unbinned_group(np.flatnonzero(keep_contig & in_none), [])
        if not selection:
            _nothing_selected(
                parsed, frames, genomes, "genome(s)", f"{mc.uri} has no contig."
            )
    return selection


def _round_groups(
    mc: MetaCooler,
    contigs: pd.DataFrame,
    round: str,
    clusters: Sequence[str] | None,
    parsed: Sequence[Filter],
    unbinned: bool,
) -> list[_Group]:
    """The clusters of *round* (those named, those the filters keep), each
    with its contigs in map order."""
    from metacooler.filters import filter_mask

    _refuse_tables(
        parsed,
        "genomes",
        f"the clusters of round {round!r} are written, not genomes: filter them "
        "with 'clusters.<column>'.",
    )
    mc._check_round(round)
    every = mc.clusters
    table = every[every["round"].to_numpy(dtype=object) == round]
    cluster_names = table["name"].to_numpy(dtype=object)
    by_cluster = clusters_filter_frame(table)
    by_contig = contigs_filter_frame(contigs, mc.genomes["name"])
    keep_cluster = filter_mask(by_cluster, parsed, table="clusters")
    if clusters:
        asked = [str(name) for name in clusters]
        _check_known(asked, cluster_names, f"clusters of round {round!r}")
        keep_cluster &= np.isin(cluster_names, asked)
    keep_contig = filter_mask(by_contig, parsed, table="contigs")

    # Cluster by cluster in the order of the round, each one's contigs in
    # map order: as the memberships are stored.
    members = mc.cluster_members(round)
    rows = pd.Index(contigs["name"]).get_indexer(members["contig"])
    names, positions = _runs(members["cluster"].to_numpy(dtype=object))
    kept = cluster_names[keep_cluster].tolist()
    wanted, held = set(kept), set(names.tolist())
    selection = []
    for name, at in zip(names, positions, strict=True):
        if name in wanted:
            part = rows[at][keep_contig[rows[at]]]
            if part.size:
                selection.append(_Group(str(name), part))
    memberless = [name for name in kept if name not in held]
    if memberless:
        logger.info(
            "%d cluster(s) of round %r hold no contig: no file for them (%s).",
            len(memberless),
            round,
            _shown(memberless),
        )
    emptied = len(kept) - len(memberless) - len(selection)
    if emptied:
        logger.info(
            "%d cluster(s) of round %r have no contig left after --filter: no file "
            "for them.",
            emptied,
            round,
        )
    restricted = bool(clusters) or any(f.table == "clusters" for f in parsed)
    frames = {"clusters": by_cluster, "contigs": by_contig}
    if not selection and (restricted or not unbinned):
        _nothing_selected(
            parsed,
            frames,
            clusters,
            "cluster(s)",
            f"Round {round!r} of {mc.uri} has no contig to write"
            + (f" in {_shown(kept)}." if clusters else "."),
        )
    if unbinned:
        in_round = np.zeros(len(contigs), dtype=bool)
        in_round[rows] = True
        rest = np.flatnonzero(keep_contig & ~in_round)
        selection += _unbinned_group(rest, cluster_names.tolist())
        if not selection:
            _nothing_selected(
                parsed,
                frames,
                clusters,
                "cluster(s)",
                f"Round {round!r} of {mc.uri} has no contig to write.",
            )
    return selection


def _given_groups(
    mc: MetaCooler,
    contigs: pd.DataFrame,
    groups: Mapping[str, Sequence[str]],
    *,
    unbinned: bool,
) -> list[_Group]:
    """A caller's groups of contigs, each in the order given."""
    if not isinstance(groups, Mapping):
        raise ValueError(
            "groups maps the name of each group to its contigs' names: got "
            f"{type(groups).__name__}."
        )
    index = pd.Index(contigs["name"])
    selection = []
    for name, members in groups.items():
        if not isinstance(name, str) or not name:
            raise ValueError(f"Every group needs a name, a non-empty text: {name!r}.")
        if isinstance(members, (str, bytes)):
            raise ValueError(
                f"Group {name!r}: give its contigs as a list of names, not one text."
            )
        asked = [str(contig) for contig in members]
        if not asked:
            raise ValueError(f"Group {name!r} has no contig.")
        rows = index.get_indexer(asked)
        missing = [contig for contig, row in zip(asked, rows, strict=True) if row < 0]
        if missing:
            raise KeyError(
                f"Group {name!r}: not contigs of {mc.uri}: {_shown(missing)}."
            )
        twice = sorted(c for c, n in Counter(asked).items() if n > 1)
        if twice:
            raise ValueError(
                f"Group {name!r} names a contig more than once: {_shown(twice)}."
            )
        selection.append(_Group(name, rows.astype(np.int64), given=True))
    if unbinned:
        grouped = np.zeros(len(contigs), dtype=bool)
        for group in selection:
            grouped[group.contig_ids] = True
        selection += _unbinned_group(np.flatnonzero(~grouped), list(groups))
    if not selection:
        raise ValueError("groups is empty: name at least one group of contigs.")
    return selection


def _check_names(selection: Sequence[_Group], *, files: bool, headers: bool) -> None:
    """Refuse a name that cannot be that of a file (*files*), or of a record
    (*headers*: stitched records, and every record with ``single_file``)."""
    names = [group.name for group in selection]
    if files:
        bad = [n for n in names if n in (".", "..") or "/" in n or "\0" in n]
        if bad:
            raise ValueError(
                f"Cannot name a file after {_shown([repr(n) for n in bad])}: write "
                "every record to one file instead (--single-file)."
            )
        folded: dict[str, list[str]] = {}
        for name in names:
            folded.setdefault(name.casefold(), []).append(name)
        same = [" and ".join(map(repr, n)) for n in folded.values() if len(n) > 1]
        if same:
            raise ValueError(
                f"{_shown(same)} would write one file where case is ignored "
                "(macOS, Windows): write every record to one file instead "
                "(--single-file)."
            )
    if headers:
        bad = [n for n in names if any(c.isspace() for c in n)]
        if bad:
            raise ValueError(
                "A record cannot be named after "
                f"{_shown([repr(n) for n in bad])}: a FASTA name ends at the "
                "first whitespace. Write one file per genome or cluster, the "
                "contigs as they are (neither --stitch nor --single-file)."
            )


def _destinations(
    selection: Sequence[_Group],
    out: Path,
    *,
    single_file: bool,
    force: bool,
    inputs: Sequence[str | Path],
) -> dict[str, Path]:
    """The file of each group, checked before the assembly is read.

    An output that is one of *inputs* (the map, the assembly) is refused even
    with *force*: a file is replaced only once its new version is written.
    """
    if single_file:
        check_output_path(out, force, inputs=inputs, defer_removal=True)
        return {group.name: out for group in selection}
    if out.exists() and not out.is_dir():
        raise NotADirectoryError(
            f"Output {str(out)!r} is a file: one FASTA file per genome or cluster "
            "is written into a directory, unless every record goes to one file "
            "(--single-file)."
        )
    paths = {group.name: out / f"{group.name}{FASTA_SUFFIX}" for group in selection}
    existing = []
    for path in paths.values():
        # Refuses an input, or a directory, in the way: replaced only once
        # written, an existing file is kept until then.
        check_output_path(path, True, inputs=inputs, defer_removal=True)
        if path.exists():
            existing.append(path.name)
    if existing and not force:
        raise FileExistsError(
            f"{len(existing)} output file(s) already exist in {str(out)!r}: "
            f"{_shown(existing)}. Use --force to overwrite."
        )
    return paths


# ----------------------------------------------------------------------
# The assembly
# ----------------------------------------------------------------------
@contextlib.contextmanager
def _assembly(
    path: Path, temp_dir: Path | None
) -> Generator[bytes | mmap.mmap, None, None]:
    """The bytes of the assembly, mapped in memory: a plain file as it is; a
    gzip-compressed one, or a pipe, copied plain into a private temporary
    directory first (in *temp_dir*, by default the system's)."""
    if stat.S_ISREG(os.stat(path).st_mode):
        with open(path, "rb") as handle:
            compressed = handle.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
        if not compressed:
            with _mapped(path) as data:
                yield data
            return
    if temp_dir is not None:
        temp_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="metacooler-fasta-", dir=temp_dir) as tmp:
        plain = Path(tmp) / "assembly.fa"
        _copy_plain(path, plain)
        with _mapped(plain) as data:
            yield data


def _copy_plain(source: Path, target: Path) -> None:
    """Copy *source*, decompressed if it is gzip, into *target*: read once, as
    a pipe can be.

    Raises:
        ValueError: If *source* is compressed but cannot be decompressed.
        OSError: If *target* cannot be written, named in the message.
    """
    with open(source, "rb") as handle:
        compressed = handle.peek(len(_GZIP_MAGIC))[: len(_GZIP_MAGIC)] == _GZIP_MAGIC
        logger.info(
            "%s %s into %s.",
            "Decompressing" if compressed else "Copying",
            source,
            target.parent,
        )
        reader = gzip.GzipFile(fileobj=handle) if compressed else handle
        try:
            with open(target, "wb") as copy:
                shutil.copyfileobj(reader, copy, _COPY_BUFFER)
        except (gzip.BadGzipFile, EOFError, zlib.error) as exc:
            raise ValueError(f"Cannot decompress {source}: {exc}") from exc
        except OSError as exc:
            raise OSError(
                exc.errno,
                f"Cannot copy {source} into {target}, its plain copy (the "
                f"temporary directory is --temp-dir's): {exc.strerror or exc}",
            ) from exc


@contextlib.contextmanager
def _mapped(path: Path) -> Generator[bytes | mmap.mmap, None, None]:
    """A file's bytes, read-only, mapped in memory (an empty file: ``b""``)."""
    with open(path, "rb") as handle:
        if os.fstat(handle.fileno()).st_size:
            with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
                yield data
        else:
            yield b""


def _find_records(
    data: bytes | mmap.mmap, wanted: set[str], fasta: Path
) -> dict[str, _Record]:
    """Where the records named in *wanted* lie in the assembly, in one pass.

    Raises:
        ValueError: If *data* is not FASTA, lacks a contig of *wanted*, or
            names one by several records.
    """
    keys = {name.encode("ascii"): name for name in wanted}  # names are ASCII
    first = _FIRST_TEXT.search(data)
    at = -1 if first is None else first.start()
    if at < 0 or data[at : at + 1] != b">" or (at and data[at - 1 : at] != b"\n"):
        raise ValueError(
            f"{fasta} is not a FASTA file: it does not start with a record, a "
            "line that starts with '>' (plain or gzip-compressed FASTA)."
        )
    records: dict[str, _Record] = {}
    twice: set[str] = set()
    pending: tuple[str, int, int] | None = None
    head = _FIRST_HEADER.match(data, at)
    assert head is not None  # data[at] is ">"
    headers = [(at, head.group(1), head.end())]
    others = ((m.start() + 1, m.group(1), m.end()) for m in _HEADER.finditer(data, at))
    for start, raw, after in chain(headers, others):
        if pending is not None:
            records[pending[0]] = _record(data, *pending[1:], start)
            pending = None
        name = keys.get(raw)
        if name is not None:
            if name in records:
                twice.add(name)
            else:
                pending = (name, start, after)
    if pending is not None:
        records[pending[0]] = _record(data, *pending[1:], len(data))
    missing = sorted(wanted - set(records))
    if missing:
        raise ValueError(
            f"{len(missing)} contig(s) to write are not in {fasta}: "
            f"{_shown(missing)}. Give the assembly the map was made from."
        )
    if twice:
        raise ValueError(
            f"{len(twice)} contig(s) to write are named by several records of "
            f"{fasta}: {_shown(sorted(twice))}."
        )
    return records


def _record(data: bytes | mmap.mmap, start: int, after_name: int, end: int) -> _Record:
    """The record from *start* to *end*, its sequence after its header line."""
    line_end = data.find(b"\n", after_name, end)
    return _Record(start, end if line_end < 0 else line_end + 1, end)


def _bases(data: bytes | mmap.mmap, record: _Record) -> bytes:
    """A record's bases: its sequence, every whitespace removed."""
    return data[record.sequence : record.end].translate(None, _WHITESPACE)


def _check_lengths(
    data: bytes | mmap.mmap,
    records: Mapping[str, _Record],
    names: np.ndarray,
    lengths: np.ndarray,
    fasta: Path,
    uri: str,
) -> None:
    """Refuse an assembly that gives a contig another length than the map.

    A map's length counts a contig's sequence as the tool that made it read
    the assembly: every byte but line breaks (pyfaidx, which ``create
    --fasta`` reads it with, counts a trailing space as a base), or its bases
    alone. Either is the map's contig.
    """
    wrong = []
    for name, length in zip(names.tolist(), lengths.tolist(), strict=True):
        record = records[name]
        sequence = data[record.sequence : record.end]
        if len(sequence) - sequence.count(b"\n") - sequence.count(b"\r") == length:
            continue
        found = len(sequence.translate(None, _WHITESPACE))
        if found != length:
            wrong.append(f"{name} ({found} bp in the FASTA, {length} in the map)")
    if wrong:
        raise ValueError(
            f"{len(wrong)} contig(s) have another length in {fasta} than in {uri}, "
            f"so it is not the assembly of this map: {_shown(wrong)}."
        )


# ----------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------
class _Wrapped:
    """Sequence lines of :data:`LINE_WIDTH` bases, from chunks of any size."""

    def __init__(self, handle: BinaryIO) -> None:
        self.handle = handle
        self.pending = b""

    def write(self, chunk: bytes) -> None:
        data = self.pending + chunk
        full = len(data) - len(data) % LINE_WIDTH
        if full:
            lines = (data[i : i + LINE_WIDTH] for i in range(0, full, LINE_WIDTH))
            self.handle.write(b"\n".join(lines) + b"\n")
        self.pending = data[full:]

    def close(self) -> None:
        if self.pending:
            self.handle.write(self.pending + b"\n")
            self.pending = b""


class _Writer:
    """Writes groups of records from the assembly, and counts them."""

    def __init__(
        self,
        data: bytes | mmap.mmap,
        records: Mapping[str, _Record],
        names: np.ndarray,
        lengths: np.ndarray,
        *,
        stitch: bool,
        spacer: int,
    ) -> None:
        self.data = data
        self._records = records
        self.names = names
        self.lengths = lengths
        self.stitch = stitch
        self.spacer = b"N" * spacer
        #: What was written: contigs, records and bases (spacers aside).
        self.contigs = 0
        self.records = 0
        self.bases = 0

    def write_all(self, handle: BinaryIO, selection: Sequence[_Group]) -> None:
        """Every group, one after the other, in one file: a contig record's
        header line prefixed with its group."""
        for group in selection:
            self.write(handle, group, prefix=group.name)

    def write(self, handle: BinaryIO, group: _Group, prefix: str | None = None) -> None:
        """One group's records: one stitched record, or its contigs' records
        as they are (their header line prefixed with *prefix*)."""
        ids = group.contig_ids
        self.contigs += len(ids)
        if self.stitch and not group.unbinned:
            if not group.given:
                # Longest first; a stable sort keeps map order among equals.
                ids = ids[np.argsort(-self.lengths[ids], kind="stable")]
            handle.write(f">{group.name}\n".encode())
            wrapped = _Wrapped(handle)
            for i, name in enumerate(self.names[ids].tolist()):
                bases = _bases(self.data, self._records[name])
                if i and self.spacer:
                    wrapped.write(self.spacer)
                wrapped.write(bases)
                self.bases += len(bases)
            wrapped.close()
            self.records += 1
            return
        head = b">" if prefix is None else f">{prefix}{SINGLE_FILE_SEPARATOR}".encode()
        for name in self.names[ids].tolist():
            record = self._records[name]
            body = self.data[record.start : record.end]
            if not body.endswith(b"\n"):
                body += b"\n"  # the last record of an assembly that ends without
            handle.write(head + body[1:])
        self.bases += int(self.lengths[ids].sum())
        self.records += len(ids)


def _write_files(plan: Sequence[tuple[Path, Callable[[BinaryIO], None]]]) -> None:
    """Write every file of *plan* beside its destination, then, once all of
    them are written, rename each over its destination.

    A failure while writing leaves no new file and every existing one as it
    was; only a failure between two renames can leave a set half updated. The
    files get the mode of any new file (``mkstemp`` makes them private).
    """
    mode = 0o666 & ~_umask()
    staged: list[tuple[Path, Path]] = []
    try:
        for destination, write in plan:
            fd, name = tempfile.mkstemp(
                dir=destination.parent, prefix=".fasta-", suffix=".tmp"
            )
            staged.append((Path(name), destination))
            with os.fdopen(fd, "wb") as handle:
                write(handle)
            os.chmod(name, mode)
        for path, destination in staged:
            os.replace(path, destination)
    except BaseException:
        for path, _ in staged:
            with contextlib.suppress(OSError):
                path.unlink()
        raise
