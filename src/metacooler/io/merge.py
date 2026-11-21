"""Merge operations for metacool files."""

from __future__ import annotations

import contextlib
import hashlib
import json
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, cast

import cooler
import h5py
import numpy as np
import pandas as pd
from cooler.api import info as cooler_info

from metacooler._version import __version__
from metacooler.io._utils import (
    check_output_path,
    cooler_writer_shims,
    store_qc,
    write_column,
    write_cooler,
)
from metacooler.io.clusters import clusters_digest, copy_clusters
from metacooler.io.staging import (
    MapChangedError,
    check_inputs_unchanged,
    input_states,
)
from metacooler.qc import (
    QC_CATEGORIES,
    QC_CONTIG_COLUMNS,
    QC_GENOME_COLUMNS,
    QC_METADATA_KEY,
    QC_PAIRS_COUNT_KEYS,
    QC_SUMMARY_KEYS,
    ContactQC,
    _contig_counts_add_up,
    _counts_problem,
    is_qc_column,
    unbinned_genome_id,
)
from metacooler.schema import (
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
    _attr,
    validate,
)
from metacooler.uri import parse_uri

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable

logger = logging.getLogger("metacooler")

# The keys of a QC block that hold a count -- of contacts, or of the pairs a
# ``create pairs`` tally read and discarded -- or a quantity derived from the
# counts: what another library changes.  The pairs counts are summed, the rest
# recomputed from the summed tally.  Every other key (``version``, ``source``,
# ``distance_threshold``, or one a later release records) says how the contacts
# were counted, and the inputs' blocks add up only where those agree --
# ``binsize`` aside (see :data:`_QC_BINSIZE_KEY`).
_QC_COUNT_KEYS = (
    QC_SUMMARY_KEYS - {"version", "source", "distance_threshold"}
) | frozenset(QC_PAIRS_COUNT_KEYS)

#: The columns of ``/contigs`` and ``/genomes`` that say *what* a map is: which
#: assembly, which genomes, which contig belongs to which, how many bins each
#: holds.  Two maps merge only when these agree.  Every other column is
#: metadata -- completeness, taxonomy, coverage, added at creation or later by
#: ``annotate`` -- and the ``qc_*`` columns describe each input's own contacts;
#: neither changes what the pixels mean, so neither may block a merge of
#: replicates.  Requiring them equal made annotating one replicate (or
#: re-running ``qc --store`` on it) enough to make it unmergeable.
_IDENTITY_COLUMNS: dict[str, tuple[str, ...]] = {
    "contigs": ("name", "length", "genome_id"),
    "genomes": ("name", "length", "n_contigs", "n_bins"),
}


def _column_bytes(dataset: h5py.Dataset) -> bytes:
    """Serialise one column for hashing, independently of its storage width.

    Integers are widened to ``int64`` and floats to ``float64`` (``from_cool``
    keeps a source's ``int32`` lengths where another writer stores ``int64``);
    text is joined on NUL, which no name can contain, so neither a fixed-length
    dtype's padding nor a vlen string's object array changes the digest.
    """
    arr = dataset[:]
    if arr.dtype.kind in "iub":
        return arr.astype(np.int64).tobytes()
    if arr.dtype.kind == "f":
        return arr.astype(np.float64).tobytes()
    return b"\0".join(
        v.encode() if isinstance(v, str) else bytes(v) for v in arr.ravel()
    )


def _hash_group(f: h5py.Group, group: str) -> str:
    """Return an MD5 hex digest of the identity columns of a table.

    Used to compare ``contigs`` and ``genomes`` tables across files for
    merge compatibility checking: only the columns in
    :data:`_IDENTITY_COLUMNS` take part (see there for why).

    Args:
        f: Open ``h5py.File`` (or any HDF5 group-like object).
        group: Table to hash: ``"contigs"`` or ``"genomes"``.

    Returns:
        Lowercase hex MD5 digest of the identity columns, in a fixed order.
    """
    h = hashlib.md5(usedforsecurity=False)
    g = f[group]
    for key in _IDENTITY_COLUMNS[group]:
        h.update(key.encode())
        h.update(_column_bytes(g[key]) if key in g else b"<absent>")
    return h.hexdigest()


def _metadata_digests(f: h5py.Group) -> dict[str, str]:
    """Digest of every metadata column, keyed ``table.column``.

    Metadata columns are those neither in :data:`_IDENTITY_COLUMNS` nor part
    of the QC block (:func:`metacooler.qc.is_qc_column`): a user's ``qc_pass``
    in ``/contigs`` is compared like any other.
    """
    digests = {}
    for table, identity in _IDENTITY_COLUMNS.items():
        for key in f[table]:
            if key in identity or is_qc_column(table, key):
                continue
            digests[f"{table}.{key}"] = hashlib.md5(
                _column_bytes(f[table][key]), usedforsecurity=False
            ).hexdigest()
    return digests


def _metadata_differences(first: dict[str, str], other: dict[str, str]) -> list[str]:
    """Describe how an input's metadata columns differ from the first input's."""
    out = []
    for name in sorted(set(first) | set(other)):
        if name not in other:
            out.append(f"{name} (absent here)")
        elif name not in first:
            out.append(f"{name} (only here)")
        elif first[name] != other[name]:
            out.append(f"{name} (other values)")
    return out


def _bins_signature(f: h5py.Group) -> tuple[int, str]:
    """Return ``(n_bins, digest)`` describing a file's bin geometry.

    Only the coordinate columns (``chrom``/``start``/``end``) take part: an
    optional balancing ``weight`` column is allowed to differ between inputs.
    Values are widened to ``int64`` first so the digest is independent of the
    on-disk integer width (cooler may store ``chrom`` as an HDF5 enum).

    Args:
        f: Open ``h5py.File`` positioned on a cooler/metacool group.

    Returns:
        The number of bins and a lowercase hex MD5 digest of the coordinates.
    """
    g = f["bins"]
    n_bins = int(g["chrom"].shape[0])
    h = hashlib.md5(usedforsecurity=False)
    for key in ("chrom", "start", "end"):
        h.update(np.asarray(g[key][:]).astype(np.int64).tobytes())
    return n_bins, h.hexdigest()


def _binsize(f: h5py.Group) -> int | None:
    """Return the ``bin-size`` attribute of *f*, or ``None`` for variable bins."""
    binsize = f.attrs.get("bin-size")
    if binsize is None or (isinstance(binsize, str) and binsize.lower() == "null"):
        return None
    return int(binsize)


@contextlib.contextmanager
def _open_uri(uri: str) -> Generator[h5py.Group, None, None]:
    """Open the map a ``file`` or ``file::group`` URI addresses, read-only.

    Args:
        uri: Path to a ``.metacool`` or an ``.mmetacool`` entry URI.

    Yields:
        The HDF5 group holding the map (the file root for a plain path).
    """
    file_path, group_path = parse_uri(uri)
    with h5py.File(file_path, "r") as f:
        yield f[group_path] if group_path != "/" else f


def _merged_count_dtype(uris: list[str]) -> np.dtype:
    """Return ``float64`` if any input stores float counts, else ``int64``.

    ``int64`` is mandatory for integer inputs -- summing replicates readily
    exceeds the ``int32`` they are usually stored in -- but forcing it on a
    float input (a balanced or rescaled map) truncates every count toward zero.

    Args:
        uris: Input URIs.

    Returns:
        The dtype to sum ``count`` into.
    """
    for uri in uris:
        with _open_uri(uri) as g:
            if g["pixels"]["count"].dtype.kind == "f":
                return np.dtype(np.float64)
    return np.dtype(np.int64)


def _refuse_repeated_inputs(inputs: list[str]) -> None:
    """Refuse a map given to ``merge`` twice, which would double its counts.

    ``merge A A`` summed every pixel of A with itself and exited 0. A map is
    identified by its file (device and inode, so a symlink, ``./`` alias or
    hard link of a listed file is caught too) and its group, so two different
    entries of one ``.mmetacool`` still merge.

    Args:
        inputs: Input paths or ``file::group`` URIs.

    Raises:
        ValueError: If two inputs name the same map.
    """
    seen: dict[tuple[int, int, str], str] = {}
    for uri in inputs:
        file_path, group = parse_uri(uri)
        st = os.stat(file_path)
        key = (st.st_dev, st.st_ino, "/" + group.strip("/"))
        if key in seen:
            raise ValueError(
                f"{uri!r} names the same map as {seen[key]!r}; merging a map "
                "with itself would count every one of its contacts twice."
            )
        seen[key] = uri


# ``binsize``: the bin size a map had when its block was stored.  The inputs of
# a merge share one bin grid -- merge checks that before anything else -- so it
# changes no count between them, and their blocks add up whether or not each
# records it.  The merged block records it when every input records the same
# value, and not at all otherwise.
_QC_BINSIZE_KEY = "binsize"

# The four ``qc_<category>`` count columns of ``/genomes``: what a tally sums.
_GENOME_COUNT_COLUMNS = QC_GENOME_COLUMNS[: len(QC_CATEGORIES)]

# What :func:`_qc_block` returns for a map whose ``metadata`` attribute is not a
# JSON object: whether it holds a QC block cannot be told.
_NOT_JSON = object()

# A key a block does not record, for :func:`_merged_extra` to compare.
_ABSENT = object()


class _SummedQC(NamedTuple):
    """The inputs' QC, summed: counts only, until the merged map can take them.

    Held through the pixel merge, so it holds no name: those are read from the
    merged map when the QC is written into it (:func:`_store_summed_qc`).

    Attributes:
        tally: The sum of the inputs' tallies.
        extra: The keys of the merged block that ``tally.summary()`` does not
            write itself (:func:`_merged_extra`).
    """

    tally: ContactQC
    extra: dict


def _qc_block(group: h5py.Group) -> object:
    """A map's ``metadata["qc"]`` value, read as cooler reads the attribute.

    ``None`` when it has none; :data:`_NOT_JSON` when its ``metadata``
    attribute is not a JSON object, which may or may not have held one.
    """
    metadata = cooler_info(group).get("metadata")
    if metadata is None:
        return None
    if not isinstance(metadata, dict):
        return _NOT_JSON
    return metadata.get(QC_METADATA_KEY)


def _columns(table: h5py.Group, columns: tuple[str, ...]) -> np.ndarray:
    """*columns* of *table* side by side; the caller has checked they are there."""
    return np.column_stack([table[column][:] for column in columns])


def _names(table: h5py.Group) -> np.ndarray:
    """The ``name`` column of ``/genomes`` or ``/contigs``, decoded."""
    return table["name"].asstr()[:]


def _input_problem(uri: str, group: h5py.Group, block: object) -> str | None:
    """Why the QC of one input cannot be summed, whatever the others hold.

    Told from its ``metadata`` attribute and the types of its ``/genomes``
    count columns: not a count is read.  metacooler reserves those names in
    ``/genomes``, as 0.2.0 already did, so a count column there that holds
    text was not written by it, and nothing of the input's QC is trusted.
    ``/contigs`` is another matter (:func:`_contig_count_problems`).
    """
    if block is _NOT_JSON:
        return f"the metadata attribute of {uri} is not a JSON object"
    if not isinstance(block, dict):
        return f"the QC block of {uri} is not a JSON object"
    genomes = group["genomes"]
    if not all(column in genomes for column in _GENOME_COUNT_COLUMNS):
        return f"{uri} lacks the qc_* count columns of /genomes that go with it"
    if any(genomes[column].dtype.kind not in "iuf" for column in _GENOME_COUNT_COLUMNS):
        return f"the qc_* columns of /genomes of {uri} are not numbers"
    return None


def _measured_differently(first: dict, other: dict) -> list[str]:
    """The keys on which two QC blocks disagree about how they were counted.

    A key recorded by one block and not the other is one of them, counts
    included: a ``create pairs`` block without its pairs counts cannot add up
    with one that has them.  ``binsize`` is not (see :data:`_QC_BINSIZE_KEY`).
    A key this release does not know is compared like ``source``: the blocks
    add up only where every input records the same value, which the merged
    block then carries.
    """
    return [
        key
        for key in sorted((set(first) | set(other)) - {_QC_BINSIZE_KEY})
        if (key in first) != (key in other)
        or (key not in _QC_COUNT_KEYS and first[key] != other[key])
    ]


def _unsummable(
    inputs: list[str], blocks: list[object], problems: list[str | None]
) -> str | None:
    """Why the inputs' QC blocks cannot be summed, or ``None`` if they can.

    They can when every input has one, with its ``/genomes`` count columns
    (:func:`_input_problem`), and every block was counted alike: the same
    keys, ``binsize`` aside, and the same value of each that is not a count
    (:func:`_measured_differently`).  What each holds is checked as it is read
    back (:func:`_summed_tally`).
    """
    missing = [uri for uri, block in zip(inputs, blocks, strict=True) if block is None]
    if missing:
        return f"{', '.join(missing)} {'has' if len(missing) == 1 else 'have'} none"
    problem = next((problem for problem in problems if problem), None)
    if problem:
        return problem
    first = blocks[0]
    assert isinstance(first, dict)
    for uri, block in zip(inputs[1:], blocks[1:], strict=True):
        assert isinstance(block, dict)
        keys = _measured_differently(first, block)
        if keys:
            return "they were not counted alike: " + "; ".join(
                f"{key} is {_described(first, key)} in {inputs[0]} but "
                f"{_described(block, key)} in {uri}"
                for key in keys
            )
    return None


def _described(block: dict, key: str) -> str:
    """The value of *key* in *block*, for a message: ``absent`` when it has none."""
    return repr(block[key]) if key in block else "absent"


class _ContigProblems(NamedTuple):
    """Why inputs' ``/contigs`` QC columns cannot be summed: empty if they can.

    Attributes:
        none: The inputs with neither column: maps written by metacooler
            0.2.0 or earlier.
        partial: One reason per other shortfall: an input with only one of
            the two, or with one that is not numbers -- a user's column of
            that name, which 0.2.0 stored as given, since it reserved no QC
            name in ``/contigs``.
        stale: The inputs whose columns no longer add up to their block, which
            a ``qc --store`` of such a release rewrote without them.
    """

    none: list[str]
    partial: list[str]
    stale: list[str]


def _contig_count_problems(inputs: list[str], blocks: list[dict]) -> _ContigProblems:
    """Why the inputs' ``/contigs`` QC columns cannot be summed, if they cannot.

    Each input's columns are read only when both are there as numbers, and let
    go before the next one's: they are as long as the map has contigs.
    """
    problems = _ContigProblems([], [], [])
    for uri, block in zip(inputs, blocks, strict=True):
        with _open_uri(uri) as group:
            contigs = group["contigs"]
            absent = [column for column in QC_CONTIG_COLUMNS if column not in contigs]
            if len(absent) == len(QC_CONTIG_COLUMNS):
                problems.none.append(uri)
                continue
            text = [
                column
                for column in QC_CONTIG_COLUMNS
                if column not in absent and contigs[column].dtype.kind not in "iuf"
            ]
            if absent:
                problems.partial.append(f"{uri} lacks {' and '.join(absent)}")
            if text:
                problems.partial.append(
                    f"the {' and '.join(text)} of {uri} "
                    f"{'is' if len(text) == 1 else 'are'} not numbers (a user's "
                    "column of that name, as metacooler 0.2.0 stored one)"
                )
            if absent or text:
                continue
            per_contig = _columns(contigs, QC_CONTIG_COLUMNS)
        if not _contig_counts_add_up(block, per_contig):
            problems.stale.append(uri)
        del per_contig
    return problems


def _summed_tally(
    inputs: list[str], blocks: list[dict], *, per_contig: bool
) -> ContactQC:
    """Rebuild every input's tally from what it stored, and add them up.

    One input at a time: its columns are read, checked and added to the sum,
    then let go, so that the sum is all that is held of them.

    Args:
        inputs: Input URIs.
        blocks: Their QC blocks, already found to be counted alike
            (:func:`_unsummable`).
        per_contig: Whether to sum the ``/contigs`` columns too.

    Raises:
        ValueError: Naming the input whose QC cannot be read back: a block of
            another version, or counts that are not those of one tally
            (:meth:`ContactQC.from_stored`), or pairs counts that are not
            whole numbers of pairs.
    """
    with _open_uri(inputs[0]) as group:
        unbinned = unbinned_genome_id(_names(group["genomes"]))
    total: ContactQC | None = None
    for uri, block in zip(inputs, blocks, strict=True):
        with _open_uri(uri) as group:
            per_genome = _columns(group["genomes"], _GENOME_COUNT_COLUMNS)
            counts = (
                _columns(group["contigs"], QC_CONTIG_COLUMNS) if per_contig else None
            )
        pairs_counts = [key for key in QC_PAIRS_COUNT_KEYS if key in block]
        try:
            tally = ContactQC.from_stored(
                block, per_genome, counts, unbinned_genome=unbinned
            )
            problem = _counts_problem(block, pairs_counts, whole=True)
            if problem:
                raise ValueError(problem)
        except ValueError as exc:
            raise ValueError(f"the QC of {uri} cannot be read back ({exc})") from exc
        del per_genome, counts
        if total is None:
            total = tally
        else:
            total.add(tally)
        del tally
    assert total is not None
    return total


def _merged_extra(blocks: list[dict]) -> dict:
    """The keys of the merged block that a tally's ``summary()`` does not write.

    The pairs counts are summed; ``binsize`` is kept when every input records
    the same value (:data:`_QC_BINSIZE_KEY`); every other key is carried as
    the inputs record it, which they all do alike (:func:`_unsummable`).
    """
    first = blocks[0]
    extra = {key: value for key, value in first.items() if key not in QC_SUMMARY_KEYS}
    for key in QC_PAIRS_COUNT_KEYS:
        if key in extra:
            extra[key] = sum(block[key] for block in blocks)
    binsizes = [block.get(_QC_BINSIZE_KEY, _ABSENT) for block in blocks]
    if any(binsize != binsizes[0] for binsize in binsizes):
        extra.pop(_QC_BINSIZE_KEY, None)
    return extra


def _warn_dropped(blocks: list[object], reason: str, output: str) -> None:
    """Say why the merged map has no QC block, and how to measure one.

    ``qc --store`` measures one from the pixels.  For blocks measured on read
    positions, that would lose what only they give -- the short/long split on
    read positions, and the pairs counts -- so ``create pairs`` with every
    library's pairs is named first.
    """
    from_pixels = f"'metacooler qc --store {output}'"
    on_pairs = any(
        isinstance(block, dict) and block.get("source") == "pairs" for block in blocks
    )
    if on_pairs:
        advice = (
            "To measure the merged map's on read positions, as blocks of "
            "source 'pairs' were, give 'metacooler create pairs' every "
            f"library's pairs file, one --pairs each; {from_pixels} measures "
            "it from the pixels instead, without the short/long split on read "
            "positions or n_pairs_read and n_pairs_discarded."
        )
    else:
        advice = f"Run {from_pixels} to measure the merged map's from its pixels."
    logger.warning(
        "Dropping the inputs' QC blocks and qc_* columns, since %s. %s", reason, advice
    )


def _report_contig_counts(inputs: list[str], problems: _ContigProblems) -> None:
    """Say why a merged map with a summed QC block lacks per-contig counts."""
    if len(problems.none) == len(inputs):
        logger.info(
            "No input has per-contig QC columns (maps written by metacooler "
            "0.2.0 and earlier lack them), so the merged map has none either."
        )
        return
    reasons = []
    if problems.none:
        has = "has" if len(problems.none) == 1 else "have"
        reasons.append(
            f"{', '.join(problems.none)} {has} none (maps written by metacooler "
            "0.2.0 and earlier lack them)"
        )
    reasons += problems.partial
    if problems.stale:
        reasons.append(
            f"those of {', '.join(problems.stale)} do not add up to the block "
            "they go with (another run stored them)"
        )
    logger.warning(
        "Summing the inputs' QC blocks and /genomes columns without the "
        "per-contig counts of /contigs, since %s. Without them, a later "
        "'reassign' of the merged map recomputes its whole QC block from the "
        "pixels: source 'pixels', with the short/long split measured between "
        "bin midpoints and no n_pairs_read or n_pairs_discarded.",
        "; ".join(reasons),
    )


def _merged_qc(inputs: list[str], output: str) -> _SummedQC | None:
    """Sum the inputs' QC blocks and ``qc_*`` count columns, if they add up.

    Every count adds up over libraries of one map: the four categories -- in
    total, per genome and per contig -- and the pairs a ``create pairs`` tally
    read and discarded.  So each input's tally is rebuilt from what it stored
    and the tallies are added (:func:`_summed_tally`); every fraction and
    ratio is recomputed from the sum when it is written, by
    :meth:`ContactQC.summary` and :meth:`ContactQC.genome_table`, the calls
    that produced each input's: the result is what one tally of every input's
    contacts gives.

    When the blocks cannot be summed, a warning says why and nothing is
    returned; nor when no input has a block, which loses nothing.  When only
    the ``/contigs`` columns cannot be (:func:`_contig_count_problems`), the
    block and the ``/genomes`` columns are summed without them, as on a map
    written by metacooler 0.2.0, with a warning when an input had them.

    Args:
        inputs: Input URIs, already checked to be libraries of one map.
        output: The output path, for the messages.

    Returns:
        The summed QC, or ``None``.
    """
    blocks: list[object] = []
    problems: list[str | None] = []
    for uri in inputs:
        with _open_uri(uri) as group:
            block = _qc_block(group)
            blocks.append(block)
            problems.append(
                None if block is None else _input_problem(uri, group, block)
            )
    if all(block is None for block in blocks):
        return None
    total = None
    contig_problems = _ContigProblems([], [], [])
    reason = _unsummable(inputs, blocks, problems)
    stored = [block for block in blocks if isinstance(block, dict)]
    if reason is None:
        contig_problems = _contig_count_problems(inputs, stored)
        per_contig = not any(contig_problems)
        try:
            total = _summed_tally(inputs, stored, per_contig=per_contig)
        except ValueError as exc:
            reason = str(exc)
    if total is None:
        _warn_dropped(blocks, str(reason), output)
        return None
    if any(contig_problems):
        _report_contig_counts(inputs, contig_problems)
    logger.info(
        "Summed the QC blocks of the %d inputs (source=%r, distance_threshold=%d).",
        len(inputs),
        total.source,
        total.distance_threshold,
    )
    return _SummedQC(total, _merged_extra(stored))


def _store_summed_qc(output: str | Path, summed: _SummedQC) -> None:
    """Write the summed QC into the merged map: its block and ``qc_*`` columns.

    The names the tables carry are read here, from the merged map, rather than
    held through the pixel merge: at millions of contigs, they outweigh the
    counts several times over.
    """
    tally = summed.tally
    with h5py.File(output, "r") as f:
        genome_names = _names(f["genomes"])
        contig_names = None if tally.per_contig is None else _names(f["contigs"])
    store_qc(
        output,
        tally.summary(**summed.extra),
        tally.genome_table(genome_names),
        None if contig_names is None else tally.contig_table(contig_names),
    )


def _n_pixels(uri: str) -> int:
    """The number of pixels of the map *uri* addresses."""
    with _open_uri(uri) as f:
        return len(f["pixels"]["bin1_id"])


def _write_empty_merge(
    first: str,
    output: str,
    binsize: int | None,
    count_dtype: np.dtype,
    storage_mode: str,
) -> None:
    """Write the merge of maps that hold no pixel: *first*'s bins, no pixel.

    ``cooler.merge_coolers`` concatenates the inputs' pixel chunks, and fails
    with "No objects to concatenate" when no input has any -- libraries with no
    contact left, say. Their sum is an empty table: it is written here as
    ``merge_coolers`` writes a merge (the bins' ``chrom``, ``start`` and
    ``end``, a ``count`` column of the summed dtype, the first input's
    assembly, the inputs' storage mode), so what :func:`merge` does next
    applies unchanged.

    Args:
        first: URI of the first input, whose bins the inputs share.
        output: Destination path.
        binsize: The inputs' bin size, for :func:`write_cooler`.
        count_dtype: dtype of the merged ``count`` column.
        storage_mode: The inputs' ``storage-mode``: ``symmetric-upper`` or
            ``square``.
    """
    logger.info("No input has a pixel: the merged map has none either.")
    clr = cooler.Cooler(first)
    pixels = pd.DataFrame(
        {
            "bin1_id": np.array([], dtype=np.int64),
            "bin2_id": np.array([], dtype=np.int64),
            "count": np.array([], dtype=count_dtype),
        }
    )
    write_cooler(
        output,
        cast(pd.DataFrame, clr.bins()[["chrom", "start", "end"]][:]),
        pixels,
        binsize=binsize,
        columns=["count"],
        dtypes={"count": count_dtype},
        assembly=clr.info.get("genome-assembly"),
        symmetric_upper=storage_mode == "symmetric-upper",
        ordered=True,
    )


def merge(
    inputs: Iterable[str | Path],
    output: str | Path,
) -> None:
    """Merge multiple ``.metacool`` files with identical axes into one.

    All input files must describe the same contigs and genomes -- names,
    lengths, contig-to-genome assignment, bins per genome -- with the same
    bins and the same ``contig-ordering`` attribute.  Their metadata columns
    (completeness, coverage, ...) may differ: the output keeps the first
    input's, and each difference is reported as a warning.  So may their
    clustering rounds (``/clusters``): the output keeps the first input's,
    with a warning about each input whose rounds differ.  Pixels from all
    inputs are concatenated, re-sorted externally by
    ``cooler.merge_coolers()``, and saved to *output* (when no input has a
    pixel, *output* is written with none, as ``merge_coolers`` would write
    it).  Metacooler-specific groups (``genomes``, ``contigs``,
    ``genome_offset``, format attributes) are then copied from the first
    input file.

    The inputs' contact-quality statistics are summed too -- the
    ``metadata["qc"]`` block and the ``qc_*`` columns of ``/genomes`` and
    ``/contigs`` -- when every input has a block, and every block was counted
    alike: the same ``version``, ``source`` and ``distance_threshold``, and the
    same keys (``binsize`` aside: it is kept when every input records the same
    value).  Every count adds up over libraries, and every fraction and ratio
    is recomputed from the sums, so the merged QC is the one a single tally of
    every input's contacts gives: for blocks measured on read positions, the
    one ``create pairs`` stores when given every input's pairs at once; for
    blocks measured on pixels, the one ``qc --store`` measures on the merged
    map.  Otherwise -- or when an input's ``/genomes`` columns do not add up
    to its block -- the output has no QC block, and a warning says why.  An
    input written by metacooler 0.2.0 has no ``/contigs`` QC columns: its
    block and ``/genomes`` columns are summed, and the output has none in
    ``/contigs`` either.  So it is when an input has only one of the two, one
    that is not numbers, or two that no longer add up to its block, and a
    warning says which.  Of the inputs' other metadata keys, this function
    carries none; ``metacooler merge`` does.

    Args:
        inputs: Iterable of paths to the input ``.metacool`` files, or
            ``file.mmetacool::entries/<name>`` URIs.  At least two are
            required.
        output: Destination path for the merged ``.metacool`` file.

    Raises:
        ValueError: If fewer than two inputs are provided, if the
            ``contigs``, ``genomes`` or ``bins`` tables differ across files,
            if the resolutions differ, or if the ``contig-ordering`` or
            ``storage-mode`` attribute is inconsistent.
        MetaCoolerFormatError: If any input file fails the metacooler
            format validation check.
    """
    inputs = [str(p) for p in inputs]
    if len(inputs) < 2:
        raise ValueError("At least two input files required for merge")
    _refuse_repeated_inputs(inputs)
    from metacooler.io.staging import writing_output

    with writing_output(output):
        _merge_into(inputs, output)


def _merge_into(inputs: list[str], output: str | Path) -> None:
    """:func:`merge`, its output checked (:func:`writing_output`)."""
    # cooler writes the output in "w" mode; refuse to truncate an input with it.
    check_output_path(output, force=True, inputs=inputs, defer_removal=True)
    # The inputs are read by path, by cooler too: checked unchanged at the end.
    states = input_states(inputs)

    # Validate and compare axes
    first = inputs[0]
    validate(first)

    with _open_uri(first) as f0:
        contigs_hash = _hash_group(f0, "contigs")
        genomes_hash = _hash_group(f0, "genomes")
        metadata_digests = _metadata_digests(f0)
        clusters_digest_first = clusters_digest(f0)
        # ``contig-ordering`` may be absent or explicitly null; normalise both
        # sides of the comparison so a missing attribute does not later become
        # ``f.attrs[...] = None``, which h5py rejects with a raw TypeError.
        # ``_attr`` decodes the attributes another writer stores as bytes, as
        # ``validate`` accepts them, whatever cooler's version: the ``info()``
        # of cooler 0.10 returns them raw.
        ordering = str(_attr(f0, "contig-ordering") or "genome-then-alpha")
        storage_mode = str(_attr(f0, "storage-mode") or "symmetric-upper")
        n_bins, bins_hash = _bins_signature(f0)
        binsize = _binsize(f0)

    for inp in inputs[1:]:
        validate(inp)
        with _open_uri(inp) as f:
            # Order matters: report the ROOT cause, not a downstream symptom.
            # ``/genomes`` carries ``n_bins``, so its hash changes whenever the
            # resolution does — checking it first made every resolution
            # mismatch report as "Genomes table mismatch", which sends the user
            # looking at the wrong table. Contigs (assembly identity) first,
            # then resolution, then the bins table, and only then genomes.
            if _hash_group(f, "contigs") != contigs_hash:
                raise ValueError(f"Contigs table mismatch: {first} vs {inp}")
            # Merging files binned at different resolutions would concatenate
            # pixels that address unrelated bin IDs, producing a valid-looking
            # but meaningless map, so refuse before any output is written.
            other_binsize = _binsize(f)
            if other_binsize != binsize:
                raise ValueError(
                    f"Resolution mismatch: {first} has bin-size {binsize!r} but "
                    f"{inp} has bin-size {other_binsize!r}. All inputs must be "
                    "binned identically."
                )
            other_n_bins, other_bins_hash = _bins_signature(f)
            if other_n_bins != n_bins:
                raise ValueError(
                    f"Bins table mismatch: {first} has {n_bins} bins but {inp} "
                    f"has {other_n_bins}."
                )
            if other_bins_hash != bins_hash:
                raise ValueError(
                    f"Bins table mismatch: {first} and {inp} have the same number "
                    "of bins but different bin coordinates."
                )
            if _hash_group(f, "genomes") != genomes_hash:
                raise ValueError(f"Genomes table mismatch: {first} vs {inp}")
            if str(_attr(f, "contig-ordering") or "genome-then-alpha") != ordering:
                raise ValueError(f"contig-ordering mismatch: {first} vs {inp}")
            other_mode = str(_attr(f, "storage-mode") or "symmetric-upper")
            if other_mode != storage_mode:
                raise ValueError(
                    f"storage-mode mismatch: {first} stores {storage_mode!r} but "
                    f"{inp} stores {other_mode!r}; a merge needs one."
                )
            if clusters_digest(f) != clusters_digest_first:
                logger.warning(
                    "%s has other clustering rounds than %s (or none): the merged "
                    "map keeps those of %s.",
                    inp,
                    first,
                    first,
                )
            differences = _metadata_differences(metadata_digests, _metadata_digests(f))
            if differences:
                logger.warning(
                    "%s has other metadata columns than %s: %s. The merged map "
                    "keeps those of %s.",
                    inp,
                    first,
                    ", ".join(differences),
                    first,
                )

    # Before any pixel is merged: reading every input's QC back is cheap, and
    # whether it adds up is better said before the long step than after it.
    # Only the summed counts are held through it.
    summed_qc = _merged_qc(inputs, str(output))

    # Use cooler's merge logic on the pixel tables.  ``dtypes`` is mandatory:
    # without it the summed counts inherit the int32 of the inputs and wrap.
    # merge_coolers calls cooler.create internally, so it needs the same shims
    # as every other writer (see write_cooler): without them its per-contig
    # costs made merge alone take minutes at 150 k contigs. They also hand
    # cooler the inputs' known bin size, which it cannot infer once every
    # contig is a single bin.
    count_dtype = _merged_count_dtype(inputs)
    if not any(_n_pixels(uri) for uri in inputs):
        _write_empty_merge(first, str(output), binsize, count_dtype, storage_mode)
    else:
        with cooler_writer_shims(binsize):
            cooler.merge_coolers(
                str(output),
                inputs,
                mergebuf=int(2e7),
                dtypes={"count": count_dtype},
            )

    # Post-process: add back metacooler-specific groups
    with _open_uri(first) as fsrc, h5py.File(output, "a") as fdst:
        # Copy genomes, minus the first input's QC columns: the merged pixel
        # table is the sum of every input, so counts copied from one of them
        # would understate it while looking authoritative. The summed ones, if
        # the inputs' add up, are stored below.
        if "genomes" in fdst:
            del fdst["genomes"]
        fdst.copy(fsrc["genomes"], "genomes")
        for key in [k for k in fdst["genomes"] if is_qc_column("genomes", k)]:
            del fdst["genomes"][key]

        # Copy contigs (overwrites the chroms-compatible one), minus their two
        # QC columns for the same reason as the genomes' -- only those: any
        # other qc_* column there is a user's contig metadata.
        if "contigs" in fdst:
            del fdst["contigs"]
        fdst.copy(fsrc["contigs"], "contigs")
        for key in [k for k in fdst["contigs"] if is_qc_column("contigs", k)]:
            del fdst["contigs"][key]

        # Soft link chroms -> contigs
        if "chroms" in fdst:
            del fdst["chroms"]
        fdst["chroms"] = h5py.SoftLink("contigs")

        # Copy genome_offset
        if "genome_offset" in fdst["indexes"]:
            del fdst["indexes"]["genome_offset"]
        fdst["indexes"].copy(fsrc["indexes"]["genome_offset"], "genome_offset")

        # The clustering rounds, as metadata: the first input's.
        copy_clusters(fsrc, fdst)

        # Strip cooler's enum dtype from bins/chrom — cooler.merge_coolers()
        # always writes bins/chrom as an HDF5 enum type, which is inconsistent
        # with the metacooler schema that stores it as plain int32.
        chrom_data = fdst["bins"]["chrom"][:]
        del fdst["bins"]["chrom"]
        write_column(fdst["bins"], "chrom", chrom_data, dtype=np.int32)

        # Set attrs
        fdst.attrs["format"] = COOLER_FORMAT
        fdst.attrs["format-url"] = COOLER_FORMAT_URL
        fdst.attrs["metacooler-format"] = METACOOLER_FORMAT
        fdst.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        fdst.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        fdst.attrs["generated-by"] = f"metacooler-{__version__}"
        fdst.attrs["ngenomes"] = len(fsrc["genomes"]["name"])
        fdst.attrs["contig-ordering"] = ordering
        # The inputs were checked to share one bin grid, so its description
        # carries over verbatim rather than being re-inferred from the table.
        for key in ("bin-size", "bin-type"):
            if key in fsrc.attrs:
                fdst.attrs[key] = fsrc.attrs[key]

        # Same reasoning as the QC columns above: whatever block an input
        # carried describes that input alone (cooler carries none today).
        raw = fdst.attrs.get("metadata")
        if raw is not None:
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
            if QC_METADATA_KEY in meta:
                del meta[QC_METADATA_KEY]
                fdst.attrs["metadata"] = json.dumps(meta)

    if summed_qc is not None:
        _store_summed_qc(output, summed_qc)

    try:
        check_inputs_unchanged(states)
    except MapChangedError:
        Path(output).unlink(missing_ok=True)
        raise
    validate(output)
