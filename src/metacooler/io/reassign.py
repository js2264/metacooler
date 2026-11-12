"""Reassign the contigs of an existing metacool file to different genomes (MAGs).

Re-binning a metagenome is iterative: a first MetaTOR/vamb/MetaBAT pass produces
one contig→MAG table, curation produces another, and the Hi-C matrix itself does
not change between them — only the labels do.  This module rewrites those labels
without re-binning a single pair.

Two things make that less trivial than "overwrite ``/contigs/genome_id``":

- ``/genomes`` is *derived*.  Its rows, their order (which defines the
  ``genome_id`` integers), and its ``length``/``n_contigs``/``n_bins`` columns
  all follow from the assignment, so the whole table has to be rebuilt and every
  genome renumbered.
- Under a genome-grouped ``contig-ordering`` the contigs are *physically* laid
  out in ascending ``genome_id`` order, which is what makes ``genome_offset``
  delimit a genome's contig block and a genome's bins a contiguous range.  Move
  one contig between two non-adjacent MAGs and that layout no longer holds:
  :func:`metacooler.schema.validate` rejects the file, and it is right to.

Restoring the layout means permuting the contigs, hence the bins, hence every
pixel's bin ids — an O(nnz) rewrite.  So :func:`reassign` runs one of two paths:

- **relabel** (O(n_contigs)): the physical order is left alone and only the
  labels, the derived ``/genomes`` table, ``genome_offset`` and the
  ``contig-ordering`` attribute are rewritten.  Cheap enough to run in place.
- **regroup** (O(nnz)): the contigs are re-ordered by the file's own ordering
  strategy, the bins are permuted with them, and the pixel table is remapped and
  re-sorted, exactly as ``create cool`` does.

``layout="auto"`` picks between them by what the source file guaranteed: a
genome-grouped file is regrouped so it stays one, and a file that never grouped
its contigs (``fasta``/``alpha``/``cool``) is only relabelled, since there is no
invariant left to preserve.  Either way the *contacts* are untouched: reassign
never changes a count, only which genome a contig is credited to.
"""

from __future__ import annotations

import json
import logging
import shutil
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple, cast

import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import (
    _finalize_metacool_file,
    check_ascii_names,
    check_output_path,
    resolve_count_dtype,
    write_column,
    write_cooler,
)
from metacooler.io.create import (
    _build_contig_order,
    _build_contigs_table,
    _build_genomes_table,
    _compute_all_offsets,
    _load_extra_genomes_table,
    _parse_contig_map,
    _remap_bins,
    _staged_output,
    _stream_remapped_pixels,
)
from metacooler.qc import (
    QC_METADATA_KEY,
    strip_qc_columns,
    strip_qc_metadata,
)
from metacooler.schema import (
    CUSTOM_ORDERING,
    DEFAULT_CHUNKSIZE,
    GENOME_GROUPED_ORDERINGS,
    LAYOUT_CHOICES,
    ContigOrderStrategy,
    Layout,
    validate,
)

if TYPE_CHECKING:
    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

UNBINNED = "__unbinned__"

# Orderings whose meaning does not depend on the contig→genome assignment, so a
# reassignment cannot invalidate the label.
_ASSIGNMENT_INDEPENDENT_ORDERINGS = frozenset(
    {"alpha", "fasta", "cool", "pairs-header"}
)

# How many names to spell out in a warning before truncating.
_MAX_REPORTED_NAMES = 10


class ReassignReport(NamedTuple):
    """What a reassignment did, as counted before anything was written.

    Attributes:
        moved: Contigs whose genome changed.
        unchanged: Contigs whose genome is the same as before.
        unknown: Names in the table that are not contigs of this file.
        created: Genome names that did not exist in the source file.
        emptied: Genome names that lost every contig they had.
    """

    moved: list[str]
    unchanged: list[str]
    unknown: list[str]
    created: list[str]
    emptied: list[str]


def _names(values: list[str]) -> str:
    """Render a name list for a log message, truncated to a readable length."""
    shown = ", ".join(values[:_MAX_REPORTED_NAMES])
    return shown + (", ..." if len(values) > _MAX_REPORTED_NAMES else "")


def _apply_reassignment(
    before: dict[str, str],
    table: dict[str, str],
    missing: Literal["keep", "unbinned"],
) -> tuple[dict[str, str], ReassignReport]:
    """Merge a contig→genome table into the file's current assignment.

    Args:
        before: Current ``{contig name: genome name}`` mapping of the file.
        table: Parsed reassignment table, genome sentinels already normalised.
        missing: What to do with a contig of the file that the table does not
            mention — ``"keep"`` leaves its current genome, ``"unbinned"``
            moves it to ``__unbinned__``.

    Returns:
        ``(after, report)``: the new mapping, over exactly the same contigs, and
        a :class:`ReassignReport` describing the change.

    Raises:
        ValueError: If the table names no contig of this file at all, which
            means the wrong table (or the wrong ``--contigs-format``) rather
            than a reassignment that happens to change nothing.
    """
    if not table:
        raise ValueError("The reassignment table is empty; there is nothing to apply.")

    known = set(before)
    unknown = sorted(name for name in table if name not in known)
    if len(unknown) == len(table):
        raise ValueError(
            f"None of the {len(table)} contig name(s) in the reassignment table "
            "match a contig of this file. Check that the table describes this "
            "assembly and that --contigs-format names the right columns "
            f"(first table names: {_names(unknown)})."
        )

    after = dict(before)
    if missing == "unbinned":
        for contig in before:
            if contig not in table:
                after[contig] = UNBINNED
    for contig, genome in table.items():
        if contig in known:
            after[contig] = genome

    moved = [c for c in before if after[c] != before[c]]
    unchanged = [c for c in before if after[c] == before[c]]
    old_genomes = set(before.values())
    new_genomes = set(after.values())
    report = ReassignReport(
        moved=moved,
        unchanged=unchanged,
        unknown=unknown,
        created=sorted(new_genomes - old_genomes),
        emptied=sorted(old_genomes - new_genomes),
    )
    return after, report


def _log_report(report: ReassignReport, n_contigs: int) -> None:
    """Log what the reassignment changed, before any file is written."""
    logger.info(
        "Reassignment: %d of %d contig(s) change genome, %d unchanged.",
        len(report.moved),
        n_contigs,
        len(report.unchanged),
    )
    if report.created:
        logger.info("%d new genome(s): %s", len(report.created), _names(report.created))
    if report.emptied:
        logger.info(
            "%d genome(s) lose their last contig and disappear from /genomes: %s",
            len(report.emptied),
            _names(report.emptied),
        )
    if report.unknown:
        logger.warning(
            "%d name(s) in the reassignment table are not contigs of this file "
            "and were ignored: %s",
            len(report.unknown),
            _names(report.unknown),
        )
    if not report.moved:
        logger.warning(
            "No contig changes genome; the output differs from the input only "
            "in its provenance metadata."
        )


def _resolve_layout(layout: str, source_ordering: str | None) -> str:
    """Return the contig-ordering strategy to lay the output out with.

    Args:
        layout: One of :data:`LAYOUT_CHOICES`.
        source_ordering: The source file's ``contig-ordering`` attribute.

    Returns:
        ``"keep"`` to leave the physical contig order alone, or the name of a
        :func:`metacooler.io.create._build_contig_order` strategy to re-derive
        it with.

    Raises:
        ValueError: If *layout* is not one of :data:`LAYOUT_CHOICES`.
    """
    if layout not in LAYOUT_CHOICES:
        raise ValueError(
            f"layout must be one of {list(LAYOUT_CHOICES)}, got {layout!r}"
        )
    if layout != "auto":
        return layout
    # Preserve whatever the source file guaranteed: regroup a genome-grouped
    # file (so genome_offset and genome_extent keep working), and leave a file
    # that never grouped its contigs exactly where it is.
    if source_ordering in GENOME_GROUPED_ORDERINGS:
        return str(source_ordering)
    return "keep"


def _scattered_genomes(
    contig_order: list[str], assignment: dict[str, str]
) -> list[str]:
    """Return the genomes whose contigs do not form one contiguous run.

    Args:
        contig_order: The physical contig order.
        assignment: Contig → genome mapping over the same contigs.

    Returns:
        Sorted names of the genomes split into two or more runs of rows.
    """
    labels = [assignment[contig] for contig in contig_order]
    runs = Counter(
        label for i, label in enumerate(labels) if i == 0 or label != labels[i - 1]
    )
    return sorted(genome for genome, n_runs in runs.items() if n_runs > 1)


def _kept_ordering(
    source_ordering: str | None,
    contig_order: list[str],
    contig_lengths: dict[str, int],
    before: dict[str, str],
    assignment: dict[str, str],
) -> str:
    """Return the ``contig-ordering`` value for a relabelled (un-permuted) file.

    A label like ``genome-then-alpha`` is a claim about the row order, and the
    reassignment may have falsified it.  ``alpha``/``fasta``/``cool`` do not
    mention genomes, so they survive untouched; a genome-grouped label survives
    only if the order it *would* produce under the new assignment is the order
    the file already has.  Otherwise the file is labelled
    :data:`CUSTOM_ORDERING` — which also stops ``validate()`` from enforcing a
    grouping the rows no longer have.

    A label that does not promise grouping can still describe a file whose
    genomes happen to be contiguous — every ``collect``ed map is, under
    ``cool`` — and relabelling may scatter them.  Nothing on disk becomes
    false, so the label stays, but the loss is warned about by name: it is
    exactly what ``layout="auto"`` does not protect against for such a file.

    Args:
        source_ordering: The source file's ``contig-ordering`` attribute.
        contig_order: The (unchanged) physical contig order.
        contig_lengths: Contig name → length, for the length-based strategy.
        before: The source's contig → genome mapping.
        assignment: The new contig → genome mapping.

    Returns:
        The ``contig-ordering`` value to write.
    """
    if source_ordering not in GENOME_GROUPED_ORDERINGS:
        scattered = _scattered_genomes(contig_order, assignment)
        if scattered and not _scattered_genomes(contig_order, before):
            logger.warning(
                "Every genome of the source occupied one contiguous block of "
                "contigs, but contig-ordering=%r does not promise it, so the "
                "rows were left where they are and %d genome(s) are now "
                "scattered: %s. Anything that needs one bin range per genome "
                "('show'/'ps' region strings, genome_extent) fails for those; "
                "selecting by contig membership is unaffected. Pass "
                "layout='genome-then-alpha' or 'genome-then-length' to re-group "
                "the file instead.",
                source_ordering,
                len(scattered),
                _names(scattered),
            )
    if source_ordering in _ASSIGNMENT_INDEPENDENT_ORDERINGS:
        return str(source_ordering)
    if source_ordering in GENOME_GROUPED_ORDERINGS:
        would_be = _build_contig_order(
            contig_lengths,
            assignment,
            ordering=cast("ContigOrderStrategy", source_ordering),
        )
        if would_be == contig_order:
            return str(source_ordering)
        logger.warning(
            "The contigs are no longer laid out as %r under the new assignment, "
            "so the file is recorded as contig-ordering=%r. A genome's contigs "
            "may now be scattered, which makes genome_offset non-delimiting and "
            "'show'/'ps' region strings that name a genome fail; selecting by "
            "contig membership ('dump --genomes', 'extract --genomes') is "
            "unaffected. Pass layout='%s' to re-group the file instead.",
            source_ordering,
            CUSTOM_ORDERING,
            source_ordering,
        )
        return CUSTOM_ORDERING
    return source_ordering or CUSTOM_ORDERING


def _carry_genome_metadata(
    old_genomes: pd.DataFrame, keep: set[str], restated: set[str]
) -> pd.DataFrame | None:
    """Extract the user metadata columns of ``/genomes`` for re-merging.

    Only the columns metacooler does not derive are carried: ``length``,
    ``n_contigs`` and ``n_bins`` are recomputed from the new assignment, and the
    ``qc_*`` columns are per-genome contact tallies that the reassignment
    invalidates outright.

    The frame is restricted to genomes that still exist, so
    ``_report_genome_table_divergence`` only reports what is genuinely new
    (a MAG with no metadata) rather than re-reporting the genomes the caller has
    already been told disappeared.

    Args:
        old_genomes: The source file's ``/genomes`` table.
        keep: Genome names present in the new assignment.
        restated: Genome names that gained or lost a contig.  Their carried
            values were computed over a different set of contigs, so they are
            named in a warning; a genome nobody touched keeps metadata that is
            still exactly true, and is not worth a line of output.

    Returns:
        A frame keyed by ``genome_id`` (the genome name) plus the carried
        columns, or ``None`` when the source carried no user metadata.
    """
    derived = {"name", "length", "n_contigs", "n_bins"}
    extra = [c for c in strip_qc_columns(old_genomes).columns if c not in derived]
    if not extra:
        return None
    frame = old_genomes.loc[
        old_genomes["name"].astype(str).isin(keep), ["name", *extra]
    ]
    logger.info(
        "Carrying %d metadata column(s) over from the source /genomes table: %s.",
        len(extra),
        ", ".join(map(str, extra)),
    )
    stale = sorted(restated & set(frame["name"].astype(str)))
    if stale:
        logger.warning(
            "%d genome(s) keep those columns although their contig set changed: "
            "%s. A value derived from a genome's contigs (completeness, "
            "contamination, …) no longer describes the genome it is attached "
            "to. Pass a --genomes table to replace them.",
            len(stale),
            _names(stale),
        )
    return frame.rename(columns={"name": "genome_id"}).reset_index(drop=True)


def _carry_contig_metadata(
    old_contigs: pd.DataFrame, assignment: dict[str, str]
) -> pd.DataFrame | None:
    """Extract the user metadata columns of ``/contigs`` for re-merging.

    Contigs themselves are untouched by a reassignment — same names, same
    lengths, same bins — so every extra column stays valid and is carried
    verbatim.

    Args:
        old_contigs: The source file's ``/contigs`` table.
        assignment: The new contig → genome mapping, used to shape the frame
            like a parsed ``--contigs`` table.

    Returns:
        A frame with ``contig_name``, ``genome_id`` and the carried columns, or
        ``None`` when the source carried no extra column.
    """
    derived = {"name", "length", "genome_id"}
    extra = [c for c in old_contigs.columns if c not in derived]
    if not extra:
        return None
    names = old_contigs["name"].astype(str)
    frame = pd.DataFrame({"contig_name": names, "genome_id": names.map(assignment)})
    for column in extra:
        frame[column] = old_contigs[column].to_numpy()
    logger.info(
        "Carrying %d metadata column(s) over from the source /contigs table: %s.",
        len(extra),
        ", ".join(map(str, extra)),
    )
    return frame


def _build_metadata(
    source_metadata: dict,
    extra: dict | None,
    *,
    source_uri: str,
    table: str,
    layout: str,
    n_moved: int,
) -> tuple[dict, bool]:
    """Assemble the output's ``metadata`` attribute.

    The source metadata is carried over minus its QC block: three of the four
    QC categories are defined by the contig→genome assignment, so a reassignment
    invalidates the block as a whole rather than only its per-genome columns.

    Args:
        source_metadata: The source file's ``metadata`` attribute.
        extra: Keys to merge in (the CLI passes ``commands`` here, which is
            *appended* to the carried provenance list rather than replacing it).
        source_uri: URI of the file being reassigned, for the history entry.
        table: Path of the reassignment table, for the history entry.
        layout: Resolved layout strategy, for the history entry.
        n_moved: Number of contigs that changed genome.

    Returns:
        ``(metadata, had_qc)`` — the metadata to write, and whether the source
        carried a QC block (so the caller can offer to recompute it).
    """
    had_qc = QC_METADATA_KEY in source_metadata
    metadata = strip_qc_metadata(source_metadata)
    commands = list(metadata.get("commands", []))
    for key, value in (extra or {}).items():
        if key == "commands":
            commands.extend(value if isinstance(value, list) else [value])
        else:
            metadata[key] = value
    history = list(metadata.get("reassign", []))
    history.append(
        {
            "source": source_uri,
            "table": table,
            "layout": layout,
            "n_contigs_reassigned": int(n_moved),
        }
    )
    metadata["reassign"] = history
    if commands:
        metadata["commands"] = commands
    return metadata, had_qc


def _materialise_copy(mc: MetaCooler, target: Path) -> None:
    """Copy the source map to *target* as a standalone ``.metacool`` file.

    A whole-file input is copied byte for byte, which is both the fastest way
    to get an editable duplicate and the only one that keeps everything the
    relabel path does not rewrite — balancing weights, extra pixel columns,
    HDF5 chunking and filters — exactly as it was.  An mmetacool *entry* is
    copied group by group into the root of a new file, so the output of
    ``reassign`` is always a single map rather than a collection.

    Args:
        mc: The open source file.
        target: Path to write; must not exist.
    """
    group = mc._group_path.strip("/")
    if not group:
        shutil.copy2(mc.store, target)
        return
    with h5py.File(mc.store, "r") as src, h5py.File(target, "w") as dst:
        node = src[group]
        for key in node:
            link = node.get(key, getlink=True)
            if isinstance(link, h5py.SoftLink):
                # Copying through the link would hard-copy /contigs a second
                # time under the name /chroms; the format wants a soft link.
                dst[key] = h5py.SoftLink(link.path)
            else:
                src.copy(node[key], dst, name=key)
        for key, value in node.attrs.items():
            dst.attrs[key] = value


def _write_relabelled(
    path: Path,
    group_path: str,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    genome_offset: np.ndarray,
    ordering: str,
    metadata: dict,
) -> None:
    """Rewrite the genome-dependent parts of a map, in place, leaving bins alone.

    Everything the reassignment cannot touch — ``/bins``, ``/pixels``,
    ``chrom_offset``, ``bin1_offset``, every extra ``/contigs`` column — is left
    exactly where it is, because the contig rows keep their positions.  Only
    ``/genomes`` (whose rows are derived), ``/contigs/genome_id``,
    ``indexes/genome_offset`` and the affected attributes are written.

    Args:
        path: File to modify.
        group_path: HDF5 group of the map (``"/"`` for a plain metacool).
        genomes: The rebuilt ``/genomes`` table.
        contigs: The rebuilt ``/contigs`` table, in the file's existing row
            order.
        genome_offset: CSR offsets over ``/contigs``, length ``n_genomes + 1``.
        ordering: Value for the ``contig-ordering`` attribute.
        metadata: Value for the ``metadata`` attribute.

    Raises:
        ValueError: If *contigs* does not describe the same contigs, in the same
            order, as the file on disk — the precondition that makes writing
            only ``genome_id`` equivalent to rewriting the table.
    """
    check_ascii_names(genomes, contigs)
    with h5py.File(path, "r+") as f:
        key = group_path.strip("/")
        grp = f[key] if key else f

        on_disk = [
            name.decode() if isinstance(name, bytes) else str(name)
            for name in grp["contigs"]["name"][:]
        ]
        if on_disk != contigs["name"].astype(str).tolist():
            raise ValueError(
                "The rebuilt contigs table does not match the file's contig "
                "rows; refusing to write genome ids that would not line up."
            )

        del grp["genomes"]
        genomes_group = grp.create_group("genomes", track_order=True)
        for column in genomes.columns:
            write_column(genomes_group, column, genomes[column].values)

        del grp["contigs"]["genome_id"]
        write_column(
            grp["contigs"], "genome_id", contigs["genome_id"].to_numpy(np.int32)
        )

        indexes = grp["indexes"]
        if "genome_offset" in indexes:
            del indexes["genome_offset"]
        write_column(indexes, "genome_offset", genome_offset, dtype=np.int64)

        grp.attrs["ngenomes"] = len(genomes)
        grp.attrs["contig-ordering"] = ordering
        grp.attrs["metadata"] = json.dumps(metadata)


def _rewrite_regrouped(
    mc: MetaCooler,
    target: Path,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    contig_order: list[str],
    ordering: str,
    metadata: dict,
    chunksize: int,
) -> None:
    """Write a new file with the contigs re-ordered and the pixels remapped.

    The permutation is a *relabelling of bin ids*, not a change of contacts:
    every pixel keeps its count and its two bins, which merely acquire new
    indices (and get swapped back into the upper triangle when the permutation
    inverts a pair).  Extra ``/bins`` columns ride along with their bins, so a
    balanced file stays balanced.

    Args:
        mc: The open source file.
        target: Path to write.
        genomes: The rebuilt ``/genomes`` table.
        contigs: The rebuilt ``/contigs`` table, in *contig_order*.
        contig_order: The new physical contig order.
        ordering: Value for the ``contig-ordering`` attribute.
        metadata: Value for the ``metadata`` attribute.
        chunksize: Source pixels read and remapped per chunk.
    """
    bins_remapped, old_to_new_bin = _remap_bins(
        mc.bins(), contig_order, keep_extra_columns=True
    )
    chrom_offset, genome_offset = _compute_all_offsets(bins_remapped, contigs, genomes)
    check_ascii_names(genomes, contigs)

    source = str(mc.store) if mc._group_path == "/" else mc.uri
    count_dtype = resolve_count_dtype(
        mc._clr.pixels()[0:1]["count"].dtype, context="reassign"
    )
    binsize = mc.binsize
    binsize = None if binsize is None else int(binsize)

    write_cooler(
        target,
        bins_remapped,
        _stream_remapped_pixels(
            source,
            int(mc.info["nnz"]),
            old_to_new_bin,
            nproc=1,
            count_dtype=count_dtype,
            chunksize=chunksize,
        ),
        binsize=binsize,
        # Neither the chunks nor their sequence are ordered: a permutation
        # moves a late bin pair anywhere. write_cooler sorts each chunk.
        ordered=False,
        symmetric_upper=True,
        assembly=mc.info.get("assembly") or "unknown",
        metadata=metadata,
        dtypes={"count": count_dtype},
    )
    _finalize_metacool_file(
        target,
        genomes,
        contigs,
        chrom_offset,
        genome_offset,
        contig_ordering=ordering,
        assembly_name=mc.info.get("assembly"),
        binsize=binsize,
    )


def _recompute_qc(uri: str, source_metadata: dict) -> None:
    """Recompute the contact-quality block of a file whose genomes just changed.

    Args:
        uri: The reassigned map to annotate, in place.
        source_metadata: The *source* file's metadata, read for the distance
            threshold its block was measured at, so the new numbers stay
            comparable with the old ones.
    """
    from metacooler.ops.qc import qc as compute_qc
    from metacooler.qc import DEFAULT_DISTANCE_THRESHOLD

    block = source_metadata.get(QC_METADATA_KEY) or {}
    threshold = block.get("distance_threshold", DEFAULT_DISTANCE_THRESHOLD)
    try:
        threshold = int(threshold)
    except (TypeError, ValueError):
        threshold = DEFAULT_DISTANCE_THRESHOLD
    logger.info("Recomputing the QC block for the new genome assignment")
    compute_qc(uri, distance_threshold=threshold, store=True)


def reassign(
    uri: str | Path,
    contig_genome_map: str | Path | dict[str, str],
    output: str | Path | None = None,
    *,
    contig_table_format: str = "contigs",
    genomes_table: str | Path | None = None,
    missing: Literal["keep", "unbinned"] = "keep",
    layout: Layout = "auto",
    inplace: bool = False,
    metadata: dict | None = None,
    qc: bool = True,
    chunksize: int = DEFAULT_CHUNKSIZE,
    force: bool = False,
) -> None:
    """Reassign the contigs of a metacool file to different genomes.

    The contact matrix is never recomputed: contigs, their lengths and their
    bins are exactly those of the source file, and every pixel keeps its count.
    What changes is which genome each contig belongs to, and everything derived
    from that — the ``/genomes`` table, the ``genome_id`` column, and (on the
    regroup path) the order the contigs are stored in.

    Args:
        uri: Source ``.metacool`` path or ``file.mmetacool::entries/name`` URI.
        contig_genome_map: Reassignment table (path), or a ready-made
            ``{contig: genome}`` mapping.  It need not cover every contig — see
            *missing*.
        output: Destination ``.metacool`` path.  Required unless *inplace*.
        contig_table_format: Layout of the reassignment table, as accepted by
            :func:`metacooler.io.contigs.parse_contigs_table` (``contigs``,
            ``contig_genome``, ``metator``, or an ``I,J[,noheader]`` column
            spec).  Extra columns of the table are ignored; only the two name
            columns are read.
        genomes_table: Optional ``--genomes`` metadata table for the *new*
            genomes.  When given it replaces the metadata columns carried over
            from the source ``/genomes`` table rather than merging with them.
        missing: What to do with a contig the table does not mention.
            ``"keep"`` (default) leaves its current genome, which makes a short
            table an incremental edit; ``"unbinned"`` moves it to
            ``__unbinned__``, which makes the table authoritative — the right
            choice when it is the output of a fresh binning run.  Neither drops
            a contig: removing contigs (and their pixels) is what
            :func:`metacooler.io.extract.extract` is for.
        layout: How to lay the output's contigs out.  ``"auto"`` (default)
            re-groups a file whose ``contig-ordering`` is genome-grouped and
            leaves any other file's order alone; ``"keep"`` never re-orders;
            ``"genome-then-alpha"``/``"genome-then-length"`` force a re-grouping
            with that strategy.  Re-ordering costs a full pixel rewrite and is
            skipped automatically when the strategy reproduces the order the
            file already has.
        inplace: Modify *uri* itself instead of writing a new file.  Only
            possible when no re-ordering is required, since permuting the bins
            means rewriting the pixel table.
        metadata: Extra keys to merge into the output's ``metadata`` attribute.
            A ``commands`` list is appended to the carried provenance rather
            than replacing it.
        qc: Recompute the contact-quality block when the source carries one.
            With ``False`` the stale block is dropped and not replaced.
        chunksize: Source pixels read and remapped per chunk, on the re-ordering
            path only.
        force: Replace *output* if it exists.  The existing file is only
            replaced once the new one has been written and validated, so a
            failed run leaves it untouched.

    Raises:
        ValueError: If the output arguments are inconsistent, if *output* is
            one of the inputs, if *inplace* is asked for a reassignment that
            requires re-ordering, if re-ordering is asked of a map whose
            storage mode is not ``symmetric-upper``, or if the table names no
            contig of this file.
        FileExistsError: If *output* exists and *force* is not set.
    """
    from metacooler.api import MetaCooler

    if inplace and output is not None:
        raise ValueError(
            "inplace rewrites the input file, so no output path is expected; "
            "pass one or the other."
        )
    if not inplace and output is None:
        raise ValueError(
            "An output path is required. Pass output=... to write a new file, "
            "or inplace=True to rewrite the input."
        )
    if output is not None:
        # Before anything is read: an output that is one of the inputs is
        # refused even under force, and an existing output is left in place
        # for _staged_output to replace once the new file has validated.
        check_output_path(
            output,
            force,
            inputs=[
                uri,
                genomes_table,
                None if isinstance(contig_genome_map, dict) else contig_genome_map,
            ],
            defer_removal=True,
        )

    mc = MetaCooler(uri)
    old_contigs = mc.contigs
    old_genomes = mc.genomes
    contig_names = old_contigs["name"].astype(str).tolist()
    genome_names = old_genomes["name"].astype(str).to_numpy()
    before = dict(
        zip(
            contig_names,
            genome_names[old_contigs["genome_id"].to_numpy(np.int64)],
            strict=True,
        )
    )

    source_table: str | Path | dict[str, str] = (
        {str(k): str(v) for k, v in contig_genome_map.items()}
        if isinstance(contig_genome_map, dict)
        else contig_genome_map
    )
    table, _ = _parse_contig_map(source_table, contig_table_format)
    assignment, report = _apply_reassignment(before, table, missing)
    _log_report(report, len(contig_names))

    contig_lengths = {
        name: int(length)
        for name, length in zip(
            contig_names, old_contigs["length"].to_numpy(np.int64), strict=True
        )
    }
    source_ordering = mc.info.get("ordering")
    strategy = _resolve_layout(layout, source_ordering)
    if strategy == "keep":
        contig_order = contig_names
    else:
        contig_order = _build_contig_order(
            contig_lengths, assignment, ordering=cast("ContigOrderStrategy", strategy)
        )
    reorder = contig_order != contig_names
    # The label follows the *strategy*, not the amount of movement: a requested
    # layout that happens to reproduce the existing order still describes the
    # rows correctly, and writing the source's old label there would understate
    # what the file now guarantees.  Only "keep" has to interrogate the rows.
    ordering = (
        _kept_ordering(
            source_ordering, contig_order, contig_lengths, before, assignment
        )
        if strategy == "keep"
        else strategy
    )

    storage_mode = mc.info.get("storage-mode", "symmetric-upper")
    if reorder and storage_mode != "symmetric-upper":
        # Remapping folds every pixel into the upper triangle, which is only
        # right for a symmetric matrix -- and square storage is chosen precisely
        # when one is not guaranteed, so a pixel and its mirror would be summed.
        # Relabelling never touches a pixel, so it stays available.
        raise ValueError(
            f"{uri} uses storage-mode {storage_mode!r}, and laying its contigs "
            f"out as {strategy!r} would fold every pixel into the upper "
            "triangle, summing each contact with its mirror. Pass layout='keep' "
            "to relabel without moving any pixel."
        )
    if reorder and inplace:
        raise ValueError(
            f"inplace cannot be used here: laying the contigs out as {strategy!r} "
            "moves them, which means permuting every bin and rewriting the whole "
            "pixel table. Write a new file with output=..., or pass "
            "layout='keep' to relabel without re-ordering (the genomes' bins "
            "then stop being contiguous)."
        )

    genomes = _build_genomes_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=assignment,
        extra_genomes_df=(
            _load_extra_genomes_table(genomes_table)
            if genomes_table is not None
            else _carry_genome_metadata(
                old_genomes,
                set(assignment.values()),
                {before[c] for c in report.moved}
                | {assignment[c] for c in report.moved},
            )
        ),
    )
    contigs = _build_contigs_table(
        contig_order=contig_order,
        contig_lengths=contig_lengths,
        contig_genome_map=assignment,
        genomes=genomes,
        contig_df=_carry_contig_metadata(old_contigs, assignment),
    )
    genomes = strip_qc_columns(genomes)
    contigs = strip_qc_columns(contigs)
    if not reorder:
        # The re-ordering path gets n_bins from _compute_all_offsets, which needs
        # the bins table; here the bins have not moved, so the per-contig counts
        # already on disk answer the same question without reading them.
        _insert_n_bins(mc, genomes, contigs)

    source_metadata = mc.metadata
    output_metadata, had_qc = _build_metadata(
        source_metadata,
        metadata,
        source_uri=str(uri),
        table=str(contig_genome_map) if not isinstance(contig_genome_map, dict) else "",
        layout=strategy,
        n_moved=len(report.moved),
    )

    if inplace:
        _relabel_in_place(mc, genomes, contigs, ordering, output_metadata)
        target_uri = mc.uri
    else:
        assert output is not None  # guaranteed by the argument check above
        target = Path(output)
        with _staged_output(target) as staged:
            if reorder:
                logger.info(
                    "Re-ordering %d contigs into %r layout: the bins move with "
                    "them, so the %s pixels are remapped and re-sorted.",
                    len(contig_order),
                    strategy,
                    f"{int(mc.info['nnz']):,}",
                )
                _rewrite_regrouped(
                    mc,
                    staged,
                    genomes,
                    contigs,
                    contig_order,
                    ordering,
                    output_metadata,
                    chunksize,
                )
            else:
                _materialise_copy(mc, staged)
                _write_relabelled(
                    staged,
                    "/",
                    genomes,
                    contigs,
                    _genome_offset(contigs, len(genomes)),
                    ordering,
                    output_metadata,
                )
            validate(staged)
        target_uri = str(target)

    if had_qc:
        # Three of the block's four categories (inter-contig/intra-genome,
        # inter-genome, and every per-genome qc_* column) are defined by the
        # assignment this command just changed, so the source block cannot be
        # carried over whatever the caller chooses here.
        if qc:
            _recompute_qc(target_uri, source_metadata)
        else:
            logger.warning(
                "The source's QC block describes the old genome assignment and "
                "was dropped; run 'metacooler qc --store %s' to measure the new "
                "one.",
                target_uri,
            )

    logger.info(
        "Reassigned %d contig(s) across %d genome(s) → %s",
        len(report.moved),
        len(genomes),
        target_uri,
    )


def _insert_n_bins(
    mc: MetaCooler, genomes: pd.DataFrame, contigs: pd.DataFrame
) -> None:
    """Insert the per-genome bin count into *genomes*, without reading ``/bins``.

    The re-ordering path gets this column from
    :func:`metacooler.io.create._compute_all_offsets`, which counts the bins of
    the table it is about to write.  A relabelling writes no bins at all, and
    ``indexes/chrom_offset`` already holds the same counts per contig, so this
    is O(n_contigs) instead of O(n_bins).

    Args:
        mc: The open source file.
        genomes: Genomes table built by
            :func:`metacooler.io.create._build_genomes_table`; modified in
            place, with ``n_bins`` inserted after ``n_contigs``.
        contigs: Contigs table in the file's own row order — the precondition
            that lets ``chrom_offset`` be read positionally.
    """
    with h5py.File(mc.store, "r") as f:
        key = mc._group_path.strip("/")
        grp = f[key] if key else f
        chrom_offset = np.asarray(grp["indexes"]["chrom_offset"][:], dtype=np.int64)
    bins_per_contig = np.diff(chrom_offset)
    counts = np.bincount(
        contigs["genome_id"].to_numpy(np.int64),
        weights=bins_per_contig.astype(np.float64),
        minlength=len(genomes),
    ).astype(np.int64)
    genomes.insert(3, "n_bins", counts)


def _genome_offset(contigs: pd.DataFrame, n_genomes: int) -> np.ndarray:
    """Return the CSR offsets of contigs per genome, length ``n_genomes + 1``."""
    counts = np.bincount(contigs["genome_id"].to_numpy(np.int64), minlength=n_genomes)
    offset = np.zeros(n_genomes + 1, dtype=np.int64)
    offset[1:] = np.cumsum(counts)
    return offset


def _relabel_in_place(
    mc: MetaCooler,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    ordering: str,
    metadata: dict,
) -> None:
    """Apply a relabelling to the source file itself, then re-validate it.

    Args:
        mc: The open source file.
        genomes: The rebuilt ``/genomes`` table.
        contigs: The rebuilt ``/contigs`` table, in the file's row order.
        ordering: Value for the ``contig-ordering`` attribute.
        metadata: Value for the ``metadata`` attribute.
    """
    logger.info("Rewriting the genome assignment of %s in place", mc.uri)
    _write_relabelled(
        mc.store,
        mc._group_path,
        genomes,
        contigs,
        _genome_offset(contigs, len(genomes)),
        ordering,
        metadata,
    )
    validate(mc.uri)
