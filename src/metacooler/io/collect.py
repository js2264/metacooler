"""Collect single-genome cool files into one metacool file.

Each input — a ``.cool`` file, or one resolution of an ``.mcool`` addressed as
``file.mcool::resolutions/<binsize>`` — becomes **one genome** of the output,
named after the input's file stem unless ``names`` says otherwise.  Every
contig is renamed ``<genome>__<contig>``, so the ``chromosome`` or
``contig_1`` that so many assemblies share cannot collide.

The layout is what makes this cheap.  Contigs are always grouped by genome, in
``genome_id`` order, so every input occupies one contiguous block of output
bins.  An input's pixels only ever address its own block — there is nothing
between genomes to address — so the output pixel table is the concatenation of
the inputs' tables, each shifted to its block.  Under the default
``contig_ordering="cool"``, which keeps every input's own contig order, each
shifted table is still sorted, and the whole collection is written in one
streaming pass with no sort at all.  Only an ordering that reorders contigs
*within* a genome (``genome-then-length``/``genome-then-alpha``) forces a
re-sort, and then only of the inputs it actually reordered, one at a time.

The corollary is worth saying out loud: the inter-genome blocks of a collected
map are **empty by construction**.  The inputs were binned separately, usually
from separate libraries, so no contact between two of them was ever observed.
"""

from __future__ import annotations

import collections
import dataclasses
import logging
import time
from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import (
    _finalize_metacool_file,
    check_ascii_names,
    check_output_not_an_input,
    chrom_codes,
    read_user_table,
    resolve_count_dtype,
    store_qc,
    to_upper_triangle,
    write_cooler,
)
from metacooler.io.create import (
    _CONTIGS_RESERVED_COLUMNS,
    _SORT_IN_MEMORY_PIXEL_THRESHOLD,
    _apply_create_filters,
    _build_contig_order,
    _build_contigs_table,
    _build_genomes_table,
    _compute_all_offsets,
    _ensure_nonempty_pixel_stream,
    _genome_id_order,
    _load_extra_genomes_table,
    _protect_reserved_columns,
    _staged_output,
    _tally_pixel_qc,
    _total_time_str,
)
from metacooler.qc import DEFAULT_DISTANCE_THRESHOLD, ContactQC, bin_maps
from metacooler.schema import (
    DEFAULT_CHUNKSIZE,
    is_metacooler,
    is_mmetacooler,
    validate,
)
from metacooler.uri import parse_uri

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence
    from typing import Literal

    CollectOrdering = Literal["cool", "genome-then-length", "genome-then-alpha"]

logger = logging.getLogger("metacooler")

#: Joins a genome name to a contig name: contig ``chr1`` of genome ``E_coli``
#: is collected as ``E_coli__chr1``.  Chosen to be safe everywhere a contig name
#: travels: not ``:`` (region strings), not ``,`` (comma-separated
#: ``--contigs`` lists), no whitespace, and nothing a shell needs quoted.
CONTIG_NAME_SEPARATOR = "__"

#: ``contig_ordering`` values ``collect`` accepts.  All three keep each genome's
#: contigs together; ``alpha`` (and ``fasta``) would interleave genomes, which
#: is never useful for a collection and would force a global sort.
COLLECT_ORDERINGS = ("cool", "genome-then-length", "genome-then-alpha")

# How many offending names to quote in a warning or an error.
_MAX_REPORTED = 10


@dataclasses.dataclass
class _Input:
    """What collecting one input needs to know about it, read up front.

    Attributes:
        source: The input exactly as the caller spelled it, for messages.
        h5_path: HDF5 file holding the input.
        h5_group: Group of the cooler inside *h5_path* (``/`` for a ``.cool``).
        genome: Name of the genome this input becomes.
        binsize: The input's ``bin-size``, or ``None`` for variable bins.
        nnz: Number of pixels.
        count_dtype: Lossless dtype for the input's ``count`` column.
        has_weight: Whether ``bins/weight`` (balancing) is present.
        chrom_names: The input's contig names, in its own ``chroms`` order.
        contig_names: The collected ``<genome>__<contig>`` names, aligned.
        chrom: Contig index (into *chrom_names*) of every bin.
        start: Start of every bin.
        end: End of every bin.
        lengths: Length of every contig: its last bin's end, as ``create cool``
            measures it, or 0 for a contig with no bin.
    """

    source: str
    h5_path: str
    h5_group: str
    genome: str
    binsize: int | None
    nnz: int
    count_dtype: np.dtype
    has_weight: bool
    chrom_names: list[str]
    contig_names: list[str]
    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    lengths: np.ndarray


@dataclasses.dataclass
class _Block:
    """One input's slice of the output bin table.

    Attributes:
        inp: The input.
        old_to_new: Output bin id of each of the input's bins, or ``-1`` for a
            bin whose contig was dropped (by ``filters``, or for having no
            length).
        monotonic: Whether *old_to_new* is increasing over the kept bins —
            i.e. the input's contigs kept their relative order, so its pixels
            stay sorted once shifted and need no re-sort.
    """

    inp: _Input
    old_to_new: np.ndarray
    monotonic: bool

    @property
    def needs_sort(self) -> bool:
        """Whether this block has pixels that must be re-sorted."""
        return not self.monotonic and self.inp.nnz > 0


def collected_contig_name(genome: str, contig: str) -> str:
    """Return the name contig *contig* of genome *genome* is collected under.

    The name is ``<genome>__<contig>``, unless *contig* already starts with
    exactly that prefix: a cool extracted from an earlier collection carries
    its genome prefix already, and collecting it again must not stack a second
    one (``E_coli__E_coli__chr1``).

    Args:
        genome: Genome name.
        contig: Contig name as it appears in the input cool.

    Returns:
        The collected contig name.
    """
    prefix = f"{genome}{CONTIG_NAME_SEPARATOR}"
    return contig if contig.startswith(prefix) else prefix + contig


def _examples(values: Iterable[object]) -> str:
    """Quote the first few *values* for a message, noting how many were cut."""
    values = list(values)
    shown = ", ".join(repr(v) for v in values[:_MAX_REPORTED])
    extra = len(values) - _MAX_REPORTED
    return shown + (f", ... (+{extra} more)" if extra > 0 else "")


def _genome_names(sources: Sequence[str], names: Sequence[str] | None) -> list[str]:
    """Name the genome each input becomes: *names*, or else the file stems.

    Args:
        sources: The inputs, as given.
        names: Explicit genome names, one per input, or ``None``.

    Returns:
        One genome name per input.

    Raises:
        ValueError: If *names* does not hold exactly one non-empty name per
            input, or if two inputs would share a genome name.
    """
    if names is not None:
        genomes = [str(n).strip() for n in names]
        if len(genomes) != len(sources):
            raise ValueError(
                f"--names gives {len(genomes)} genome name(s) for "
                f"{len(sources)} input(s); give exactly one per input, in the "
                "same order."
            )
        if not all(genomes):
            raise ValueError("--names contains an empty genome name.")
        remedy = "Genome names must be unique."
    else:
        genomes = [parse_uri(src)[0].stem for src in sources]
        remedy = (
            "Genome names default to the input file names; use --names to give "
            "each input a distinct one."
        )
    counts = collections.Counter(genomes)
    clashes = sorted(g for g, n in counts.items() if n > 1)
    if clashes:
        detail = "\n".join(
            f"  {g!r} <- "
            + ", ".join(
                s for s, name in zip(sources, genomes, strict=True) if name == g
            )
            for g in clashes[:_MAX_REPORTED]
        )
        raise ValueError(
            f"{len(clashes)} genome name(s) would be used by more than one "
            f"input:\n{detail}\n{remedy}"
        )
    return genomes


def _open_input(source: str) -> cooler.Cooler:
    """Open one input, refusing anything that is not a single-genome cooler.

    Args:
        source: A ``.cool`` path or a cooler URI (``file.mcool::resolutions/N``).

    Returns:
        The open ``cooler.Cooler``.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If *source* is a metacool or mmetacool file, the root of a
            multi-resolution ``.mcool``, not a cooler at all, or a cooler that
            stores a full square matrix.
    """
    path, group = parse_uri(source)
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")
    if group == "/" and is_mmetacooler(path):
        raise ValueError(
            f"{source} is an mmetacool collection, not a cool file. collect "
            "turns each input into one genome; to collect a map stored in a "
            "collection, extract its genomes to cool files first ('metacooler "
            "extract --split')."
        )
    if is_metacooler(source):
        raise ValueError(
            f"{source} is already a metacool file. collect turns each input into "
            "ONE genome, which would erase the genome hierarchy it already has. "
            "To collect its contigs as a single genome anyway, flatten it first: "
            f"'metacooler flatten {source} -o flat.cool'. To collect its genomes "
            "separately, extract them first: 'metacooler extract --split'."
        )
    uri = str(path) if group == "/" else f"{path}::{group}"
    try:
        is_cool = cooler.fileops.is_cooler(uri)
    except (OSError, KeyError):
        is_cool = False
    if not is_cool:
        if group == "/" and cooler.fileops.is_multires_file(str(path)):
            resolutions = [
                p.rsplit("/", 1)[-1]
                for p in cooler.fileops.list_coolers(str(path))
                if p.startswith("/resolutions/")
            ]
            example = resolutions[0] if resolutions else "<binsize>"
            raise ValueError(
                f"{source} is a multi-resolution .mcool file; pick one resolution "
                f"per input, e.g. {path}::resolutions/{example} (available: "
                f"{', '.join(resolutions) or 'none'})."
            )
        raise ValueError(f"{source} is not a cooler file.")
    clr = cooler.Cooler(uri)
    if clr.storage_mode != "symmetric-upper":
        raise ValueError(
            f"{source} stores a full square matrix (storage-mode "
            f"{clr.storage_mode!r}); metacooler only reads symmetric-upper maps."
        )
    return clr


def _read_input(source: str, genome: str) -> _Input:
    """Read everything about *source* that is needed before its pixels.

    Args:
        source: The input, as given.
        genome: Name of the genome it becomes.

    Returns:
        The input's :class:`_Input`.

    Raises:
        ValueError: If the input's bins are not sorted by contig.
    """
    clr = _open_input(source)
    chrom_names = [str(name) for name in clr.chromnames]
    bins = pd.DataFrame(clr.bins()[["chrom", "start", "end"]][:])
    chrom = chrom_codes(bins)
    start = bins["start"].to_numpy(dtype=np.int64)
    end = bins["end"].to_numpy(dtype=np.int64)
    if chrom.size and np.any(np.diff(chrom) < 0):
        raise ValueError(f"{source}: the bins are not sorted by contig.")

    # A contig's length is where its last bin ends — the same measure
    # ``create cool`` uses, so collecting one cool agrees with importing it.
    # One ``reduceat`` over the non-empty contigs, not one groupby per contig.
    counts = np.bincount(chrom, minlength=len(chrom_names))
    lengths = np.zeros(len(chrom_names), dtype=np.int64)
    nonempty = np.flatnonzero(counts)
    if nonempty.size:
        first_bin = (np.cumsum(counts) - counts)[nonempty]
        lengths[nonempty] = np.maximum.reduceat(end, first_bin)

    with h5py.File(clr.filename, "r") as f:
        grp = f[clr.root]
        source_count_dtype = grp["pixels"]["count"].dtype
        nnz = int(grp["pixels"]["bin1_id"].shape[0])
        has_weight = "weight" in grp["bins"]
    binsize = clr.binsize
    return _Input(
        source=source,
        h5_path=clr.filename,
        h5_group=clr.root,
        genome=genome,
        binsize=None if binsize is None else int(binsize),
        nnz=nnz,
        count_dtype=resolve_count_dtype(source_count_dtype, context=source),
        has_weight=has_weight,
        chrom_names=chrom_names,
        contig_names=[collected_contig_name(genome, c) for c in chrom_names],
        chrom=chrom,
        start=start,
        end=end,
        lengths=lengths,
    )


def _check_resolutions(inputs: Sequence[_Input]) -> int | None:
    """Return the bin size every input shares, refusing a mixed collection.

    A metacool file has one bin table.  Collecting a 1 kb map with a 5 kb one
    would produce a file that is ``variable``-binned in name and neither
    resolution in practice, and every later ``binnify``/``zoomify`` would then
    reason about it as a fragment-level map.

    Args:
        inputs: Every input.

    Returns:
        The common bin size, or ``None`` when every input is variable-binned.

    Raises:
        ValueError: If the inputs do not all share one bin size (or all have
            variable bins).
    """
    binsizes = {inp.binsize for inp in inputs}
    if len(binsizes) == 1:
        return next(iter(binsizes))
    detail = "\n".join(
        f"  {inp.source}: "
        + (f"{inp.binsize} bp" if inp.binsize is not None else "variable bins")
        for inp in inputs
    )
    raise ValueError(
        "The inputs are not binned alike, but a metacool file has a single bin "
        f"table:\n{detail}\nCollect maps of one resolution — pick the same "
        "resolution from each .mcool (file.mcool::resolutions/<binsize>), or "
        "coarsen the finer ones first ('cooler coarsen')."
    )


def _contig_catalog(inputs: Sequence[_Input]) -> pd.DataFrame:
    """Tabulate every collectable contig of every input.

    Args:
        inputs: Every input.

    Returns:
        One row per contig that holds at least one bin: ``name`` (collected),
        ``original`` (as in its input), ``genome``, ``length``, ``n_bins`` and
        ``source``, in input order.

    Raises:
        ValueError: If two contigs would be collected under the same name.
    """
    catalog = pd.concat(
        [
            pd.DataFrame(
                {
                    "name": inp.contig_names,
                    "original": inp.chrom_names,
                    "genome": inp.genome,
                    "length": inp.lengths,
                    "n_bins": np.bincount(inp.chrom, minlength=len(inp.chrom_names)),
                    "source": inp.source,
                }
            )
            for inp in inputs
        ],
        ignore_index=True,
    )

    # A contig with no bin has nothing to collect, and a zero-length contig
    # would fail validation; ``create cool`` drops both the same way.
    unusable = catalog["length"] <= 0
    if unusable.any():
        logger.warning(
            "Left out %d contig(s) that have no bins in their input: %s",
            int(unusable.sum()),
            _examples(
                f"{o} ({s})"
                for o, s in zip(
                    catalog.loc[unusable, "original"],
                    catalog.loc[unusable, "source"],
                    strict=True,
                )
            ),
        )
        catalog = catalog.loc[~unusable].reset_index(drop=True)
    if catalog.empty:
        raise ValueError("Nothing to collect: no input has a contig with a bin.")

    # The genome prefix makes names from different genomes distinct, but not
    # quite always: genome 'A' with contig 'B__x' and genome 'A__B' with contig
    # 'x' both give 'A__B__x', and a contig already carrying its genome's
    # prefix meets its unprefixed twin. Rare, so refuse rather than guess.
    clashing = catalog["name"].duplicated(keep=False)
    if clashing.any():
        clashes = catalog.loc[clashing].groupby("name", sort=True)
        detail = "\n".join(
            f"  {name!r} <- "
            + ", ".join(
                f"{o!r} in {s}"
                for o, s in zip(group["original"], group["source"], strict=True)
            )
            for name, group in list(clashes)[:_MAX_REPORTED]
        )
        raise ValueError(
            f"{clashes.ngroups} collected contig name(s) would not be unique:\n"
            f"{detail}\nUse --names to rename the genomes involved."
        )
    return catalog


def _report_genome_divergence(genomes: Sequence[str], extra: pd.DataFrame) -> None:
    """Warn when the genomes table and the inputs name different genomes.

    ``collect``'s counterpart to ``create``'s divergence report, in terms of
    inputs rather than of a contigs table.  Run once, before any ``filters``,
    so a genome a filter removes is not reported as missing from the inputs.

    Args:
        genomes: The collected genome names, one per input.
        extra: Parsed ``--genomes`` table, keyed by ``genome_id``.
    """
    listed = set(extra["genome_id"])
    unlisted = sorted(set(genomes) - listed - {"__unbinned__"})
    unused = sorted(listed - set(genomes))
    if unlisted:
        logger.warning(
            "%d collected genome(s) have no row in the genomes table: %s. They "
            "are kept, with empty metadata columns.",
            len(unlisted),
            _examples(unlisted),
        )
    if unused:
        logger.warning(
            "%d genome(s) in the genomes table match no input: %s. They are "
            "absent from /genomes. Genome names come from the input file names "
            "unless --names is given; the collected genomes are: %s.",
            len(unused),
            _examples(unused),
            _examples(sorted(genomes)),
        )


def _load_contig_metadata(
    path: str | Path, catalog: pd.DataFrame
) -> pd.DataFrame | None:
    """Read a per-contig metadata table and resolve its rows to collected contigs.

    The table is a TSV with a header whose first column is ``contig_name``.
    A row may name a contig as it appears in its input cool (``chromosome``)
    or as it is collected (``E_coli__chromosome``).  A ``genome_name`` column
    is optional; when present, each row is looked up inside that genome, which
    is what disambiguates a name several inputs share, and a row that puts a
    contig under a genome it does not come from is an error — collect assigns
    every contig to the genome of its input, and never re-assigns one.  Every
    other column is metadata, stored in ``/contigs``.

    Args:
        path: The table.
        catalog: Output of :func:`_contig_catalog`.

    Returns:
        A frame shaped like :func:`metacooler.io.contigs.parse_contigs_table`'s:
        ``contig_name`` (collected name), ``genome_id`` (genome name) and the
        metadata columns, one row per described contig.  ``None`` when the
        table has no metadata column.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If the table is empty or malformed, names a contig
            ambiguously, or contradicts the genome a contig comes from.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Contigs table not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Contigs table is empty: {path}")
    table = read_user_table(path, kind="Contigs table")
    if str(table.columns[0]) != "contig_name":
        raise ValueError(
            f"Contigs table expects its first column to be 'contig_name', got "
            f"{str(table.columns[0])!r}: {path}"
        )
    has_genome = "genome_name" in table.columns
    names = table["contig_name"]
    genome_of = pd.Series(catalog["genome"].to_numpy(), index=catalog["name"])

    if has_genome:
        genomes = table["genome_name"]
        # Inside the genome first: that lookup is the one a shared name needs.
        by_key = pd.Series(
            catalog["name"].to_numpy(),
            index=(catalog["genome"] + "\t" + catalog["original"]).to_numpy(),
        )
        resolved = (genomes + "\t" + names).map(by_key)
        direct = resolved.isna() & names.isin(genome_of.index)
        wrong = direct & (names.map(genome_of) != genomes)
        resolved = resolved.fillna(names.where(direct & ~wrong))
        # A contig that does exist, but under another collected genome.
        elsewhere = (
            resolved.isna()
            & genomes.isin(set(catalog["genome"]))
            & names.isin(set(catalog["original"]))
        )
        misplaced = wrong | elsewhere
        if misplaced.any():
            origin = catalog["genome"].map(repr) + " (" + catalog["source"] + ")"
            origin_of_collected = pd.Series(origin.to_numpy(), index=catalog["name"])
            origin_of_original = origin.groupby(catalog["original"]).agg(" and ".join)
            detail = "\n".join(
                f"  {c!r} is listed under {g!r} but comes from genome "
                + (
                    origin_of_collected[c]
                    if c in origin_of_collected.index
                    else origin_of_original[c]
                )
                for c, g in list(
                    zip(names[misplaced], genomes[misplaced], strict=True)
                )[:_MAX_REPORTED]
            )
            raise ValueError(
                f"Contigs table ({path}): {int(misplaced.sum())} row(s) put a "
                f"contig under a genome it does not come from:\n{detail}\n"
                "collect assigns every contig to the genome of its input file; "
                "a genome_name column only says which input a row refers to."
            )
    else:
        originals = catalog["original"]
        shared = originals.value_counts()
        shared = shared[shared > 1]
        unique_original = pd.Series(
            catalog["name"].to_numpy(), index=originals.to_numpy()
        ).drop(labels=shared.index)
        # An original name that is also *another* contig's collected name.
        crossed = set(originals[originals != catalog["name"]]) & set(catalog["name"])
        ambiguous = names.isin(shared.index) | names.isin(crossed)
        if ambiguous.any():
            raise ValueError(
                f"Contigs table ({path}): {int(ambiguous.sum())} row(s) name a "
                "contig ambiguously — the name belongs to more than one input, "
                f"e.g. {_examples(pd.unique(names[ambiguous]))}. Add a "
                "genome_name column saying which genome each row describes, or "
                "use the collected names ('<genome>__<contig>')."
            )
        resolved = names.map(unique_original)
        resolved = resolved.fillna(names.where(names.isin(genome_of.index)))

    matched = resolved.notna()
    if not matched.all():
        logger.warning(
            "Contigs table (%s): %d of %d row(s) match no collected contig and "
            "are ignored, e.g. %s. A row names a contig as it appears in its "
            "input cool, or as collected ('<genome>__<contig>')%s.",
            path,
            int((~matched).sum()),
            len(table),
            _examples(pd.unique(names[~matched])),
            ", in the genome its genome_name column gives" if has_genome else "",
        )
    repeated = resolved[matched].duplicated(keep=False)
    if repeated.any():
        logger.warning(
            "Contigs table (%s): %d contig(s) are described by more than one "
            "row; the last row wins: %s",
            path,
            int(resolved[matched][repeated].nunique()),
            _examples(pd.unique(resolved[matched][repeated])),
        )

    metadata_columns = [
        c for c in table.columns if c not in ("contig_name", "genome_name")
    ]
    if not metadata_columns:
        logger.info(
            "Contigs table (%s) has no column beyond contig_name/genome_name; "
            "nothing to add to /contigs.",
            path,
        )
        return None
    # Rename a user column called name/length/genome_id before adding the
    # genome_id column this frame needs for itself.
    extras = _protect_reserved_columns(
        table.loc[matched, metadata_columns].reset_index(drop=True),
        reserved=_CONTIGS_RESERVED_COLUMNS,
        keys=frozenset(),
        context="Contigs table",
    )
    keys = resolved[matched].to_numpy()
    frame = pd.concat(
        [
            pd.DataFrame(
                {"contig_name": keys, "genome_id": genome_of.loc[keys].to_numpy()}
            ),
            extras,
        ],
        axis=1,
    )
    return frame.drop_duplicates(subset="contig_name", keep="last")


def _contig_order(
    inputs: Sequence[_Input],
    contig_lengths: dict[str, int],
    contig_genome: dict[str, str],
    ordering: CollectOrdering,
) -> list[str]:
    """Order the surviving contigs: genome by genome, each as *ordering* says.

    Args:
        inputs: Every input.
        contig_lengths: Collected name to length, for the contigs to keep.
        contig_genome: Collected name to genome name.
        ordering: One of :data:`COLLECT_ORDERINGS`.

    Returns:
        The collected contig names, grouped by genome in ``genome_id`` order.
    """
    if ordering != "cool":
        return _build_contig_order(contig_lengths, contig_genome, ordering=ordering)
    # Genomes in genome_id order, each input's contigs in its own order.
    by_genome = {inp.genome: inp for inp in inputs}
    present = {contig_genome[c] for c in contig_lengths}
    return [
        name
        for genome in _genome_id_order(present)
        for name in by_genome[genome].contig_names
        if name in contig_lengths
    ]


def _layout_bins(
    inputs: Sequence[_Input], contig_id: dict[str, int]
) -> tuple[pd.DataFrame, list[_Block]]:
    """Build the output bin table, one contiguous block per input.

    Vectorised per input — one ``np.fromiter`` over its contig names and array
    operations over its bins — never per contig.

    Args:
        inputs: The inputs that still have a contig, in ``genome_id`` order,
            which is the order their blocks take in the output.
        contig_id: Collected contig name to output contig id, for every
            surviving contig.

    Returns:
        ``(bins, blocks)``: the output bins (int32 ``chrom``/``start``/``end``,
        sorted by ``chrom`` then ``start``), and one :class:`_Block` per input.
    """
    frames: list[pd.DataFrame] = []
    blocks: list[_Block] = []
    offset = 0
    for inp in inputs:
        new_cid = np.fromiter(
            (contig_id.get(name, -1) for name in inp.contig_names),
            dtype=np.int64,
            count=len(inp.contig_names),
        )
        new_chrom = new_cid[inp.chrom]
        kept = np.flatnonzero(new_chrom >= 0)
        # A contig's bins are contiguous and sorted in its input, so the
        # stable sort below keeps them so: only the contigs change places.
        monotonic = bool(np.all(np.diff(new_chrom[kept]) >= 0))
        if not monotonic:
            kept = kept[np.argsort(new_chrom[kept], kind="stable")]
        old_to_new = np.full(inp.chrom.size, -1, dtype=np.int64)
        old_to_new[kept] = offset + np.arange(kept.size, dtype=np.int64)
        offset += kept.size
        frames.append(
            pd.DataFrame(
                {
                    "chrom": new_chrom[kept],
                    "start": inp.start[kept],
                    "end": inp.end[kept],
                }
            )
        )
        blocks.append(_Block(inp, old_to_new, monotonic))
    bins = pd.concat(frames, ignore_index=True).astype(
        {"chrom": np.int32, "start": np.int32, "end": np.int32}
    )
    # Holds by construction (genome-grouped order, blocks in genome_id order);
    # everything downstream — the offsets, cooler's index — depends on it.
    if np.any(np.diff(bins["chrom"].to_numpy()) < 0):
        raise RuntimeError("collect laid out bins out of contig order")
    return bins, blocks


def _block_pixels(
    block: _Block, count_dtype: np.dtype, chunksize: int
) -> Iterator[pd.DataFrame]:
    """Yield one input's pixels, shifted into its block, in the input's order.

    Pixels touching a dropped contig are discarded.  When the block is
    :attr:`_Block.monotonic` the chunks come out sorted and upper-triangular
    as they are; otherwise the upper triangle is restored here and the order
    is someone else's job.

    The three columns are sliced straight out of HDF5 rather than through
    ``cooler.Cooler(...).pixels()``: the input was already vetted by
    :func:`_open_input`, and a second ``Cooler`` per input — plus a DataFrame
    per chunk — was the largest cost left on a collection of many small maps.

    Args:
        block: The input's block.
        count_dtype: Output ``count`` dtype.
        chunksize: Pixels read per chunk.

    Yields:
        Pixel frames (``bin1_id``, ``bin2_id``, ``count``).
    """
    inp = block.inp
    with h5py.File(inp.h5_path, "r") as f:
        pixels = f[inp.h5_group]["pixels"]
        for lo in range(0, inp.nnz, chunksize):
            hi = min(lo + chunksize, inp.nnz)
            bin1 = block.old_to_new[pixels["bin1_id"][lo:hi]]
            bin2 = block.old_to_new[pixels["bin2_id"][lo:hi]]
            count = pixels["count"][lo:hi]
            keep = (bin1 >= 0) & (bin2 >= 0)
            if not keep.all():
                bin1, bin2, count = bin1[keep], bin2[keep], count[keep]
            if not block.monotonic:
                bin1, bin2 = to_upper_triangle(bin1, bin2)
            yield pd.DataFrame(
                {
                    "bin1_id": bin1,
                    "bin2_id": bin2,
                    "count": count.astype(count_dtype, copy=False),
                }
            )


def _sorted_block_pixels(
    block: _Block, count_dtype: np.dtype, chunksize: int
) -> Iterator[pd.DataFrame]:
    """Yield one input's pixels, shifted into its block and re-sorted in memory.

    A block is a bijective relabelling of one input's bins, so its pixels hold
    no duplicate to aggregate: a ``lexsort`` is all it takes.

    Args:
        block: The input's block.
        count_dtype: Output ``count`` dtype.
        chunksize: Pixels read, and yielded, per chunk.

    Yields:
        Sorted pixel frames.
    """
    frames = list(_block_pixels(block, count_dtype, chunksize))
    if not frames:
        return
    bin1 = np.concatenate([f["bin1_id"].to_numpy() for f in frames])
    bin2 = np.concatenate([f["bin2_id"].to_numpy() for f in frames])
    count = np.concatenate([f["count"].to_numpy() for f in frames])
    del frames
    order = np.lexsort((bin2, bin1))
    bin1, bin2, count = bin1[order], bin2[order], count[order]
    for lo in range(0, bin1.size, chunksize):
        hi = lo + chunksize
        yield pd.DataFrame(
            {"bin1_id": bin1[lo:hi], "bin2_id": bin2[lo:hi], "count": count[lo:hi]}
        )


def _coalesce(frames: Iterable[pd.DataFrame], rows: int) -> Iterator[pd.DataFrame]:
    """Merge consecutive frames until each holds at least *rows* rows.

    Inputs are often small — a bacterial map at 1 kb holds a few million
    pixels at most — and under an unordered write every frame handed to cooler
    becomes one ``cooler.create`` call with a fixed per-call cost, so a
    collection of hundreds of genomes must not become hundreds of calls.
    Concatenating *consecutive* frames keeps an ordered stream ordered.

    Args:
        frames: Pixel frames.
        rows: Target rows per yielded frame.

    Yields:
        Pixel frames of at least *rows* rows, except possibly the last.
    """
    pending: list[pd.DataFrame] = []
    n = 0
    for frame in frames:
        if frame.empty:
            continue
        pending.append(frame)
        n += len(frame)
        if n >= rows:
            yield (
                pending[0]
                if len(pending) == 1
                else pd.concat(pending, ignore_index=True)
            )
            pending, n = [], 0
    if pending:
        yield pending[0] if len(pending) == 1 else pd.concat(pending, ignore_index=True)


def _choose_sort(blocks: Sequence[_Block], sort_in_memory: bool | None) -> str:
    """Decide how the pixel stream gets into ``(bin1_id, bin2_id)`` order.

    Args:
        blocks: Every block.
        sort_in_memory: ``True``/``False`` to force re-sorting reordered
            inputs in memory or through cooler's external sort; ``None`` to
            decide from the largest input that needs it.

    Returns:
        ``"none"`` when no input needs re-sorting, else ``"memory"`` or
        ``"external"``.
    """
    to_sort = [b for b in blocks if b.needs_sort]
    if not to_sort:
        logger.info(
            "Every input keeps its own contig order, so the pixels are streamed "
            "into place without sorting."
        )
        return "none"
    if sort_in_memory is None:
        largest = max(to_sort, key=lambda b: b.inp.nnz)
        if largest.inp.nnz > _SORT_IN_MEMORY_PIXEL_THRESHOLD:
            logger.warning(
                "%s holds %s pixels, above the %s-pixel in-memory budget, and "
                "its contigs were reordered — switching to cooler's external "
                "sort. Pass sort_in_memory=True to override, or keep each "
                "input's own order (contig ordering 'cool'), which needs no "
                "sort at all.",
                largest.inp.source,
                f"{largest.inp.nnz:,}",
                f"{_SORT_IN_MEMORY_PIXEL_THRESHOLD:,}",
            )
            return "external"
        sort_in_memory = True
    logger.info(
        "Re-sorting the pixels of %d input(s) whose contigs were reordered, %s.",
        len(to_sort),
        "in memory, one input at a time"
        if sort_in_memory
        else "through cooler's external sort",
    )
    return "memory" if sort_in_memory else "external"


def collect(
    inputs: Sequence[str | Path],
    output: str | Path,
    names: Sequence[str] | None = None,
    genomes_table: str | Path | None = None,
    contigs_table: str | Path | None = None,
    contig_ordering: CollectOrdering = "cool",
    assembly_name: str | None = None,
    filters: Sequence[str] | None = None,
    metadata: dict | None = None,
    sort_in_memory: bool | None = None,
    chunksize: int = DEFAULT_CHUNKSIZE,
    qc: bool = True,
    qc_distance_threshold: int = DEFAULT_DISTANCE_THRESHOLD,
) -> None:
    """Collect single-genome cool files into one metacool file, one genome each.

    Every input becomes one genome of the output, named after its file stem
    (``E_coli.cool`` → ``E_coli``) or by *names*.  Every contig is renamed
    ``<genome>__<contig>`` (see :func:`collected_contig_name`), so contig names
    shared between assemblies cannot collide.  Bins and pixels are carried
    over unchanged; contacts *between* two inputs were never observed, so the
    output's inter-genome blocks are empty.

    Args:
        inputs: Cool files, or cooler URIs such as
            ``file.mcool::resolutions/1000``.  All must share one bin size, or
            all have variable (e.g. restriction-fragment) bins.  Metacool files
            are refused: each input must be a single genome.
        output: Destination ``.metacool`` path.  Written to a ``.partial``
            sibling first and moved into place once it validates.
        names: Genome name for each input, in order.  Defaults to the file
            stems, which must then be unique.
        genomes_table: Optional genomes metadata TSV with a header whose first
            column is ``genome_name``, matched against the genome names; any
            further column is stored in ``/genomes``.  A mismatch in either
            direction is reported as a warning.
        contigs_table: Optional per-contig metadata TSV with a header whose
            first column is ``contig_name``; any further column is stored in
            ``/contigs``.  A row names a contig as in its input cool or as
            collected.  An optional ``genome_name`` column says which genome
            each row describes — it disambiguates names several inputs share,
            and is checked, never used to move a contig to another genome.
        contig_ordering: ``"cool"`` (default) keeps every input's own contig
            order; ``"genome-then-length"`` puts each genome's longest contig
            first; ``"genome-then-alpha"`` sorts each genome's contigs by name.
            Genomes are always laid out in ``genome_id`` (alphabetical) order.
            ``"cool"`` never needs a pixel sort.
        assembly_name: Stored in the ``genome-assembly`` attribute
            (``"unknown"`` when not given).
        filters: ``--filter`` expressions (see :mod:`metacooler.filters`)
            applied before anything is written, e.g.
            ``["contigs.length >= 5kb"]``.  A dropped contig keeps none of its
            bins and every pixel touching it is discarded; a genome left with
            no contig is left out.  Unlike ``create``, ``genomes.n_bins`` is
            available, since the bins come with the inputs.
        metadata: Metadata dict stored in the ``metadata`` attribute.
        sort_in_memory: Only matters when *contig_ordering* reorders an
            input's contigs.  ``True`` re-sorts such inputs in memory, one at a
            time; ``False`` hands the stream to cooler's external sort;
            ``None`` (default) sorts in memory unless the largest such input
            exceeds 100 M pixels.
        chunksize: Pixels read per chunk, and the minimum size of a chunk
            handed to cooler.
        qc: Tally contact-quality statistics off the pixel stream and store
            them in ``metadata["qc"]`` and the ``/genomes`` ``qc_*`` columns,
            as ``create cool`` does (``source="pixels"``).  The
            ``inter_genome`` category is 0 by construction.
        qc_distance_threshold: Intra-contig short/long split for *qc*, in bp.

    Raises:
        ValueError: On an empty or ill-formed input list, inputs of different
            resolutions, genome or contig names that would not be unique, a
            malformed metadata table, or an unknown *contig_ordering*.
        FileNotFoundError: If an input or a table does not exist.
        metacooler.filters.FilterError: If a filter is malformed or keeps
            nothing.
        MetaCoolerFormatError: If a genome or contig name is not ASCII.
    """
    timer = time.perf_counter()
    # A lone path is one input, not a sequence of one-character paths.
    if isinstance(inputs, str | Path):
        inputs = [inputs]
    if isinstance(names, str):
        # Iterating it would silently name two inputs 'A' and 'B' from "AB".
        raise TypeError(
            "names must be a sequence of genome names, one per input, not a "
            f"single string ({names!r})."
        )
    sources = [str(src) for src in inputs]
    output = Path(output)
    if not sources:
        raise ValueError("collect needs at least one input cool file.")
    if contig_ordering not in COLLECT_ORDERINGS:
        raise ValueError(
            f"Unknown contig ordering {contig_ordering!r} for collect; expected "
            f"one of {', '.join(COLLECT_ORDERINGS)}."
        )
    if chunksize < 1:
        raise ValueError(f"chunksize must be a positive integer, got {chunksize}.")
    check_output_not_an_input(output, sources)
    genome_names = _genome_names(sources, names)

    # --- Read every input's bins before a single pixel ---
    inputs_ = [
        _read_input(src, genome)
        for src, genome in zip(sources, genome_names, strict=True)
    ]
    binsize = _check_resolutions(inputs_)
    count_dtype = np.result_type(*(inp.count_dtype for inp in inputs_))
    for inp in inputs_:
        logger.debug(
            "Genome %r <- %s: %d contig(s), %d bin(s), %d pixel(s)",
            inp.genome,
            inp.source,
            len(inp.chrom_names),
            inp.chrom.size,
            inp.nnz,
        )
    weighted = [inp.source for inp in inputs_ if inp.has_weight]
    if weighted:
        logger.info(
            "Balancing weights (bins/weight) are not carried over from %d "
            "input(s) (%s); re-balance the collected map with 'cooler balance' "
            "if you need them.",
            len(weighted),
            _examples(weighted),
        )

    # --- Name every contig, then read the optional tables against the names ---
    catalog = _contig_catalog(inputs_)
    contig_lengths: dict[str, int] = dict(
        zip(catalog["name"], catalog["length"].tolist(), strict=True)
    )
    contig_genome: dict[str, str] = dict(
        zip(catalog["name"], catalog["genome"], strict=True)
    )
    extra_genomes_df = _load_extra_genomes_table(genomes_table)
    if extra_genomes_df is not None:
        _report_genome_divergence(genome_names, extra_genomes_df)
    contig_df = (
        _load_contig_metadata(contigs_table, catalog)
        if contigs_table is not None
        else None
    )

    # --- Apply --filter predicates before anything is written ---
    contig_lengths, _ = _apply_create_filters(
        contig_lengths,
        contig_genome,
        contig_df,
        extra_genomes_df,
        filters,
        genome_n_bins=catalog.groupby("genome")["n_bins"].sum().to_dict(),
        report_divergence=False,
        dropped_effect="keep none of their bins, and every pixel touching one "
        "is discarded",
    )

    # --- Order the contigs and build the /genomes and /contigs tables ---
    contig_order = _contig_order(
        inputs_, contig_lengths, contig_genome, contig_ordering
    )
    genomes_df = _build_genomes_table(
        contig_order,
        contig_lengths,
        contig_genome,
        extra_genomes_df,
        report_divergence=False,
    )
    contigs_df = _build_contigs_table(
        contig_order, contig_lengths, contig_genome, genomes_df, contig_df
    )
    # Before the pixels are streamed: a name cooler cannot read costs nothing
    # here and a whole pass over every input if it is only caught at write time.
    check_ascii_names(genomes_df, contigs_df)

    # --- Lay out the bins, one block per surviving genome ---
    by_genome = {inp.genome: inp for inp in inputs_}
    contig_id = {name: i for i, name in enumerate(contig_order)}
    bins_df, blocks = _layout_bins(
        [by_genome[g] for g in genomes_df["name"]], contig_id
    )
    chrom_offset, genome_offset = _compute_all_offsets(bins_df, contigs_df, genomes_df)
    del catalog, contig_df, extra_genomes_df, contig_genome, contig_lengths
    logger.info(
        "Collecting %d genome(s) from %d input(s): %d contig(s), %d bin(s).",
        len(genomes_df),
        len(inputs_),
        len(contigs_df),
        len(bins_df),
    )

    # --- Stream the pixels, block after block ---
    sort = _choose_sort(blocks, sort_in_memory)
    ordered = sort != "external"

    def _pixels() -> Iterator[pd.DataFrame]:
        for block in blocks:
            if sort == "memory" and block.needs_sort:
                yield from _sorted_block_pixels(block, count_dtype, chunksize)
            else:
                yield from _block_pixels(block, count_dtype, chunksize)

    pixel_source: Iterable[pd.DataFrame] = _coalesce(_pixels(), chunksize)
    qc_tally = (
        ContactQC(len(genomes_df), qc_distance_threshold, source="pixels")
        if qc
        else None
    )
    if qc_tally is not None:
        pixel_source = _tally_pixel_qc(
            pixel_source, qc_tally, bin_maps(bins_df, contigs_df)
        )
    pixel_source = _ensure_nonempty_pixel_stream(pixel_source)
    with _staged_output(output) as staged:
        write_cooler(
            staged,
            bins_df,
            pixel_source,
            binsize=binsize,
            ordered=ordered,
            symmetric_upper=True,
            assembly=assembly_name,
            metadata=metadata or {},
            dtypes={"count": count_dtype},
        )

        logger.info("Post-processing cooler file to add metacool structure")
        _finalize_metacool_file(
            staged,
            genomes_df,
            contigs_df,
            chrom_offset,
            genome_offset,
            contig_ordering=contig_ordering,
            assembly_name=assembly_name,
            binsize=binsize,
        )

        # Complete only now: the tally rides on the stream write_cooler consumed.
        if qc_tally is not None:
            store_qc(
                staged,
                qc_tally.summary(binsize=binsize),
                qc_tally.genome_table(genomes_df["name"].to_numpy()),
            )

        logger.info("Validating collected metacool file")
        validate(staged)

    with h5py.File(output, "r") as f:
        nnz = int(f.attrs["nnz"])
        total = f.attrs["sum"]
    logger.info(
        "Collected %d genome(s), %d contig(s), %d bin(s). NNZ=%d, sum=%s.",
        len(genomes_df),
        len(contigs_df),
        len(bins_df),
        nnz,
        total,
    )
    logger.info(
        "Created metacool file: %s (%s elapsed)", output, _total_time_str(timer)
    )
