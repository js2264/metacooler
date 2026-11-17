"""HDF5 format constants, dtypes, and validators for metacooler."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal, NamedTuple

import h5py
import numpy as np
import pandas as pd

logger = logging.getLogger("metacooler")

METACOOLER_FORMAT = "HDF5::MetaCooler"
METACOOLER_FORMAT_VERSION = 1
METACOOLER_FORMAT_URL = "https://github.com/js2264/metacooler"
MMETACOOLER_FORMAT = "HDF5::MMetaCooler"
MMETACOOLER_FORMAT_VERSION = 1
# These two mirror ``cooler.create._constants.MAGIC`` / ``.URL``.  They are
# vendored as literals rather than imported: that module is private, so an
# upstream refactor (rename, relocation, lazy submodule) would otherwise break
# every metacooler import.  The values are part of the on-disk cooler format
# and therefore cannot change without a format version bump.
COOLER_FORMAT: Literal["HDF5::Cooler"] = "HDF5::Cooler"
COOLER_FORMAT_URL: Literal["https://github.com/open2c/cooler"] = (
    "https://github.com/open2c/cooler"
)

# Dense matrix size guard: refuse to allocate matrices larger than this
DENSE_MATRIX_BIN_LIMIT = 50_000

# Rows sampled from each probed window of the pixel table (see ``validate``).
PIXEL_SAMPLE_SIZE = 1000

# Storage modes understood by the validator.  metacooler only ever writes
# ``symmetric-upper``; ``square`` is accepted (cooler writes it) but disables
# the upper-triangle check, which is meaningless for a full matrix.
STORAGE_MODES = frozenset({"symmetric-upper", "square"})

# ``contig-ordering`` values that lay contigs out grouped by ascending
# ``genome_id``.  Only for these is ``indexes/genome_offset`` a meaningful
# offset array (see ``validate``).
GENOME_GROUPED_ORDERINGS = frozenset({"genome-then-length", "genome-then-alpha"})

# ``contig-ordering`` value for a file whose row order is inherited from an
# earlier layout and is no longer described by any creation strategy.  Written
# by ``reassign --layout keep`` when the new assignment breaks genome grouping:
# keeping ``genome-then-alpha`` there would claim an invariant the file no
# longer has (and ``validate`` below would refuse the file for it), while
# ``alpha``/``fasta`` would claim a different one it never had.
CUSTOM_ORDERING = "custom"

# The ``--layout`` vocabulary of ``reassign``: the two orderings it can lay a
# file out with, plus the two policy words that choose between re-ordering and
# leaving the rows alone.  Here rather than in ``io.reassign`` so the CLI can
# build its ``click.Choice`` at import time without pulling the io layer (and
# with it cooler and pyfaidx) into every ``metacooler --help``.
Layout = Literal["auto", "keep", "genome-then-alpha", "genome-then-length"]
LAYOUT_CHOICES: tuple[str, ...] = (
    "auto",
    "keep",
    "genome-then-alpha",
    "genome-then-length",
)
# The subset of ``io.create._build_contig_order`` strategies a reassignment can
# use: ``fasta`` and ``cool`` are creation-time orders that need a source file
# to read them from, so they are not on offer.
ContigOrderStrategy = Literal["genome-then-alpha", "genome-then-length"]

# Column dtypes for metacooler-specific tables.
# ``n_bins`` is deliberately absent from GENOMES_DTYPE: the column is inserted
# into the genomes table by ``io.create._compute_all_offsets`` *after* the
# dtype cast in ``_build_genomes_table``, so declaring it here had no effect.
GENOMES_DTYPE = {
    "name": "object",
    "length": "i8",
    "n_contigs": "i4",
}

CONTIGS_DTYPE = {
    "name": "object",
    "length": "i4",
    "genome_id": "i4",
}

# On-disk dtype *kinds* the core columns must have, as ``numpy.dtype.kind``
# characters.  Deliberately a kind and not an exact dtype: a legitimately
# written file may hold ``length`` as int32 or int64 (``from_cool`` copies the
# source width) and ``name`` as fixed-length bytes or a vlen str, and rejecting
# either would fail files metacooler itself wrote.  What is *not* legitimate is
# a float bin coordinate or a string offset, which is what this catches.
CORE_COLUMN_KINDS: dict[str, str] = {
    "genomes/name": "SOU",
    "genomes/length": "iu",
    "genomes/n_contigs": "iu",
    "genomes/n_bins": "iu",
    "contigs/name": "SOU",
    "contigs/length": "iu",
    "contigs/genome_id": "iu",
    "bins/chrom": "iu",
    "bins/start": "iu",
    "bins/end": "iu",
    "pixels/bin1_id": "iu",
    "pixels/bin2_id": "iu",
    "pixels/count": "iuf",
    "indexes/chrom_offset": "iu",
    "indexes/bin1_offset": "iu",
    "indexes/genome_offset": "iu",
}

# Pixel rows read per iteration by the ``--deep`` scan.
DEEP_PIXEL_CHUNK = 1_000_000

# Rows (pixels, or pairs on the ``create`` path) handled per chunk by every
# streaming operation.  High on purpose: any writer that hands cooler an
# unordered chunk stream pays a fresh ``cooler.create`` call *per chunk* plus
# one per merge pass, and that call's cost scales with contig count -- on a
# metagenome it dominates everything else (see CLAUDE.md).  Chunk count, not
# chunk size, is what makes a metagenomic run intractable, so the default
# favours few large chunks.  The cost is peak memory: roughly 24 bytes per
# buffered pixel, so ~1.2 GB at this size.  Lower it on a memory-capped host.
DEFAULT_CHUNKSIZE = 50_000_000

# How many offending values to name per column when refusing non-ASCII text.
_MAX_REPORTED_NON_ASCII = 5


class MetaCoolerFormatError(ValueError):
    """Raised when a file does not conform to the metacooler format."""


def _attr(node: h5py.Group, key: str) -> object:
    """Return an HDF5 attribute, decoding byte strings to ``str``.

    Args:
        node: HDF5 file or group carrying the attribute.
        key: Attribute name.

    Returns:
        The attribute value, or ``None`` if it is absent.
    """
    value = node.attrs.get(key)
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def _int_attr(node: h5py.Group, key: str) -> int | None:
    """Return an integer-valued HDF5 attribute, or ``None`` if unset.

    ``bin-size`` is written as the string ``"null"`` for variable-width bins,
    so string sentinels are treated as "absent" rather than as an error.

    Args:
        node: HDF5 file or group carrying the attribute.
        key: Attribute name.

    Returns:
        The attribute as an ``int``, or ``None`` if absent or explicitly null.

    Raises:
        MetaCoolerFormatError: If the attribute is present but not an integer.
    """
    value = _attr(node, key)
    if value is None:
        return None
    if isinstance(value, str):
        if value.strip().lower() in ("", "null", "none", "nan"):
            return None
        try:
            return int(value)
        except ValueError:
            raise MetaCoolerFormatError(
                f"attrs[{key!r}] must be an integer, got {value!r}"
            ) from None
    if isinstance(value, (bool, np.bool_)):
        raise MetaCoolerFormatError(f"attrs[{key!r}] must be an integer, got {value!r}")
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return int(value)
    raise MetaCoolerFormatError(f"attrs[{key!r}] must be an integer, got {value!r}")


def _table_length(group: h5py.Group, name: str) -> int:
    """Return the shared row count of an HDF5 table group.

    Every dataset under a table group is a column of the same table, so all of
    them must have the same length.  Without this check a truncated column
    (e.g. a ``/pixels/count`` that is shorter than ``/pixels/bin1_id``) passes
    validation and silently mis-pairs counts with coordinates.

    Only dataset *shapes* are read, so the cost is O(n_columns).

    Args:
        group: HDF5 group holding one dataset per column.
        name: Group name, used in the error message.

    Returns:
        The common column length, or ``0`` if the group holds no datasets.

    Raises:
        MetaCoolerFormatError: If two columns have different lengths.
    """
    for key in group:
        if not isinstance(group[key], h5py.Dataset):
            # A column name containing '/' becomes a nested group when written,
            # and cooler -- which reads every member of /chroms, i.e. /contigs,
            # as a column -- then fails with an AttributeError on open.
            raise MetaCoolerFormatError(
                f"/{name}/{key} is a group, not a column; every member of a "
                "table group must be a dataset (a column name containing '/' "
                "is written as a nested group)"
            )
    lengths = {key: len(group[key]) for key in group if group[key].ndim > 0}
    if not lengths:
        return 0
    if len(set(lengths.values())) > 1:
        detail = ", ".join(f"{k}={v}" for k, v in sorted(lengths.items()))
        raise MetaCoolerFormatError(
            f"/{name} columns have unequal lengths ({detail}); every column of a "
            "table must have one row per record"
        )
    return next(iter(lengths.values()))


def _check_offset_array(
    values: np.ndarray, name: str, expected_len: int, total: int
) -> None:
    """Validate the universal invariants of a CSR-style offset array.

    An offset array must have exactly ``expected_len`` entries, start at 0, end
    at *total*, and never decrease.  ``bin1_offset`` in particular is the most
    load-bearing index in the format — every range query walks it — yet it was
    previously only checked at its two endpoints.  ``np.diff`` makes the full
    check O(len(values)).

    Args:
        values: The offset array, as int64.
        name: Dataset name, used in error messages.
        expected_len: Required number of entries.
        total: Required value of the final entry.

    Raises:
        MetaCoolerFormatError: If any invariant is violated.
    """
    if len(values) != expected_len:
        raise MetaCoolerFormatError(f"{name} length ({len(values)}) != {expected_len}")
    if len(values) == 0:
        return
    if int(values[0]) != 0 or int(values[-1]) != total:
        raise MetaCoolerFormatError(
            f"{name} endpoints ({int(values[0])}, {int(values[-1])}) "
            f"do not bracket the table (0, {total})"
        )
    deltas = np.diff(values)
    if np.any(deltas < 0):
        bad = int(np.argmax(deltas < 0))
        raise MetaCoolerFormatError(
            f"{name} is not non-decreasing: {name}[{bad}]={int(values[bad])} > "
            f"{name}[{bad + 1}]={int(values[bad + 1])}"
        )


def _check_offset_ends(
    ds: h5py.Dataset, name: str, expected_len: int, total: int
) -> None:
    """The O(1) part of :func:`_check_offset_array`: length and both endpoints.

    Reads two elements, not the array, so it is what an open can afford on a
    file whose ``bin1_offset`` has millions of entries.

    Raises:
        MetaCoolerFormatError: If the length or an endpoint is wrong.
    """
    if len(ds) != expected_len:
        raise MetaCoolerFormatError(f"{name} length ({len(ds)}) != {expected_len}")
    if len(ds) == 0:
        return
    first, last = int(ds[0]), int(ds[len(ds) - 1])
    if first != 0 or last != total:
        raise MetaCoolerFormatError(
            f"{name} endpoints ({first}, {last}) do not bracket the table (0, {total})"
        )


def _has_duplicates(values: np.ndarray) -> bool:
    """Whether *values* repeats an entry: a sort and one adjacent comparison.

    Four times faster than ``len(np.unique(values)) != len(values)`` on two
    million fixed-width byte names, whose hash-based ``unique`` is the slow
    path in NumPy 2.
    """
    if len(values) < 2:
        return False
    ordered = np.sort(values)
    return bool((ordered[1:] == ordered[:-1]).any())


def _is_text_dtype(dtype: object) -> bool:
    """Return ``True`` for a dtype that can hold text (NumPy, h5py, or pandas)."""
    return getattr(dtype, "kind", "") in ("O", "S", "U") or str(dtype) == "string"


def non_ascii_values(values: object) -> list[str]:
    """Return the distinct values of *values* that are not pure ASCII, in order.

    Sized for the hot path: :func:`validate` calls this on every string column
    of ``/contigs`` on every ``MetaCooler`` construction, so an all-ASCII
    fixed-width byte column is settled by a single vectorised pass over the raw
    buffer and is never decoded value by value.  Only once a high byte turns up
    do we pay for the per-value decode needed to name the offenders.
    """
    arr = np.asarray(values)
    if arr.size == 0:
        return []
    if (
        arr.dtype.kind == "S"
        and arr.dtype.itemsize > 0
        and int(np.max(np.ascontiguousarray(arr).view(np.uint8))) < 0x80
    ):
        return []
    if arr.dtype.kind in ("O", "U") and _all_ascii_text(arr):
        return []
    out: list[str] = []
    seen: set[str] = set()
    for value in arr.tolist():
        if value is None:
            continue
        # "replace" rather than "strict": this runs to *report* bad text, so a
        # column that is not even valid UTF-8 must still produce a message.
        text = (
            value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)
        )
        if not text.isascii() and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _all_ascii_text(arr: np.ndarray) -> bool:
    """Vectorised "is every value plain ASCII" for a ``str``/``object`` column.

    As fixed-width UTF-32 (NumPy ``U``), a value is ASCII exactly when every
    code point is below 0x80, which one ``max`` over the buffer settles.  Bytes
    or mixed objects take the per-value path instead (``False`` here).
    """
    if arr.dtype.kind == "O":
        # Only text (or missing values, which the per-value path skips or
        # renders as ASCII anyway): str() of a bytes value is its ASCII repr.
        if not set(map(type, arr.tolist())) <= {str, type(None), float}:
            return False
        arr = arr.astype("U")
    return int(np.max(np.ascontiguousarray(arr).view(np.uint32))) < 0x80


def ascii_violations(genomes: Any, contigs: Any) -> list[tuple[str, str, list[str]]]:
    """Find text that would break cooler compatibility in the two name tables.

    *genomes* and *contigs* are anything that maps a column name to a sequence
    of values — an ``h5py.Group`` or a ``pandas.DataFrame`` both qualify — so
    the same rule serves the writers (which hold DataFrames) and
    :func:`validate` (which holds the on-disk groups).

    The two tables are checked to different depths on purpose:

    - **Every** string column of ``/contigs``, because ``/chroms`` is a soft
      link to it and cooler reads the whole group.
    - Only ``name`` in ``/genomes``, because cooler never sees that group; the
      names escape into cooler only as chrom names, when ``aggregate --level
      genome`` writes them into a ``.cool``.  An accented ``phylum`` column
      alongside them is metacooler-only metadata and stays legal.

    Returns:
        ``(group, column, offending values)`` triples, empty if the tables are
        clean.
    """

    def _column(table: Any, name: str) -> Any | None:
        try:
            column = table[name]
        except (KeyError, IndexError):
            return None
        return column if _is_text_dtype(getattr(column, "dtype", None)) else None

    found: list[tuple[str, str, list[str]]] = []
    for column_name in list(contigs):
        column = _column(contigs, str(column_name))
        if column is None:
            continue
        bad = non_ascii_values(column)
        if bad:
            found.append(("contigs", str(column_name), bad))
    name_column = _column(genomes, "name")
    if name_column is not None:
        bad = non_ascii_values(name_column)
        if bad:
            found.append(("genomes", "name", bad))
    return found


def format_ascii_violations(violations: list[tuple[str, str, list[str]]]) -> str:
    """Render :func:`ascii_violations` output as the body of an error message."""
    lines = []
    for group, column, values in violations:
        shown = ", ".join(repr(v) for v in values[:_MAX_REPORTED_NON_ASCII])
        extra = len(values) - _MAX_REPORTED_NON_ASCII
        if extra > 0:
            shown += f", ... (+{extra} more)"
        lines.append(f"  /{group}.{column}: {shown}")
    lines.append(
        "cooler decodes /chroms — a soft link to /contigs — as ASCII "
        '(`.astype("U")` in cooler/core/_tableops.py), so a non-ASCII byte in '
        "any /contigs column makes cooler.Cooler(), and therefore MetaCooler() "
        "which opens through it, raise UnicodeDecodeError on a file that is "
        "otherwise perfectly well formed. Genome names meet the same assumption "
        "on the way out, where `aggregate --level genome` writes them into a "
        ".cool as chrom names. Rename to ASCII, or drop the offending column."
    )
    return "\n".join(lines)


def is_metacooler(uri: str | Path) -> bool:
    """Return ``True`` if *uri* points to a valid ``.metacool`` file.

    Checks for the ``metacooler-format`` HDF5 attribute equal to
    ``HDF5::MetaCooler`` and the presence of the ``/genomes`` group.
    Note that ``cooler.fileops.is_cooler()`` will *also* return ``True``
    for a metacool file because of the soft-link ``chroms → contigs``
    and the ``format = "HDF5::Cooler"`` attribute.

    Accepts the same ``file::group`` URI syntax as :func:`validate`, so
    ``is_metacooler("reps.mmetacool::entries/rep1")`` tests the entry rather
    than the collection root.

    Args:
        uri: Path to the file to test, optionally followed by ``::group``.

    Returns:
        ``True`` if the target is a valid metacool file or entry, ``False``
        otherwise (including if the file does not exist or cannot be opened).
    """
    from metacooler.uri import parse_uri

    try:
        path, group_path = parse_uri(str(uri))
    except ValueError:
        return False
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as froot:
            key = group_path.lstrip("/")
            if key and key not in froot:
                return False
            node = froot[key] if key else froot
            if not isinstance(node, h5py.Group):
                return False
            if _attr(node, "metacooler-format") != METACOOLER_FORMAT:
                return False
            return all(
                g in node for g in ("genomes", "contigs", "bins", "pixels", "indexes")
            )
    except Exception:
        return False


def validate(uri: str | Path, *, full: bool = True) -> None:
    """Validate that *uri* is a well-formed ``.metacool`` file.

    Performs a structural integrity check: required groups, columns,
    indexes, and cross-count consistency.

    Every check is O(n_bins) / O(n_contigs) at worst — the pixel table is
    never scanned end to end.  That is still seconds on a multi-million-contig
    metagenome, so ``MetaCooler()`` only runs the *structural* subset
    (``full=False``) on open; ``metacooler check`` (:func:`check`) and every
    writer run the whole thing.

    Args:
        uri: Path to the ``.metacool`` file, or a URI of the form
            ``file.mmetacool::entries/entry`` to validate a sub-group inside
            an mmetacool file.
        full: If ``False``, stop after the checks that cost O(n_columns)
            reads: format attributes, required groups and columns, equal
            column lengths, the ``/chroms`` link, the count attributes, and
            the length and endpoints of every offset index.  Those are what a
            reader needs before it can trust an index lookup; the O(n_bins)
            content checks (bin geometry, offset monotonicity, genome tables,
            names, pixel samples) are skipped.

    Raises:
        MetaCoolerFormatError: If any required group, dataset, or attribute
            is missing or malformed, or if the URI names a group that does
            not exist.
    """
    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(str(uri))
    if not file_path.exists():
        raise MetaCoolerFormatError(f"File not found: {file_path}")

    with h5py.File(file_path, "r") as _froot:
        # Resolve the group named by the URI.  A typo used to surface as a bare
        # h5py KeyError, contradicting the documented contract that every
        # failure of this function is a MetaCoolerFormatError.
        key = group_path.lstrip("/")
        if key:
            if key not in _froot:
                available = sorted(_froot.keys())
                raise MetaCoolerFormatError(
                    f"Group {group_path!r} not found in {file_path}. "
                    f"Top-level keys: {available[:20]}"
                )
            node = _froot[key]
            if not isinstance(node, h5py.Group):
                raise MetaCoolerFormatError(
                    f"{group_path!r} in {file_path} is a "
                    f"{type(node).__name__}, not an HDF5 group"
                )
            f = node
        else:
            f = _froot
        # --- format attributes ---
        fmt = f.attrs.get("format")
        if fmt != COOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected format='{COOLER_FORMAT}', got {fmt!r}"
            )

        mfmt = f.attrs.get("metacooler-format")
        if mfmt != METACOOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected metacooler-format='{METACOOLER_FORMAT}', got {mfmt!r}"
            )

        # --- required groups ---
        for grp in ("genomes", "contigs", "bins", "pixels", "indexes"):
            if grp not in f:
                raise MetaCoolerFormatError(f"Missing required group: /{grp}")

        # --- genomes table ---
        g = f["genomes"]
        for col in ("name", "length", "n_contigs", "n_bins"):
            if col not in g:
                raise MetaCoolerFormatError(f"Missing column in /genomes: {col}")

        # --- contigs table ---
        c = f["contigs"]
        for col in ("name", "length", "genome_id"):
            if col not in c:
                raise MetaCoolerFormatError(f"Missing column in /contigs: {col}")

        # --- bins table ---
        b = f["bins"]
        for col in ("chrom", "start", "end"):
            if col not in b:
                raise MetaCoolerFormatError(f"Missing column in /bins: {col}")

        # --- pixels table ---
        p = f["pixels"]
        for col in ("bin1_id", "bin2_id", "count"):
            if col not in p:
                raise MetaCoolerFormatError(f"Missing column in /pixels: {col}")

        # --- indexes ---
        idx = f["indexes"]
        for ds in ("chrom_offset", "bin1_offset"):
            if ds not in idx:
                raise MetaCoolerFormatError(f"Missing index dataset: /indexes/{ds}")

        # --- per-group column lengths must agree, then cross-check counts ---
        n_genomes = _table_length(g, "genomes")
        n_contigs = _table_length(c, "contigs")
        n_bins = _table_length(b, "bins")
        n_pixels = _table_length(p, "pixels")
        if n_pixels == 0:
            logger.warning("Empty pixels table")
        if n_contigs == 0:
            raise MetaCoolerFormatError("/contigs table is empty")
        if n_genomes == 0:
            raise MetaCoolerFormatError("/genomes table is empty")

        # --- /chroms must be a soft link to /contigs (cooler compatibility) ---
        link = f.get("chroms", getlink=True)
        if not isinstance(link, h5py.SoftLink):
            raise MetaCoolerFormatError(
                "/chroms must be a soft link to /contigs for cooler compatibility"
            )
        # Resolve the link rather than compare its last path component: an
        # absolute "/contigs" is dangling once the map is copied into an
        # .mmetacool entry, and "/nowhere/contigs" or another entry's contigs
        # used to pass, after which MetaCooler() raised an AttributeError.
        target = f.get("chroms")
        if target is None or target != c:
            where = "resolves to nothing" if target is None else "resolves elsewhere"
            raise MetaCoolerFormatError(
                f"/chroms soft link must resolve to this map's /contigs, but "
                f"{link.path!r} {where}"
            )

        # --- attribute consistency ---
        # Every one of these is optional (docs/format.md allows a file written
        # by an older writer to omit them), but must agree with the tables when
        # present.  ``ngenomes`` used to be the only documented count with no
        # check, which made it silently drift from /genomes.
        for attr_name, actual, label in (
            ("nchroms", n_contigs, "len(/contigs)"),
            ("ngenomes", n_genomes, "len(/genomes)"),
            ("nbins", n_bins, "len(/bins)"),
            ("nnz", n_pixels, "len(/pixels)"),
        ):
            declared = _int_attr(f, attr_name)
            if declared is not None and declared != actual:
                raise MetaCoolerFormatError(
                    f"attrs[{attr_name!r}] ({declared}) != {label} ({actual})"
                )

        mfmt_version = _int_attr(f, "metacooler-format-version")
        if mfmt_version is not None and mfmt_version != METACOOLER_FORMAT_VERSION:
            raise MetaCoolerFormatError(
                f"attrs['metacooler-format-version'] is {mfmt_version}; this "
                f"metacooler only understands version {METACOOLER_FORMAT_VERSION}"
            )

        storage_mode = _attr(f, "storage-mode")
        if storage_mode is not None and storage_mode not in STORAGE_MODES:
            raise MetaCoolerFormatError(
                f"attrs['storage-mode'] must be one of {sorted(STORAGE_MODES)}, "
                f"got {storage_mode!r}"
            )
        # --- offset indexes: length and endpoints, read without the body ---
        for ds_name, expected_len, total in (
            ("chrom_offset", n_contigs + 1, n_bins),
            ("bin1_offset", n_bins + 1, n_pixels),
            ("genome_offset", n_genomes + 1, n_contigs),
        ):
            if ds_name not in idx:
                raise MetaCoolerFormatError(f"indexes/{ds_name} is missing")
            _check_offset_ends(idx[ds_name], ds_name, expected_len, total)

        if not full:
            return

        if storage_mode == "square":
            logger.warning(
                "storage-mode='square': skipping the upper-triangle check "
                "(metacooler writes 'symmetric-upper')"
            )

        # --- name uniqueness ---
        if _has_duplicates(np.array(c["name"])):
            raise MetaCoolerFormatError("contig names are not unique")
        if _has_duplicates(np.array(g["name"])):
            raise MetaCoolerFormatError("genome names are not unique")

        # --- names must be ASCII (cooler compatibility) ---
        # Without this, a file with an accented contig name passes every other
        # check here and then raises UnicodeDecodeError the moment cooler opens
        # it — validate() vouching for a file that cannot be read.
        violations = ascii_violations(g, c)
        if violations:
            raise MetaCoolerFormatError(
                f"{file_path} holds names cooler cannot decode:\n"
                + format_ascii_violations(violations)
            )

        # --- contig lengths must be positive ---
        contig_lengths = np.array(c["length"]).astype(np.int64)
        if (contig_lengths <= 0).any():
            raise MetaCoolerFormatError("contig lengths must be positive (> 0)")

        # --- bins/chrom range and ordering ---
        bins_chrom = np.array(b["chrom"]).astype(np.int64)
        if n_bins > 0:
            if np.min(bins_chrom) < 0 or np.max(bins_chrom) >= n_contigs:
                raise MetaCoolerFormatError(
                    f"bins/chrom values out of range [0, {n_contigs})"
                )
            if np.any(np.diff(bins_chrom) < 0):
                raise MetaCoolerFormatError("bins are not sorted by contig (chrom)")

        # --- bins geometry ---
        # Fully vectorised, O(n_bins).  Without this an int32 overflow in a bin
        # tiler can write negative starts, or a stale /contigs/length can leave
        # bins hanging off the end of their contig, and the file still validates.
        bins_start = np.array(b["start"]).astype(np.int64)
        bins_end = np.array(b["end"]).astype(np.int64)
        if n_bins > 0:
            if (bins_start < 0).any():
                bad = int(np.argmax(bins_start < 0))
                raise MetaCoolerFormatError(
                    f"bins/start must be non-negative (bin {bad} starts at "
                    f"{int(bins_start[bad])}); this is the signature of an "
                    "integer overflow in the bin tiler"
                )
            empty_or_inverted = bins_start >= bins_end
            if empty_or_inverted.any():
                bad = int(np.argmax(empty_or_inverted))
                raise MetaCoolerFormatError(
                    f"bins must satisfy start < end (bin {bad}: "
                    f"start={int(bins_start[bad])}, end={int(bins_end[bad])})"
                )
            past_end = bins_end > contig_lengths[bins_chrom]
            if past_end.any():
                bad = int(np.argmax(past_end))
                raise MetaCoolerFormatError(
                    f"bin {bad} ends at {int(bins_end[bad])} but its contig "
                    f"({int(bins_chrom[bad])}) is only "
                    f"{int(contig_lengths[bins_chrom[bad]])} bp long"
                )
            same_contig = bins_chrom[1:] == bins_chrom[:-1]
            overlapping = same_contig & (bins_start[1:] < bins_end[:-1])
            if overlapping.any():
                bad = int(np.argmax(overlapping))
                raise MetaCoolerFormatError(
                    f"bins {bad} and {bad + 1} overlap or are out of order within "
                    f"contig {int(bins_chrom[bad])}: "
                    f"[{int(bins_start[bad])}, {int(bins_end[bad])}) then "
                    f"[{int(bins_start[bad + 1])}, {int(bins_end[bad + 1])})"
                )

        # --- bin-size must describe the bins it claims to describe ---
        bin_type = _attr(f, "bin-type")
        if bin_type is not None and bin_type not in ("fixed", "variable"):
            raise MetaCoolerFormatError(
                f"attrs['bin-type'] must be 'fixed' or 'variable', got {bin_type!r}"
            )
        bin_size = _int_attr(f, "bin-size")
        if bin_type == "fixed" and bin_size is None:
            raise MetaCoolerFormatError(
                "attrs['bin-type']='fixed' but attrs['bin-size'] is missing or null"
            )
        if bin_type == "variable" and bin_size is not None:
            logger.warning(
                "attrs['bin-type']='variable' but attrs['bin-size']=%d; "
                "downstream tools will disagree about the resolution",
                bin_size,
            )
        if bin_size is not None and bin_size <= 0:
            raise MetaCoolerFormatError(
                f"attrs['bin-size'] must be positive, got {bin_size}"
            )
        if bin_size is not None and n_bins > 0:
            # cooler's own rule (``cooler.util.get_binsize``): the last bin of
            # each contig may be short, every other bin is exactly binsize.
            is_last_in_contig = np.append(bins_chrom[1:] != bins_chrom[:-1], True)
            widths = bins_end - bins_start
            wrong = ~is_last_in_contig & (widths != bin_size)
            if wrong.any():
                bad = int(np.argmax(wrong))
                raise MetaCoolerFormatError(
                    f"attrs['bin-size'] is {bin_size} but bin {bad} is "
                    f"{int(widths[bad])} bp wide"
                )

        # --- index datasets: shape, endpoints and monotonicity ---
        chrom_offset = np.array(idx["chrom_offset"]).astype(np.int64)
        _check_offset_array(chrom_offset, "chrom_offset", n_contigs + 1, n_bins)
        bin1_offset = np.array(idx["bin1_offset"]).astype(np.int64)
        _check_offset_array(bin1_offset, "bin1_offset", n_bins + 1, n_pixels)
        genome_offset = np.array(idx["genome_offset"]).astype(np.int64)
        _check_offset_array(genome_offset, "genome_offset", n_genomes + 1, n_contigs)

        # --- contig genome_ids must be in range ---
        genome_ids = np.array(c["genome_id"]).astype(np.int64)
        if np.min(genome_ids) < 0 or np.max(genome_ids) >= n_genomes:
            raise MetaCoolerFormatError(
                f"genome_id values out of range [0, {n_genomes})"
            )

        # --- chrom_offset must be the bin-count CSR over the (sorted) bins ---
        expected_chrom = np.zeros(n_contigs + 1, dtype=np.int64)
        if n_bins > 0:
            expected_chrom[1:] = np.cumsum(np.bincount(bins_chrom, minlength=n_contigs))
        if not np.array_equal(chrom_offset, expected_chrom):
            raise MetaCoolerFormatError(
                "chrom_offset is inconsistent with the bins table"
            )

        genome_contig_counts = np.bincount(genome_ids, minlength=n_genomes)
        # genome_offset only *means* anything when the contigs are physically
        # laid out grouped by ascending genome_id: for `alpha`, `fasta` and
        # `cool` orderings the array is still written (docs/format.md:111
        # describes it unconditionally) but it cannot delimit a genome's contig
        # block, so only its shape/monotonicity — checked above — is validated.
        ordering = _attr(f, "contig-ordering")
        if ordering in GENOME_GROUPED_ORDERINGS:
            if np.any(np.diff(genome_ids) < 0):
                raise MetaCoolerFormatError(
                    f"contig-ordering={ordering!r} requires contigs grouped by "
                    "ascending genome_id, but the genome_id sequence is not "
                    "non-decreasing, so genome_offset does not match the physical "
                    "contig layout"
                )
            expected_genome = np.zeros(n_genomes + 1, dtype=np.int64)
            expected_genome[1:] = np.cumsum(genome_contig_counts)
            if not np.array_equal(genome_offset, expected_genome):
                raise MetaCoolerFormatError(
                    f"genome_offset {genome_offset.tolist()} is inconsistent with "
                    f"the contig genome_id counts; expected "
                    f"{expected_genome.tolist()}"
                )

        # --- genomes/n_contigs, n_bins and length must match reality ---
        if not np.array_equal(
            np.array(g["n_contigs"]).astype(np.int64), genome_contig_counts
        ):
            raise MetaCoolerFormatError(
                "genomes/n_contigs does not match the contigs table"
            )
        if n_bins > 0:
            genome_bin_counts = np.bincount(genome_ids[bins_chrom], minlength=n_genomes)
        else:
            genome_bin_counts = np.zeros(n_genomes, dtype=np.int64)
        if not np.array_equal(
            np.array(g["n_bins"]).astype(np.int64), genome_bin_counts
        ):
            raise MetaCoolerFormatError("genomes/n_bins does not match the bins table")
        genome_lengths = np.array(g["length"]).astype(np.int64)
        expected_lengths = np.bincount(
            genome_ids, weights=contig_lengths.astype(np.float64), minlength=n_genomes
        ).astype(np.int64)
        if not np.array_equal(genome_lengths, expected_lengths):
            bad = int(np.argmax(genome_lengths != expected_lengths))
            raise MetaCoolerFormatError(
                f"genomes/length does not match the contigs table: genome {bad} "
                f"declares {int(genome_lengths[bad])} bp but its contigs sum to "
                f"{int(expected_lengths[bad])} bp"
            )

        # --- pixel ordering, range, and upper-triangle ---
        # Global bin1_id ordering is proven by the monotonicity of bin1_offset
        # (checked above, O(n_bins)); a full pass over the pixel table would be
        # O(n_pixels) and this function runs on every file open.  What is left
        # is checked on three sampled windows — head, middle and tail — each of
        # which is additionally reconciled against bin1_offset, so an interior
        # scramble no longer slips through the way it did with head/tail-only
        # sampling.
        if n_pixels > 0:
            window = min(PIXEL_SAMPLE_SIZE, n_pixels)
            starts = sorted(
                {0, max(0, (n_pixels - window) // 2), max(0, n_pixels - window)}
            )
            for start in starts:
                stop = min(start + window, n_pixels)
                positions = np.arange(start, stop)
                b1 = np.array(p["bin1_id"][start:stop]).astype(np.int64)
                b2 = np.array(p["bin2_id"][start:stop]).astype(np.int64)
                if np.min(b1) < 0 or np.max(b2) >= n_bins:
                    raise MetaCoolerFormatError(
                        f"pixel bin ids out of range [0, {n_bins})"
                    )
                if storage_mode != "square" and np.any(b1 > b2):
                    raise MetaCoolerFormatError(
                        "pixels are not upper-triangular (bin1_id <= bin2_id)"
                    )
                # Strictly increasing: a repeated (bin1, bin2) pair is a
                # duplicate pixel, which cooler sums on some reads and not on
                # others, so the stored counts no longer define one matrix.
                if len(b1) > 1 and not np.all(
                    (b1[:-1] < b1[1:]) | ((b1[:-1] == b1[1:]) & (b2[:-1] < b2[1:]))
                ):
                    raise MetaCoolerFormatError(
                        "pixels are not sorted by (bin1_id, bin2_id), or a "
                        "(bin1_id, bin2_id) pair is repeated"
                    )
                # Each sampled row must live in the bin1_offset slot that claims
                # it; this ties the (cheap) index check to the actual pixels.
                owner = np.searchsorted(bin1_offset, positions, side="right") - 1
                if not np.array_equal(owner, b1):
                    bad = int(np.argmax(owner != b1))
                    raise MetaCoolerFormatError(
                        f"bin1_offset disagrees with pixels/bin1_id at row "
                        f"{start + bad}: the index places it in bin "
                        f"{int(owner[bad])}, the pixel says {int(b1[bad])}"
                    )


# ------------------------------------------------------------------
# MMetaCooler helpers
# ------------------------------------------------------------------


def is_mmetacooler(uri: str | Path) -> bool:
    """Return ``True`` if *uri* points to a valid ``.mmetacool`` file.

    Checks for the ``format`` attribute equal to ``HDF5::MMetaCooler`` and
    the presence of the top-level ``entries/`` group.

    Args:
        uri: Path to the file to test.

    Returns:
        ``True`` if the file is a valid mmetacool collection, ``False``
        otherwise (including if the file does not exist or cannot be
        opened).
    """
    path = Path(uri)
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as f:
            fmt = f.attrs.get("format")
            if fmt != MMETACOOLER_FORMAT:
                return False
            return "entries" in f
    except Exception:
        return False


def validate_mmetacooler(uri: str | Path) -> None:
    """Validate that *uri* is a well-formed ``.mmetacool`` file.

    Checks top-level attributes and the ``entries/`` group, then runs the
    full :func:`validate` on every entry.  Entry failures are collected and
    reported together rather than stopping at the first one, so a
    multi-resolution file is diagnosed in a single pass.

    Args:
        uri: Path to the ``.mmetacool`` file to validate.

    Raises:
        MetaCoolerFormatError: If the structure is invalid, or if any entry
            is not a well-formed metacooler.
    """
    path = Path(uri)
    if not path.exists():
        raise MetaCoolerFormatError(f"File not found: {path}")

    errors: list[str] = []
    to_validate: list[str] = []

    with h5py.File(path, "r") as f:
        fmt = f.attrs.get("format")
        if fmt != MMETACOOLER_FORMAT:
            raise MetaCoolerFormatError(
                f"Expected format='{MMETACOOLER_FORMAT}', got {fmt!r}"
            )

        fv = f.attrs.get("format-version")
        if fv is None:
            raise MetaCoolerFormatError("Missing 'format-version' attribute.")

        if "entries" not in f:
            raise MetaCoolerFormatError("Missing required key: /entries")

        grp = f["entries"]
        for name in grp:
            entry = grp[name]
            # Each entry must have the core metacooler tables
            missing = [
                tbl
                for tbl in ("bins", "pixels", "contigs", "genomes", "indexes")
                if tbl not in entry
            ]
            if missing:
                errors.append(
                    f"Entry {name!r} missing required group: {', '.join(missing)}"
                )
            else:
                to_validate.append(name)

    # Run the per-entry validation with the collection closed: validate()
    # reopens the file itself.
    for name in to_validate:
        try:
            validate(f"{path}::entries/{name}")
        except MetaCoolerFormatError as exc:
            errors.append(f"Entry {name!r}: {exc}")

    if errors:
        raise MetaCoolerFormatError(
            f"{path} is not a well-formed .mmetacool file:\n  - "
            + "\n  - ".join(errors)
        )


def list_entry_names(uri: str | Path) -> list[str]:
    """Return the bare names of every entry in an ``.mmetacool`` file.

    Single implementation of the ``/entries`` scan.  The two public spellings
    differ only in their return convention and both derive from this:
    :func:`list_metacoolers` yields URI-ready group paths, while
    :func:`metacooler.io.mmetacool.list_entries` yields bare names.

    Args:
        uri: Path to the ``.mmetacool`` file.

    Returns:
        Entry names in HDF5 iteration order, or ``[]`` if the file has no
        ``entries`` group.
    """
    with h5py.File(Path(uri), "r") as f:
        if "entries" not in f:
            return []
        return list(f["entries"].keys())


def list_metacoolers(uri: str | Path) -> list[str]:
    """Return HDF5 group paths of all metacooler entries in an mmetacool file.

    Returns paths of the form ``/entries/<name>`` for each entry in the
    ``/entries`` group.

    Args:
        uri: Path to the ``.mmetacool`` file.

    Returns:
        List of HDF5 group paths (e.g. ``["/entries/sample1",
        "/entries/sample2"]``).  Returns an empty list if the file has no
        ``entries`` group.
    """
    return [f"/entries/{name}" for name in list_entry_names(uri)]


# ------------------------------------------------------------------
# Integrity report
# ------------------------------------------------------------------


class CheckResult(NamedTuple):
    """One line of a :func:`check` report.

    Attributes:
        name: Short name of the check (e.g. ``"structure"``).
        ok: Whether the check passed.
        detail: What was found — the reason on failure, the evidence on success.
    """

    name: str
    ok: bool
    detail: str


def _resolve_node(froot: h5py.File, group_path: str) -> h5py.Group:
    """Return the group named by a URI group path inside an open file."""
    key = group_path.lstrip("/")
    return froot[key] if key else froot


def _core_dtype_problems(path: Path, group_path: str) -> list[str]:
    """Report core columns whose on-disk dtype kind is unusable.

    :func:`validate` checks that the columns exist and that their *values* are
    consistent, but never looks at their storage type, so a ``bins/start``
    written as float64 or a ``genome_id`` written as a string passes it and then
    fails much later, inside cooler or numpy indexing.

    Args:
        path: Path to the HDF5 file.
        group_path: URI group path of the map to inspect (``"/"`` for a
            top-level ``.metacool``).

    Returns:
        One message per offending column; empty when every core column is
        stored in an acceptable type.
    """
    problems: list[str] = []
    with h5py.File(path, "r") as froot:
        node = _resolve_node(froot, group_path)
        for column, kinds in CORE_COLUMN_KINDS.items():
            if column not in node:
                continue  # absence is validate()'s business, not this check's
            dataset = node[column]
            if not isinstance(dataset, h5py.Dataset):
                problems.append(
                    f"/{column} is a {type(dataset).__name__}, not a table column"
                )
                continue
            if dataset.dtype.kind not in kinds:
                problems.append(
                    f"/{column} is {dataset.dtype} (kind {dataset.dtype.kind!r}); "
                    f"expected one of the kinds {list(kinds)}"
                )
    return problems


def _deep_pixel_problems(path: Path, group_path: str) -> list[str]:
    """Scan the *whole* pixel table and report every invariant it breaks.

    :func:`validate` samples three windows of the pixel table because it runs on
    every ``MetaCooler`` construction; this is the O(nnz) version behind
    ``check --deep``, and it additionally reconciles ``attrs['sum']`` with the
    counts actually stored — the one declared total nothing else checks.

    Args:
        path: Path to the HDF5 file.
        group_path: URI group path of the map to scan.

    Returns:
        At most one message per broken invariant, in the order encountered.
    """
    problems: list[str] = []
    with h5py.File(path, "r") as froot:
        node = _resolve_node(froot, group_path)
        pixels = node["pixels"]
        n_bins = _table_length(node["bins"], "bins")
        n_pixels = _table_length(pixels, "pixels")
        bin1_offset = np.array(node["indexes"]["bin1_offset"]).astype(np.int64)
        square = _attr(node, "storage-mode") == "square"
        integer_counts = np.issubdtype(pixels["count"].dtype, np.integer)

        reported: set[str] = set()

        def once(kind: str, message: str) -> None:
            """Record the first occurrence of each kind of breakage."""
            if kind not in reported:
                reported.add(kind)
                problems.append(message)

        # Only integer counts are totalled: attrs['sum'] is an integer contact
        # count, so there is nothing to reconcile a float weight column against.
        total = 0
        previous: tuple[int, int] | None = None
        for lo in range(0, n_pixels, DEEP_PIXEL_CHUNK):
            hi = min(lo + DEEP_PIXEL_CHUNK, n_pixels)
            b1 = np.array(pixels["bin1_id"][lo:hi]).astype(np.int64)
            b2 = np.array(pixels["bin2_id"][lo:hi]).astype(np.int64)
            if integer_counts:
                total += int(np.array(pixels["count"][lo:hi]).astype(np.int64).sum())

            if (
                np.min(b1) < 0
                or np.min(b2) < 0
                or np.max(b1) >= n_bins
                or np.max(b2) >= n_bins
            ):
                once("range", f"pixel bin ids out of range [0, {n_bins})")
            if not square and np.any(b1 > b2):
                bad = int(np.argmax(b1 > b2))
                once(
                    "triangle",
                    f"pixels are not upper-triangular (bin1_id <= bin2_id): row "
                    f"{lo + bad} is ({int(b1[bad])}, {int(b2[bad])})",
                )
            # Prepend the last row of the previous chunk so the ordering is
            # checked *across* chunk boundaries too, not just within them.
            j1 = np.concatenate(([previous[0]], b1)) if previous else b1
            j2 = np.concatenate(([previous[1]], b2)) if previous else b2
            ordered = (j1[:-1] < j1[1:]) | ((j1[:-1] == j1[1:]) & (j2[:-1] < j2[1:]))
            if len(j1) > 1 and not np.all(ordered):
                k = int(np.argmax(~ordered))  # compares rows k and k + 1 of j1/j2
                bad = k + (0 if previous else 1)
                if j1[k] == j1[k + 1] and j2[k] == j2[k + 1]:
                    message = (
                        f"pixels repeat the (bin1_id, bin2_id) pair of the row "
                        f"before at row {lo + bad} (a duplicate pixel)"
                    )
                else:
                    message = (
                        f"pixels are not sorted by (bin1_id, bin2_id) at row {lo + bad}"
                    )
                once("order", message)
            owner = np.searchsorted(bin1_offset, np.arange(lo, hi), side="right") - 1
            if not np.array_equal(owner, b1):
                bad = int(np.argmax(owner != b1))
                once(
                    "index",
                    f"bin1_offset disagrees with pixels/bin1_id at row {lo + bad}: "
                    f"the index places it in bin {int(owner[bad])}, the pixel says "
                    f"{int(b1[bad])}",
                )
            if hi > lo:
                previous = (int(b1[-1]), int(b2[-1]))

        declared = _int_attr(node, "sum") if integer_counts else None
        if declared is not None and declared != total:
            problems.append(
                f"attrs['sum'] ({declared}) != the sum of pixels/count ({total})"
            )
    return problems


def check(uri: str | Path, deep: bool = False) -> list[CheckResult]:
    """Check the integrity of a ``.metacool`` or ``.mmetacool`` file.

    This is the reporting front end to the validators: it runs
    :func:`validate` / :func:`validate_mmetacooler` and turns their
    all-or-nothing exception into one result row per check, so a caller sees
    every problem at once instead of only the first.  Two checks go beyond
    ``validate``, which is deliberately cheap because it runs on every file
    open: the on-disk dtypes of the core columns, and — under *deep* — a full
    scan of the pixel table.

    Checks are ordered cheapest first and stop early only when a later one
    could not be meaningful (an unreadable file, or one that is not a
    metacooler at all).  A collection is checked entry by entry.

    Args:
        uri: Path to the file, or a ``file.mmetacool::entries/entry`` URI.  A
            bare ``.mmetacool`` path checks every entry.
        deep: Also scan the whole pixel table (O(nnz)) rather than the three
            sampled windows ``validate`` uses, and reconcile ``attrs['sum']``
            with the stored counts.

    Returns:
        The report, one :class:`CheckResult` per check.  The file is sound iff
        every result has ``ok=True``; the function itself does not raise for a
        malformed file.
    """
    from metacooler.uri import parse_uri

    results: list[CheckResult] = []

    try:
        path, group_path = parse_uri(str(uri))
    except ValueError as exc:
        return [CheckResult("readable", False, str(exc))]
    if not path.exists():
        return [CheckResult("readable", False, f"File not found: {path}")]
    try:
        with h5py.File(path, "r") as froot:
            if group_path.lstrip("/") and group_path.lstrip("/") not in froot:
                return [
                    CheckResult(
                        "readable",
                        False,
                        f"Group {group_path!r} not found in {path}; top-level keys: "
                        f"{sorted(froot.keys())[:20]}",
                    )
                ]
    except OSError as exc:
        return [CheckResult("readable", False, f"Cannot open {path} as HDF5: {exc}")]
    results.append(CheckResult("readable", True, f"{path} opens as HDF5"))

    # --- what kind of file is this? ---
    collection = group_path == "/" and is_mmetacooler(path)
    if collection:
        entries = list_metacoolers(path)
        results.append(
            CheckResult(
                "format",
                True,
                f"{MMETACOOLER_FORMAT} collection with {len(entries)} "
                f"{'entry' if len(entries) == 1 else 'entries'}",
            )
        )
        targets = [(name.rsplit("/", 1)[-1], name) for name in entries]
    elif is_metacooler(uri):
        results.append(CheckResult("format", True, METACOOLER_FORMAT))
        targets = [("", group_path)]
    else:
        import cooler.fileops

        cooler_uri = str(path) if group_path == "/" else f"{path}::{group_path}"
        extra = (
            " It is a plain cooler file: convert it with `metacooler create cool`."
            if cooler.fileops.is_cooler(cooler_uri)
            else ""
        )
        results.append(
            CheckResult(
                "format",
                False,
                f"Not a metacool or mmetacool file: no {METACOOLER_FORMAT!r} "
                f"metacooler-format attribute with the expected groups.{extra}",
            )
        )
        return results

    # --- structure: the existing validators, reported instead of raised ---
    def _structure() -> str:
        if collection:
            validate_mmetacooler(path)
            return f"every entry conforms to {METACOOLER_FORMAT}"
        validate(uri)
        return "required groups, columns, indexes and counts are consistent"

    try:
        results.append(CheckResult("structure", True, _structure()))
    except (MetaCoolerFormatError, OSError, KeyError) as exc:
        results.append(CheckResult("structure", False, str(exc)))

    # --- per-map checks ---
    import cooler

    for label, target_group in targets:
        suffix = f" [{label}]" if label else ""
        problems = _core_dtype_problems(path, target_group)
        results.append(
            CheckResult(
                f"dtypes{suffix}",
                not problems,
                "; ".join(problems)
                if problems
                else f"{len(CORE_COLUMN_KINDS)} core columns are stored in usable types",
            )
        )

        target_uri = str(path) if target_group == "/" else f"{path}::{target_group}"
        try:
            clr = cooler.Cooler(target_uri)
            n_bins, _ = clr.shape
            results.append(
                CheckResult(
                    f"cooler-compat{suffix}",
                    True,
                    f"cooler.Cooler() opens it: {len(clr.chromnames)} chroms "
                    f"(= contigs), {n_bins} bins, {clr.info.get('nnz')} pixels",
                )
            )
        except Exception as exc:
            # Deliberately broad: *any* way in which cooler fails to open the
            # file breaks the dual-identity invariant, and the point of this
            # check is to report that rather than propagate it.
            results.append(
                CheckResult(
                    f"cooler-compat{suffix}",
                    False,
                    f"cooler.Cooler() cannot open it — the central metacooler "
                    f"invariant is broken: {type(exc).__name__}: {exc}",
                )
            )

        if deep:
            try:
                problems = _deep_pixel_problems(path, target_group)
            except (OSError, KeyError, ValueError, MetaCoolerFormatError) as exc:
                problems = [f"the pixel scan could not complete: {exc}"]
            results.append(
                CheckResult(
                    f"pixels{suffix}",
                    not problems,
                    "; ".join(problems)
                    if problems
                    else "every pixel is in range, upper-triangular, sorted, "
                    "indexed by bin1_offset, and summed by attrs['sum']",
                )
            )

    return results


def read_table(uri: str | Path, group: str) -> pd.DataFrame:
    """Read an HDF5 table group as a pandas DataFrame.

    Matches cooler behaviour: fixed-length byte strings are decoded to
    unicode, and column order follows the HDF5 dataset creation order.

    Args:
        uri: Path to the HDF5 file.
        group: HDF5 group path containing the table datasets (e.g.
            ``"genomes"`` or ``"entries/sample1/contigs"``).

    Returns:
        DataFrame where each HDF5 dataset in *group* becomes a column.
        Byte-string columns are decoded to ``str``; all other dtypes are
        returned as-is.
    """
    path = Path(uri)
    with h5py.File(path, "r", track_order=True) as f:
        g = f[group]
        columns = {}
        for key in g:
            ds = g[key]
            data = ds[:]
            # Decode fixed-length byte strings to str (like cooler). Use an
            # explicit UTF-8 decode: ``astype("U")`` decodes as ASCII and raises
            # UnicodeDecodeError on any non-ASCII byte (e.g. a "é" in a name),
            # which would break every read of a file that metacooler itself
            # wrote.
            if data.dtype.kind == "S":
                data = np.char.decode(data, "utf-8")
            elif data.dtype == object and len(data) > 0 and isinstance(data[0], bytes):
                data = np.array([x.decode("utf-8") for x in data])
            columns[key] = data
    return pd.DataFrame(columns)
