"""Shared I/O utilities used across multiple metacooler modules."""

from __future__ import annotations

import contextlib
import gzip
import importlib
import itertools
import logging
import os
import re
import types
from pathlib import Path
from typing import TYPE_CHECKING, Literal, cast

import h5py
import numpy as np
import pandas as pd

from metacooler._version import __version__
from metacooler.schema import (
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
    MetaCoolerFormatError,
    ascii_violations,
    format_ascii_violations,
)

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Iterator, Mapping, Sequence

    import numpy.typing as npt

    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

# pandas' placeholder name for a header cell that is empty — what a trailing tab
# on the header line produces.
_UNNAMED_COLUMN_RE = re.compile(r"Unnamed: \d+")
# How many ragged line numbers to name in the warning.
_MAX_REPORTED_RAGGED_LINES = 5


def _raw_field_counts(path: Path) -> list[int]:
    """Count the tab-separated fields on every non-blank line of *path*.

    Counted from the raw bytes rather than read off the parsed frame, because
    the parsed frame does not remember: ``pandas.read_csv`` pads a row that has
    fewer fields than the header, and under ``keep_default_na=False`` (which
    every metacooler table reader uses, so that a genome literally called
    ``"NA"`` survives) the pad is an empty string — indistinguishable from a
    field that was present and empty.

    Args:
        path: Table to scan.  Transparently decompressed when it ends in ``.gz``.

    Returns:
        One field count per non-blank line, in file order.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:
        return [line.count(b"\t") + 1 for line in fh if line.strip()]


def _drop_trailing_unnamed_columns(df: pd.DataFrame, kind: str) -> pd.DataFrame:
    """Drop trailing header-less columns that no row fills.

    A trailing tab on the header line is an artifact of whatever wrote the file,
    not a column: pandas turns it into an ``Unnamed: N`` column of empty
    strings, which metacooler would then write out as a junk ``/genomes`` or
    ``/contigs`` HDF5 dataset.  Only *trailing* ones are dropped — an empty
    column in the middle is still a positional placeholder, and removing it
    would shift every column after it.

    Args:
        df: Frame just read with ``header=0``.
        kind: Human-readable table name used in the log message.

    Returns:
        *df*, or a narrower view of it.
    """
    n = df.shape[1]
    while n > 1 and _UNNAMED_COLUMN_RE.fullmatch(str(df.columns[n - 1])):
        if not (df.iloc[:, n - 1] == "").all():
            break
        n -= 1
    if n == df.shape[1]:
        return df
    logger.debug(
        "%s: dropped %d unnamed, empty trailing column(s) (trailing tab on the "
        "header line).",
        kind,
        df.shape[1] - n,
    )
    return df.iloc[:, :n]


def _warn_ragged_rows(
    counts: list[int], expected: int, kind: str, path: Path, first_line: int
) -> None:
    """Warn about individual lines narrower than *expected*.

    Args:
        counts: Field count of each data line.
        expected: Field count the table should have.
        kind: Human-readable table name.
        path: Table path, for the message.
        first_line: 1-based file line number of ``counts[0]``.
    """
    short = [i + first_line for i, c in enumerate(counts) if c < expected]
    if not short:
        return
    shown = ", ".join(str(n) for n in short[:_MAX_REPORTED_RAGGED_LINES])
    logger.warning(
        "%s (%s): %d line(s) have fewer than %d tab-separated field(s) "
        "(line(s) %s%s); their missing trailing fields were read as empty.",
        kind,
        path,
        len(short),
        expected,
        shown,
        ", ..." if len(short) > _MAX_REPORTED_RAGGED_LINES else "",
    )


def read_user_table(
    path: str | Path,
    *,
    kind: str,
    header: bool = True,
) -> pd.DataFrame:
    """Read a user-supplied TSV, refusing to silently misalign a ragged one.

    Single entry point for every user-facing text table metacooler parses
    (genomes tables, contig tables), because ``pandas.read_csv`` resolves a
    header/data field-count mismatch silently and in a different direction
    depending on which side is wider:

    - Header **wider** than the rows: the trailing columns are padded, so every
      value lands under a column name to the left of the one it belongs to.  A
      real MetaTOR genome table shipped an 11-name header over 8-field rows,
      and its GC content was stored as ``Size``, its completeness as
      ``n_contigs`` and its taxonomy as ``GC`` — all without a single warning.
    - Rows **wider** than the header: pandas promotes the surplus leading
      field(s) into the DataFrame index instead, shifting the values the other
      way.  This is the same trap already documented for short ``names=`` lists
      in :func:`metacooler.io.contigs.parse_contigs_table`, reached through the
      header rather than through ``names``.

    Either way a mislabelled table parses cleanly, so the misalignment only
    surfaces much later as impossible metadata.  Both shapes are rejected here
    instead.

    Args:
        path: Path to the table.  ``.gz`` is handled.
        kind: Human-readable table name used in errors and warnings, e.g.
            ``"Genomes table"``.
        header: Whether the first line names the columns.  When False the frame
            keeps pandas' positional integer column labels.

    Returns:
        DataFrame of ``str`` values, with NA-like strings preserved verbatim.

    Raises:
        ValueError: If *path* has no parseable line, or if its header and data
            rows systematically disagree on the number of fields.
    """
    path = Path(path)
    counts = _raw_field_counts(path)
    if not counts:
        raise ValueError(f"{kind} has no non-blank lines: {path}")

    # Every read below passes keep_default_na=False so that a genome or contig
    # literally named "NA", "nan" or "None" is not turned into a missing value.
    if not header:
        widest = max(counts)
        _warn_ragged_rows(counts, widest, kind, path, first_line=1)
        return pd.read_csv(
            path, sep="\t", header=None, dtype=str, keep_default_na=False
        )

    n_header, data_counts = counts[0], counts[1:]
    widest = max(data_counts, default=0)
    if widest > n_header:
        raise ValueError(
            f"{kind} is misaligned: {path}\n"
            f"  The header line has {n_header} tab-separated field(s), but some "
            f"data row(s) have up to {widest}. The surplus leading field(s) "
            f"would be absorbed into the row index, storing every value under "
            f"the name of a column to its left. Give the header one name per "
            f"column."
        )

    df = _drop_trailing_unnamed_columns(
        pd.read_csv(path, sep="\t", header=0, dtype=str, keep_default_na=False),
        kind,
    )
    n_named = df.shape[1]
    if data_counts and widest < n_named:
        unfilled = ", ".join(repr(str(c)) for c in df.columns[widest:])
        raise ValueError(
            f"{kind} is misaligned: {path}\n"
            f"  The header names {n_named} column(s), but no data row has more "
            f"than {widest} field(s), leaving {unfilled} unfilled on every "
            f"row. metacooler will not guess whether the absent fields belong "
            f"at the start or at the end of each row: assuming the end (what "
            f"pandas does) would store every value under the name of the "
            f"column {n_named - widest} position(s) to its left. Fix the table "
            f"so that the header and the data rows agree."
        )
    _warn_ragged_rows(data_counts, n_named, kind, path, first_line=2)
    return df


def rename_key_columns(
    df: pd.DataFrame, keys: Mapping[int, str], *, kind: str
) -> pd.DataFrame:
    """Give a user table's key columns their internal names, by position.

    Every user table has one or two *key* columns that metacooler renames to
    its own internal labels (``contig_name``, ``genome_id``); every other
    column is metadata that ends up as an HDF5 dataset of the same name.  This
    is the single place that relabelling happens, because each way of doing it
    by *name* went wrong:

    - ``df.rename(columns={"genome_name": "genome_id"})`` on a table that also
      carries its own ``genome_id`` column yields two columns of that label,
      and ``df["genome_id"]`` then returns a frame -- the contig->genome map
      came out as ``{contig: "genome_id"}`` and every contig fell into
      ``__unbinned__`` without an error.
    - A headerless table keeps pandas' integer labels on its extra columns,
      which ``h5py`` cannot use as dataset names (a ``TypeError`` raised only
      after the whole pairs pass).
    - A column name containing ``/`` is an HDF5 *path*: ``create_dataset``
      quietly builds a nested group, and the file then opens in neither
      cooler nor metacooler.

    So key columns are relabelled by position, and every other label is made
    a string, freed of ``/``, and suffixed if it would collide -- each change
    logged, never silent.

    Args:
        df: Table as read by :func:`read_user_table`.
        keys: Mapping of 0-based column position to internal label.
        kind: Human-readable table name for log messages.

    Returns:
        A copy of *df* with the new column labels.
    """
    labels: list[str] = []
    taken = set(keys.values())
    for i, raw in enumerate(df.columns):
        if i in keys:
            labels.append(keys[i])
            continue
        name = f"col{i + 1}" if isinstance(raw, int | np.integer) else str(raw)
        if "/" in name:
            fixed = name.replace("/", "_")
            logger.warning(
                "%s column %r contains '/', which HDF5 reads as a group path; "
                "stored as %r.",
                kind,
                name,
                fixed,
            )
            name = fixed
        if name in taken:
            fixed = f"{name}_col{i + 1}"
            logger.warning(
                "%s column %d is named %r, which collides with a column "
                "metacooler uses internally; renamed it to %r.",
                kind,
                i + 1,
                name,
                fixed,
            )
            name = fixed
        taken.add(name)
        labels.append(name)
    out = df.copy()
    out.columns = labels
    return out


def _str_to_hdf5(arr: np.ndarray) -> np.ndarray:
    """Convert a string/object NumPy array to fixed-length bytes for HDF5 storage.

    h5py cannot write arrays with ``dtype=object`` or pandas ``StringDtype``
    directly.  This function converts them to NumPy fixed-length byte-string
    arrays (``dtype="S<n>"``) suitable for ``h5py.Dataset.create_dataset()``.

    Args:
        arr: Input array.  Arrays with dtype ``"O"``, ``"U"``, ``"S"``, or
            pandas ``"string"`` are converted; all other dtypes are returned
            unchanged via ``numpy.asarray()``.

    Returns:
        A NumPy array with dtype ``S<n>`` where *n* is the length of the
        longest encoded string.  Non-string inputs are returned as-is.
    """
    if arr.dtype.kind == "S":
        return arr
    if arr.dtype.kind in ("O", "U") or str(arr.dtype) == "string":
        fast = _ascii_to_hdf5(arr)
        if fast is not None:
            return fast
        str_arr = ["" if pd.isna(v) else str(v) for v in arr]
        encoded = [s.encode("utf-8") for s in str_arr]
        # Size the fixed-length dtype by ENCODED BYTE length, not character
        # count: a multi-byte character (e.g. "é") is longer as UTF-8 than as a
        # Python str, so sizing by len(str) would silently truncate the bytes
        # (and could collide two distinct names).
        max_len = max((len(b) for b in encoded), default=1)
        return np.array(encoded, dtype=f"S{max_len}")
    return np.asarray(arr)


def _ascii_to_hdf5(arr: np.ndarray) -> np.ndarray | None:
    """:func:`_str_to_hdf5` of an all-ASCII column, vectorised; ``None`` otherwise.

    NumPy's object-to-bytes cast does what the per-value loop does -- ``str()``
    of each value, ASCII-encoded -- in C, three times faster on two million
    contig names; it raises on the first non-ASCII value, which is left to the
    UTF-8 loop.  The width is then trimmed to the longest value, as the loop
    sizes it (NumPy widens to 32 bytes for some scalar types).
    """
    values = np.asarray(arr, dtype=object)
    missing = pd.isna(values)
    if missing.any():
        values = values.copy()
        values[missing] = ""
    try:
        out = values.astype("S")
    except UnicodeEncodeError:
        return None
    width = max(int(np.char.str_len(out).max(initial=0)), 1)
    return out if out.dtype.itemsize == width else out.astype(f"S{width}")


def to_upper_triangle(
    bin1: np.ndarray, bin2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Reorder each ``(bin1, bin2)`` pair so that ``bin1 <= bin2``.

    metacooler stores pixels in ``symmetric-upper`` form, so every writer that
    remaps bin ids has to restore that invariant afterwards.  This used to be
    open-coded five times in two different spellings — a NumPy
    ``a[swap], b[swap] = b[swap], a[swap]`` dance (which relies on fancy
    indexing returning copies) and a pandas ``.loc[...] = ....values`` swap
    (which is only correct because of the trailing ``.values``: without it,
    pandas aligns on column *labels* and the assignment is a silent no-op).
    Both are correct today and both are one careless edit away from not being.
    ``minimum``/``maximum`` says the same thing with nothing to get wrong.

    Args:
        bin1: Row bin ids.
        bin2: Column bin ids, same shape as *bin1*.

    Returns:
        ``(lo, hi)`` — new arrays with ``lo <= hi`` elementwise.  Inputs are
        not modified.
    """
    return np.minimum(bin1, bin2), np.maximum(bin1, bin2)


#: Storage filters for every dataset metacooler writes itself (``/genomes``,
#: ``/contigs``, ``bins/chrom``, ``chrom_offset``, ``genome_offset``, the QC
#: columns): cooler's own defaults -- gzip level 6 behind the shuffle filter,
#: chunked by h5py -- so these tables compress like the ``/bins`` and
#: ``/pixels`` cooler writes beside them.  Filters are transparent to readers.
H5_STORAGE: dict[str, object] = {
    "compression": "gzip",
    "compression_opts": 6,
    "shuffle": True,
}


def write_column(
    group: h5py.Group, name: str, data: object, dtype: npt.DTypeLike | None = None
) -> h5py.Dataset:
    """Create the dataset *name* in *group* from *data*, compressed.

    The single way metacooler writes a table column of its own: with
    :data:`H5_STORAGE`, and with text converted to fixed-length bytes (h5py
    cannot write ``object`` or ``StringDtype`` arrays).

    Args:
        group: Table group to write into.
        name: Column name.
        data: Column values.
        dtype: Optional dtype to cast to.

    Returns:
        The new dataset.
    """
    arr = _str_to_hdf5(np.asarray(data))
    if dtype is not None:
        arr = arr.astype(dtype)
    return group.create_dataset(name, data=arr, **H5_STORAGE)


#: A pixel's sort key is ``bin1 << 32 | bin2``: one int64 whose integer order is
#: ``(bin1, bin2)`` order, as long as every bin id is below ``2**31``.
_PIXEL_KEY_SHIFT = 32
_PIXEL_KEY_MASK = (1 << _PIXEL_KEY_SHIFT) - 1
_MAX_KEYED_BIN = 1 << 31


def pixel_keys(bin1: np.ndarray, bin2: np.ndarray) -> np.ndarray:
    """Fold ``(bin1, bin2)`` into one ``int64`` that sorts as the pair does.

    Sorting, deduplicating and merging on this one key is three to four times
    faster, and needs half the memory, of a pandas ``groupby`` over two columns.

    Args:
        bin1: Row bin ids.
        bin2: Column bin ids, same length.

    Returns:
        ``int64`` keys.

    Raises:
        OverflowError: If a bin id is negative or does not fit in 31 bits, where
            the key would no longer order (or decode) the pair.
    """
    b1 = np.asarray(bin1).astype(np.int64, copy=False)
    b2 = np.asarray(bin2).astype(np.int64, copy=False)
    if len(b1):
        lo = min(int(np.min(b1)), int(np.min(b2)))
        hi = max(int(np.max(b1)), int(np.max(b2)))
        if lo < 0 or hi >= _MAX_KEYED_BIN:
            raise OverflowError(
                f"bin ids must lie in [0, 2**31) to be keyed; got [{lo}, {hi}]"
            )
    return (b1 << _PIXEL_KEY_SHIFT) | b2


def pixels_from_keys(keys: np.ndarray, counts: np.ndarray) -> pd.DataFrame:
    """Unfold :func:`pixel_keys` back into a ``bin1_id, bin2_id, count`` frame.

    ``copy=False`` keeps the three columns as separate blocks: by default pandas
    consolidates same-dtype columns into one 2-D block, a copy of every pixel.
    """
    return pd.DataFrame(
        {
            "bin1_id": keys >> _PIXEL_KEY_SHIFT,
            "bin2_id": keys & _PIXEL_KEY_MASK,
            "count": counts,
        },
        copy=False,
    )


def _is_sorted(values: np.ndarray) -> bool:
    return bool(len(values) < 2 or (values[1:] >= values[:-1]).all())


#: Mean length of the already-sorted runs in a key array above which a stable
#: sort (timsort, which merges existing runs) beats quicksort.  Measured on a
#: 200 k-contig map: a genome-grouped remap has runs of ~40-60 keys and sorts
#: 1.5-4x faster stably; scrambled keys have runs of ~2 and sort 2x faster
#: with quicksort.
_STABLE_SORT_MIN_RUN = 8


def _sort_kind(
    keys: np.ndarray, runs_known_sorted: bool = False
) -> Literal["stable", "quicksort"]:
    """Pick the faster ``np.argsort`` kind for *keys* (see above)."""
    if runs_known_sorted:
        return "stable"
    descents = int(np.count_nonzero(keys[1:] < keys[:-1]))
    return (
        "stable" if len(keys) >= _STABLE_SORT_MIN_RUN * (descents + 1) else "quicksort"
    )


def reduce_sorted_keys(
    keys: np.ndarray, counts: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray]:
    """Sum the counts of equal, adjacent keys (*keys* must be sorted).

    Args:
        keys: Sorted pixel keys.
        counts: Count of each key, or ``None`` for one each (raw contacts),
            in which case the sums are run lengths, as ``int32``.

    Returns:
        ``(unique_keys, sums)``.  *sums* keeps the dtype of *counts*: never
        narrowed, as a pandas ``groupby().sum()`` did not narrow it either.
    """
    if not len(keys):
        dtype = np.int32 if counts is None else counts.dtype
        return keys, np.empty(0, dtype=dtype)
    starts = np.flatnonzero(np.concatenate(([True], keys[1:] != keys[:-1])))
    if counts is None:
        sums = np.diff(np.append(starts, len(keys))).astype(np.int32)
    elif len(starts) == len(keys):
        sums = counts
    else:
        sums = np.add.reduceat(counts, starts, dtype=counts.dtype)
    return keys[starts], sums


def aggregate_pixels(
    frames: Sequence[pd.DataFrame], *, consume: bool = False
) -> pd.DataFrame:
    """Merge pixel frames into one, summing repeated ``(bin1, bin2)``, sorted.

    The single place metacooler aggregates pixels it built itself: the frames
    are folded onto :func:`pixel_keys`, sorted (a stable sort, which runs in
    linear time over the already-sorted frames this produces, when every input
    is one), and reduced with ``np.add.reduceat``.

    Args:
        frames: Non-empty sequence of frames with ``bin1_id``, ``bin2_id`` and
            ``count`` columns (other columns are dropped).
        consume: Empty *frames* (a list) as it goes, so each frame is freed as
            soon as it has been keyed; for the large in-memory sorts.

    Returns:
        ``bin1_id``, ``bin2_id`` (``int64``) and ``count`` (the input dtype),
        sorted by ``(bin1_id, bin2_id)`` with no repeated pair.
    """
    n_frames = len(frames)
    keys_list: list[np.ndarray] = []
    counts_list: list[np.ndarray] = []
    queue = frames if consume else list(frames)
    while queue:
        f = cast("list[pd.DataFrame]", queue).pop(0)
        keys_list.append(pixel_keys(f["bin1_id"].to_numpy(), f["bin2_id"].to_numpy()))
        counts_list.append(f["count"].to_numpy())
        del f
    all_sorted = all(_is_sorted(k) for k in keys_list)
    keys = np.concatenate(keys_list) if n_frames > 1 else keys_list[0]
    counts = np.concatenate(counts_list) if n_frames > 1 else counts_list[0]
    del keys_list, counts_list
    if not (all_sorted and n_frames == 1):
        order = np.argsort(keys, kind=_sort_kind(keys, all_sorted))
        keys = keys[order]
        counts = counts[order]
        del order
    return pixels_from_keys(*reduce_sorted_keys(keys, counts))


def enforce_upper_triangle(
    pixels: pd.DataFrame, col1: str = "bin1_id", col2: str = "bin2_id"
) -> pd.DataFrame:
    """Reorder a pixel frame's bin-id columns in place so ``col1 <= col2``.

    DataFrame-shaped wrapper over :func:`to_upper_triangle`.

    Args:
        pixels: Pixel frame to fix up, modified in place.
        col1: Name of the row bin-id column.
        col2: Name of the column bin-id column.

    Returns:
        The same *pixels* object, for chaining.
    """
    lo, hi = to_upper_triangle(pixels[col1].to_numpy(), pixels[col2].to_numpy())
    pixels[col1] = lo
    pixels[col2] = hi
    return pixels


def chrom_codes(bins: pd.DataFrame) -> np.ndarray:
    """Return the integer contig id of every row of a *bins* frame.

    Single spelling for an idiom that was repeated ~12 times across the
    codebase as ``bins["chrom"].cat.codes.values``.  That form assumes cooler
    handed back a *categorical* ``chrom`` column and raises ``AttributeError``
    the moment it did not (a plain integer column read straight from HDF5, or a
    frame that has been through a groupby); this handles both.

    Args:
        bins: A bins DataFrame with a ``chrom`` column — either cooler's
            categorical form or plain integer contig ids.

    Returns:
        A fresh ``int64`` array of contig ids, safe for the caller to mutate.

    Raises:
        KeyError: If *bins* has no ``chrom`` column.
    """
    col = bins["chrom"]
    if isinstance(col.dtype, pd.CategoricalDtype):
        return col.cat.codes.to_numpy(dtype=np.int64)
    return col.to_numpy(dtype=np.int64)


def check_output_path(
    path: str | Path,
    force: bool = False,
    *,
    inputs: Iterable[str | Path | None] = (),
    defer_removal: bool = False,
) -> None:
    """Refuse to overwrite an existing output file unless *force* is set.

    Single home for the "output exists" guard shared by the CLI commands and
    the ops layer, so every command reports the same message and every command
    that accepts ``--force`` actually honours it.

    Args:
        path: Output path about to be written.
        force: If True, an existing file may be replaced instead of raising.
        inputs: Paths or ``file::group`` URIs the command is about to read.
            An output that is the same file as any of them is refused even
            under *force*: removing it would delete the input before it is
            read (``binnify a.metacool -o a.metacool -f`` used to do exactly
            that, then fail with "File not found").
        defer_removal: Leave an existing file in place for the caller to
            replace on success, instead of removing it now.  For writers that
            stage their output and rename it over *path* at the end, so a run
            that fails does not also destroy the previous result.

    Raises:
        ValueError: If *path* is the same file as one of *inputs*.
        FileExistsError: If *path* exists and *force* is False.
        IsADirectoryError: If *path* exists and is a directory.
    """
    path = Path(path)
    _refuse_input_as_output(path, inputs)
    if not path.exists():
        return
    if path.is_dir():
        raise IsADirectoryError(f"Output path {str(path)!r} is a directory.")
    if not force:
        raise FileExistsError(
            f"Output file {str(path)!r} already exists. Use --force to overwrite."
        )
    if not defer_removal:
        path.unlink()


def _refuse_input_as_output(path: Path, inputs: Iterable[str | Path | None]) -> None:
    """Raise if *path* names the same file as any of *inputs*.

    Args:
        path: Output path about to be written.
        inputs: Input paths or URIs; ``None`` and ``"-"`` (stdin) are skipped.

    Raises:
        ValueError: If *path* and an input resolve to the same file.
    """
    from metacooler.uri import parse_uri

    if not path.exists():
        return
    for raw in inputs:
        if raw is None or str(raw) == "-":
            continue
        file_path, _ = parse_uri(str(raw))
        if file_path.exists() and os.path.samefile(file_path, path):
            raise ValueError(
                f"Output {str(path)!r} is the same file as the input "
                f"{str(raw)!r}; writing it would destroy the input before it "
                "is read. Choose a different output path."
            )


def check_output_not_an_input(output: str | Path, inputs: Sequence[str | Path]) -> None:
    """Refuse an output path that names one of the command's own inputs.

    Must run *before* :func:`check_output_path`: under ``--force`` that guard
    deletes an existing output up front, so ``-o`` pointing at an input would
    destroy the input before a single byte of it had been read, and the command
    would then fail on the file it had just removed.

    Args:
        output: Output path about to be written.
        inputs: Input paths or ``file::group`` URIs; only the file part counts.

    Raises:
        ValueError: If *output* resolves to the same file as any input.
    """
    from metacooler.uri import parse_uri

    target = Path(output).resolve()
    for uri in inputs:
        if parse_uri(str(uri))[0].resolve() == target:
            raise ValueError(
                f"Output {str(output)!r} is also the input {str(uri)!r}; choose "
                "another output path. An input is never overwritten."
            )


def resolve_count_dtype(dtype: np.dtype | str, *, context: str = "pixels") -> np.dtype:
    """Pick a lossless on-disk dtype for a ``pixels/count`` column being copied.

    Blindly casting to ``int32`` silently wraps counts written by
    :mod:`metacooler.ops.binnify` (which deliberately writes ``int64``) and
    floors float counts towards zero.  This preserves the source width instead.

    Args:
        dtype: dtype of the source ``count`` column.
        context: Short description of the caller, used in the float warning.

    Returns:
        ``int32`` for narrow integer sources, ``int64`` for wide or unsigned
        ones, and ``float64`` for float sources (with a warning).

    Raises:
        ValueError: If *dtype* is neither integer, boolean, nor float.
    """
    dtype = np.dtype(dtype)
    if dtype.kind == "f":
        logger.warning(
            "%s: source 'count' column is %s; keeping float64 to avoid "
            "truncating fractional counts.",
            context,
            dtype,
        )
        return np.dtype(np.float64)
    if dtype.kind == "b":
        return np.dtype(np.int32)
    if dtype.kind not in ("i", "u"):
        raise ValueError(f"Unsupported 'count' dtype {dtype!r} in {context}.")
    if dtype.itemsize >= 8 or (dtype.kind == "u" and dtype.itemsize >= 4):
        return np.dtype(np.int64)
    return np.dtype(np.int32)


def check_ascii_names(genomes: pd.DataFrame, contigs: pd.DataFrame) -> None:
    """Refuse to write name tables that cooler would not be able to read.

    The one hard rule metacooler places on the text it stores.  Everything
    else about a ``.metacool`` file is metacooler's business, but ``/contigs``
    doubles as cooler's ``/chroms``, and cooler assumes ASCII there — see
    :func:`metacooler.schema.ascii_violations` for exactly which columns that
    reaches and why.

    Callers run this twice: once as soon as the tables exist, so a bad name
    costs seconds rather than a whole binning pass, and once inside
    :func:`_finalize_metacool_file`, so that no writer — present or future —
    can put an unreadable file on disk.

    Raises:
        MetaCoolerFormatError: If any checked column holds non-ASCII text.
    """
    violations = ascii_violations(genomes, contigs)
    if violations:
        raise MetaCoolerFormatError(
            "Refusing to write a metacool file that cooler could not open — "
            "these names are not ASCII:\n" + format_ascii_violations(violations)
        )


def infer_binsize(bins: pd.DataFrame) -> int | None:
    """Infer a uniform bin width, or ``None`` if the bins are not uniform.

    Reproduces ``cooler.util.get_binsize`` -- the width shared by every bin
    except the last of each contig, or ``None`` when more than one width is in
    play -- without its per-contig ``groupby`` loop, which costs ~50 us per
    contig (measured) and so minutes on a metagenome.

    Assumes each contig's rows are contiguous, which holds for every bins table
    metacooler hands to cooler (they are sorted by ``chrom`` then ``start``).

    Args:
        bins: Bins table with ``chrom``, ``start`` and ``end`` columns.

    Returns:
        The uniform bin width, or ``None``.
    """
    chrom = bins["chrom"].to_numpy()
    n = len(chrom)
    if n == 0:
        return None
    # A contig's last bin is the row before a change of chrom, plus the final row.
    is_last = np.empty(n, dtype=bool)
    is_last[:-1] = chrom[1:] != chrom[:-1]
    is_last[-1] = True
    widths = (bins["end"].to_numpy() - bins["start"].to_numpy())[~is_last]
    if widths.size == 0:
        # Every contig holds a single bin, so no width is ever compared —
        # the same "no uniform size" answer cooler gives from an empty set.
        return None
    first = widths[0]
    return int(first) if bool((widths == first).all()) else None


# Modules in which a supported cooler keeps ``write_bins`` and ``get_binsize`` as
# module globals that its writer calls by name, newest layout first.  cooler
# 0.11 split the writer out of ``create/_create.py`` into ``create/_writer.py``;
# patching the old location there is a silent no-op, which is exactly how the
# quadratic enum build came back unnoticed on cooler 0.12.
_COOLER_WRITER_MODULES = ("cooler.create._writer", "cooler.create._create")

_warned_no_writer_module = False


def cooler_writer_module() -> types.ModuleType | None:
    """Return the cooler module whose ``write_bins``/``get_binsize`` to shim.

    Returns:
        The first module of :data:`_COOLER_WRITER_MODULES` that imports and
        defines both functions, or ``None`` if the installed cooler keeps them
        somewhere metacooler does not know about.
    """
    for name in _COOLER_WRITER_MODULES:
        try:
            module = importlib.import_module(name)
        except ImportError:
            continue
        if callable(getattr(module, "write_bins", None)) and callable(
            getattr(module, "get_binsize", None)
        ):
            return module
    return None


@contextlib.contextmanager
def cooler_writer_shims(binsize: int | None) -> Generator[None, None, None]:
    """Neutralise two ``cooler.create`` costs that scale with *contig* count.

    Both are invisible on a genome assembly and dominate everything else on a
    metagenome, and neither is reachable through ``create_cooler``'s signature:

    - ``write_bins`` types ``/bins/chrom`` as an HDF5 enum with one member per
      contig.  ``H5Tenum_insert`` rescans every member already inserted, so the
      build is quadratic: 2.3 s at 32 k contigs, extrapolating to **~2 h at
      1.9 M** (measured, x3.8 per doubling).  Worse, it is doomed -- past
      ~4 k contigs the enum header exceeds HDF5's object-header limit, the
      ``create_dataset`` raises ``ValueError``, and cooler silently retries as
      a plain int.  The whole cost buys nothing, and
      :func:`_finalize_metacool_file` rewrites that dataset as ``int32``
      regardless.
    - ``get_binsize`` re-derives the bin width with a per-contig ``groupby``
      loop over a table metacooler tiled itself and already knows the width of.
      It does not short-circuit on a uniform grid, so it runs to completion:
      ~1.6 min at 1.9 M contigs.

    Patching a dependency's module globals is not free of risk, so this only
    installs shims it can find, and restores them unconditionally.  Neither
    changes a byte of the output file.

    Args:
        binsize: The bin width to hand back in place of cooler's inference,
            typically from :func:`infer_binsize`.  ``None`` means non-uniform,
            which is what cooler would have concluded.

    Yields:
        ``None``, with the shims installed for the duration of the block.
    """
    global _warned_no_writer_module

    cooler_create = cooler_writer_module()
    if cooler_create is None:
        # A cooler that has moved these on: leave it alone and pay the cost —
        # but say so, once, because on a metagenome that cost is hours.
        if not _warned_no_writer_module:
            _warned_no_writer_module = True
            logger.warning(
                "This cooler version keeps its bin writer somewhere metacooler "
                "does not recognise, so writes will scale quadratically with "
                "the number of contigs. Please report this at "
                "https://github.com/js2264/metacooler/issues."
            )
        yield
        return
    original_write_bins = getattr(cooler_create, "write_bins")  # noqa: B009
    original_get_binsize = getattr(cooler_create, "get_binsize")  # noqa: B009

    # Signatures mirror cooler's exactly: these are installed as its own module
    # globals, and everything calling them is cooler's code, not ours.
    def write_bins_without_enum(
        grp: h5py.Group,
        bins: pd.DataFrame,
        chromnames: list[str],
        h5opts: dict,
        chrom_as_enum: bool = True,
    ) -> None:
        # cooler zips the names into a dict: over an object array, not the
        # Arrow-backed column pandas 3 builds, which iterates 4x slower.
        names = np.asarray(chromnames, dtype=object)
        original_write_bins(grp, bins, names, h5opts, chrom_as_enum=False)

    def get_binsize_precomputed(bins: pd.DataFrame) -> int | None:
        return binsize

    # Rebinding another module's globals is the whole point here; setattr
    # because the module is only known at run time.
    setattr(cooler_create, "write_bins", write_bins_without_enum)  # noqa: B010
    setattr(cooler_create, "get_binsize", get_binsize_precomputed)  # noqa: B010
    try:
        yield
    finally:
        setattr(cooler_create, "write_bins", original_write_bins)  # noqa: B010
        setattr(cooler_create, "get_binsize", original_get_binsize)  # noqa: B010


def write_cooler(
    path: str | Path,
    bins: pd.DataFrame,
    pixels: object,
    *,
    binsize: int | None = None,
    **kwargs: object,
) -> None:
    """``cooler.create_cooler`` with :func:`cooler_writer_shims` installed.

    **Every** metacooler writer must go through this rather than calling
    ``cooler.create_cooler`` directly, because the per-call cost the shims
    remove is not paid once per file.  With ``ordered=False`` -- which
    ``binnify``, ``aggregate`` and ``extract`` all use, since their chunks are
    individually sorted but not globally ordered -- cooler's
    ``create_from_unordered`` calls ``create()`` once per pixel chunk plus once
    per merge pass, and each of those calls rebuilds the doomed contig enum and
    re-derives the bin size from scratch.  At ``binnify``'s default 1 M-pixel
    chunk that is tens of calls, so a cost of ~2 h per call on a 1.9 M-contig
    metagenome turns a re-binning run into one that never finishes.

    An unordered stream also has its chunks sorted on the way in (see
    :func:`_sorted_pixel_chunks`): cooler's unordered writer is only correct
    when each of them is, and several writers stream theirs in source order.

    Args:
        path: Destination ``.cool``/``.metacool`` path.
        bins: Bins table, sorted by ``chrom`` then ``start``.
        pixels: Anything ``cooler.create_cooler`` accepts as *pixels*.
        binsize: Bin width to record, when the caller means to assert one
            regardless of what the table looks like (``create`` and ``binnify``
            do: cooler infers ``None`` as soon as every contig holds a single
            bin).  Defaults to :func:`infer_binsize` of *bins*, which is what
            cooler would have concluded.
        **kwargs: Passed through to ``cooler.create_cooler`` untouched.
    """
    import cooler

    if kwargs.get("ordered") is False and not isinstance(pixels, pd.DataFrame | dict):
        stream, was_empty = _skip_empty_frames(cast("Iterable[pd.DataFrame]", pixels))
        if was_empty:
            # cooler's unordered writer concatenates its per-chunk files and
            # fails with "No objects to concatenate" when there are none; an
            # empty table is trivially sorted, so write it on the ordered path.
            kwargs["ordered"] = True
            pixels = stream
        else:
            pixels = _sorted_pixel_chunks(
                cast("Iterable[pd.DataFrame]", stream), len(bins)
            )
    resolved = binsize if binsize is not None else infer_binsize(bins)
    with cooler_writer_shims(resolved):
        cooler.create_cooler(str(path), bins=bins, pixels=pixels, **kwargs)


def _skip_empty_frames(
    pixels: Iterable[pd.DataFrame],
) -> tuple[pd.DataFrame | Iterator[pd.DataFrame], bool]:
    """Drop empty frames from a pixel stream, reporting if nothing was left.

    Only the first non-empty frame is pulled to decide, so a lazy stream stays
    lazy.  A zero-pixel source (a valid, if unusual, metacool) reached
    ``binnify``, ``zoomify`` and ``aggregate`` as either an empty stream or one
    of empty frames, and both crashed cooler's unordered writer.

    Args:
        pixels: An iterable of pixel DataFrames.

    Returns:
        ``(stream, was_empty)``; when *was_empty*, *stream* is a single empty
        frame that ``cooler.create_cooler`` accepts.
    """
    frames = (df for df in pixels if len(df))
    first = next(frames, None)
    if first is None:
        empty = pd.DataFrame(
            {
                "bin1_id": np.array([], dtype=np.int64),
                "bin2_id": np.array([], dtype=np.int64),
                "count": np.array([], dtype=np.int32),
            }
        )
        return empty, True
    return itertools.chain([first], frames), False


def _sorted_pixel_chunks(
    pixels: Iterable[pd.DataFrame], n_bins: int
) -> Iterator[pd.DataFrame]:
    """Yield each chunk of an unordered pixel stream sorted by ``(bin1, bin2)``.

    cooler's unordered writer writes every chunk to a temporary cooler and
    merges those, reading each one through its ``bin1_offset`` index, one
    ``mergebuf`` (20 M records) epoch at a time.  Before cooler 0.12 a chunk is
    written as given and indexed as if it were sorted.  While a single epoch
    holds every record the merge reads every row regardless, so an unsorted
    chunk goes unnoticed; past that, rows are read into the wrong epoch, and the
    output has out-of-order and *duplicate* pixels.  That is how
    ``create pairs --sort-strategy external`` failed ``validate()`` (``bin1_offset
    disagrees with pixels/bin1_id``) on 30 M pairs over 200 k contigs, after the
    whole run.  ``create pairs`` (aggregated with ``sort=False``), ``create
    cool`` and ``reassign`` (re-mapped in source order) all stream such chunks.

    A chunk that is already sorted, as ``binnify``'s, ``aggregate``'s and
    ``extract``'s are, costs one vectorised comparison and is passed through
    as is.

    Args:
        pixels: Pixel chunks with ``bin1_id`` and ``bin2_id`` columns.
        n_bins: Number of bins, so that ``bin1_id * n_bins + bin2_id`` orders
            pixels as ``(bin1_id, bin2_id)`` does.

    Yields:
        Each chunk, sorted.
    """
    for chunk in pixels:
        key = chunk["bin1_id"].to_numpy(np.int64) * n_bins + chunk["bin2_id"].to_numpy(
            np.int64
        )
        if not (key[1:] >= key[:-1]).all():
            chunk = chunk.iloc[np.argsort(key)]
        yield chunk


def _finalize_metacool_file(
    path: Path,
    genomes: pd.DataFrame,
    contigs: pd.DataFrame,
    chrom_offset: np.ndarray | None,
    genome_offset: np.ndarray,
    *,
    contig_ordering: str = "genome-then-alpha",
    assembly_name: str | None = None,
    binsize: int | None = None,
) -> None:
    """Add metacooler-specific groups and attributes to a cooler file.

    Opens *path* in HDF5 ``"a"`` (append) mode and performs the following
    in-place modifications:

    - Replaces ``/chroms`` with a soft-link to ``/contigs`` for cooler
      compatibility.
    - Writes the ``/genomes`` and ``/contigs`` HDF5 groups (with
      fixed-length byte-string encoding).
    - Strips the enum dtype from ``/bins/chrom`` and rewrites it as plain
      ``int32``.
    - Writes ``indexes/chrom_offset`` and ``indexes/genome_offset``.
    - Sets all format attributes (``format``, ``metacooler-format``,
      ``ngenomes``, ``contig-ordering``, etc.).

    Args:
        path: Path to the cooler file to finalize (opened in append mode).
        genomes: ``/genomes`` DataFrame with columns ``name``, ``length``,
            ``n_contigs``, and ``n_bins``.
        contigs: ``/contigs`` DataFrame with columns ``name``, ``length``,
            and ``genome_id``.
        chrom_offset: Pre-computed ``chrom_offset`` array of length
            ``n_contigs + 1``.  If ``None``, it is computed on the fly
            from the ``bins/chrom`` dataset.
        genome_offset: Pre-computed ``genome_offset`` array of length
            ``n_genomes + 1``.
        contig_ordering: Value written to the ``contig-ordering`` HDF5
            attribute: one of the ``--contig-order`` strategies
            (``"genome-then-length"``, ``"genome-then-alpha"``, ``"fasta"``,
            ``"alpha"``) or ``"cool"``.
        assembly_name: Optional assembly name.  Written to cooler's
            standard ``genome-assembly`` attribute (the only attribute
            metacooler uses for this), and only when ``cooler.create``
            has not already stored a real value there.  Ignored if
            ``None``.
        binsize: The bin size that was *requested*, when the file was tiled at
            a fixed resolution.  Written explicitly because cooler infers
            ``bin-size`` from the bin table and returns ``None`` whenever every
            contig holds a single bin — which would leave a genuinely
            fixed-resolution file recorded as ``bin-type="variable"``,
            ``bin-size=null``, and let a later ``binnify`` to a *finer* size
            slip past the "target must be coarser" guard.  Ignored if ``None``.

    Raises:
        MetaCoolerFormatError: If *genomes* or *contigs* hold non-ASCII names
            (see :func:`check_ascii_names`).  Checked before the first write,
            so a refused file is never left half-finalized.
    """
    check_ascii_names(genomes, contigs)
    with h5py.File(path, "a") as f:
        # Recreate chroms as a soft link to contigs (cooler compat)
        if "chroms" in f:
            del f["chroms"]
        f["chroms"] = h5py.SoftLink("contigs")

        # Write /genomes
        if "genomes" in f:
            del f["genomes"]
        g = f.create_group("genomes", track_order=True)
        for col in genomes.columns:
            write_column(g, col, genomes[col].values)

        # Write /contigs
        if "contigs" in f:
            del f["contigs"]
        c = f.create_group("contigs", track_order=True)
        for col in contigs.columns:
            write_column(c, col, contigs[col].values)

        # Strip cooler's enum dtype from bins/chrom so it stores plain int32
        chrom_data = f["bins"]["chrom"][:]
        del f["bins"]["chrom"]
        write_column(f["bins"], "chrom", chrom_data, dtype=np.int32)

        # Compute chrom_offset from bins if not provided
        if chrom_offset is None:
            n_contigs = len(contigs)
            contig_counts = np.bincount(
                chrom_data.astype(np.int64), minlength=n_contigs
            )
            chrom_offset = np.zeros(n_contigs + 1, dtype=np.int64)
            chrom_offset[1:] = np.cumsum(contig_counts)

        # Write index offsets
        idx = f.require_group("indexes")
        for key, data in (
            ("chrom_offset", chrom_offset),
            ("genome_offset", genome_offset),
        ):
            if key in idx:
                del idx[key]
            write_column(idx, key, data, dtype=np.int64)

        # Format attributes
        f.attrs["format"] = COOLER_FORMAT
        f.attrs["format-url"] = COOLER_FORMAT_URL
        f.attrs["metacooler-format"] = METACOOLER_FORMAT
        f.attrs["metacooler-format-url"] = METACOOLER_FORMAT_URL
        f.attrs["metacooler-format-version"] = METACOOLER_FORMAT_VERSION
        f.attrs["generated-by"] = f"metacooler-{__version__}"
        f.attrs["ngenomes"] = len(genomes)
        f.attrs["contig-ordering"] = contig_ordering
        if binsize is not None:
            # Assert the requested resolution rather than trusting cooler's
            # inference; see the *binsize* docstring entry.
            f.attrs["bin-size"] = int(binsize)
            f.attrs["bin-type"] = "fixed"
        # One attribute for the assembly name, and it is cooler's:
        # "genome-assembly".  cooler.create already writes it (defaulting to
        # "unknown"), so only fill it in when it carries no real value —
        # writing a second "assembly-name" attribute here is what made
        # `metacooler info` disagree with `cooler info` on the same file.
        if assembly_name is not None:
            existing = f.attrs.get("genome-assembly")
            if isinstance(existing, bytes):
                existing = existing.decode("utf-8")
            if existing in (None, "", "unknown"):
                f.attrs["genome-assembly"] = assembly_name


def store_qc(
    path: str | Path,
    summary: dict,
    per_genome: pd.DataFrame | None = None,
    *,
    group_path: str = "/",
) -> None:
    """Write a QC block into an existing metacool file.

    The single writer for both producers — ``create`` (right after
    :func:`_finalize_metacool_file`) and ``metacooler qc --store`` on a file
    that already exists — so the two cannot disagree about where the numbers
    live or what they are called.

    Two things are written:

    - ``metadata["qc"]``, merged into whatever the ``metadata`` attribute
      already holds, so a user's own keys survive;
    - one ``/genomes`` column per :data:`metacooler.qc.QC_GENOME_COLUMNS`,
      replacing any column of that name left by an earlier run.

    Args:
        path: Metacool file to annotate, in place.
        summary: The block from :meth:`metacooler.qc.ContactQC.summary`.
        per_genome: Per-genome table from
            :meth:`metacooler.qc.ContactQC.genome_table`, or ``None`` to write
            the file-level block only.
        group_path: HDF5 group of the map, for an mmetacool entry.

    Raises:
        ValueError: If *per_genome* does not have one row per ``/genomes`` row.
    """
    import json

    from metacooler.qc import QC_GENOME_COLUMNS, QC_METADATA_KEY

    with h5py.File(path, "r+") as f:
        grp = f[group_path.lstrip("/")] if group_path.strip("/") else f

        raw = grp.attrs.get("metadata", "{}")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        metadata = json.loads(raw) if isinstance(raw, str) else dict(raw)
        metadata[QC_METADATA_KEY] = summary
        grp.attrs["metadata"] = json.dumps(metadata)

        if per_genome is None:
            return
        genomes = grp["genomes"]
        n_genomes = len(genomes["name"])
        if len(per_genome) != n_genomes:
            raise ValueError(
                f"QC table has {len(per_genome)} row(s) but /genomes has "
                f"{n_genomes}; refusing to write columns that would not line up."
            )
        for column in QC_GENOME_COLUMNS:
            if column in genomes:
                del genomes[column]
            write_column(genomes, column, per_genome[column].to_numpy())


def _filter_by_genomes_and_contigs(
    mc: MetaCooler,
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
) -> set[int]:
    """Return the set of contig IDs matching the genome/contig filters.

    Args:
        mc: Open ``MetaCooler`` instance.
        genomes: Genome names to include (all their contigs are selected).
        contigs: Contig names to include directly.

    Returns:
        Set of integer contig IDs matching the combined filter.

    Raises:
        KeyError: If any genome or contig name is not found in *mc*.
    """
    all_genomes = mc.genomes
    all_contigs = mc.contigs
    selected: set[int] = set()
    known_genomes = set(all_genomes["name"])
    known_contigs = set(all_contigs["name"])

    if genomes:
        for g in genomes:
            if g not in known_genomes:
                available = all_genomes["name"].tolist()
                raise KeyError(
                    f"Genome {g!r} not found. Available ({len(available)} total): {available[:20]}..."
                )
            gid = int(all_genomes.loc[all_genomes["name"] == g].index[0])
            cids = all_contigs.index[all_contigs["genome_id"] == gid].tolist()
            selected.update(cids)

    if contigs:
        for c in contigs:
            if c not in known_contigs:
                available = all_contigs["name"].tolist()
                raise KeyError(
                    f"Contig {c!r} not found. Available ({len(available)} total): {available[:20]}..."
                )
            cid = int(all_contigs.loc[all_contigs["name"] == c].index[0])
            selected.add(cid)

    return selected


def contigs_with_genome_names(mc: MetaCooler) -> pd.DataFrame:
    """The ``/contigs`` table with a ``genome`` column holding the genome *name*.

    ``/contigs`` stores only the integer ``genome_id``, but a filter is written
    against names (``contigs.genome == MAG_01``), and ``create`` evaluates the
    same expressions before any id exists.  Deriving the column here is what
    makes the filter language identical on both sides.

    A contig table that already carries its own ``genome`` column keeps it:
    ``genome`` is not one of ``_CONTIGS_RESERVED_COLUMNS``, so a user column of
    that name survives creation verbatim, and overwriting it here would filter
    on something other than what ``dump`` prints.

    Args:
        mc: Open ``MetaCooler`` instance.

    Returns:
        The contigs table, with a ``genome`` column.
    """
    contigs = mc.contigs
    if "genome" in contigs.columns:
        return contigs
    return contigs.assign(genome=contigs["genome_id"].map(mc.genomes["name"]))


def resolve_contig_selection(
    mc: MetaCooler,
    *,
    genomes: Sequence[str] | None = None,
    contigs: Sequence[str] | None = None,
    filters: Sequence[str] | None = None,
) -> set[int] | None:
    """Resolve name selectors plus ``--filter`` predicates to a set of contig IDs.

    Names choose the candidate set; filters narrow it.  ``--genomes MAG_01
    --filter 'contigs.length > 50kb'`` is therefore *MAG_01's long contigs*,
    not "MAG_01, plus every long contig in the file".  With no names given the
    candidate set is the whole file, so ``--filter`` alone selects from
    everything.

    Args:
        mc: Open ``MetaCooler`` instance.
        genomes: Genome names to include (all their contigs).
        contigs: Contig names to include directly.
        filters: Raw ``--filter`` expressions (see :mod:`metacooler.filters`).

    Returns:
        The selected contig IDs, or ``None`` when no selector was given at all
        — which callers read as "no filtering, use the whole file".

    Raises:
        KeyError: If a genome or contig name is not found.
        metacooler.filters.FilterError: If a filter is malformed or names a
            column the tables do not have, or if nothing survives it.
    """
    from metacooler.filters import (
        FilterError,
        explain_empty,
        filter_mask,
        parse_filters,
    )

    parsed = parse_filters(filters)
    named = bool(genomes) or bool(contigs)
    if not named and not parsed:
        return None

    all_contigs = contigs_with_genome_names(mc)
    if named:
        selected = _filter_by_genomes_and_contigs(mc, genomes=genomes, contigs=contigs)
    else:
        selected = {int(i) for i in all_contigs.index}
    if not parsed:
        return selected

    all_genomes = mc.genomes
    keep_genome_ids = set(
        all_genomes.index[filter_mask(all_genomes, parsed, table="genomes")]
    )
    keep = set(
        all_contigs.index[
            filter_mask(all_contigs, parsed, table="contigs")
            & all_contigs["genome_id"].isin(keep_genome_ids).to_numpy()
        ]
    )
    before = len(selected)
    selected &= {int(i) for i in keep}
    logger.info(
        "--filter kept %d of %d contig(s) (%d of %d genome(s) passed the "
        "genome-level filters).",
        len(selected),
        before,
        len(keep_genome_ids),
        len(all_genomes),
    )
    if not selected:
        context = ""
        if named:
            context = (
                f"The {before} contig(s) named by --genomes/--contigs were then "
                f"intersected with the above, leaving none."
            )
        raise FilterError(
            explain_empty(
                {"genomes": all_genomes, "contigs": all_contigs},
                parsed,
                context=context,
            )
        )
    return selected
