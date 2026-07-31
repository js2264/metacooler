"""Shared I/O utilities used across multiple metacooler modules."""

from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING

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
    if arr.dtype.kind in ("O", "U", "S") or str(arr.dtype) == "string":
        str_arr = ["" if pd.isna(v) else str(v) for v in arr]
        encoded = [s.encode("utf-8") for s in str_arr]
        # Size the fixed-length dtype by ENCODED BYTE length, not character
        # count: a multi-byte character (e.g. "é") is longer as UTF-8 than as a
        # Python str, so sizing by len(str) would silently truncate the bytes
        # (and could collide two distinct names).
        max_len = max((len(b) for b in encoded), default=1)
        return np.array(encoded, dtype=f"S{max_len}")
    return np.asarray(arr)


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


def check_output_path(path: str | Path, force: bool = False) -> None:
    """Refuse to overwrite an existing output file unless *force* is set.

    Single home for the "output exists" guard shared by the CLI commands and
    the ops layer, so every command reports the same message and every command
    that accepts ``--force`` actually honours it.

    Args:
        path: Output path about to be written.
        force: If True, remove an existing file instead of raising.

    Raises:
        FileExistsError: If *path* exists and *force* is False.
        IsADirectoryError: If *path* exists and is a directory.
    """
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        raise IsADirectoryError(f"Output path {str(path)!r} is a directory.")
    if not force:
        raise FileExistsError(
            f"Output file {str(path)!r} already exists. Use --force to overwrite."
        )
    path.unlink()


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
            attribute.  Typically ``"genome-then-alpha"``,
            ``"pairs-header"``, or ``"genomeid-then-alpha"``.
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
            g.create_dataset(col, data=_str_to_hdf5(genomes[col].values))

        # Write /contigs
        if "contigs" in f:
            del f["contigs"]
        c = f.create_group("contigs", track_order=True)
        for col in contigs.columns:
            c.create_dataset(col, data=_str_to_hdf5(contigs[col].values))

        # Strip cooler's enum dtype from bins/chrom so it stores plain int32
        chrom_data = f["bins"]["chrom"][:]
        del f["bins"]["chrom"]
        f["bins"].create_dataset("chrom", data=chrom_data, dtype=np.int32)

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
            idx.create_dataset(key, data=data.astype(np.int64))

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


def _filter_by_genomes_and_contigs(
    mc: MetaCooler,
    genomes: list[str] | None = None,
    contigs: list[str] | None = None,
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
