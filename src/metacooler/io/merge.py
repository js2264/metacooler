"""Merge operations for metacool files."""

from __future__ import annotations

import contextlib
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING

import cooler
import h5py
import numpy as np

from metacooler._version import __version__
from metacooler.io._utils import check_output_path, cooler_writer_shims, write_column
from metacooler.qc import QC_COLUMN_PREFIX, QC_METADATA_KEY
from metacooler.schema import (
    COOLER_FORMAT,
    COOLER_FORMAT_URL,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_URL,
    METACOOLER_FORMAT_VERSION,
    validate,
)
from metacooler.uri import parse_uri

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable


def _hash_group(f: h5py.Group, group: str) -> str:
    """Return an MD5 hex digest of all datasets in a group.

    Used to compare ``contigs`` and ``genomes`` tables across files for
    merge compatibility checking.

    Args:
        f: Open ``h5py.File`` (or any HDF5 group-like object).
        group: HDF5 group path to hash (e.g. ``"contigs"`` or ``"genomes"``).

    Returns:
        Lowercase hex MD5 digest string of all dataset values, sorted by
        key name for determinism.
    """
    h = hashlib.md5(usedforsecurity=False)
    g = f[group]
    for key in sorted(g.keys()):
        # QC columns describe each input's own contacts, not the assembly they
        # share, so two libraries of the same metagenome legitimately differ
        # there. Hashing them would reject exactly the merge this command
        # exists for — replicates over one assembly.
        if str(key).startswith(QC_COLUMN_PREFIX):
            continue
        arr = g[key][:]
        # Variable-length string datasets read back as object arrays of Python
        # strings; encode them to bytes for a stable hash.
        if arr.dtype == object or str(arr.dtype).startswith("S"):
            arr = np.array(
                [v.encode() if isinstance(v, str) else bytes(v) for v in arr.ravel()]
            )
        h.update(arr.tobytes())
    return h.hexdigest()


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


def merge(
    inputs: Iterable[str | Path],
    output: str | Path,
) -> None:
    """Merge multiple ``.metacool`` files with identical axes into one.

    All input files must have identical ``contigs`` and ``genomes`` tables
    and the same ``contig-ordering`` attribute.  Pixels from all inputs are
    concatenated, re-sorted externally by ``cooler.merge_coolers()``, and
    saved to *output*.  Metacooler-specific groups (``genomes``, ``contigs``,
    ``genome_offset``, format attributes) are then copied from the first
    input file.

    Args:
        inputs: Iterable of paths to the input ``.metacool`` files, or
            ``file.mmetacool::entries/<name>`` URIs.  At least two are
            required.
        output: Destination path for the merged ``.metacool`` file.

    Raises:
        ValueError: If fewer than two inputs are provided, if the
            ``contigs``, ``genomes`` or ``bins`` tables differ across files,
            if the resolutions differ, or if the ``contig-ordering``
            attribute is inconsistent.
        MetaCoolerFormatError: If any input file fails the metacooler
            format validation check.
    """
    inputs = [str(p) for p in inputs]
    if len(inputs) < 2:
        raise ValueError("At least two input files required for merge")
    # cooler writes the output in "w" mode; refuse to truncate an input with it.
    check_output_path(output, force=True, inputs=inputs, defer_removal=True)

    # Validate and compare axes
    first = inputs[0]
    validate(first)

    with _open_uri(first) as f0:
        contigs_hash = _hash_group(f0, "contigs")
        genomes_hash = _hash_group(f0, "genomes")
        # ``contig-ordering`` may be absent or explicitly null; normalise both
        # sides of the comparison so a missing attribute does not later become
        # ``f.attrs[...] = None``, which h5py rejects with a raw TypeError.
        ordering = f0.attrs.get("contig-ordering") or "genome-then-alpha"
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
            if (f.attrs.get("contig-ordering") or "genome-then-alpha") != ordering:
                raise ValueError(f"contig-ordering mismatch: {first} vs {inp}")

    # Use cooler's merge logic on the pixel tables.  ``dtypes`` is mandatory:
    # without it the summed counts inherit the int32 of the inputs and wrap.
    # merge_coolers calls cooler.create internally, so it needs the same shims
    # as every other writer (see write_cooler): without them its per-contig
    # costs made merge alone take minutes at 150 k contigs. They also hand
    # cooler the inputs' known bin size, which it cannot infer once every
    # contig is a single bin.
    with cooler_writer_shims(binsize):
        cooler.merge_coolers(
            str(output),
            inputs,
            mergebuf=int(2e7),
            dtypes={"count": _merged_count_dtype(inputs)},
        )

    # Post-process: add back metacooler-specific groups
    with _open_uri(first) as fsrc, h5py.File(output, "a") as fdst:
        # Copy genomes, minus the first input's QC columns: the merged pixel
        # table is the sum of every input, so counts copied from one of them
        # would understate it while looking authoritative. Re-run
        # ``metacooler qc --store`` on the output to get honest ones.
        if "genomes" in fdst:
            del fdst["genomes"]
        fdst.copy(fsrc["genomes"], "genomes")
        for key in [k for k in fdst["genomes"] if str(k).startswith(QC_COLUMN_PREFIX)]:
            del fdst["genomes"][key]

        # Copy contigs (overwrites the chroms-compatible one)
        if "contigs" in fdst:
            del fdst["contigs"]
        fdst.copy(fsrc["contigs"], "contigs")

        # Soft link chroms -> contigs
        if "chroms" in fdst:
            del fdst["chroms"]
        fdst["chroms"] = h5py.SoftLink("contigs")

        # Copy genome_offset
        if "genome_offset" in fdst["indexes"]:
            del fdst["indexes"]["genome_offset"]
        fdst["indexes"].copy(fsrc["indexes"]["genome_offset"], "genome_offset")

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
        # carried describes that input alone.
        raw = fdst.attrs.get("metadata")
        if raw is not None:
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
            if QC_METADATA_KEY in meta:
                del meta[QC_METADATA_KEY]
                fdst.attrs["metadata"] = json.dumps(meta)

    validate(output)
