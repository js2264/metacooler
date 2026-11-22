"""Normalization weights: per-bin columns of ``/bins``, from any tool.

A normalization that *factorises* -- the normalized value of a pixel is its
raw count times one factor for each of its two bins, ``raw[i, j] * w[i] *
w[j]`` -- is stored the way cooler stores a balancing vector: one ``float64``
column of ``/bins``, ``NaN`` for a masked bin, multiplicative.  A map keeps
working in every tool that reads cooler weights
(``cooler.Cooler.matrix(balance=...)``, ``cooler dump --balanced``,
cooltools, HiGlass), several normalizations sit side by side under different
names, and the pixels are never rewritten.

``metacooler normalize ice`` stores its weights here, and so can any other
tool: :func:`store_weights` writes a weight column, :func:`read_weights`
reads one back -- what ``show``, ``dump`` and ``network --balance`` do.

A weight column carries, as HDF5 attributes, the method that produced it
(:data:`METHOD_ATTR`), its parameters, cooler's ``divisive_weights=False``,
the metacooler version that stored it, and what it was computed from
(:data:`DEPENDS_ATTR`): the pixels, the genome assignment, named ``/contigs``
columns.  A weight is *derived data*, like the QC block: the writers that
change what it depends on drop it (:func:`drop_stale_weights`) rather than
leave a map whose normalized values no longer mean anything.

This module is the single place that reads and writes those columns.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np

from metacooler._version import __version__
from metacooler.io._utils import write_column

if TYPE_CHECKING:
    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

#: Columns of ``/bins`` that hold the bin coordinates, never a weight.
BIN_COORDINATE_COLUMNS = ("chrom", "start", "end")

#: Column names cooler reads as *divisive* weights (``normalized = raw / (w_i
#: w_j)``): the vectors the 4DN portal and hic2cool store.  metacooler writes
#: multiplicative weights only, so it refuses these names rather than write a
#: column cooler would silently invert.
DIVISIVE_WEIGHT_NAMES = frozenset({"KR", "VC", "VC_SQRT"})

#: Attribute naming the method that produced a column: ``ice``, or the name
#: the tool that stored it gave.  A weight column without it was written by
#: another route, typically ``cooler balance``.
METHOD_ATTR = "normalization"

#: Attribute listing what a weight was computed from, so that a later edit of
#: the map can tell whether it still holds (see :func:`stale_weights`).
DEPENDS_ATTR = "depends_on"

#: Attribute recording the metacooler version that stored a column.
VERSION_ATTR = "metacooler_version"

#: :data:`DEPENDS_ATTR` entry for the pixel table.
DEPENDS_ON_PIXELS = "pixels"

#: :data:`DEPENDS_ATTR` entry for the contig-to-genome assignment
#: (``/contigs/genome_id``): a weight computed per genome, or fitted on the
#: pairs of contigs one genome holds.
DEPENDS_ON_ASSIGNMENT = "genome_id"

#: Prefix of a :data:`DEPENDS_ATTR` entry naming a ``/contigs`` column a weight
#: read, e.g. ``contigs/Shotgun_coverage``.
CONTIG_COLUMN_PREFIX = "contigs/"

#: Attributes recording a spurious-contact cutoff (NormCC, HiCzin): contig
#: pairs whose normalized contact is below it -- or at it, under the ``<=``
#: rule -- are deemed spurious.  Any tool may store them, as parameters of its
#: weight column, for ``network --drop-spurious`` to apply.
SPURIOUS_CUTOFF_ATTR = "spurious_cutoff"
SPURIOUS_RULE_ATTR = "spurious_rule"

#: Attributes :func:`store_weights` writes itself, which a method's parameters
#: cannot take over.
RESERVED_WEIGHT_ATTRS = frozenset(
    {METHOD_ATTR, DEPENDS_ATTR, VERSION_ATTR, "divisive_weights"}
)


def check_weight_name(name: str) -> None:
    """Refuse a name that cannot hold a metacooler weight column.

    Raises:
        ValueError: If *name* is empty, contains ``/`` (HDF5 would make a
            nested group, which cooler cannot read as a column), is one of the
            bin coordinate columns, or is a name cooler reads as divisive.
    """
    if not name or not name.strip():
        raise ValueError("A weight column needs a non-empty name.")
    if "/" in name:
        raise ValueError(
            f"Weight column name {name!r} contains '/', which HDF5 would turn "
            "into a nested group that cooler cannot read as a column."
        )
    if name in BIN_COORDINATE_COLUMNS:
        raise ValueError(
            f"{name!r} is a bin coordinate column of /bins, not a weight name."
        )
    if name in DIVISIVE_WEIGHT_NAMES:
        raise ValueError(
            f"cooler reads a weight column named {name!r} as divisive "
            f"(raw / (w_i * w_j)), but metacooler writes multiplicative weights; "
            "pick another name."
        )


def check_depends_on(
    depends_on: str | Iterable[str], contig_columns: Iterable[str]
) -> list[str]:
    """Validate the :data:`DEPENDS_ATTR` entries of a weight about to be stored.

    An entry is ``pixels``, ``genome_id`` (the contig-to-genome assignment) or
    ``contigs/<column>``, naming a column stored in the map's ``/contigs``.  A
    typo would leave a weight no edit can invalidate, so anything else is
    refused, as is a column the map does not hold -- a user's column named
    ``genome``, say, is stored as ``user_genome``.  ``pixels`` is not
    required: a weight computed from a coverage column alone depends on that
    column only.

    Args:
        depends_on: The entries, in any order (repeats are dropped), or a
            single entry.
        contig_columns: The columns of the map's ``/contigs``.

    Returns:
        The entries, without repeats, in the order given.

    Raises:
        ValueError: If there is no entry, or an entry is none of the above.
    """
    if isinstance(depends_on, str):
        depends_on = [depends_on]
    entries = list(dict.fromkeys(str(entry) for entry in depends_on))
    if not entries:
        raise ValueError(
            "depends_on must name what the weights were computed from: "
            "'pixels', 'genome_id' or 'contigs/<column>'."
        )
    columns = [c for c in contig_columns if c != "genome_id"]
    bad = []
    for entry in entries:
        if entry in (DEPENDS_ON_PIXELS, DEPENDS_ON_ASSIGNMENT):
            continue
        column = entry.removeprefix(CONTIG_COLUMN_PREFIX)
        if column == "genome_id":
            raise ValueError(
                "Name the genome assignment as 'genome_id' in depends_on, not "
                "'contigs/genome_id': that is the entry a reassignment drops a "
                "weight for."
            )
        if column == entry or column not in columns:
            bad.append(entry)
    if bad:
        raise ValueError(
            f"depends_on entries must be 'pixels', 'genome_id' or "
            f"'contigs/<column>' naming a stored /contigs column; got "
            f"{', '.join(repr(b) for b in bad)}. Stored /contigs columns: "
            f"{', '.join(columns) or 'none'}."
        )
    return entries


def weight_attributes(
    method: str,
    depends_on: Sequence[str],
    parameters: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """The attributes :func:`store_weights` writes with a weight column.

    ``None`` parameters are left out.
    """
    attrs: dict[str, object] = {
        METHOD_ATTR: method,
        DEPENDS_ATTR: list(depends_on),
        "divisive_weights": False,
    }
    attrs.update({k: v for k, v in (parameters or {}).items() if v is not None})
    attrs[VERSION_ATTR] = __version__
    return attrs


def weight_columns(bins: h5py.Group) -> list[str]:
    """Return the names of the weight-like columns of a ``/bins`` group.

    A weight-like column is any floating-point dataset beside the coordinate
    columns, whoever wrote it.
    """
    names = []
    for key in bins:
        if key in BIN_COORDINATE_COLUMNS:
            continue
        node = bins.get(key)
        if isinstance(node, h5py.Dataset) and node.dtype.kind == "f":
            names.append(key)
    return names


def _attr_value(value: object) -> object:
    """Convert *value* into something h5py stores as an attribute."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return json.dumps(value, default=str)
    if isinstance(value, (list, tuple)):
        if all(isinstance(v, str) for v in value):
            return np.array(value, dtype=h5py.string_dtype())
        return np.asarray(value)
    return value


def _decoded(value: object) -> object:
    """Return an attribute value with HDF5 byte strings decoded."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.dtype.kind in "OS":
        return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in value]
    return value


def _open_map(source: MetaCooler | str | Path) -> MetaCooler:
    """*source* as an open map: a ``MetaCooler`` as it is, a URI opened."""
    from metacooler.api import MetaCooler

    return source if isinstance(source, MetaCooler) else MetaCooler(source)


def depends_on(dataset: h5py.Dataset) -> set[str]:
    """Return the :data:`DEPENDS_ATTR` entries of a weight column.

    A column without the attribute -- one ``cooler balance`` wrote -- is
    treated as depending on the pixels alone, what ``cooler balance``
    computes from.
    """
    raw = _decoded(dataset.attrs.get(DEPENDS_ATTR))
    if raw is None:
        return {DEPENDS_ON_PIXELS}
    if isinstance(raw, str):
        return {raw}
    return {str(v) for v in np.atleast_1d(np.asarray(raw, dtype=object))}


def stale_weights(
    bins: h5py.Group,
    *,
    assignment: bool = False,
    contig_columns: Iterable[str] = (),
) -> list[str]:
    """Return the weight columns an edit of the file invalidates.

    Args:
        bins: The ``/bins`` group of the map being edited.
        assignment: The edit changes the contig-to-genome assignment.
        contig_columns: ``/contigs`` columns the edit adds, replaces or drops.

    Returns:
        Names of the weight columns that depend on what the edit changes, in
        ``/bins`` order.
    """
    columns = {f"{CONTIG_COLUMN_PREFIX}{c}" for c in contig_columns}
    stale = []
    for name in weight_columns(bins):
        needs = depends_on(bins[name])
        if (assignment and DEPENDS_ON_ASSIGNMENT in needs) or (needs & columns):
            stale.append(name)
    return stale


def drop_stale_weights(
    bins: h5py.Group,
    *,
    assignment: bool = False,
    contig_columns: Iterable[str] = (),
    why: str,
) -> list[str]:
    """Delete the weight columns an edit invalidates, and say so.

    Args:
        bins: The ``/bins`` group of the map being edited, open for writing.
        assignment: The edit changes the contig-to-genome assignment.
        contig_columns: ``/contigs`` columns the edit adds, replaces or drops.
        why: What the edit changed, completing "Dropped ... computed from".

    Returns:
        The names of the columns deleted.
    """
    stale = stale_weights(bins, assignment=assignment, contig_columns=contig_columns)
    for name in stale:
        del bins[name]
    warn_dropped_weights(stale, why)
    return stale


def warn_dropped_weights(names: Sequence[str], why: str) -> None:
    """Say that the weight columns *names* were not carried over, and why."""
    if names:
        logger.warning(
            "Dropped normalization weight(s) %s: they were computed from %s, "
            "which this edit changed. Recompute them: 'metacooler normalize "
            "ice' for an ICE weight, or the tool that stored them.",
            ", ".join(repr(n) for n in names),
            why,
        )


def weight_attrs(bins: h5py.Group, names: Iterable[str]) -> dict[str, dict]:
    """Return the HDF5 attributes of the named ``/bins`` columns.

    Used by writers that rebuild ``/bins`` through cooler -- which writes the
    values of an extra column but not its attributes -- to restore them.
    """
    return {name: dict(bins[name].attrs) for name in names if name in bins}


def restore_weight_attrs(bins: h5py.Group, attrs: Mapping[str, Mapping]) -> None:
    """Write back attributes captured with :func:`weight_attrs`."""
    for name, values in attrs.items():
        if name not in bins:
            continue
        for key, value in values.items():
            bins[name].attrs[key] = value


def store_weights(
    uri: str | Path,
    name: str,
    weights: np.ndarray,
    *,
    method: str,
    depends_on: str | Sequence[str],
    parameters: Mapping[str, object] | None = None,
    force: bool = False,
    commands: Sequence[str] = (),
    source: MetaCooler | None = None,
) -> None:
    """Store a normalization as a weight column of a map's ``/bins``, in place.

    Through a staged copy of the map, under its lock, which replaces the map
    once it validates (see :mod:`metacooler.io.staging`).  Weights computed
    from the map belong to the map as it was read: hold its lock
    (:func:`metacooler.lock`) from before that read, and pass the map read as
    *source*, so that a map changed meanwhile is refused rather than given
    weights that do not line up with its bins.

    The way any tool -- a binner, a normalization of its own -- puts a
    normalization in a map: ``show``, ``dump`` and ``network --balance NAME``
    then read it, as cooler, cooltools and HiGlass read ``/bins/NAME``.  The
    normalized value of a pixel is ``raw * weights[bin1] * weights[bin2]``.

    *depends_on* is what lets the map keep the column honest: ``reassign``
    drops a weight that lists ``genome_id`` when it moves a contig, and
    ``annotate`` one that lists a ``/contigs`` column it changes.  A weight
    computed from a coverage column, say, must name that column, or it would
    outlive the values it was computed from.

    Args:
        uri: ``file.metacool`` or ``file.mmetacool::entries/<name>``.
        name: The ``/bins`` column (see :func:`check_weight_name`).
        weights: One multiplicative weight per bin, ``NaN`` for a masked bin.
        method: The normalization's name, stored as the ``normalization``
            attribute.
        depends_on: What the weights were computed from: ``pixels``,
            ``genome_id`` (the contig-to-genome assignment) and
            ``contigs/<column>`` for each ``/contigs`` column read (see
            :func:`check_depends_on`).
        parameters: The method's own values -- a mode, the tool's version --
            stored as attributes of the column; ``None`` values are left out.
        force: Replace a column of that name if there is one.
        commands: Command lines to append to ``metadata["commands"]``.
        source: The map the weights were computed from, as opened for it:
            refused (``MapChangedError``) if its file changed since -- or, in
            an edit session, if a write of the session changed what the
            weight depends on (see
            :func:`metacooler.io.staging.check_source_depends`); one opened
            on the session's staged map is judged by that alone, so several
            weights may be computed from one read of it.

    Raises:
        ValueError: If *name* cannot hold a weight column, or is taken and
            *force* is not set; if *method* is empty; if *depends_on* is not
            valid; if a parameter takes the name of an attribute this function
            writes; or if *weights* does not hold one value per bin.

    Example:
        A tool normalizing each contig by its Hi-C hits, as MetaTOR's
        ``empirical_hit`` does, stores ``1 / sqrt(hits)`` for every bin of
        the contig::

            metacooler.store_weights(
                "map.metacool",
                "empirical_hit",
                weights,
                method="metator",
                depends_on=["pixels"],
                parameters={"mode": "empirical_hit"},
            )

        and ``metacooler network map.metacool --balance empirical_hit``
        gives its normalized contig network.
    """
    check_weight_name(name)
    if not isinstance(method, str) or not method.strip():
        raise ValueError("method must name the normalization the weights are.")
    params = dict(parameters or {})
    taken = sorted(RESERVED_WEIGHT_ATTRS & params.keys())
    if taken:
        raise ValueError(
            f"Parameter name(s) {', '.join(repr(k) for k in taken)} are attributes "
            "store_weights writes itself; rename them."
        )
    values = np.asarray(weights, dtype=np.float64)
    from metacooler.api import MetaCooler
    from metacooler.io.staging import check_source_depends, lock, staged_inplace
    from metacooler.uri import parse_uri

    file_path, group = parse_uri(str(uri))
    with lock(file_path):
        mc = MetaCooler(uri)
        # Checked on the map before it is copied.
        with mc.open() as grp:
            entries = check_depends_on(depends_on, grp["contigs"].keys())
            n_bins = int(grp["bins"]["start"].shape[0])
            taken = name in grp["bins"]
        if values.shape != (n_bins,):
            raise ValueError(
                f"Got {values.size} weight(s) for a map of {n_bins} bin(s); "
                "refusing to write a column that would not line up."
            )
        if taken and not force:
            raise ValueError(
                f"/bins/{name} already exists in {mc.uri}: store these "
                "weights under another name, or replace that column with "
                "force=True."
            )
        # In an edit session: a source that holds what a write of the session
        # changed of what the weight depends on is not current.  One on the
        # staged map is judged by the session's log of its writes alone, as
        # the writes of the session change its file.
        judged = check_source_depends(uri, source, entries)
        planned = [mc] if source is None or judged else [mc, source]
        with (
            staged_inplace(uri, planned=planned, touches=()) as staged,
            h5py.File(staged, "r+") as f,
        ):
            key = group.strip("/")
            _store_column(
                f[key] if key else f, name, values, method, entries, params, commands
            )


def _store_column(
    grp: h5py.Group,
    name: str,
    values: np.ndarray,
    method: str,
    entries: list[str],
    params: dict,
    commands: Sequence[str],
) -> None:
    """Write the weight column *name* of a map's group, and the commands.

    The new column is written whole, under a scratch name, before it takes
    the place of a column of that name: a write that fails leaves the column
    it would replace as it was, and no scratch column.
    """
    bins = grp["bins"]
    scratch = f".{name}.new"
    if scratch in bins:
        del bins[scratch]
    try:
        dataset = write_column(bins, scratch, values)
        for key, value in weight_attributes(method, entries, params).items():
            dataset.attrs[key] = _attr_value(value)
    except BaseException:
        if scratch in bins:
            del bins[scratch]
        raise
    if name in bins:
        del bins[name]
    bins.move(scratch, name)
    if commands:
        raw = grp.attrs.get("metadata", "{}")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        metadata = json.loads(raw) if isinstance(raw, str) else dict(raw)
        metadata["commands"] = [*metadata.get("commands", []), *commands]
        grp.attrs["metadata"] = json.dumps(metadata)


def read_weights(source: MetaCooler | str | Path, name: str) -> np.ndarray:
    """Return a map's weight column as multiplicative weights, one per bin.

    A column cooler reads as divisive (:data:`DIVISIVE_WEIGHT_NAMES`, or one
    whose ``divisive_weights`` attribute is true) is inverted, so the result
    always means ``normalized = raw * w[bin1] * w[bin2]``.  Masked bins are
    ``NaN``.

    Args:
        source: ``file.metacool``, ``file.mmetacool::entries/<name>``, or an
            open ``MetaCooler``.
        name: The ``/bins`` column.

    Raises:
        KeyError: If the map has no such weight column.
        ValueError: If the column does not hold one number per bin.
    """
    mc = _open_map(source)
    with mc.open() as grp:
        bins = grp["bins"]
        node = bins.get(name) if name not in BIN_COORDINATE_COLUMNS else None
        if not isinstance(node, h5py.Dataset):
            available = weight_columns(bins)
            listing = ", ".join(available) if available else "none"
            raise KeyError(
                f"No weight column {name!r} in /bins of {mc.uri} (weight columns: "
                f"{listing}). Compute one with 'metacooler normalize ice', or "
                "store another tool's with metacooler.store_weights."
            )
        n_bins = int(bins["start"].shape[0])
        if node.dtype.kind not in "fiu" or node.shape != (n_bins,):
            raise ValueError(
                f"/bins/{name} does not hold one number per bin "
                f"(dtype {node.dtype}, shape {node.shape})."
            )
        values = node[()].astype(np.float64)
        divisive = name in DIVISIVE_WEIGHT_NAMES or bool(
            node.attrs.get("divisive_weights", False)
        )
    if divisive:
        with np.errstate(divide="ignore"):
            values = 1.0 / values
    values[~np.isfinite(values)] = np.nan
    return values


def spurious_cutoff(
    source: MetaCooler | str | Path, name: str
) -> tuple[float, str] | None:
    """Return the spurious-contact cutoff stored with a weight column, if any.

    Args:
        source: ``file.metacool``, ``file.mmetacool::entries/<name>``, or an
            open ``MetaCooler``.
        name: The ``/bins`` column.

    Returns:
        ``(cutoff, rule)`` -- a contig pair whose normalized contact is below
        the cutoff (rule ``"<"``) or at most the cutoff (rule ``"<="``) is
        spurious -- or ``None`` if the column records no cutoff.
    """
    with _open_map(source).open() as grp:
        node = grp["bins"].get(name)
        if not isinstance(node, h5py.Dataset):
            return None
        cutoff = node.attrs.get(SPURIOUS_CUTOFF_ATTR)
        rule = _decoded(node.attrs.get(SPURIOUS_RULE_ATTR, "<"))
    if cutoff is None:
        return None
    return float(cutoff), str(rule)


def describe_weights(source: MetaCooler | str | Path) -> list[dict[str, object]]:
    """Describe every weight column of a map: name, method, and attributes.

    Args:
        source: ``file.metacool``, ``file.mmetacool::entries/<name>``, or an
            open ``MetaCooler``.

    Returns:
        One dict per weight column, in ``/bins`` order, with ``name``,
        ``method`` (``None`` for a column without a ``normalization``
        attribute), ``n_masked`` (bins whose weight is ``NaN``) and ``attrs``
        (every attribute, decoded).
    """
    rows: list[dict[str, object]] = []
    with _open_map(source).open() as grp:
        bins = grp["bins"]
        for name in weight_columns(bins):
            node = bins[name]
            attrs = {key: _decoded(value) for key, value in node.attrs.items()}
            values = node[()]
            rows.append(
                {
                    "name": name,
                    "method": attrs.get(METHOD_ATTR),
                    "n_masked": int(np.count_nonzero(~np.isfinite(values))),
                    "attrs": attrs,
                }
            )
    return rows
