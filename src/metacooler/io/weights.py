"""Normalization weights: the per-bin columns ``metacooler normalize`` stores.

Every normalization metacooler offers *factorises*: the normalized value of a
pixel is its raw count times one factor for each of its two bins,
``raw[i, j] * w[i] * w[j]``.  That is exactly how cooler stores a balancing
vector, so a normalization is stored the same way -- one ``float64`` column of
``/bins``, ``NaN`` for a masked bin, multiplicative -- and a file keeps working
in every tool that reads cooler weights (``cooler.Cooler.matrix(balance=...)``,
``cooler dump --balanced``, cooltools, HiGlass).  Several normalizations can
sit side by side under different names.  The pixels are never rewritten.

A weight column written by metacooler carries, as HDF5 attributes, the method
that produced it (:data:`METHOD_ATTR`), its parameters, cooler's
``divisive_weights=False``, and what it was computed from
(:data:`DEPENDS_ATTR`): the pixels, and possibly the genome assignment and
named ``/contigs`` columns.  A weight is *derived data*, like the QC block:
the writers that change what it depends on drop it (:func:`stale_weights`)
rather than leave a file whose normalized values no longer mean anything.

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

#: Attribute naming the metacooler method that wrote a column (``ice``,
#: ``metator``, ``mge``, ``normcc``, ``hiczin``).  A weight column without it
#: was written by another tool, typically ``cooler balance``.
METHOD_ATTR = "normalization"

#: Attribute listing what a weight was computed from, so that a later edit of
#: the file can tell whether it still holds (see :func:`stale_weights`).
DEPENDS_ATTR = "depends_on"

#: :data:`DEPENDS_ATTR` entry for the pixel table.  Every weight has it.
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
#: rule -- are deemed spurious.
SPURIOUS_CUTOFF_ATTR = "spurious_cutoff"
SPURIOUS_RULE_ATTR = "spurious_rule"


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


def depends_on(dataset: h5py.Dataset) -> set[str]:
    """Return the :data:`DEPENDS_ATTR` entries of a weight column.

    A column metacooler did not write has none, and is treated as depending on
    the pixels alone -- what ``cooler balance`` computes from.
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
            "which this edit changed. Re-run 'metacooler normalize' to "
            "recompute them.",
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
    path: str | Path,
    name: str,
    weights: np.ndarray,
    attrs: Mapping[str, object],
    *,
    group_path: str = "/",
    force: bool = False,
    commands: Sequence[str] = (),
) -> None:
    """Write one weight column into ``/bins`` of an existing map, in place.

    Args:
        path: The metacool file.
        name: Column name (see :func:`check_weight_name`).
        weights: One multiplicative weight per bin, ``NaN`` for a masked bin.
        attrs: Attributes for the column; ``None`` values are left out.
        group_path: HDF5 group of the map, for an mmetacool entry.
        force: Replace a column of that name if there is one.
        commands: Command lines to append to ``metadata["commands"]``.

    Raises:
        ValueError: If *name* is not a valid weight name, *weights* does not
            have one value per bin, or the column exists and *force* is not set.
    """
    check_weight_name(name)
    values = np.asarray(weights, dtype=np.float64)
    with h5py.File(path, "r+") as f:
        grp = f[group_path.strip("/")] if group_path.strip("/") else f
        bins = grp["bins"]
        n_bins = int(bins["start"].shape[0])
        if values.shape != (n_bins,):
            raise ValueError(
                f"Got {values.size} weight(s) for a map of {n_bins} bin(s); "
                "refusing to write a column that would not line up."
            )
        if name in bins:
            if not force:
                raise ValueError(
                    f"/bins/{name} already exists in {path}. Pass --force to "
                    "replace it, or --name to store this normalization beside it."
                )
            del bins[name]
        dataset = write_column(bins, name, values)
        dataset.attrs["divisive_weights"] = False
        for key, value in attrs.items():
            if value is not None:
                dataset.attrs[key] = _attr_value(value)
        if commands:
            raw = grp.attrs.get("metadata", "{}")
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            metadata = json.loads(raw) if isinstance(raw, str) else dict(raw)
            metadata["commands"] = [*metadata.get("commands", []), *commands]
            grp.attrs["metadata"] = json.dumps(metadata)


def read_weights(mc: MetaCooler, name: str) -> np.ndarray:
    """Return a map's weight column as multiplicative weights, one per bin.

    A column cooler reads as divisive (:data:`DIVISIVE_WEIGHT_NAMES`, or one
    whose ``divisive_weights`` attribute is true) is inverted, so the result
    always means ``normalized = raw * w[i] * w[j]``.  Masked bins are ``NaN``.

    Args:
        mc: The open map.
        name: The ``/bins`` column.

    Raises:
        KeyError: If the map has no such weight column.
        ValueError: If the column does not hold one number per bin.
    """
    with mc.open() as grp:
        bins = grp["bins"]
        node = bins.get(name) if name not in BIN_COORDINATE_COLUMNS else None
        if not isinstance(node, h5py.Dataset):
            available = weight_columns(bins)
            listing = ", ".join(available) if available else "none"
            raise KeyError(
                f"No weight column {name!r} in /bins of {mc.uri} (weight columns: "
                f"{listing}). Compute one with 'metacooler normalize'."
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


def spurious_cutoff(mc: MetaCooler, name: str) -> tuple[float, str] | None:
    """Return the spurious-contact cutoff stored with a weight column, if any.

    Returns:
        ``(cutoff, rule)`` -- a contig pair whose normalized contact is below
        the cutoff (rule ``"<"``) or at most the cutoff (rule ``"<="``) is
        spurious -- or ``None`` if the method records no cutoff.
    """
    with mc.open() as grp:
        node = grp["bins"].get(name)
        if not isinstance(node, h5py.Dataset):
            return None
        cutoff = node.attrs.get(SPURIOUS_CUTOFF_ATTR)
        rule = _decoded(node.attrs.get(SPURIOUS_RULE_ATTR, "<"))
    if cutoff is None:
        return None
    return float(cutoff), str(rule)


def describe_weights(mc: MetaCooler) -> list[dict[str, object]]:
    """Describe every weight column of a map: name, method, and attributes.

    Returns:
        One dict per weight column, in ``/bins`` order, with ``name``,
        ``method`` (``None`` for a column metacooler did not write),
        ``n_masked`` (bins whose weight is ``NaN``) and ``attrs`` (every
        attribute, decoded).
    """
    rows: list[dict[str, object]] = []
    with mc.open() as grp:
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
