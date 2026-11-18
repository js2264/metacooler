"""The ``/clusters`` group: the clustering rounds that led to a map's genomes.

A tool that groups contigs into genomes writes the rounds of its latest run
here (:func:`write_clusters`): metaTOR's overlapping, recursive and final
clusters, say. ``/genomes`` and ``/contigs/genome_id`` keep the current
assignment. :class:`metacooler.MetaCooler` reads the rounds back
(``cluster_rounds``, ``clusters``, ``cluster_members``), and the writers that
make a map from a map carry them (:func:`copy_clusters`). See
``docs/format.md#clusters`` for the layout.
"""

from __future__ import annotations

import json
import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from metacooler.io._utils import write_column
from metacooler.schema import (
    CLUSTERS_GROUP,
    MetaCoolerFormatError,
    _check_clusters,
    _check_clusters_structure,
    _decode_table_group,
    cluster_rounds,
    is_metacooler,
    validate,
)
from metacooler.uri import parse_uri

logger = logging.getLogger("metacooler")

#: Where :func:`write_clusters` builds the new group before swapping it in.
_SCRATCH = ".clusters-scratch"


@dataclass(frozen=True)
class _Clusters:
    """A ``/clusters`` group in memory, by row ids, as it is stored."""

    rounds: list[dict]
    #: ``table/``: ``round``, ``name``, ``parent``, ``n_contigs``, ``length``,
    #: then the tool's own columns.
    table: pd.DataFrame
    cluster_id: np.ndarray
    contig_id: np.ndarray


def _rounds(rounds: Sequence[str | Mapping[str, object]]) -> list[dict]:
    """The rounds as stored: one dict per round, each with a ``name``."""
    out = []
    for item in rounds:
        if not isinstance(item, (str, Mapping)):
            raise ValueError(
                "Every round is a name, or a mapping with a name and what its "
                f"tool stores with it: got {item!r}."
            )
        entry = {"name": item} if isinstance(item, str) else dict(item)
        name = entry.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(
                f"Every round needs a name, a non-empty string: got {item!r}."
            )
        if not name.isascii():
            raise ValueError(f"Round names must be ASCII: {name!r}.")
        out.append(entry)
    names = [str(entry["name"]) for entry in out]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ValueError(f"Round names must be unique: {', '.join(repeated)}.")
    try:
        json.dumps(out)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"The rounds must be JSON-serialisable: {exc}") from exc
    return out


def _first(values: Sequence[object], n: int = 5) -> str:
    """The first *n* of *values*, for an error message."""
    shown = ", ".join(str(value) for value in list(values)[:n])
    return shown + (", ..." if len(values) > n else "")


def _build(
    rounds: Sequence[str | Mapping[str, object]],
    clusters: pd.DataFrame,
    members: pd.DataFrame,
    contig_names: np.ndarray,
    contig_lengths: np.ndarray,
) -> _Clusters:
    """Resolve the frames :func:`write_clusters` takes to the stored row ids."""
    stored_rounds = _rounds(rounds)
    round_index = {entry["name"]: i for i, entry in enumerate(stored_rounds)}

    for column in ("round", "name"):
        if column not in clusters.columns:
            raise ValueError(f"The clusters table needs a {column!r} column.")
    has_parent = {"parent", "parent_round"} & set(clusters.columns)
    if has_parent and has_parent != {"parent", "parent_round"}:
        raise ValueError(
            "A cluster's parent is given by two columns, 'parent_round' and "
            "'parent': the round and the name of the cluster it was split from."
        )
    reserved = {"n_contigs", "length"}
    extra = [
        c
        for c in clusters.columns
        if c not in {"round", "name", "parent", "parent_round"}
    ]
    if reserved & set(extra):
        raise ValueError(
            "These columns of the clusters table are computed from the members, "
            f"and cannot be given: {', '.join(sorted(reserved & set(extra)))}."
        )
    unnamed = [repr(str(c)) for c in extra if str(c) in ("", ".") or "/" in str(c)]
    if unnamed:
        raise ValueError(
            "Column names must be HDF5 dataset names: not empty, not '.', and "
            f"without '/' (HDF5 nests them): {', '.join(unnamed)}."
        )
    unstored = [
        f"{c} ({clusters[c].dtype})"
        for c in extra
        if clusters[c].dtype.kind in "mMc"
        or isinstance(
            clusters[c].dtype, (pd.CategoricalDtype, pd.PeriodDtype, pd.IntervalDtype)
        )
    ]
    if unstored:
        raise ValueError(
            "Columns are stored as numbers or as text: convert these (to text, "
            f"say) first: {', '.join(unstored)}."
        )
    gaps = [
        f"{c} ({clusters[c].dtype})"
        for c in extra
        if isinstance(clusters[c].dtype, pd.api.extensions.ExtensionDtype)
        and clusters[c].dtype.kind in "biu"
        and clusters[c].isna().any()
    ]
    if gaps:
        raise ValueError(
            "Nullable boolean or integer columns hold missing values, which "
            "HDF5 stores in neither: fill them, or convert the columns to "
            f"floats or text: {', '.join(gaps)}."
        )

    _refuse_missing(clusters, "clusters", ("round", "name"))
    frame = clusters.reset_index(drop=True)
    round_names = frame["round"].astype(str)
    unknown = sorted(set(round_names) - set(round_index))
    if unknown:
        raise ValueError(f"Clusters of undeclared rounds: {_first(unknown)}.")
    names = frame["name"].astype(str)
    if not all(name.isascii() and name for name in names):
        raise ValueError("Cluster names must be non-empty ASCII strings.")
    # Rows come round by round, in the order of the rounds; a round's clusters
    # keep the order they were given in.
    order = np.argsort(round_names.map(round_index).to_numpy(), kind="stable")
    frame = frame.iloc[order].reset_index(drop=True)
    round_names = round_names.iloc[order].reset_index(drop=True)
    names = names.iloc[order].reset_index(drop=True)
    keys = list(zip(round_names, names, strict=True))
    row = {key: i for i, key in enumerate(keys)}
    if len(row) != len(keys):
        repeated = sorted(k for k, n in Counter(keys).items() if n > 1)
        raise ValueError(
            "Cluster names must be unique within their round: "
            f"{_first([f'{r}/{n}' for r, n in repeated])}."
        )
    round_ids = round_names.map(round_index).to_numpy(np.int32)

    parents = np.full(len(frame), -1, dtype=np.int32)
    if has_parent:
        for i, (parent_round, parent) in enumerate(
            zip(frame["parent_round"], frame["parent"], strict=True)
        ):
            if _missing(parent_round) and _missing(parent):
                continue
            if _missing(parent_round) or _missing(parent):
                raise ValueError(
                    f"Cluster {keys[i][0]}/{keys[i][1]} gives a parent with "
                    "only one of 'parent_round' and 'parent': a parent needs "
                    "both, and a cluster without one neither."
                )
            key = (str(parent_round), str(parent))
            if key not in row:
                raise ValueError(
                    f"Cluster {keys[i][0]}/{keys[i][1]} names a parent that is not "
                    f"a cluster: {key[0]}/{key[1]}.{_numbers_hint(frame)}"
                )
            if round_index[key[0]] >= round_ids[i]:
                raise ValueError(
                    f"Cluster {keys[i][0]}/{keys[i][1]}: its parent, "
                    f"{key[0]}/{key[1]}, must be of an earlier round."
                )
            parents[i] = row[key]

    for column in ("round", "cluster", "contig"):
        if column not in members.columns:
            raise ValueError(f"The members table needs a {column!r} column.")
    _refuse_missing(members, "members", ("round", "cluster", "contig"))
    member_keys = list(
        zip(members["round"].astype(str), members["cluster"].astype(str), strict=True)
    )
    unknown_clusters = sorted({key for key in member_keys if key not in row})
    if unknown_clusters:
        raise ValueError(
            "Members of clusters that are not in the clusters table: "
            f"{_first([f'{r}/{n}' for r, n in unknown_clusters])}."
        )
    cluster_ids = np.fromiter(
        (row[key] for key in member_keys), dtype=np.int64, count=len(member_keys)
    )
    contig_index = pd.Index(contig_names)
    contig_ids = contig_index.get_indexer(members["contig"].astype(str))
    if np.any(contig_ids < 0):
        missing = sorted(set(members["contig"].astype(str)[contig_ids < 0]))
        raise ValueError(f"Members that are not contigs of the map: {_first(missing)}.")
    key = cluster_ids * max(len(contig_names), 1) + contig_ids
    repeated = np.flatnonzero(pd.Series(key).duplicated().to_numpy())
    if repeated.size:
        first = int(repeated[0])
        raise ValueError(
            f"{repeated.size} membership(s) given more than once, e.g. contig "
            f"{contig_names[contig_ids[first]]} of cluster "
            f"{member_keys[first][0]}/{member_keys[first][1]}."
        )
    table = pd.DataFrame(
        {
            "round": round_ids,
            "name": names.to_numpy(),
            "parent": parents,
        }
    )
    for column in extra:
        table[column] = frame[column].to_numpy()
    return _sorted(stored_rounds, table, cluster_ids, contig_ids, contig_lengths)


def _missing(value: object) -> bool:
    """Whether *value* says "no parent": a missing value (``None``, ``NaN``,
    ``pd.NA``, ``NaT``) or an empty string."""
    if isinstance(value, str):
        return value == ""
    return bool(pd.isna(value))


def _refuse_missing(table: pd.DataFrame, what: str, columns: Sequence[str]) -> None:
    """Refuse a missing value in a key column of *table*, which ``astype(str)``
    would otherwise turn into a name (``"nan"``, ``"None"``)."""
    for column in columns:
        values = table[column]
        missing = values.isna().to_numpy(bool, copy=True)
        missing[~missing] = values[~missing].eq("").to_numpy(bool)
        if missing.any():
            rows = np.flatnonzero(missing)
            raise ValueError(
                f"The {what} table has {len(rows)} row(s) without a {column!r}, "
                f"e.g. row {int(rows[0])}: every row needs one."
            )


def _numbers_hint(frame: pd.DataFrame) -> str:
    """A hint for parents read as numbers: ``read_csv`` makes ``1.0`` of a
    cluster called ``1`` in a column that also holds missing values."""
    if any(frame[c].dtype.kind in "fiu" for c in ("parent_round", "parent")):
        return (
            " The parent columns hold numbers: read them as text (for "
            "pandas.read_csv, dtype=str), as cluster and round names are."
        )
    return ""


def _sorted(
    rounds: list[dict],
    table: pd.DataFrame,
    cluster_ids: np.ndarray,
    contig_ids: np.ndarray,
    contig_lengths: np.ndarray,
) -> _Clusters:
    """The members sorted, each pair once, and the derived columns set."""
    n_clusters = len(table)
    n_contigs = len(contig_lengths)
    key = cluster_ids.astype(np.int64) * max(n_contigs, 1) + contig_ids.astype(np.int64)
    order = np.argsort(key, kind="stable")
    key = key[order]
    if np.any(np.diff(key) == 0):
        repeated = key[1:][np.diff(key) == 0]
        raise ValueError(
            f"{len(repeated)} membership(s) given more than once, e.g. contig "
            f"{int(repeated[0] % max(n_contigs, 1))} of cluster row "
            f"{int(repeated[0] // max(n_contigs, 1))}."
        )
    cluster_ids = cluster_ids[order].astype(np.int32)
    contig_ids = contig_ids[order].astype(np.int32)
    derived = pd.DataFrame(
        {
            "n_contigs": np.bincount(cluster_ids, minlength=n_clusters).astype(
                np.int32
            ),
            "length": np.bincount(
                cluster_ids,
                weights=contig_lengths[contig_ids].astype(np.float64),
                minlength=n_clusters,
            ).astype(np.int64),
        }
    )
    head = table[["round", "name", "parent"]].reset_index(drop=True)
    rest = table.drop(
        columns=["round", "name", "parent", "n_contigs", "length"], errors="ignore"
    ).reset_index(drop=True)
    table = pd.concat([head, derived, rest], axis=1)
    return _Clusters(rounds, table, cluster_ids, contig_ids)


def _write_group(group: h5py.Group, clusters: _Clusters) -> None:
    """Write *clusters* into the empty *group* (see docs/format.md#clusters)."""
    group.attrs["rounds"] = json.dumps(clusters.rounds)
    table = group.create_group("table", track_order=True)
    dtypes = {"round": np.int32, "parent": np.int32, "n_contigs": np.int32}
    for column in clusters.table.columns:
        write_column(
            table, str(column), clusters.table[column].to_numpy(), dtypes.get(column)
        )
    members = group.create_group("members", track_order=True)
    write_column(members, "cluster_id", clusters.cluster_id, np.int32)
    write_column(members, "contig_id", clusters.contig_id, np.int32)
    offsets = np.zeros(len(clusters.table) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(
        np.bincount(clusters.cluster_id, minlength=len(clusters.table))
    )
    write_column(group.create_group("indexes"), "cluster_offset", offsets, np.int64)


def _read_group(group: h5py.Group) -> _Clusters:
    """Read a stored ``/clusters`` group, by row ids."""
    return _Clusters(
        cluster_rounds(group),
        _decode_table_group(group["table"]),
        np.asarray(group["members"]["cluster_id"][()]).astype(np.int32),
        np.asarray(group["members"]["contig_id"][()]).astype(np.int32),
    )


def _map_group(uri: str | Path) -> tuple[Path, str]:
    """The file and the HDF5 group of the map *uri* addresses."""
    path, group = parse_uri(str(uri))
    return Path(path), group.strip("/")


def write_clusters(
    uri: str | Path,
    rounds: Sequence[str | Mapping[str, object]],
    clusters: pd.DataFrame,
    members: pd.DataFrame,
) -> None:
    """Replace the clustering rounds of a map, in place.

    The map keeps one run: its previous ``/clusters``, if any, is replaced
    whole. ``/genomes`` and ``/contigs`` are left alone; to make a round the
    map's genomes, reassign the contigs to it (:func:`metacooler.reassign`).
    The map is checked in full first, and the new group is built aside,
    checked, and swapped in once complete, so a failure leaves the map as it
    was. As every edit of an HDF5 file in place, it leaves the space of the
    group it replaces in the file, which a copy of the map reclaims.

    Args:
        uri: A ``.metacool`` path, or an entry of an ``.mmetacool``
            (``file.mmetacool::entries/<name>``).
        rounds: The rounds, in the order they were made: a name each, or a
            mapping with a ``name`` and any JSON-serialisable keys (by
            convention ``tool``, ``version``, ``parameters`` and ``date``).
        clusters: One row per cluster: ``round`` (a round's name) and
            ``name`` (unique within its round), optionally ``parent_round`` and
            ``parent`` (the round and name of the cluster it was split from,
            of an earlier round; empty for none), and any other column, stored
            as given (completeness, contamination, ...).
        members: One row per membership: ``round``, ``cluster`` (its name)
            and ``contig`` (a contig name of the map). A contig may be a
            member of several clusters.

    Raises:
        ValueError: If a round, a cluster, a parent or a member is not what
            the description above says, or a column of *clusters* is one
            metacooler computes (``n_contigs``, ``length``).
    """
    from metacooler.api import MetaCooler

    # The map is checked in full before anything is written, so that the new
    # group, itself checked before the swap, is the only change, and a
    # failure leaves the map as it was.
    try:
        mc = MetaCooler(uri)
        validate(uri)
    except MetaCoolerFormatError as exc:
        _hint_broken_clusters(uri, exc)
        raise
    contigs = mc.contigs
    contig_lengths = contigs["length"].to_numpy(np.int64)
    built = _build(
        rounds,
        clusters,
        members,
        contigs["name"].astype(str).to_numpy(),
        contig_lengths,
    )
    path, key = _map_group(uri)
    with h5py.File(path, "r+") as f:
        grp = f[key] if key else f
        if _SCRATCH in grp:
            del grp[_SCRATCH]
        try:
            scratch = grp.create_group(_SCRATCH)
            _write_group(scratch, built)
            # Checked before it replaces the map's group: a group validate()
            # would refuse never reaches the map.
            _check_clusters_structure(scratch)
            _check_clusters(scratch, len(contig_lengths), contig_lengths)
        except BaseException:
            if _SCRATCH in grp:
                del grp[_SCRATCH]
            raise
        if CLUSTERS_GROUP in grp:
            del grp[CLUSTERS_GROUP]
        grp.move(_SCRATCH, CLUSTERS_GROUP)
    logger.info(
        "Wrote %d cluster(s) in %d round(s) to %s",
        len(built.table),
        len(built.rounds),
        uri,
    )


def drop_clusters(uri: str | Path) -> bool:
    """Remove the clustering rounds of a map, in place.

    Args:
        uri: A ``.metacool`` path, or an ``.mmetacool`` entry URI.

    Returns:
        Whether the map had any.

    Raises:
        KeyError: If *uri* is a ``.mmetacool`` collection or names no entry
            of one, as :class:`~metacooler.api.MetaCooler` raises.
        ValueError: If *uri* is not a ``.metacool`` map (a plain cooler).
    """
    from metacooler.api import MetaCooler

    # Resolves the URI without checking the map: drop_clusters is what
    # removes a /clusters group that keeps it from opening.
    MetaCooler(uri, validate=False)
    if not is_metacooler(uri):
        raise ValueError(
            f"{uri} is not a .metacool map: its 'clusters' group, if any, is "
            "not metacooler's, and is left alone."
        )
    path, key = _map_group(uri)
    with h5py.File(path, "r+") as f:
        grp = f[key] if key else f
        if CLUSTERS_GROUP not in grp:
            return False
        del grp[CLUSTERS_GROUP]
    return True


def _hint_broken_clusters(uri: str | Path, error: MetaCoolerFormatError) -> None:
    """Name ``drop_clusters`` when the map's own ``/clusters`` fails the check,
    since writing new rounds would replace it anyway."""
    path, key = _map_group(uri)
    with h5py.File(path, "r") as f:
        grp = f[key] if key else f
        node = grp.get(CLUSTERS_GROUP)
        lengths = grp.get("contigs/length")
        # Without /contigs/length, the group cannot be checked, and the error
        # is the map's: no hint.
        if node is None or not isinstance(lengths, h5py.Dataset):
            return
        lengths = np.asarray(lengths[()]).astype(np.int64)
        try:
            if not isinstance(node, h5py.Group):
                raise MetaCoolerFormatError("/clusters is not a group")
            _check_clusters_structure(node)
            _check_clusters(node, len(lengths), lengths)
        except MetaCoolerFormatError:
            raise MetaCoolerFormatError(
                f"{error}. The /clusters group of {uri} fails the check too: "
                f"drop_clusters({str(uri)!r}) removes it."
            ) from error


def copy_clusters(
    source: h5py.Group,
    target: h5py.Group,
    contig_map: np.ndarray | None = None,
    contig_lengths: np.ndarray | None = None,
) -> bool:
    """Carry a map's ``/clusters`` into the map a writer makes from it.

    Args:
        source: The source map's group (the file root for a ``.metacool``).
        target: The new map's group, which must not hold ``/clusters`` yet.
        contig_map: The new row of each source contig, when the writer moved
            them (``reassign`` re-grouping the contigs). ``None`` copies the
            group as stored.
        contig_lengths: The new map's ``/contigs/length``; required with
            *contig_map*.

    Returns:
        Whether the source had clusters.
    """
    if CLUSTERS_GROUP not in source:
        return False
    if contig_map is None:
        source.copy(source[CLUSTERS_GROUP], target, name=CLUSTERS_GROUP)
        return True
    if contig_lengths is None:
        raise ValueError(
            "copy_clusters needs the new contig lengths with a contig map."
        )
    stored = _read_group(source[CLUSTERS_GROUP])
    remapped = _sorted(
        stored.rounds,
        stored.table,
        stored.cluster_id,
        np.asarray(contig_map, dtype=np.int64)[stored.contig_id],
        np.asarray(contig_lengths, dtype=np.int64),
    )
    _write_group(target.create_group(CLUSTERS_GROUP), remapped)
    return True


def clusters_digest(group: h5py.Group) -> str | None:
    """A digest of a map's ``/clusters``, to compare two maps' (``merge``).

    Returns:
        ``None`` for a map without clusters.
    """
    import hashlib

    if CLUSTERS_GROUP not in group:
        return None
    # Of the values, not of how they are stored: text of fixed or variable
    # length, integers of any width, columns in any order.
    stored = _read_group(group[CLUSTERS_GROUP])
    digest = hashlib.sha256(json.dumps(stored.rounds, sort_keys=True).encode())
    for column in sorted(stored.table.columns, key=str):
        values = json.dumps(stored.table[column].tolist(), default=str)
        digest.update(f"{column}:{values}".encode())
    for ids in (stored.cluster_id, stored.contig_id):
        digest.update(np.asarray(ids, dtype=np.int64).tobytes())
    return digest.hexdigest()
