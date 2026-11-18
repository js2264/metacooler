"""The /clusters group: a map's clustering rounds (docs/format.md#clusters).

A tool that groups contigs into genomes (metaTOR) writes the rounds of its
latest run with ``write_clusters``; ``MetaCooler`` reads them back; ``validate``
checks them; and the writers that make a map from a map carry them.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import metacooler
from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io import clusters as clusters_module
from metacooler.io.annotate import annotate
from metacooler.io.clusters import copy_clusters, drop_clusters, write_clusters
from metacooler.io.create import from_pairs
from metacooler.io.merge import merge
from metacooler.io.mmetacool import add_entry, create_mmetacool
from metacooler.io.reassign import reassign
from metacooler.ops.binnify import binnify
from metacooler.ops.zoomify import zoomify
from metacooler.schema import MetaCoolerFormatError, validate

LENGTHS = {"c1": 4_000, "c2": 3_000, "c3": 2_000, "c4": 1_000}
ASSIGNMENT = {"c1": "MAG_A", "c2": "MAG_A", "c3": "MAG_B", "c4": "MAG_B"}

ROUNDS = [
    {"name": "overlapping", "tool": "metaTOR", "parameters": {"iterations": 100}},
    "recursive",
    {"name": "final", "tool": "metaTOR"},
]
CLUSTERS = pd.DataFrame(
    {
        "round": ["overlapping", "overlapping", "recursive", "recursive", "final"],
        "name": ["O1", "O2", "O1_1", "O1_2", "F1"],
        "parent_round": [None, None, "overlapping", "overlapping", "recursive"],
        "parent": [None, None, "O1", "O1", "O1_1"],
        "completeness": [90.5, 10.0, 60.0, 30.0, 60.0],
        "lineage": ["Bacteria", "Archaea", "Bacteria", "Bacteria", "Bacteria"],
    }
)
MEMBERS = pd.DataFrame(
    {
        "round": ["overlapping"] * 4 + ["recursive"] * 3 + ["final"] * 2,
        "cluster": ["O1", "O1", "O1", "O2", "O1_2", "O1_1", "O1_1", "F1", "F1"],
        "contig": ["c3", "c1", "c2", "c4", "c3", "c2", "c1", "c1", "c2"],
    }
)


def _map(path: Path, binsize: int | None = None):
    pairs = path.with_suffix(".pairs")
    rng = np.random.default_rng(1)
    names = list(LENGTHS)
    rows = []
    for i in range(400):
        a, b = names[rng.integers(4)], names[rng.integers(4)]
        rows.append(
            f"r{i}\t{a}\t{rng.integers(1, LENGTHS[a])}\t{b}\t{rng.integers(1, LENGTHS[b])}\n"
        )
    pairs.write_text(
        "## pairs format v1.0\n"
        + "".join(f"#chromsize: {c} {n}\n" for c, n in LENGTHS.items())
        + "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(rows)
    )
    bins = None
    if binsize is None:
        bins = pd.DataFrame(
            {"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())}
        )
    from_pairs(
        path,
        bins=bins,
        binsize=binsize,
        pairs=pairs,
        contig_genome_map=ASSIGNMENT,
        ps=False,
    )
    return path


@pytest.fixture
def clustered(tmp_path: Path) -> Path:
    path = _map(tmp_path / "map.metacool")
    write_clusters(path, ROUNDS, CLUSTERS, MEMBERS)
    return path


def _by_name(mc: MetaCooler) -> pd.DataFrame:
    """The memberships of a map, as names, in a stable order."""
    return (
        mc.cluster_members()
        .sort_values(["round", "cluster", "contig"])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Writing and reading
# ---------------------------------------------------------------------------


def test_clusters_written_are_read_back(clustered: Path):
    validate(clustered, full=True)
    mc = MetaCooler(clustered)

    assert mc.cluster_rounds == [
        {"name": "overlapping", "tool": "metaTOR", "parameters": {"iterations": 100}},
        {"name": "recursive"},
        {"name": "final", "tool": "metaTOR"},
    ]
    table = mc.clusters
    assert list(table.columns) == [
        "round",
        "name",
        "parent_round",
        "parent",
        "n_contigs",
        "length",
        "completeness",
        "lineage",
    ]
    assert table["name"].tolist() == ["O1", "O2", "O1_1", "O1_2", "F1"]
    assert table["n_contigs"].tolist() == [3, 1, 2, 1, 2]
    assert table["length"].tolist() == [9_000, 1_000, 7_000, 2_000, 7_000]
    assert table["parent"].tolist()[2:] == ["O1", "O1", "O1_1"]
    assert table["parent_round"].tolist()[2:] == [
        "overlapping",
        "overlapping",
        "recursive",
    ]
    assert table[["parent", "parent_round"]].iloc[:2].isna().all().all()
    assert table["completeness"].tolist() == [90.5, 10.0, 60.0, 30.0, 60.0]
    assert table["lineage"].tolist()[1] == "Archaea"

    members = mc.cluster_members("recursive")
    assert members.to_dict("list") == {
        "round": ["recursive"] * 3,
        "cluster": ["O1_1", "O1_1", "O1_2"],
        "contig": ["c1", "c2", "c3"],
    }
    assert len(mc.cluster_members()) == len(MEMBERS)
    # The map stays a cooler, which never looks at /clusters.
    assert cooler.Cooler(str(clustered)).info["nchroms"] == 4


def test_clusters_are_kept_in_round_order_whatever_the_rows_order(tmp_path: Path):
    path = _map(tmp_path / "map.metacool")
    write_clusters(path, ROUNDS, CLUSTERS.iloc[::-1], MEMBERS.iloc[::-1])
    table = MetaCooler(path).clusters
    assert table["round"].tolist() == [
        "overlapping",
        "overlapping",
        "recursive",
        "recursive",
        "final",
    ]
    # Within a round, the clusters keep the order they were given in.
    assert table["name"].tolist() == ["O2", "O1", "O1_2", "O1_1", "F1"]
    validate(path, full=True)


def test_rounds_stored_as_bytes_read_as_text(clustered: Path):
    """Another writer may store the attribute as bytes rather than a string."""
    rounds = MetaCooler(clustered).cluster_rounds
    with h5py.File(clustered, "r+") as f:
        f["clusters"].attrs["rounds"] = np.bytes_(json.dumps(rounds))
    assert MetaCooler(clustered).cluster_rounds == rounds
    validate(clustered, full=True)


def test_a_tools_text_column_of_variable_length_reads_as_text(clustered: Path):
    """h5py reads variable-length strings as bytes: they are decoded like the
    fixed-length ones."""
    with h5py.File(clustered, "r+") as f:
        f["clusters/table"].create_dataset(
            "note",
            data=["a", "b", "c", "d", "e"],
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    assert MetaCooler(clustered).clusters["note"].tolist() == ["a", "b", "c", "d", "e"]
    validate(clustered, full=True)


def test_a_map_without_clusters_reads_as_none(tmp_path: Path):
    mc = MetaCooler(_map(tmp_path / "map.metacool"))
    assert mc.cluster_rounds == []
    assert mc.clusters.empty
    assert list(mc.clusters.columns)[:2] == ["round", "name"]
    assert mc.cluster_members().empty
    with pytest.raises(KeyError, match="its rounds: none"):
        mc.cluster_members("final")


def test_an_unknown_round_is_refused_by_name(clustered: Path):
    with pytest.raises(KeyError, match="'partition' is not a clustering round"):
        MetaCooler(clustered).cluster_members("partition")


def test_a_round_without_clusters_reads_as_empty(tmp_path: Path):
    path = _map(tmp_path / "map.metacool")
    write_clusters(path, [*ROUNDS, "empty"], CLUSTERS, MEMBERS)
    mc = MetaCooler(path)
    assert mc.cluster_members("empty").empty
    assert [r["name"] for r in mc.cluster_rounds][-1] == "empty"


def test_writing_clusters_replaces_the_previous_run(clustered: Path):
    write_clusters(
        clustered,
        ["only"],
        pd.DataFrame({"round": ["only"], "name": ["X"]}),
        pd.DataFrame({"round": ["only"], "cluster": ["X"], "contig": ["c4"]}),
    )
    mc = MetaCooler(clustered)
    assert mc.cluster_rounds == [{"name": "only"}]
    assert mc.clusters["name"].tolist() == ["X"]
    validate(clustered, full=True)


def test_clusters_of_an_entry_of_a_collection(tmp_path: Path):
    collection = tmp_path / "maps.mmetacool"
    create_mmetacool(collection, {"rep1": _map(tmp_path / "map.metacool")})
    uri = f"{collection}::entries/rep1"
    write_clusters(uri, ROUNDS, CLUSTERS, MEMBERS)
    assert MetaCooler(uri).clusters["name"].tolist() == [
        "O1",
        "O2",
        "O1_1",
        "O1_2",
        "F1",
    ]
    assert drop_clusters(uri)
    assert MetaCooler(uri).cluster_rounds == []


def test_drop_clusters(clustered: Path):
    assert drop_clusters(clustered) is True
    assert MetaCooler(clustered).cluster_rounds == []
    assert drop_clusters(clustered) is False
    validate(clustered, full=True)


def test_the_package_exports_the_writers():
    assert metacooler.write_clusters is write_clusters
    assert metacooler.drop_clusters is drop_clusters


def test_a_failed_write_leaves_the_map_as_it_was(
    clustered: Path, monkeypatch: pytest.MonkeyPatch
):
    before = _by_name(MetaCooler(clustered))

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(clusters_module, "_write_group", interrupted)
    with pytest.raises(KeyboardInterrupt):
        write_clusters(clustered, ["only"], CLUSTERS.iloc[:0], MEMBERS.iloc[:0])

    pd.testing.assert_frame_equal(_by_name(MetaCooler(clustered)), before)
    with h5py.File(clustered, "r") as f:
        assert ".clusters-scratch" not in f


def test_a_scratch_group_left_behind_is_cleared(clustered: Path):
    with h5py.File(clustered, "r+") as f:
        f.create_group(".clusters-scratch/table")
    write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)
    with h5py.File(clustered, "r") as f:
        assert ".clusters-scratch" not in f


@pytest.mark.parametrize(
    ("rounds", "clusters", "members", "message"),
    [
        ([{"tool": "x"}], CLUSTERS, MEMBERS, "needs a name"),
        (["über"], CLUSTERS, MEMBERS, "must be ASCII"),
        (["final", "final"], CLUSTERS, MEMBERS, "must be unique: final"),
        ([{"name": "a", "when": object()}], CLUSTERS, MEMBERS, "JSON-serialisable"),
        (ROUNDS, CLUSTERS.drop(columns="round"), MEMBERS, "needs a 'round' column"),
        (ROUNDS, CLUSTERS.drop(columns="parent"), MEMBERS, "given by two columns"),
        (ROUNDS, CLUSTERS.assign(n_contigs=1), MEMBERS, "computed from the members"),
        (
            ROUNDS,
            CLUSTERS.rename(columns={"lineage": "a/b"}),
            MEMBERS,
            "without '/'",
        ),
        (["final"], CLUSTERS, MEMBERS, "undeclared rounds: overlapping"),
        (ROUNDS, CLUSTERS.assign(name=["O1", "O1", "a", "b", "c"]), MEMBERS, "unique"),
        (ROUNDS, CLUSTERS.assign(name=["O1", "é", "a", "b", "c"]), MEMBERS, "ASCII"),
        (
            ROUNDS,
            CLUSTERS.assign(parent=[None, None, "O9", "O1", "O1_1"]),
            MEMBERS,
            "names a parent that is not a cluster: overlapping/O9",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(parent_round=[None, None, "o", "o", "r"]),
            MEMBERS,
            "names a parent that is not a cluster: o/O1",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(
                parent_round=[None, None, "final", "overlapping", "recursive"],
                parent=[None, None, "F1", "O1", "O1_1"],
            ),
            MEMBERS,
            "must be of an earlier round",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(name=["O1", None, "a", "b", "c"]),
            MEMBERS,
            "clusters table has 1 row\\(s\\) without a 'name', e.g. row 1",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(round=[np.nan, *CLUSTERS["round"][1:]]),
            MEMBERS,
            "clusters table has 1 row\\(s\\) without a 'round', e.g. row 0",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(
                parent_round=[None, None, "overlapping", None, "recursive"]
            ),
            MEMBERS,
            "recursive/O1_2 gives a parent with only one of",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(parent=[None, None, "O1", "O1", ""]),
            MEMBERS,
            "final/F1 gives a parent with only one of",
        ),
        (
            ROUNDS,
            CLUSTERS,
            MEMBERS.assign(contig=[*MEMBERS["contig"][:-1], None]),
            "members table has 1 row\\(s\\) without a 'contig', e.g. row 8",
        ),
        (
            ROUNDS,
            CLUSTERS,
            MEMBERS.assign(cluster=pd.array([pd.NA] * 9, dtype="string")),
            "members table has 9 row\\(s\\) without a 'cluster', e.g. row 0",
        ),
        (ROUNDS, CLUSTERS, MEMBERS.drop(columns="contig"), "needs a 'contig' column"),
        (
            ROUNDS,
            CLUSTERS,
            MEMBERS.assign(cluster=[*MEMBERS["cluster"][:-1], "Z9"]),
            "not in the clusters table: final/Z9",
        ),
        (
            ROUNDS,
            CLUSTERS,
            MEMBERS.assign(contig=[*MEMBERS["contig"][:-1], "c9"]),
            "not contigs of the map: c9",
        ),
        (
            ROUNDS,
            CLUSTERS,
            pd.concat([MEMBERS, MEMBERS.iloc[:2]]),
            "2 membership\\(s\\) given more than once, e.g. contig c3 of cluster "
            "overlapping/O1",
        ),
        ([1, 2], CLUSTERS, MEMBERS, "Every round is a name, or a mapping"),
        (ROUNDS, CLUSTERS.rename(columns={"lineage": ""}), MEMBERS, "not empty"),
        (ROUNDS, CLUSTERS.rename(columns={"lineage": "."}), MEMBERS, "not '.'"),
        (
            ROUNDS,
            CLUSTERS.assign(ok=pd.array([True, None, False, True, True], "boolean")),
            MEMBERS,
            "Nullable boolean or integer columns hold missing values",
        ),
        (
            ROUNDS,
            CLUSTERS.assign(when=pd.Timestamp("2026-01-01")),
            MEMBERS,
            "stored as numbers or as text: convert these \\(to text, say\\) "
            "first: when \\(datetime64",
        ),
    ],
)
def test_bad_clusters_are_refused_before_anything_is_written(
    tmp_path: Path, rounds, clusters, members, message
):
    path = _map(tmp_path / "map.metacool")
    with pytest.raises(ValueError, match=message):
        write_clusters(path, rounds, clusters, members)
    assert MetaCooler(path).cluster_rounds == []


def test_parents_missing_as_pandas_na_are_no_parents(tmp_path: Path):
    """Guard: ``pd.NA`` in the parent columns of a ``string`` frame means no
    parent, as ``None`` does (it raised "boolean value of NA is ambiguous")."""
    path = _map(tmp_path / "map.metacool")
    clusters = CLUSTERS.astype({"parent_round": "string", "parent": "string"})
    write_clusters(path, ROUNDS, clusters, MEMBERS)
    stored = MetaCooler(path).clusters
    assert stored["parent"].isna().tolist() == [True, True, False, False, False]


def test_parents_read_as_numbers_get_a_hint(tmp_path: Path):
    """``read_csv`` reads clusters called 1 and 2, in a parent column that
    also holds missing values, as 1.0 and 2.0: the refusal says why."""
    path = _map(tmp_path / "map.metacool")
    clusters = pd.DataFrame(
        {
            "round": ["a", "b"],
            "name": ["1", "2"],
            "parent_round": [None, "a"],
            "parent": [np.nan, 1.0],
        }
    )
    members = pd.DataFrame(
        {"round": ["a", "b"], "cluster": ["1", "2"], "contig": ["c1", "c1"]}
    )
    with pytest.raises(ValueError, match="a/1.0. The parent columns hold numbers"):
        write_clusters(path, ["a", "b"], clusters, members)
    write_clusters(path, ["a", "b"], clusters.astype({"parent": "Int64"}), members)


def test_member_columns_may_name_columns_of_the_clusters_table(tmp_path: Path):
    """``cluster_id`` and ``contig_id`` are columns of ``members/``, not of
    ``table/``: a tool may use either name for a column of its own."""
    path = _map(tmp_path / "map.metacool")
    write_clusters(path, ROUNDS, CLUSTERS.assign(cluster_id=range(5)), MEMBERS)
    assert MetaCooler(path).clusters["cluster_id"].tolist() == [0, 1, 2, 3, 4]


def test_names_that_join_alike_in_two_rounds_are_distinct(tmp_path: Path):
    """Guard: cluster 0x of round 1 and cluster x of round 10 are two clusters
    (validate() joined round and name into one byte string, and both made
    b"10x")."""
    path = _map(tmp_path / "map.metacool")
    rounds = [f"r{i}" for i in range(11)]
    clusters = pd.DataFrame({"round": ["r1", "r10"], "name": ["0x", "x"]})
    members = pd.DataFrame(
        {"round": ["r1", "r10"], "cluster": ["0x", "x"], "contig": ["c1", "c2"]}
    )
    write_clusters(path, rounds, clusters, members)
    validate(path, full=True)
    assert MetaCooler(path).clusters["name"].tolist() == ["0x", "x"]


def test_a_group_validate_would_refuse_never_reaches_the_map(
    clustered: Path, monkeypatch: pytest.MonkeyPatch
):
    """Guard: the new group is checked before it replaces the map's, so a
    wrong one leaves the map as it was (it used to be swapped in, then
    refused by validate(), leaving a map that does not open)."""
    before = _by_name(MetaCooler(clustered))
    write_group = clusters_module._write_group

    def wrong(group, built):
        write_group(group, built)
        _replace(group["table"], "n_contigs", np.zeros(len(built.table), np.int32))

    monkeypatch.setattr(clusters_module, "_write_group", wrong)
    with pytest.raises(MetaCoolerFormatError, match="n_contigs does not match"):
        write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)
    pd.testing.assert_frame_equal(_by_name(MetaCooler(clustered)), before)
    validate(clustered, full=True)
    with h5py.File(clustered, "r") as f:
        assert ".clusters-scratch" not in f


@pytest.mark.parametrize(
    "edit",
    [
        lambda f: f["clusters"].__delitem__("indexes"),
        lambda f: (f.__delitem__("clusters"), f.create_dataset("clusters", data=[1])),
    ],
    ids=["broken-group", "not-a-group"],
)
def test_a_broken_clusters_group_is_named_with_its_remedy(clustered: Path, edit):
    with h5py.File(clustered, "r+") as f:
        edit(f)
    with pytest.raises(MetaCoolerFormatError, match="drop_clusters\\("):
        write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)
    assert drop_clusters(clustered)
    write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)
    validate(clustered, full=True)


@pytest.mark.parametrize("with_clusters", [True, False])
def test_a_map_broken_elsewhere_gets_no_clusters_hint(
    tmp_path: Path, with_clusters: bool
):
    path = _map(tmp_path / "map.metacool")
    if with_clusters:
        write_clusters(path, ROUNDS, CLUSTERS, MEMBERS)
    with h5py.File(path, "r+") as f:
        del f["indexes/genome_offset"]
    with pytest.raises(MetaCoolerFormatError) as caught:
        write_clusters(path, ROUNDS, CLUSTERS, MEMBERS)
    assert "drop_clusters" not in str(caught.value)


def test_drop_clusters_needs_one_map(tmp_path: Path):
    """Guard: a bare collection names its entries, and a missing entry is
    named (drop_clusters returned False for the first, and raised h5py's
    KeyError for the second)."""
    path = _map(tmp_path / "map.metacool")
    collection = tmp_path / "maps.mmetacool"
    create_mmetacool(collection, {"a": path})
    with pytest.raises(KeyError, match="Metacool entries found"):
        drop_clusters(collection)
    with pytest.raises(KeyError, match="entries/b"):
        drop_clusters(f"{collection}::entries/b")
    assert not drop_clusters(f"{collection}::entries/a")


def test_a_map_that_fails_the_check_is_left_as_it_was(clustered: Path):
    """The map is checked in full before the new group replaces its own (the
    group was swapped in, then the check of the rest of the map raised)."""
    with h5py.File(clustered, "r+") as f:
        lengths = f["genomes/length"][()]
        _replace(f["genomes"], "length", lengths + 1)
    with pytest.raises(MetaCoolerFormatError, match="length"):
        write_clusters(clustered, ["only"], CLUSTERS.iloc[:0], MEMBERS.iloc[:0])
    assert [entry["name"] for entry in MetaCooler(clustered).cluster_rounds] == [
        "overlapping",
        "recursive",
        "final",
    ]


def test_a_clusters_group_that_fails_the_full_check_is_named(clustered: Path):
    """A /clusters that opens (the structural check passes) but fails the full
    check is named with its remedy too."""
    _corrupt(clustered, lambda g: _replace(g["table"], "length", [1, 1, 1, 1, 1]))
    with pytest.raises(MetaCoolerFormatError, match="fails the check too: drop_clus"):
        write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)


def test_a_map_broken_without_its_contig_lengths_gets_no_clusters_hint(
    clustered: Path,
):
    """Without /contigs/length the group cannot be checked: the error is the
    map's, and names no remedy that would not help (it blamed /clusters,
    which dropping did not fix)."""
    with h5py.File(clustered, "r+") as f:
        del f["contigs/length"]
    with pytest.raises(MetaCoolerFormatError) as caught:
        write_clusters(clustered, ROUNDS, CLUSTERS, MEMBERS)
    assert "drop_clusters" not in str(caught.value)


def test_drop_clusters_leaves_a_plain_cooler_alone(tmp_path: Path):
    """A ``clusters`` group of a plain cooler is not metacooler's."""
    from _cools import ECOLI, make_cool

    cool = make_cool(tmp_path / "plain.cool", ECOLI)
    with h5py.File(cool, "r+") as f:
        f.create_group("clusters")
    with pytest.raises(ValueError, match="is not a .metacool map"):
        drop_clusters(cool)
    with h5py.File(cool, "r") as f:
        assert "clusters" in f


# ---------------------------------------------------------------------------
# validate()
# ---------------------------------------------------------------------------


def _corrupt(path: Path, edit) -> None:
    with h5py.File(path, "r+") as f:
        edit(f["clusters"])


def _replace(group: h5py.Group, column: str, values) -> None:
    del group[column]
    group.create_dataset(column, data=np.asarray(values))


CORRUPTIONS = {
    "no rounds attribute": (
        lambda g: g.attrs.__delitem__("rounds"),
        "no 'rounds' attribute",
    ),
    "rounds not JSON": (
        lambda g: g.attrs.__setitem__("rounds", "{not json"),
        "is not JSON",
    ),
    "rounds not a list": (
        lambda g: g.attrs.__setitem__("rounds", json.dumps({"name": "a"})),
        "must be a list of objects",
    ),
    "repeated round": (
        lambda g: g.attrs.__setitem__("rounds", json.dumps([{"name": "a"}] * 3)),
        "round names are not unique",
    ),
    "non-ASCII round": (
        lambda g: g.attrs.__setitem__(
            "rounds", json.dumps([{"name": "a"}, {"name": "b"}, {"name": "ü"}])
        ),
        "round names must be ASCII",
    ),
    "no members": (
        lambda g: g.__delitem__("members"),
        "Missing group: /clusters/members",
    ),
    "no column": (
        lambda g: g["table"].__delitem__("parent"),
        "Missing column in /clusters/table: parent",
    ),
    "no index": (
        lambda g: g.__delitem__("indexes"),
        "Missing index: /clusters/indexes/cluster_offset",
    ),
    "unequal columns": (
        lambda g: _replace(g["table"], "parent", [-1, -1]),
        "columns have unequal lengths",
    ),
    "offset ends": (
        lambda g: _replace(g["indexes"], "cluster_offset", [0, 3, 4, 6, 7, 8]),
        "cluster_offset endpoints",
    ),
    "round out of range": (
        lambda g: _replace(g["table"], "round", [0, 0, 1, 1, 7]),
        "round values out of range",
    ),
    "rounds not grouped": (
        lambda g: _replace(g["table"], "round", [0, 1, 0, 1, 2]),
        "not grouped by round",
    ),
    "non-ASCII name": (
        lambda g: _replace(
            g["table"], "name", np.array([b"O1", b"O2", "\xe9".encode(), b"x", b"y"])
        ),
        "name must be non-empty ASCII",
    ),
    "empty name": (
        lambda g: _replace(
            g["table"], "name", np.array(["O1", "O2", "", "b", "c"], dtype="S")
        ),
        "name must be non-empty ASCII",
    ),
    "non-ASCII name of variable length": (
        lambda g: (
            g["table"].__delitem__("name"),
            g["table"].create_dataset(
                "name",
                data=["O1", "O2", "\xe9", "b", "c"],
                dtype=h5py.string_dtype(),
            ),
        ),
        "name must be non-empty ASCII",
    ),
    "names that are not text": (
        lambda g: _replace(g["table"], "name", [1, 2, 3, 4, 5]),
        "/clusters/table/name must hold text",
    ),
    "rounds that are not integers": (
        lambda g: _replace(g["table"], "round", [0.0, 0.0, 1.0, 1.0, 2.0]),
        "/clusters/table/round must hold integers",
    ),
    "contig ids that are not integers": (
        lambda g: _replace(g["members"], "contig_id", np.linspace(0, 1, 9)),
        "/clusters/members/contig_id must hold integers",
    ),
    "names repeated in a round": (
        lambda g: _replace(
            g["table"], "name", np.array(["O1", "O1", "a", "b", "c"], dtype="S")
        ),
        "not unique within their round",
    ),
    "parent out of range": (
        lambda g: _replace(g["table"], "parent", [-1, -1, 0, 0, 9]),
        "parent values out of range",
    ),
    "parent of a later round": (
        lambda g: _replace(g["table"], "parent", [4, -1, 0, 0, 2]),
        "must be a cluster of an earlier round",
    ),
    "cluster id out of range": (
        lambda g: _replace(g["members"], "cluster_id", [0, 0, 0, 1, 2, 2, 3, 4, 9]),
        "cluster_id values out of range",
    ),
    "contig id out of range": (
        lambda g: _replace(g["members"], "contig_id", [0, 1, 2, 3, 0, 1, 2, 0, 9]),
        "contig_id values out of range",
    ),
    "members unsorted": (
        lambda g: _replace(g["members"], "contig_id", [1, 0, 2, 3, 0, 1, 2, 0, 1]),
        "not sorted by",
    ),
    "offsets inconsistent": (
        lambda g: _replace(g["indexes"], "cluster_offset", [0, 2, 4, 6, 7, 9]),
        "cluster_offset is inconsistent",
    ),
    "n_contigs wrong": (
        lambda g: _replace(g["table"], "n_contigs", [3, 1, 2, 1, 1]),
        "n_contigs does not match",
    ),
    "length wrong": (
        lambda g: _replace(g["table"], "length", [9_000, 1_000, 7_000, 2_000, 1]),
        "length does not match",
    ),
}


@pytest.mark.parametrize("name", sorted(CORRUPTIONS))
def test_validate_refuses_a_broken_clusters_group(clustered: Path, name: str):
    edit, message = CORRUPTIONS[name]
    _corrupt(clustered, edit)
    with pytest.raises(MetaCoolerFormatError, match=message):
        validate(clustered, full=True)


def test_validate_refuses_clusters_that_are_not_a_group(tmp_path: Path):
    path = _map(tmp_path / "map.metacool")
    with h5py.File(path, "r+") as f:
        f.create_dataset("clusters", data=[1])
    with pytest.raises(MetaCoolerFormatError, match="/clusters is not a group"):
        validate(path, full=False)


def test_opening_a_map_checks_the_structure_only(clustered: Path):
    """``MetaCooler()`` runs the O(n_columns) checks: a wrong length passes it,
    and the full check names it."""
    _corrupt(clustered, lambda g: _replace(g["table"], "length", [1, 1, 1, 1, 1]))
    MetaCooler(clustered)
    with pytest.raises(MetaCoolerFormatError, match="length does not match"):
        validate(clustered, full=True)


# ---------------------------------------------------------------------------
# The writers that make a map from a map carry the clusters
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("inplace", [False, True], ids=["new-file", "in-place"])
def test_a_relabelling_keeps_the_clusters(clustered: Path, inplace: bool):
    """A guard: the relabel path copies the map as stored, or edits it in
    place, so it kept the group before it was part of the format."""
    before = _by_name(MetaCooler(clustered))
    out = clustered if inplace else clustered.with_name("out.metacool")
    reassign(
        clustered,
        {"c1": "MAG_A", "c2": "MAG_A", "c3": "MAG_B", "c4": "MAG_C"},
        None if inplace else out,
        inplace=inplace,
        layout="keep",
    )
    pd.testing.assert_frame_equal(_by_name(MetaCooler(out)), before)
    validate(out, full=True)


def test_a_regrouping_renumbers_the_members(clustered: Path):
    """The contigs move, and the members with them: by name, nothing changes."""
    before = _by_name(MetaCooler(clustered))
    out = clustered.with_name("out.metacool")
    reassign(
        clustered,
        {"c1": "MAG_Z", "c2": "MAG_B", "c3": "MAG_A", "c4": "MAG_Z"},
        out,
        layout="genome-then-length",
    )
    mc = MetaCooler(out)
    assert mc.contigs["name"].tolist() != list(LENGTHS)  # the rows did move
    pd.testing.assert_frame_equal(_by_name(mc), before)
    pd.testing.assert_frame_equal(mc.clusters, MetaCooler(clustered).clusters)
    validate(out, full=True)


def test_an_annotation_keeps_the_clusters(clustered: Path, tmp_path: Path):
    """A guard: annotate edits /genomes and /contigs only, in place or in a
    copy of the map."""
    table = tmp_path / "gc.tsv"
    table.write_text("contig_name\tgc\nc1\t0.4\nc2\t0.5\nc3\t0.6\nc4\t0.7\n")
    before = _by_name(MetaCooler(clustered))
    annotate(clustered, contigs_table=table, inplace=True)
    pd.testing.assert_frame_equal(_by_name(MetaCooler(clustered)), before)


def test_a_merge_keeps_the_first_inputs_clusters(
    clustered: Path, caplog: pytest.LogCaptureFixture
):
    replica = clustered.with_name("replica.metacool")
    shutil.copy(clustered, replica)
    out = clustered.with_name("merged.metacool")
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([clustered, replica], out)
    assert "other clustering rounds" not in caplog.text
    pd.testing.assert_frame_equal(
        _by_name(MetaCooler(out)), _by_name(MetaCooler(clustered))
    )
    validate(out, full=True)


def test_a_merge_compares_the_clusters_not_their_storage(
    clustered: Path, caplog: pytest.LogCaptureFixture
):
    """The same clusters, names stored as text of variable length: not other
    clustering rounds (the digest hashed the stored bytes and dtypes)."""
    other = clustered.with_name("other.metacool")
    shutil.copy(clustered, other)
    with h5py.File(other, "r+") as f:
        names = [name.decode() for name in f["clusters/table/name"][()]]
        del f["clusters/table/name"]
        f["clusters/table"].create_dataset(
            "name", data=names, dtype=h5py.string_dtype()
        )
    validate(other, full=True)
    out = clustered.with_name("merged.metacool")
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([clustered, other], out)
    assert "other clustering rounds" not in caplog.text


def test_a_merge_warns_about_an_input_with_other_clusters(
    clustered: Path, caplog: pytest.LogCaptureFixture
):
    other = clustered.with_name("other.metacool")
    shutil.copy(clustered, other)
    drop_clusters(other)
    out = clustered.with_name("merged.metacool")
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([clustered, other], out)
    assert f"{other} has other clustering rounds than {clustered}" in caplog.text
    assert MetaCooler(out).clusters["name"].tolist() == [
        "O1",
        "O2",
        "O1_1",
        "O1_2",
        "F1",
    ]


def test_a_rebinning_keeps_the_clusters(tmp_path: Path):
    source = _map(tmp_path / "map.metacool", binsize=500)
    write_clusters(source, ROUNDS, CLUSTERS, MEMBERS)
    out = tmp_path / "coarse.metacool"
    binnify(source, out, 1_000)
    pd.testing.assert_frame_equal(
        _by_name(MetaCooler(out)), _by_name(MetaCooler(source))
    )
    validate(out, full=True)


def test_a_rebinning_by_genome_drops_the_clusters(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    """Guard: binned by genome, the map's contigs are its genomes, so the
    rounds, which group the source's contigs, are not carried (they were,
    with contig ids out of range: a map validate() refused)."""
    source = _map(tmp_path / "map.metacool", binsize=500)
    write_clusters(source, ROUNDS, CLUSTERS, MEMBERS)
    out = tmp_path / "by-genome.metacool"
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        binnify(source, out, 1_000, level="genome")
    assert MetaCooler(out).cluster_rounds == []
    validate(out, full=True)
    assert "binned by genome, has not: they are not carried" in caplog.text


def test_a_rebinning_by_genome_of_a_map_without_clusters_says_nothing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    source = _map(tmp_path / "map.metacool", binsize=500)
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        binnify(source, tmp_path / "by-genome.metacool", 1_000, level="genome")
    assert "clustering rounds" not in caplog.text


def test_every_level_of_a_pyramid_keeps_the_clusters(tmp_path: Path):
    source = _map(tmp_path / "map.metacool", binsize=500)
    write_clusters(source, ROUNDS, CLUSTERS, MEMBERS)
    out = tmp_path / "pyramid.mmetacool"
    zoomify(source, out, resolutions="500,1000,2000")
    expected = _by_name(MetaCooler(source))
    for level in (500, 1_000, 2_000):
        uri = f"{out}::entries/{level}"
        pd.testing.assert_frame_equal(_by_name(MetaCooler(uri)), expected)
        validate(uri, full=True)


def test_an_entry_added_to_a_collection_keeps_the_clusters(
    clustered: Path, tmp_path: Path
):
    """A guard: an entry is a copy of the map's group, whatever it holds."""
    collection = tmp_path / "maps.mmetacool"
    create_mmetacool(collection, {"rep1": clustered})
    add_entry(collection, "rep2", clustered)
    for entry in ("rep1", "rep2"):
        uri = f"{collection}::entries/{entry}"
        pd.testing.assert_frame_equal(
            _by_name(MetaCooler(uri)), _by_name(MetaCooler(clustered))
        )


def test_copy_clusters_needs_the_lengths_with_a_contig_map(clustered: Path):
    with (
        h5py.File(clustered, "r") as src,
        h5py.File(clustered.with_name("x.h5"), "w") as dst,
    ):
        assert not copy_clusters(dst, src)  # a map without clusters: nothing
        with pytest.raises(ValueError, match="needs the new contig lengths"):
            copy_clusters(src, dst, np.arange(4))


# ---------------------------------------------------------------------------
# dump
# ---------------------------------------------------------------------------


def test_dump_the_clusters_and_their_members(clustered: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["dump", str(clustered), "-t", "clusters", "-H"])
    assert result.exit_code == 0, result.output
    lines = result.output.strip().splitlines()
    assert lines[0].split("\t")[:6] == [
        "round",
        "name",
        "parent_round",
        "parent",
        "n_contigs",
        "length",
    ]
    assert len(lines) == 6

    result = runner.invoke(
        cli, ["dump", str(clustered), "-t", "members", "-H", "--contigs", "c3"]
    )
    assert result.exit_code == 0, result.output
    assert result.output.strip().splitlines()[1:] == [
        "overlapping\tO1\tc3",
        "recursive\tO1_2\tc3",
    ]

    result = runner.invoke(
        cli, ["dump", str(clustered), "-t", "clusters", "--contigs", "c4"]
    )
    assert result.exit_code == 0, result.output
    assert [line.split("\t")[1] for line in result.output.strip().splitlines()] == [
        "O2"
    ]


def test_dump_all_writes_the_clusters_of_a_map_that_has_them(
    clustered: Path, tmp_path: Path
):
    prefix = tmp_path / "dumped"
    result = CliRunner().invoke(
        cli, ["dump", str(clustered), "-t", "all", "--prefix", str(prefix)]
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "dumped.clusters.tsv").exists()
    assert (tmp_path / "dumped.members.tsv").exists()

    drop_clusters(clustered)
    other = tmp_path / "plain"
    result = CliRunner().invoke(
        cli, ["dump", str(clustered), "-t", "all", "--prefix", str(other)]
    )
    assert result.exit_code == 0, result.output
    assert not (tmp_path / "plain.clusters.tsv").exists()


@pytest.mark.parametrize("collection", [False, True])
def test_tree_shows_the_clusters(clustered: Path, tmp_path: Path, collection: bool):
    uri = clustered
    if collection:
        uri = tmp_path / "maps.mmetacool"
        create_mmetacool(uri, {"a": clustered})
    for target in [uri] + ([f"{uri}::entries/a"] if collection else []):
        result = CliRunner().invoke(cli, ["tree", str(target)])
        assert result.exit_code == 0, result.output
        assert "clusters" in result.output


def test_a_membership_given_twice_is_refused_wherever_it_comes_from():
    """The members are sorted and checked as they are stored, whoever hands
    them over: copy_clusters, renumbering the contigs of a regrouped map,
    as well as write_clusters, which names the cluster and contig first. A
    guard: no writer hands over a repeated pair today."""
    table = pd.DataFrame({"round": [0], "name": ["a"], "parent": [-1]})
    with pytest.raises(ValueError, match="1 membership\\(s\\) given more than once"):
        clusters_module._sorted(
            [{"name": "r"}],
            table,
            np.array([0, 0]),
            np.array([1, 1]),
            np.array([10, 20]),
        )
