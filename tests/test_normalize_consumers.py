"""Stored normalization weights: read by network/dump/show, dropped by edits."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from _metahic import MetaHiC, make_metahic
from numpy.testing import assert_allclose, assert_array_equal

from metacooler.api import MetaCooler
from metacooler.io import weights as W
from metacooler.io.annotate import annotate
from metacooler.io.clusters import write_clusters
from metacooler.io.dump import dump_table
from metacooler.io.network import network
from metacooler.io.reassign import reassign
from metacooler.ops import normalize as N
from metacooler.ops.aggregate import _aggregate_sparse
from metacooler.ops.show import _weigh_block, show


@pytest.fixture(scope="module")
def community(tmp_path_factory) -> MetaHiC:
    return make_metahic(tmp_path_factory.mktemp("community"), seed=11)


def _hits_from_pairs(community: MetaHiC, names: list[str]) -> np.ndarray:
    """Hi-C hits per contig as ``metator network`` counts them: one per read
    end, so a pair inside a contig counts twice for it."""
    pairs = pd.read_csv(community.pairs, sep="\t", comment="#", header=None)
    counts = pd.concat([pairs[1], pairs[3]]).value_counts()
    return counts.reindex(names, fill_value=0).to_numpy(dtype=float)


@pytest.fixture
def mapfile(tmp_path: Path, community: MetaHiC) -> Path:
    """The community map with four weights.

    ``weight`` is ``normalize ice`` (per genome) and ``normcc`` is
    ``normalize normcc``.  The other two are MetaTOR's ``empirical_hit`` and
    ``abundance`` normalizations, computed here, outside metacooler, and stored
    as any tool would store its own: ``1 / sqrt(hits)`` and
    ``1 / sqrt(coverage)`` for every bin of a contig.
    """
    path = tmp_path / "map.metacool"
    shutil.copy2(community.metacool, path)
    N.normalize_ice(path)
    N.normalize_normcc(path)
    mc = MetaCooler(path)
    contigs = mc.contigs
    hits = _hits_from_pairs(community, contigs["name"].tolist())
    contig_of_bin = mc._bin_contig_ids
    W.store_weights(
        path,
        "empirical_hit",
        (1 / np.sqrt(hits))[contig_of_bin],
        method="metator",
        depends_on=["pixels"],
        parameters={"mode": "empirical_hit"},
    )
    W.store_weights(
        path,
        "abundance",
        (1 / np.sqrt(contigs["cov"].to_numpy(float)))[contig_of_bin],
        method="metator",
        depends_on=["pixels", "contigs/cov"],
        parameters={"mode": "abundance"},
    )
    return path


def _read_tsv(path: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kw)


def _weighted_pixels(mc: MetaCooler, name: str) -> pd.DataFrame:
    px = mc.pixels()
    w = W.read_weights(mc, name)
    px["value"] = (
        px["count"].to_numpy()
        * w[px["bin1_id"].to_numpy()]
        * w[px["bin2_id"].to_numpy()]
    )
    return px


# ----------------------------------------------------------------------------
# network
# ----------------------------------------------------------------------------


def test_network_balances_bin_edges(tmp_path: Path, mapfile: Path):
    out = tmp_path / "edges.tsv"
    network(mapfile, level="bin", balance="weight", out=out, header=True)
    edges = _read_tsv(out)
    expected = _weighted_pixels(MetaCooler(mapfile), "weight")
    expected = expected[np.isfinite(expected["value"])]
    assert_array_equal(edges["node1"], expected["bin1_id"])
    assert_allclose(edges["value"], expected["value"])


def test_network_contig_edges_are_metators_network(
    tmp_path: Path, mapfile: Path, community: MetaHiC
):
    """A contig edge normalized by stored ``1 / sqrt(hits)`` weights is MetaTOR's
    ``empirical_hit`` edge: its contacts over ``sqrt(hits_i * hits_j)``."""
    out = tmp_path / "edges.tsv"
    network(mapfile, level="contig", balance="empirical_hit", out=out, header=True)
    edges = _read_tsv(out)
    names = MetaCooler(mapfile).contigs["name"].tolist()
    hits = pd.Series(_hits_from_pairs(community, names), index=names)
    raw = _read_tsv(_raw_contig_edges(tmp_path, mapfile))
    merged = edges.merge(raw, on=["node1", "node2"], suffixes=("", "_raw"))
    assert len(merged) == len(raw) == len(edges)
    hits_1 = hits[merged["node1"]].to_numpy()
    hits_2 = hits[merged["node2"]].to_numpy()
    assert_allclose(merged["value"], merged["value_raw"] / np.sqrt(hits_1 * hits_2))


def _raw_contig_edges(tmp_path: Path, mapfile: Path) -> Path:
    out = tmp_path / "raw.tsv"
    network(mapfile, level="contig", out=out, header=True, force=True)
    return out


def test_network_genome_edges_sum_the_normalized_contacts(
    tmp_path: Path, mapfile: Path
):
    out = tmp_path / "edges.tsv"
    network(mapfile, level="genome", balance="abundance", out=out, header=True)
    edges = _read_tsv(out)
    mc = MetaCooler(mapfile)
    px = _weighted_pixels(mc, "abundance")
    genome_of_bin = mc.contigs["genome_id"].to_numpy()[mc._bin_contig_ids]
    g1 = genome_of_bin[px["bin1_id"].to_numpy()]
    g2 = genome_of_bin[px["bin2_id"].to_numpy()]
    lo, hi = np.minimum(g1, g2), np.maximum(g1, g2)
    names = mc.genomes["name"].to_numpy()
    expected = (
        pd.DataFrame({"node1": names[lo], "node2": names[hi], "value": px["value"]})
        .groupby(["node1", "node2"], as_index=False)["value"]
        .sum()
    )
    merged = edges.merge(expected, on=["node1", "node2"], suffixes=("", "_ref"))
    assert len(merged) == len(edges) == len(expected)
    assert_allclose(merged["value"], merged["value_ref"])


def test_network_drops_spurious_contig_pairs(tmp_path: Path, mapfile: Path):
    mc = MetaCooler(mapfile)
    found = W.spurious_cutoff(mc, "normcc")
    assert found is not None
    cutoff, rule = found
    assert rule == "<"
    kept = tmp_path / "kept.tsv"
    every = tmp_path / "every.tsv"
    network(mapfile, balance="normcc", drop_spurious=True, out=kept, header=True)
    network(mapfile, balance="normcc", out=every, header=True)
    kept_edges, all_edges = _read_tsv(kept), _read_tsv(every)
    inter = all_edges["node1"] != all_edges["node2"]
    spurious = inter & (all_edges["value"] < cutoff)
    assert spurious.any()
    assert len(kept_edges) == len(all_edges) - int(spurious.sum())
    kept_inter = kept_edges[kept_edges["node1"] != kept_edges["node2"]]
    assert (kept_inter["value"] >= cutoff).all()

    genome_level = tmp_path / "genomes.tsv"
    network(
        mapfile,
        level="genome",
        balance="normcc",
        drop_spurious=True,
        out=genome_level,
        header=True,
    )
    # Genome pairs are ordered by genome id, as network writes them.
    genome_id = dict(zip(mc.contigs["name"], mc.contigs["genome_id"], strict=True))
    names = mc.genomes["name"].to_numpy()
    g1 = kept_edges["node1"].map(genome_id).to_numpy()
    g2 = kept_edges["node2"].map(genome_id).to_numpy()
    expected = (
        pd.DataFrame(
            {
                "node1": names[np.minimum(g1, g2)],
                "node2": names[np.maximum(g1, g2)],
                "value": kept_edges["value"],
            }
        )
        .groupby(["node1", "node2"], as_index=False)["value"]
        .sum()
    )
    merged = _read_tsv(genome_level).merge(expected, on=["node1", "node2"])
    assert len(merged) == len(expected)
    assert_allclose(merged["value_x"], merged["value_y"])


def test_network_spurious_rule_at_the_cutoff(tmp_path: Path, mapfile: Path):
    n = len(MetaCooler(mapfile).bins())
    W.store_weights(
        mapfile,
        "flat",
        np.ones(n),
        method="x",
        depends_on=["pixels"],
        parameters={"spurious_cutoff": 1.0, "spurious_rule": "<="},
    )
    out = tmp_path / "e.tsv"
    network(mapfile, balance="flat", drop_spurious=True, out=out, header=True)
    edges = _read_tsv(out)
    inter = edges[edges["node1"] != edges["node2"]]
    assert (inter["value"] > 1.0).all()


def test_network_with_every_bin_masked(tmp_path: Path, mapfile: Path):
    n = len(MetaCooler(mapfile).bins())
    W.store_weights(
        mapfile,
        "masked",
        np.full(n, np.nan),
        method="x",
        depends_on=["pixels"],
        parameters={"spurious_cutoff": 1.0},
    )
    for level in ("bin", "contig", "genome"):
        out = tmp_path / f"{level}.tsv"
        network(mapfile, level=level, balance="masked", out=out)
        assert out.read_text() == ""
    out = tmp_path / "dropped.tsv"
    network(mapfile, level="genome", balance="masked", drop_spurious=True, out=out)
    assert out.read_text() == ""


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"drop_spurious": True}, ValueError, "pass --balance"),
        (
            {"balance": "normcc", "drop_spurious": True, "level": "bin"},
            ValueError,
            "contig pairs",
        ),
        (
            {"balance": "weight", "drop_spurious": True},
            ValueError,
            "no spurious-contact cutoff",
        ),
        ({"balance": "nope"}, KeyError, "No weight column"),
    ],
)
def test_network_refuses(mapfile: Path, kwargs, error, match):
    with pytest.raises(error, match=match):
        network(mapfile, **kwargs)


# ----------------------------------------------------------------------------
# dump
# ----------------------------------------------------------------------------


def test_dump_long_adds_a_balanced_column(tmp_path: Path, mapfile: Path):
    out = tmp_path / "px.tsv"
    dump_table(mapfile, "pixels", out=out, header=True, balance="weight")
    dumped = _read_tsv(out)
    expected = _weighted_pixels(MetaCooler(mapfile), "weight")
    assert list(dumped.columns) == ["bin1_id", "bin2_id", "count", "balanced"]
    assert_allclose(dumped["balanced"], expected["value"], equal_nan=True)
    assert dumped["balanced"].isna().any()  # masked bins come out empty

    joined = tmp_path / "joined.tsv"
    dump_table(mapfile, "pixels", out=joined, header=True, balance="weight", join=True)
    assert _read_tsv(joined).columns[-1] == "balanced"


def test_dump_wide_is_the_balanced_matrix(tmp_path: Path, mapfile: Path):
    out = tmp_path / "wide.tsv"
    dump_table(
        mapfile,
        "pixels",
        out=out,
        pixels_format="wide",
        balance="weight",
        genomes=["G1"],
    )
    mat = _read_tsv(out, header=None).to_numpy()
    mc = MetaCooler(mapfile)
    lo, hi = mc.genome_extent("G1")
    w = W.read_weights(mc, "weight")[lo:hi]
    raw = mc.matrix(balance=False)[lo:hi, lo:hi].astype(float)
    expected = raw * w[:, None] * w[None, :]
    masked = np.isnan(w)
    expected[masked, :] = np.nan
    expected[:, masked] = np.nan
    assert_allclose(mat, expected, equal_nan=True)


def test_dump_balance_only_applies_to_pixels(tmp_path: Path, mapfile: Path):
    with pytest.raises(ValueError, match="only applies to the pixels table"):
        dump_table(mapfile, "bins", balance="weight")
    dump_table(
        mapfile, "all", prefix=str(tmp_path / "all"), header=True, balance="weight"
    )
    assert "balanced" in _read_tsv(tmp_path / "all.pixels.tsv").columns
    assert "balanced" not in _read_tsv(tmp_path / "all.bins.tsv").columns


def test_dump_header_of_an_empty_selection_names_the_balanced_column(
    tmp_path: Path, community: MetaHiC
):
    empty = tmp_path / "empty.pairs"
    empty.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    from metacooler.io.create import from_pairs

    path = tmp_path / "empty.metacool"
    from_pairs(
        uri=path,
        bins=None,
        binsize=1000,
        pairs=empty,
        contig_genome_map=community.contigs,
        fasta=community.fasta,
        qc=False,
        ps=False,
    )
    n = len(MetaCooler(path).bins())
    W.store_weights(path, "w", np.ones(n), method="x", depends_on=["pixels"])
    out = tmp_path / "px.tsv"
    dump_table(path, "pixels", out=out, header=True, balance="w")
    assert out.read_text().strip().split("\t") == [
        "bin1_id",
        "bin2_id",
        "count",
        "balanced",
    ]


# ----------------------------------------------------------------------------
# show
# ----------------------------------------------------------------------------


def test_weigh_block_multiplies_both_axes():
    mat = np.arange(1, 7, dtype=float).reshape(2, 3)
    weights = np.array([1.0, 2.0, 3.0, np.nan, 10.0])
    out = _weigh_block(mat, [(0, 1), (4, 5)], [(1, 3), (3, 4)], weights)
    expected = (
        mat * np.array([1.0, 10.0])[:, None] * np.array([2.0, 3.0, np.nan])[None, :]
    )
    assert_allclose(out, expected, equal_nan=True)
    assert _weigh_block(mat, [(0, 2)], [(0, 3)], None) is mat


def test_aggregate_sparse_sums_normalized_pixels(mapfile: Path):
    mc = MetaCooler(mapfile)
    w = W.read_weights(mc, "empirical_hit")
    agg = _aggregate_sparse(mc, level="genome", weights=w).toarray()
    raw = _aggregate_sparse(mc, level="genome").toarray()
    assert agg.dtype == np.float64
    assert not np.allclose(agg, raw)
    px = _weighted_pixels(mc, "empirical_hit")
    assert np.trace(agg) + np.triu(agg, 1).sum() == pytest.approx(px["value"].sum())


@pytest.mark.parametrize(
    "kwargs",
    [
        {"level": "bin", "genomes": ["G0", "G2"]},
        {"level": "bin", "range1": "c0", "range2": "c4"},
        {"level": "contig"},
        {"level": "genome"},
    ],
)
def test_show_renders_a_normalized_map(
    tmp_path: Path, mapfile: Path, kwargs, monkeypatch
):
    titles = []
    import matplotlib.axes

    original = matplotlib.axes.Axes.set_title

    def spy(self, label, *args, **kw):
        titles.append(label)
        return original(self, label, *args, **kw)

    monkeypatch.setattr(matplotlib.axes.Axes, "set_title", spy)
    out = tmp_path / "map.png"
    show(mapfile, out=out, dpi=50, balance="weight", transform="none", **kwargs)
    assert out.stat().st_size > 0
    assert any("normalized: weight" in t for t in titles)


# ----------------------------------------------------------------------------
# Edits that invalidate weights
# ----------------------------------------------------------------------------


def _table(path: Path, rows: dict[str, str]) -> Path:
    path.write_text(
        "contig_name\tgenome_name\n" + "".join(f"{c}\t{g}\n" for c, g in rows.items())
    )
    return path


def _bins_columns(path: Path) -> set[str]:
    with h5py.File(path) as f:
        return set(f["bins"])


def test_reassign_relabel_drops_what_the_assignment_defined(
    tmp_path: Path, mapfile: Path, caplog
):
    out = tmp_path / "moved.metacool"
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        reassign(mapfile, _table(tmp_path / "t.tsv", {"c0": "G1"}), out, layout="keep")
    columns = _bins_columns(out)
    assert "weight" not in columns  # balanced per genome
    assert {"empirical_hit", "abundance", "normcc"} <= columns
    assert "Dropped normalization weight(s) 'weight'" in caplog.text
    with h5py.File(out) as f:
        assert f["bins/abundance"].attrs["normalization"] == "metator"
        assert f["bins/normcc"].attrs["normalization"] == "normcc"


def test_reassign_regroup_carries_weights_with_their_attributes(
    tmp_path: Path, mapfile: Path
):
    source = MetaCooler(mapfile)
    before = W.read_weights(source, "empirical_hit")
    out = tmp_path / "regrouped.metacool"
    reassign(mapfile, _table(tmp_path / "t.tsv", {"c0": "G2"}), out)
    result = MetaCooler(out)
    assert result.info.get("contig-ordering") != "custom"
    columns = _bins_columns(out)
    assert "weight" not in columns
    with h5py.File(out) as f:
        attrs = dict(f["bins/empirical_hit"].attrs)
    assert attrs["normalization"] == "metator"
    assert attrs["mode"] == "empirical_hit"
    assert list(attrs["depends_on"]) == ["pixels"]
    # The weights moved with their bins: each bin keeps its contig's weight.
    after = W.read_weights(result, "empirical_hit")
    old = dict(zip(source.contigs["name"], before[_first_bins(source)], strict=True))
    new = dict(zip(result.contigs["name"], after[_first_bins(result)], strict=True))
    assert new == pytest.approx(old)


def _first_bins(mc: MetaCooler) -> np.ndarray:
    """The first bin of every contig, in ``/contigs`` order."""
    with mc.open() as grp:
        return grp["indexes/chrom_offset"][:-1]


def test_reassign_that_moves_nothing_keeps_every_weight(tmp_path: Path, mapfile: Path):
    out = tmp_path / "same.metacool"
    reassign(mapfile, _table(tmp_path / "t.tsv", {"c0": "G0"}), out)
    assert {"weight", "empirical_hit", "abundance"} <= _bins_columns(out)


def test_reassign_in_place_drops_what_the_assignment_defined(
    tmp_path: Path, mapfile: Path
):
    reassign(
        mapfile, _table(tmp_path / "t.tsv", {"c0": "G1"}), layout="keep", inplace=True
    )
    columns = _bins_columns(mapfile)
    assert "weight" not in columns
    assert {"empirical_hit", "abundance"} <= columns


def test_reassign_to_a_round_drops_a_genome_scoped_ice_weight(
    tmp_path: Path, mapfile: Path, caplog
):
    """A clustering round is a new assignment like any table: the per-genome
    ICE weight goes, the weights that never read the assignment stay."""
    rounds = {"A": ["c0", "c3", "c4", "c5"], "B": ["c1", "c2"], "C": ["c6", "c7", "c8"]}
    write_clusters(
        mapfile,
        ["r"],
        pd.DataFrame({"round": "r", "name": list(rounds)}),
        pd.DataFrame(
            [("r", name, c) for name, contigs in rounds.items() for c in contigs],
            columns=["round", "cluster", "contig"],
        ),
    )
    out = tmp_path / "round.metacool"
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        reassign(mapfile, output=out, round="r")
    columns = _bins_columns(out)
    assert "weight" not in columns
    assert {"empirical_hit", "abundance"} <= columns
    assert "Dropped normalization weight(s) 'weight'" in caplog.text
    assert set(MetaCooler(out).genomes["name"]) >= {"A", "B", "C"}


def test_annotate_drops_weights_that_read_a_changed_column(
    tmp_path: Path, mapfile: Path, caplog
):
    table = tmp_path / "cov.tsv"
    names = MetaCooler(mapfile).contigs["name"].tolist()
    table.write_text(
        "contig_name\tcov\n" + "".join(f"{n}\t{i + 1}\n" for i, n in enumerate(names))
    )
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        annotate(mapfile, contigs_table=table, replace=True, inplace=True)
    columns = _bins_columns(mapfile)
    assert "abundance" not in columns  # stored with depends_on contigs/cov
    assert {"weight", "empirical_hit", "normcc"} <= columns
    assert "Dropped normalization weight(s) 'abundance'" in caplog.text


def test_annotate_of_another_column_keeps_every_weight(tmp_path: Path, mapfile: Path):
    table = tmp_path / "gc.tsv"
    names = MetaCooler(mapfile).contigs["name"].tolist()
    table.write_text("contig_name\tgc\n" + "".join(f"{n}\t0.5\n" for n in names))
    out = tmp_path / "annotated.metacool"
    annotate(mapfile, out, contigs_table=table)
    assert {"weight", "empirical_hit", "abundance"} <= _bins_columns(out)


def test_describe_weights_lists_every_method(mapfile: Path):
    rows = {r["name"]: r["method"] for r in W.describe_weights(MetaCooler(mapfile))}
    assert rows == {
        "weight": "ice",
        "normcc": "normcc",
        "empirical_hit": "metator",
        "abundance": "metator",
    }
