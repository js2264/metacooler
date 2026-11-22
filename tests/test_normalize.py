"""Tests of ``metacooler normalize``: ICE, and the weight columns any tool stores."""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

import cooler
import h5py
import numpy as np
import pytest
from _metahic import MetaHiC, make_metahic
from numpy.testing import assert_allclose, assert_array_equal

import metacooler
from metacooler import normalize as engine
from metacooler.api import MetaCooler
from metacooler.io import weights as W
from metacooler.io.create import from_pairs
from metacooler.ops import normalize as N
from metacooler.schema import validate

# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def community(tmp_path_factory) -> MetaHiC:
    """Three genomes of three contigs and one unbinned contig, 1 kb bins."""
    return make_metahic(tmp_path_factory.mktemp("community"))


@pytest.fixture
def mapfile(tmp_path: Path, community: MetaHiC) -> Path:
    """A fresh, writable copy of the community map."""
    path = tmp_path / "map.metacool"
    shutil.copy2(community.metacool, path)
    return path


def _rebuild(tmp_path: Path, community: MetaHiC, genome_of, name: str) -> Path:
    """The community's pairs binned again under another contig assignment."""
    names = [
        line[1:].strip()
        for line in community.fasta.read_text().splitlines()
        if line.startswith(">")
    ]
    table = tmp_path / f"{name}.tsv"
    table.write_text(
        "contig_name\tgenome_name\n" + "".join(f"{n}\t{genome_of(n)}\n" for n in names)
    )
    out = tmp_path / f"{name}.metacool"
    from_pairs(
        uri=out,
        bins=None,
        binsize=1000,
        pairs=community.pairs,
        contig_genome_map=table,
        fasta=community.fasta,
        qc=False,
        ps=False,
    )
    return out


def _n_bins(uri: Path | str) -> int:
    return len(MetaCooler(uri).bins())


def _store(uri: Path | str, name: str, values, **kwargs) -> None:
    """``store_weights`` with the two required keywords filled in."""
    kwargs.setdefault("method", "test")
    kwargs.setdefault("depends_on", ["pixels"])
    W.store_weights(uri, name, values, **kwargs)


# ----------------------------------------------------------------------------
# Engine: matrix balancing
# ----------------------------------------------------------------------------


def test_ice_marginals_count_the_diagonal_once():
    marg = engine.ice_marginals(
        np.array([0, 0, 1]), np.array([0, 1, 1]), np.array([5.0, 2.0, 3.0]), 3
    )
    assert_array_equal(marg, [7.0, 5.0, 0.0])


def test_mad_filter_compares_a_bin_with_its_own_genome():
    # Genome 1 is a hundred times rarer than genome 0: compared with the whole
    # map its bins all look like outliers, compared with its own median none
    # does -- except the one far below the rest of its genome.
    rng = np.random.default_rng(0)
    abundant = 100 * (1 + 0.05 * rng.standard_normal(20))
    marg = np.concatenate([abundant, [1.0, 1.1, 0.9, 1.05, 1e-6]])
    groups = np.repeat([0, 1], [20, 5])
    dropped = engine.ice_mad_drop(marg, groups, mad_max=5)
    assert_array_equal(dropped, [False] * 24 + [True])
    one_group = engine.ice_mad_drop(marg, np.zeros(25, dtype=int), mad_max=5)
    assert one_group[20:].all()
    assert not one_group[:20].any()


def test_mad_filter_of_an_empty_map_drops_nothing():
    assert not engine.ice_mad_drop(np.zeros(3), np.zeros(3, dtype=int), 5).any()


def test_ice_iterate_handles_an_empty_map():
    stats = engine.ice_iterate(
        np.ones(0),
        np.zeros(0, dtype=np.int64),
        lambda b: np.zeros(0),
        tol=1e-5,
        max_iters=10,
    )
    assert stats.scale.size == 0


def test_compact_groups():
    assert_array_equal(engine.compact_groups(np.array([7, 3, 7, 10])), [1, 0, 1, 2])


# ----------------------------------------------------------------------------
# Weight storage: what any tool stores
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "match"),
    [
        ("", "non-empty"),
        ("a/b", "nested group"),
        ("chrom", "coordinate"),
        ("VC_SQRT", "divisive"),
    ],
)
def test_weight_names_that_cannot_be_stored(name, match):
    with pytest.raises(ValueError, match=match):
        W.check_weight_name(name)


def test_store_and_read_weights(mapfile: Path):
    n = _n_bins(mapfile)
    values = np.linspace(0.5, 2.0, n)
    values[3] = np.nan
    W.store_weights(
        mapfile,
        "w",
        values,
        method="test",
        depends_on=["pixels", "genome_id"],
        parameters={"mode": "a", "skip": None},
        commands=["sometool normalize map.metacool"],
    )
    mc = MetaCooler(mapfile)
    assert_array_equal(W.read_weights(mc, "w"), values)
    assert_array_equal(W.read_weights(mapfile, "w"), values)  # a URI too
    with h5py.File(mapfile) as f:
        attrs = dict(f["bins/w"].attrs)
        assert W.depends_on(f["bins/w"]) == {"pixels", "genome_id"}
    assert attrs["normalization"] == "test"
    assert attrs["mode"] == "a"
    assert not attrs["divisive_weights"]
    assert attrs["metacooler_version"] == metacooler.__version__
    assert "skip" not in attrs
    assert mc.metadata["commands"][-1] == "sometool normalize map.metacool"


def test_store_weights_keeps_parameters_of_every_kind(mapfile: Path):
    _store(
        mapfile,
        "w",
        np.ones(_n_bins(mapfile)),
        parameters={
            "tool_version": "1.3.4",
            "inputs": ["a.bam", "b.bam"],
            "cutoffs": [0.5, 2.0],
            "config": {"k": 1},
            "fasta": Path("assembly.fa"),
            "converged": True,
        },
    )
    attrs = W.describe_weights(mapfile)[0]["attrs"]
    assert isinstance(attrs, dict)
    assert attrs["tool_version"] == "1.3.4"
    assert attrs["inputs"] == ["a.bam", "b.bam"]
    assert_array_equal(attrs["cutoffs"], [0.5, 2.0])
    assert json.loads(str(attrs["config"])) == {"k": 1}
    assert attrs["fasta"] == "assembly.fa"
    assert attrs["converged"]


def test_store_weights_into_an_mmetacool_entry(tiny_mmetacool: Path):
    uri = f"{tiny_mmetacool}::entries/entry2"
    n = _n_bins(uri)
    _store(uri, "w", np.full(n, 3.0), depends_on=["contigs/length"])
    assert_array_equal(W.read_weights(uri, "w"), np.full(n, 3.0))
    assert [row["name"] for row in W.describe_weights(uri)] == ["w"]
    assert W.describe_weights(f"{tiny_mmetacool}::entries/entry1") == []


def test_store_weights_refuses_a_taken_name_and_a_wrong_length(mapfile: Path):
    n = _n_bins(mapfile)
    _store(mapfile, "w", np.ones(n))
    with pytest.raises(ValueError, match="already exists"):
        _store(mapfile, "w", np.ones(n))
    _store(mapfile, "w", np.full(n, 2.0), force=True)
    assert_array_equal(W.read_weights(mapfile, "w"), np.full(n, 2.0))
    with pytest.raises(ValueError, match="would not line up"):
        _store(mapfile, "w2", np.ones(n + 1))
    with pytest.raises(ValueError, match="divisive"):
        _store(mapfile, "KR", np.ones(n))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"method": ""}, "method must name"),
        ({"parameters": {"depends_on": "x"}}, "'depends_on' are attributes"),
        (
            {"parameters": {"normalization": "x", "divisive_weights": True}},
            "'divisive_weights', 'normalization' are attributes",
        ),
        ({"depends_on": []}, "must name what the weights were computed from"),
        ({"depends_on": ["pixel"]}, "got 'pixel'"),
        (
            {"depends_on": ["contigs/coverage"]},
            "Stored /contigs columns: name, length, cov, RS, qc_",
        ),
        ({"depends_on": ["contigs/genome_id"]}, "not 'contigs/genome_id'"),
    ],
)
def test_store_weights_refuses(mapfile: Path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        _store(mapfile, "w", np.ones(_n_bins(mapfile)), **kwargs)
    with h5py.File(mapfile) as f:
        assert "w" not in f["bins"]


def test_check_depends_on():
    columns = ["name", "length", "genome_id", "cov"]
    # Repeats dropped, order kept; the pixels are not required.
    assert W.check_depends_on(["contigs/cov", "pixels", "contigs/cov"], columns) == [
        "contigs/cov",
        "pixels",
    ]
    assert W.check_depends_on(["contigs/cov"], columns) == ["contigs/cov"]
    # One entry given as a string is one entry, not its characters.
    assert W.check_depends_on("genome_id", columns) == ["genome_id"]
    with pytest.raises(ValueError, match="Stored /contigs columns: none"):
        W.check_depends_on(["contigs/cov"], [])


def test_read_weights_inverts_a_divisive_column(mapfile: Path):
    n = _n_bins(mapfile)
    with h5py.File(mapfile, "r+") as f:
        f["bins"].create_dataset("KR", data=np.full(n, 4.0))
        ds = f["bins"].create_dataset("vc", data=np.full(n, 2.0))
        ds.attrs["divisive_weights"] = True
    mc = MetaCooler(mapfile)
    assert_array_equal(W.read_weights(mc, "KR"), np.full(n, 0.25))
    assert_array_equal(W.read_weights(mc, "vc"), np.full(n, 0.5))


def test_read_weights_names_the_available_columns(mapfile: Path):
    _store(mapfile, "w", np.ones(_n_bins(mapfile)))
    with pytest.raises(KeyError, match="weight columns: w"):
        W.read_weights(mapfile, "nope")
    with pytest.raises(KeyError, match="store_weights"):
        W.read_weights(MetaCooler(mapfile), "start")


def test_read_weights_refuses_a_column_of_another_shape(mapfile: Path):
    n = _n_bins(mapfile)
    with h5py.File(mapfile, "r+") as f:
        f["bins"].create_dataset("bad", data=np.array([b"x"] * n))
    with pytest.raises(ValueError, match="one number per bin"):
        W.read_weights(MetaCooler(mapfile), "bad")


def test_stale_weights(mapfile: Path):
    n = _n_bins(mapfile)
    _store(mapfile, "by_genome", np.ones(n), depends_on=["pixels", "genome_id"])
    _store(mapfile, "by_cov", np.ones(n), depends_on=["pixels", "contigs/cov"])
    _store(mapfile, "single", np.ones(n), depends_on="pixels")
    with h5py.File(mapfile, "r+") as f:
        f["bins"].create_dataset("foreign", data=np.ones(n))  # cooler balance
        bins = f["bins"]
        assert W.stale_weights(bins, assignment=True) == ["by_genome"]
        assert W.stale_weights(bins, contig_columns=["cov", "x"]) == ["by_cov"]
        assert W.stale_weights(bins) == []
        assert W.depends_on(bins["single"]) == {"pixels"}
        assert W.depends_on(bins["foreign"]) == {"pixels"}
        assert W.drop_stale_weights(bins, assignment=True, why="it") == ["by_genome"]
        assert "by_genome" not in bins


def test_depends_on_of_a_column_another_writer_stored(mapfile: Path):
    """A single entry, stored as a fixed-length byte string, is one entry."""
    n = _n_bins(mapfile)
    with h5py.File(mapfile, "r+") as f:
        ds = f["bins"].create_dataset("other", data=np.ones(n))
        ds.attrs["depends_on"] = np.bytes_(b"contigs/cov")
        assert W.depends_on(ds) == {"contigs/cov"}
        assert W.stale_weights(f["bins"], contig_columns=["cov"]) == ["other"]


def test_restore_weight_attrs_skips_a_column_that_was_not_written(mapfile: Path):
    _store(mapfile, "w", np.ones(_n_bins(mapfile)))
    with h5py.File(mapfile, "r+") as f:
        W.restore_weight_attrs(f["bins"], {"gone": {"a": 1}, "w": {"a": 2}})
        assert f["bins/w"].attrs["a"] == 2
        assert "gone" not in f["bins"]


def test_describe_weights(mapfile: Path):
    n = _n_bins(mapfile)
    values = np.ones(n)
    values[:2] = np.nan
    _store(mapfile, "w", values, method="x")
    with h5py.File(mapfile, "r+") as f:
        f["bins"].create_dataset("foreign", data=np.ones(n))
    rows = {row["name"]: row for row in W.describe_weights(MetaCooler(mapfile))}
    assert rows["w"]["method"] == "x"
    assert rows["w"]["n_masked"] == 2
    assert rows["foreign"]["method"] is None


def test_the_package_exports_the_weight_api():
    assert metacooler.store_weights is W.store_weights
    assert metacooler.read_weights is W.read_weights


# ----------------------------------------------------------------------------
# ICE
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "filters",
    [{"min_nnz": 2}, {"min_nnz": 0, "min_count": 50}],
    ids=["min_nnz", "min_count"],
)
def test_ice_by_genome_is_cooler_cis_only_when_every_contig_is_a_genome(
    tmp_path: Path, community: MetaHiC, filters
):
    # Contigs of 2 to 8 bins: at cooler's default min_nnz of 10, every bin of
    # a contig balanced on its own is masked, on both sides, and the
    # comparison would hold vacuously.
    path = _rebuild(tmp_path, community, lambda n: f"G_{n}", "per_contig")
    ours = N.normalize_ice(path, store=False, tol=1e-12, max_iters=2000, **filters)
    ref, _ = cooler.balance_cooler(
        cooler.Cooler(str(path)), cis_only=True, tol=1e-12, max_iters=2000, **filters
    )
    assert ours.attrs["converged"]
    assert np.isfinite(ours.weights).sum() >= len(ours.weights) // 3
    assert_array_equal(np.isnan(ours.weights), np.isnan(ref))
    assert_allclose(ours.weights, ref, rtol=1e-9)


def test_ice_by_genome_is_cooler_global_when_there_is_one_genome(
    tmp_path: Path, community: MetaHiC
):
    path = _rebuild(tmp_path, community, lambda n: "G", "one_genome")
    # cooler's MAD filter takes its medians per contig, ours per genome, so it
    # is off.  And one diagonal: cooler before 0.12 also dropped the pixel
    # between the last bin of a contig and the first of the next, which is
    # |bin2 - bin1| = 1 but not a diagonal of either (cooler #211); with a
    # single diagonal ignored, every cooler version drops the same pixels.
    ours = N.normalize_ice(
        path, store=False, mad_max=0, ignore_diags=1, tol=1e-12, max_iters=2000
    )
    ref, _ = cooler.balance_cooler(
        cooler.Cooler(str(path)), mad_max=0, ignore_diags=1, tol=1e-12, max_iters=2000
    )
    assert ours.attrs["converged"]
    assert np.isfinite(ours.weights).sum() >= len(ours.weights) // 2
    assert_allclose(ours.weights, ref, rtol=1e-9)


def test_ice_by_genome_rereads_a_map_too_large_to_cache(mapfile: Path, monkeypatch):
    """Past the cache, every iteration re-reads the pixels: same weights."""
    cached = N.normalize_ice(mapfile, store=False, chunksize=97)
    monkeypatch.setattr(N, "_ICE_CACHE_PIXELS", 50)
    reread = N.normalize_ice(mapfile, store=False, chunksize=97)
    assert np.isfinite(cached.weights).any()
    assert_array_equal(reread.weights, cached.weights)


def test_ice_ignore_diags_stops_at_a_contig_boundary(community: MetaHiC):
    # As cooler does since 0.12: two adjacent bins of two contigs are not on
    # a diagonal of the matrix of either contig, so their pixel is kept.
    mc = MetaCooler(community.metacool)
    contig = mc._bin_contig_ids
    groups = N._contig_groups(mc)[contig]
    kept = N._GroupPixels(mc, groups, contig, ignore_diags=2, chunksize=97)
    parts = list(kept._chunks())
    b1 = np.concatenate([part[0] for part in parts])
    b2 = np.concatenate([part[1] for part in parts])
    px = mc.pixels()
    p1, p2 = px["bin1_id"].to_numpy(), px["bin2_id"].to_numpy()
    same_group = groups[p1] == groups[p2]
    near = (np.abs(p2 - p1) < 2) & (contig[p1] == contig[p2])
    expected = same_group & ~near
    assert_array_equal(b1, p1[expected])
    assert_array_equal(b2, p2[expected])
    across = (np.abs(b2 - b1) == 1) & (contig[b1] != contig[b2])
    assert across.any()


def test_ice_by_genome_ignores_contacts_between_genomes(
    tmp_path: Path, community: MetaHiC
):
    noisy = tmp_path / "noisy.pairs"
    lines = community.pairs.read_text().splitlines()
    extra = [
        f"x{i}\tc0\t{100 + i}\tc4\t{200 + i}\t+\t+" for i in range(300)
    ]  # G0 <-> G1
    noisy.write_text("\n".join(lines + extra) + "\n")
    out = tmp_path / "noisy.metacool"
    from_pairs(
        uri=out,
        bins=None,
        binsize=1000,
        pairs=noisy,
        contig_genome_map=community.contigs,
        fasta=community.fasta,
        qc=False,
        ps=False,
    )
    clean = N.normalize_ice(community.metacool, store=False).weights
    assert_array_equal(N.normalize_ice(out, store=False).weights, clean)
    # cooler's whole-matrix balancing does see them.
    assert not np.allclose(
        N.normalize_ice(out, scope="global", store=False).weights,
        N.normalize_ice(community.metacool, scope="global", store=False).weights,
        equal_nan=True,
    )


def test_ice_balances_each_genome_block(mapfile: Path):
    result = N.normalize_ice(mapfile)
    mc = MetaCooler(mapfile)
    assert result.attrs["converged"]
    assert result.attrs["depends_on"] == ["pixels", "genome_id"]
    w = W.read_weights(mc, "weight")
    dense = mc.matrix(balance=False)[:, :].astype(float)
    contig = mc._bin_contig_ids
    near = np.abs(np.subtract.outer(np.arange(len(w)), np.arange(len(w)))) < 2
    dense[near & (contig[:, None] == contig[None, :])] = 0
    balanced = dense * w[:, None] * w[None, :]
    groups = N._contig_groups(mc)[contig]
    for g in np.unique(groups):
        rows = np.flatnonzero((groups == g) & np.isfinite(w))
        if rows.size:
            sums = np.nansum(balanced[np.ix_(rows, rows)], axis=1)
            assert_allclose(sums, 1.0, rtol=1e-2)


@pytest.mark.parametrize("scope", ["contig", "global"])
def test_ice_with_coolers_scopes(mapfile: Path, scope: str):
    result = N.normalize_ice(mapfile, scope=scope)
    ref, _ = cooler.balance_cooler(
        cooler.Cooler(str(mapfile)), cis_only=scope == "contig"
    )
    assert_allclose(result.weights, ref, rtol=1e-12)
    assert result.attrs["depends_on"] == ["pixels"]
    with h5py.File(mapfile) as f:
        attrs = f["bins/weight"].attrs
        # Per-contig statistics are summarised to scalars.
        assert np.ndim(attrs["converged"]) == 0
        assert np.ndim(attrs["var"]) == 0


def test_summarize_stats():
    out = N._summarize_stats(
        {
            "converged": np.array([True, False, True]),
            "var": np.array([1e-6, np.nan, 3e-6]),
            "scale": np.array([1.0, 2.0, 3.0]),
            "tol": np.float64(1e-5),
        }
    )
    assert out == {
        "converged": False,
        "n_unconverged": 1,
        "n_groups": 3,
        "var": 3e-6,
        "tol": 1e-5,
    }
    assert np.isnan(N._summarize_stats({"var": np.array([np.nan])})["var"])


def test_ice_warns_when_it_does_not_converge(mapfile: Path, caplog):
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        result = N.normalize_ice(mapfile, max_iters=1, store=False)
    assert not result.attrs["converged"]
    assert "without convergence" in caplog.text


def test_ice_passes_coolers_warnings_on(mapfile: Path, caplog):
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        result = N.normalize_ice(mapfile, scope="global", max_iters=1, store=False)
    assert not result.attrs["converged"]
    assert "cooler: " in caplog.text


def test_ice_on_a_map_without_pixels(tmp_path: Path, community: MetaHiC, caplog):
    empty = tmp_path / "empty.pairs"
    empty.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    out = tmp_path / "empty.metacool"
    from_pairs(
        uri=out,
        bins=None,
        binsize=1000,
        pairs=empty,
        contig_genome_map=community.contigs,
        fasta=community.fasta,
        qc=False,
        ps=False,
    )
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        for scope in ("genome", "global"):
            result = N.normalize_ice(out, scope=scope, name=f"w_{scope}")
            assert np.isnan(result.weights).all()
    assert "Every bin is masked" in caplog.text


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"scope": "chromosome"}, "scope must be"),
        ({"mad_max": -1}, "mad_max must be"),
        ({"tol": 0}, "tol must be"),
    ],
)
def test_ice_refuses_bad_parameters(mapfile: Path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        N.normalize_ice(mapfile, **kwargs)


def test_ice_stays_off_a_name_that_is_taken(mapfile: Path):
    N.normalize_ice(mapfile)
    with pytest.raises(ValueError, match="already exists"):
        N.normalize_ice(mapfile)
    with pytest.raises(ValueError, match="divisive"):
        N.normalize_ice(mapfile, name="KR")
    N.normalize_ice(mapfile, scope="global", force=True)
    with h5py.File(mapfile) as f:
        assert f["bins/weight"].attrs["scope"] == "global"


def test_stored_weights_are_read_by_cooler(mapfile: Path):
    N.normalize_ice(mapfile, commands=["metacooler normalize ice map"])
    validate(mapfile)
    clr = cooler.Cooler(str(mapfile))
    w = clr.bins()["weight"][:].to_numpy()
    raw = clr.matrix(balance=False)[:20, :20].astype(float)
    assert_allclose(
        clr.matrix(balance="weight")[:20, :20],
        raw * w[:20, None] * w[None, :20],
        equal_nan=True,
    )
    assert (
        MetaCooler(mapfile).metadata["commands"][-1] == "metacooler normalize ice map"
    )


def test_store_false_writes_nothing(mapfile: Path):
    N.normalize_ice(mapfile, store=False)
    with h5py.File(mapfile) as f:
        assert "weight" not in f["bins"]
