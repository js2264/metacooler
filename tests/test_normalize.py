"""Tests of ``metacooler normalize``: ICE, NormCC, HiCzin, and the weight columns."""

from __future__ import annotations

import gzip
import json
import logging
import shutil
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _metahic import MetaHiC, make_metahic
from numpy.testing import assert_allclose, assert_array_equal

import metacooler
from metacooler import normalize as engine
from metacooler.api import MetaCooler
from metacooler.io import weights as W
from metacooler.io.annotate import annotate
from metacooler.io.create import from_pairs
from metacooler.io.reassign import reassign
from metacooler.io.sites import count_sites, enzyme_sites, iter_fasta, site_pattern
from metacooler.ops import normalize as N
from metacooler.schema import validate

# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def community(tmp_path_factory) -> MetaHiC:
    """Three genomes of three contigs and one unbinned contig, 1 kb bins."""
    return make_metahic(tmp_path_factory.mktemp("community"))


@pytest.fixture(scope="module")
def rich_community(tmp_path_factory) -> MetaHiC:
    """Genomes of many contigs, so that HiCzin has pairs to fit on."""
    return make_metahic(
        tmp_path_factory.mktemp("rich"),
        genomes=3,
        contigs_per_genome=9,
        unbinned=2,
        n_pairs=20_000,
        seed=5,
    )


@pytest.fixture
def mapfile(tmp_path: Path, community: MetaHiC) -> Path:
    """A fresh, writable copy of the community map."""
    path = tmp_path / "map.metacool"
    shutil.copy2(community.metacool, path)
    return path


@pytest.fixture
def rich_mapfile(tmp_path: Path, rich_community: MetaHiC) -> Path:
    path = tmp_path / "rich.metacool"
    shutil.copy2(rich_community.metacool, path)
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


def _floats(value: object) -> np.ndarray:
    """A stored attribute (typed ``object``) as a float array."""
    return np.asarray(value, dtype=np.float64)


def _contig_matrix(mc: MetaCooler) -> np.ndarray:
    """Dense contig-by-contig contacts, intra-contig contacts on the diagonal."""
    contig_of_bin = mc._bin_contig_ids
    n = len(mc.contigs)
    out = np.zeros((n, n))
    pixels = mc.pixels()
    c1 = contig_of_bin[pixels["bin1_id"].to_numpy()]
    c2 = contig_of_bin[pixels["bin2_id"].to_numpy()]
    np.add.at(out, (c1, c2), pixels["count"].to_numpy())
    return out + out.T - np.diag(np.diag(out))


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
# Engine: regressions
# ----------------------------------------------------------------------------


def _nb_negloglik(beta, y, X, alpha=1.0):
    from scipy.special import gammaln

    mu = np.exp(X @ beta)
    r = 1.0 / alpha
    return -np.sum(
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1)
        + r * np.log(r / (r + mu))
        + y * np.log(mu / (r + mu))
    )


def _nb_mle(y, X, alpha=1.0):
    """Maximum likelihood of a log-link NB model of known dispersion, by Newton.

    The log-likelihood is strictly concave in the coefficients (its second
    derivative in the linear predictor, -mu (1 + alpha y) / (1 + alpha mu)**2,
    is negative), so the one stationary point Newton's method converges to,
    quadratically, is the maximum -- to machine precision, where a quasi-Newton
    optimiser stops "due to precision loss" wherever its scipy version does.
    """
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(100):
        mu = np.exp(X @ beta)
        score = X.T @ ((y - mu) / (1 + alpha * mu))
        info = (X.T * (mu * (1 + alpha * y) / (1 + alpha * mu) ** 2)) @ X
        step = np.linalg.solve(info, score)
        t = 1.0
        while (
            _nb_negloglik(beta + t * step, y, X, alpha)
            > _nb_negloglik(beta, y, X, alpha)
            and t > 1e-8
        ):
            t /= 2
        beta = beta + t * step
        if np.max(np.abs(t * step)) < 1e-13 * (1 + np.max(np.abs(beta))):
            break
    return beta


def test_nb_glm_reaches_the_maximum_likelihood():
    rng = np.random.default_rng(2)
    X = np.column_stack([np.ones(2000), rng.normal(size=(2000, 2))])
    mu = np.exp(X @ np.array([1.5, 0.4, -0.7]))
    y = rng.negative_binomial(1.0, 1.0 / (1.0 + mu)).astype(float)
    assert_allclose(engine.fit_nb_glm(y, X), _nb_mle(y, X), rtol=1e-7)


def test_nb_glm_refuses_an_all_zero_response():
    with pytest.raises(ValueError, match="no signal"):
        engine.fit_nb_glm(np.zeros(5), np.ones((5, 1)))


def test_nb_glm_warns_when_it_stops_early(caplog):
    rng = np.random.default_rng(0)
    X = np.column_stack([np.ones(50), rng.normal(size=50)])
    y = rng.poisson(5, 50).astype(float)
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        engine.fit_nb_glm(y, X, max_iter=1)
    assert "did not converge" in caplog.text


def _simulate_zinb(n, seed):
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
    beta = np.array([1.0, 0.5, -0.3, 0.8])
    gamma = np.array([-0.5, 0.4, 0.0, -0.6])
    theta = 2.5
    mu = np.exp(X @ beta)
    pi = 1 / (1 + np.exp(-(X @ gamma)))
    counts = rng.negative_binomial(theta, theta / (theta + mu))
    y = np.where(rng.random(n) < pi, 0, counts).astype(float)
    return y, X, beta, gamma, theta


def test_zinb_gradient_matches_finite_differences():
    from scipy.optimize import check_grad

    y, X, beta, gamma, theta = _simulate_zinb(400, 4)
    params = np.concatenate([beta, gamma, [np.log(theta)]]) + 0.1

    def f(p):
        return -engine.zinb_loglik(p, y, X, X)[0]

    def g(p):
        return -engine.zinb_loglik(p, y, X, X)[1]

    assert check_grad(f, g, params) / np.linalg.norm(g(params)) < 1e-5


def test_zinb_fit_recovers_the_simulated_model():
    y, X, beta, gamma, theta = _simulate_zinb(8000, 3)
    fit = engine.fit_zinb(y, X)
    assert fit.converged
    assert_allclose(fit.beta, beta, atol=0.06)
    assert_allclose(fit.gamma, gamma, atol=0.15)
    assert fit.theta == pytest.approx(theta, rel=0.15)


def test_zinb_fit_without_zeros_settles_on_no_inflation(caplog):
    rng = np.random.default_rng(3)
    X = np.column_stack([np.ones(300), rng.normal(size=(300, 3))])
    y = rng.poisson(5 + np.exp(X @ np.array([1.0, 0.5, -0.3, 0.8]))) + 1.0
    with caplog.at_level(logging.INFO, logger="metacooler"):
        fit = engine.fit_zinb(y, X)
    assert fit.converged
    assert fit.gamma[0] < -10  # a zero-inflation probability of ~0
    assert "No zero inflation in these counts" in caplog.text


def test_zinb_fit_of_underdispersed_counts_reduces_to_a_poisson(caplog):
    # Binomial counts vary less than their mean: no negative binomial fits
    # them better than its Poisson limit, so the dispersion runs to its bound.
    rng = np.random.default_rng(9)
    X = np.column_stack([np.ones(2000), rng.normal(size=(2000, 2))])
    mu = np.exp(X @ np.array([2.0, 0.3, -0.2]))
    y = rng.binomial(200, mu / 200).astype(float)
    with caplog.at_level(logging.INFO, logger="metacooler"):
        fit = engine.fit_zinb(y, X)
    assert fit.converged
    assert fit.theta == pytest.approx(np.exp(12.0))
    assert_allclose(fit.beta, [2.0, 0.3, -0.2], atol=0.03)
    assert "not overdispersed" in caplog.text


def test_zinb_fit_warns_when_the_optimiser_gives_up(caplog):
    y, X, *_ = _simulate_zinb(500, 6)
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        fit = engine.fit_zinb(y, X, max_iter=1)
    assert not fit.converged
    assert "stopped before converging" in caplog.text


def test_hiczin_weights_factorise_the_fitted_mean():
    rng = np.random.default_rng(7)
    logs = rng.normal(5, 2, size=(20, 3))
    beta = np.array([0.3, 0.8, -0.2, 0.5])
    i, j = rng.integers(0, 20, 50), rng.integers(0, 20, 50)
    pair = logs[i] + logs[j]
    means, sds = pair.mean(axis=0), pair.std(axis=0, ddof=1)
    mu = np.exp(beta[0] + ((pair - means) / sds) @ beta[1:])
    w = engine.hiczin_contig_weights(beta, means, sds, logs)
    assert_allclose(w[i] * w[j], 1.0 / mu)


def test_pairs_within_groups_decodes_every_pair():
    sizes = np.array([0, 1, 4, 2, 7, 3])
    expected = [
        (g, a, b)
        for g, n in enumerate(sizes)
        for a in range(n)
        for b in range(a + 1, n)
    ]
    group, a, b = engine.pairs_within_groups(sizes, np.arange(len(expected)))
    assert list(zip(group.tolist(), a.tolist(), b.tolist(), strict=True)) == expected


def test_pairs_within_groups_decodes_a_large_group():
    n = 100_000
    total = n * (n - 1) // 2
    index = np.array([0, 1, n - 2, n - 1, total // 2, total - 2, total - 1])
    _, a, b = engine.pairs_within_groups(np.array([n]), index)
    rank = a * n - a * (a + 1) // 2 + (b - a - 1)
    assert_array_equal(rank, index)
    assert (a < b).all() and (b < n).all()


# ----------------------------------------------------------------------------
# Restriction sites
# ----------------------------------------------------------------------------


def test_enzyme_sites_by_name_and_iupac():
    assert enzyme_sites(["dpnii", "HinfI"]) == ["GATC", "GANTC"]
    assert enzyme_sites(["Arima"]) == ["GATC", "GANTC"]
    assert enzyme_sites(["ggatcc", " "]) == ["GGATCC"]


def test_enzyme_sites_refuse_unknown_and_empty():
    with pytest.raises(ValueError, match="Unknown enzyme"):
        enzyme_sites(["NotAnEnzyme"])
    with pytest.raises(ValueError, match="No enzyme"):
        enzyme_sites([""])


def test_site_pattern_joins_the_sites_of_every_enzyme():
    # One alternation, so a position two enzymes cut counts once.
    assert site_pattern(["GATC", "GANTC"]).pattern == "GA.TC|GATC"


def test_site_pattern_searches_both_strands_of_a_non_palindrome():
    pattern = site_pattern(["GAAG"])
    assert pattern.pattern == "CTTC|GAAG"
    assert site_pattern(["RAAT"]).pattern == "ATT[CT]|[AG]AAT"


def test_iter_fasta_reads_plain_and_gzip(tmp_path: Path):
    text = ">a desc\nACGT\nAC\n>b\nGG\n"
    plain = tmp_path / "x.fa"
    plain.write_text(text)
    zipped = tmp_path / "x.fa.gz"
    with gzip.open(zipped, "wt") as f:
        f.write(text)
    for path in (plain, zipped):
        assert list(iter_fasta(path)) == [("a", "ACGTAC"), ("b", "GG")]


def test_count_sites(tmp_path: Path):
    fasta = tmp_path / "x.fa"
    fasta.write_text(">a\nGATCgatc\nGAATC\n>b\nTTTT\n>c\nGATC\n")
    counts = count_sites(fasta, ["b", "a"], [4, 13], ["DpnII", "HinfI"])
    assert_array_equal(counts, [0, 3])


def test_count_sites_refuses_another_assembly(tmp_path: Path):
    fasta = tmp_path / "x.fa"
    fasta.write_text(">a\nGATC\n")
    with pytest.raises(ValueError, match="another length"):
        count_sites(fasta, ["a"], [5], ["DpnII"])
    with pytest.raises(ValueError, match="not in"):
        count_sites(fasta, ["a", "zz"], [4, 4], ["DpnII"])


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


def test_spurious_cutoff(mapfile: Path):
    n = _n_bins(mapfile)
    _store(
        mapfile,
        "w",
        np.ones(n),
        parameters={"spurious_cutoff": 1.5, "spurious_rule": "<="},
    )
    _store(mapfile, "plain", np.ones(n))
    assert W.spurious_cutoff(mapfile, "w") == (1.5, "<=")
    assert W.spurious_cutoff(MetaCooler(mapfile), "plain") is None
    assert W.spurious_cutoff(mapfile, "missing") is None


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


# ----------------------------------------------------------------------------
# Restriction sites and coverages, from /contigs or the assembly
# ----------------------------------------------------------------------------


def test_sites_from_a_column_or_the_fasta(mapfile: Path, community: MetaHiC):
    by_fasta = N.normalize_normcc(
        mapfile, enzymes=["DpnII"], fasta=community.fasta, store=False
    )
    assert by_fasta.attrs["enzymes"] == ["DpnII"]
    assert by_fasta.attrs["depends_on"] == ["pixels"]
    by_column = N.normalize_normcc(mapfile, sites="RS", store=False)
    assert by_column.attrs["sites"] == "contigs/RS"
    assert by_column.attrs["depends_on"] == ["pixels", "contigs/RS"]


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"sites": "RS", "enzymes": ["DpnII"]}, "not both"),
        ({"fasta": "x.fa"}, "pass the enzyme too"),
        ({"enzymes": ["DpnII"]}, "--fasta"),
        ({"sites": "nope"}, "no column 'nope'"),
    ],
)
def test_site_sources_that_are_refused(mapfile: Path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        N.normalize_normcc(mapfile, **kwargs)


def test_unusable_column_values_are_masked(caplog):
    contigs = pd.DataFrame(
        {"name": ["a", "b", "c", "d"], "cov": ["1.5", "NA", "-2", "4"]}
    )
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        values = N._numeric_contig_column(contigs, "cov", what="coverage")
    assert_array_equal(np.isnan(values), [False, True, True, False])
    assert "2 contig(s) have no usable coverage" in caplog.text


# ----------------------------------------------------------------------------
# NormCC
# ----------------------------------------------------------------------------


def _normcc_reference(mc: MetaCooler, sites, min_length=1000, min_signal=2):
    """MetaCC's NormCC, written out on a dense contig matrix."""
    dense = _contig_matrix(mc)
    lengths = mc.contigs["length"].to_numpy(float)
    intra = np.diag(dense).copy()
    off = dense - np.diag(intra)
    long_enough = lengths >= min_length
    keep = long_enough & ((off * long_enough[None, :]).max(axis=1) >= min_signal)
    keep &= intra > 0
    columns = [np.ones(len(lengths))]
    if sites is not None:
        keep &= sites > 0
        columns.append(np.log(np.where(keep, sites, 1)))
    columns += [np.log(lengths), np.log(np.where(keep, intra, 1))]
    X = np.column_stack(columns)[keep]
    signal = (off[:, keep]).sum(axis=1)[keep]
    beta = _nb_mle(signal, X)
    mu = np.exp(X @ beta)
    w = np.full(len(lengths), np.nan)
    w[keep] = np.sqrt(mu.max() / mu)
    return beta, w, off, keep


def test_normcc_matches_metaccs_model(mapfile: Path, community: MetaHiC):
    result = N.normalize_normcc(
        mapfile, enzymes=["DpnII", "HinfI"], fasta=community.fasta
    )
    mc = MetaCooler(mapfile)
    sites = count_sites(
        community.fasta,
        mc.contigs["name"].tolist(),
        mc.contigs["length"],
        ["DpnII", "HinfI"],
    ).astype(float)
    beta, w, off, keep = _normcc_reference(mc, sites)
    assert_allclose(_floats(result.attrs["coefficients"]), beta, rtol=1e-7)
    assert_allclose(result.weights, w[mc._bin_contig_ids], rtol=1e-7)
    assert result.attrs["covariates"] == [
        "intercept",
        "log_sites",
        "log_length",
        "log_intra_contacts",
    ]
    i, j = np.triu_indices(len(w), k=1)
    sel = keep[i] & keep[j] & (off[i, j] > 0)
    normalized = off[i, j][sel] * w[i][sel] * w[j][sel]
    assert result.attrs["spurious_cutoff"] == pytest.approx(
        np.percentile(normalized, 5), rel=1e-7
    )
    assert result.attrs["spurious_rule"] == "<"


def test_normcc_without_sites(mapfile: Path):
    result = N.normalize_normcc(mapfile)
    beta, w, _, _ = _normcc_reference(MetaCooler(mapfile), None)
    assert result.attrs["covariates"] == [
        "intercept",
        "log_length",
        "log_intra_contacts",
    ]
    assert_allclose(_floats(result.attrs["coefficients"]), beta, rtol=1e-7)


def test_normcc_needs_contigs_to_fit(mapfile: Path):
    with pytest.raises(ValueError, match="pass NormCC's filters"):
        N.normalize_normcc(mapfile, min_length=10**9)
    with pytest.raises(ValueError, match="must be a fraction"):
        N.normalize_normcc(mapfile, spurious=1.5)


# ----------------------------------------------------------------------------
# HiCzin
# ----------------------------------------------------------------------------


def test_hiczin_normalizes_by_the_fitted_count_model(rich_mapfile: Path):
    result = N.normalize_hiczin(rich_mapfile, coverage="cov", sites="RS")
    mc = MetaCooler(rich_mapfile)
    attrs = result.attrs
    assert attrs["converged"]
    assert attrs["covariates"] == [
        "intercept",
        "log_sites",
        "log_length",
        "log_coverage",
    ]
    assert attrs["n_pairs_fitted"] == attrs["n_pairs_total"] == 3 * 36
    contigs = mc.contigs
    logs = np.column_stack(
        [
            np.log(contigs["RS"].to_numpy(float)),
            np.log(contigs["length"].to_numpy(float)),
            np.log(contigs["cov"].to_numpy(float)),
        ]
    )
    w_contig = np.full(len(contigs), np.nan)
    first_bin = mc._chrom_offset[:-1]
    w_contig[:] = result.weights[first_bin]
    # Every pair of contigs of one genome: c * w_i * w_j == c / mu_ij.
    genome = contigs["genome_id"].to_numpy()
    binned = mc.genomes["name"].to_numpy()[genome] != "__unbinned__"
    i, j = np.triu_indices(len(contigs), k=1)
    sel = (genome[i] == genome[j]) & binned[i]
    pair = logs[i[sel]] + logs[j[sel]]
    means = _floats(attrs["covariate_means"])
    coefficients = _floats(attrs["coefficients"])
    z = (pair - means) / _floats(attrs["covariate_sds"])
    mu = np.exp(coefficients[0] + z @ coefficients[1:])
    assert_allclose(w_contig[i[sel]] * w_contig[j[sel]], 1 / mu)
    assert_allclose(means, pair.mean(axis=0))
    assert attrs["spurious_rule"] == "<="
    assert attrs["depends_on"] == ["pixels", "genome_id", "contigs/cov", "contigs/RS"]


def test_hiczin_samples_pairs_reproducibly(rich_mapfile: Path):
    a = N.normalize_hiczin(
        rich_mapfile, coverage="cov", max_pairs=50, seed=3, store=False
    )
    b = N.normalize_hiczin(
        rich_mapfile, coverage="cov", max_pairs=50, seed=3, store=False
    )
    assert a.attrs["n_pairs_fitted"] == 50
    assert a.attrs["seed"] == 3
    assert_array_equal(a.weights, b.weights)
    assert a.attrs["covariates"] == ["intercept", "log_length", "log_coverage"]


def test_hiczin_raises_a_zero_coverage(rich_mapfile: Path):
    with h5py.File(rich_mapfile, "r+") as f:
        cov = f["contigs/cov"][()].astype(float)
        cov[0] = 0.0
        del f["contigs/cov"]
        f["contigs"].create_dataset("cov", data=cov)
    result = N.normalize_hiczin(rich_mapfile, coverage="cov", store=False)
    assert np.isfinite(result.weights[0])


def test_hiczin_needs_genomes_with_two_contigs(tmp_path: Path):
    lonely = make_metahic(tmp_path / "lonely", contigs_per_genome=1, n_pairs=2000)
    with pytest.raises(ValueError, match="no genome has two contigs"):
        N.normalize_hiczin(lonely.metacool, coverage="cov")


def test_hiczin_refuses_a_constant_covariate(rich_mapfile: Path):
    with h5py.File(rich_mapfile, "r+") as f:
        f["contigs"].create_dataset("flat", data=np.full(len(f["contigs/name"]), 7.0))
    with pytest.raises(ValueError, match="log_sites is constant"):
        N.normalize_hiczin(rich_mapfile, coverage="cov", sites="flat")


def test_hiczin_refuses_bad_parameters(rich_mapfile: Path):
    with pytest.raises(ValueError, match="max_pairs must be"):
        N.normalize_hiczin(rich_mapfile, coverage="cov", max_pairs=-1)
    with h5py.File(rich_mapfile, "r+") as f:
        f["contigs"].create_dataset("zero", data=np.zeros(len(f["contigs/name"])))
    with pytest.raises(ValueError, match="positive coverage"):
        N.normalize_hiczin(rich_mapfile, coverage="zero")


def test_hiczin_goes_stale_with_the_assignment_and_the_coverage(
    tmp_path: Path, rich_mapfile: Path
):
    """HiCzin is fitted on the pairs of contigs of one genome, on a coverage:
    a reassignment drops it, and so does an edit of that coverage column."""
    N.normalize_hiczin(rich_mapfile, coverage="cov")
    N.normalize_normcc(rich_mapfile)
    moved = tmp_path / "moved.metacool"
    reassign(rich_mapfile, {"c0": "G1"}, moved, layout="keep")
    with h5py.File(moved) as f:
        assert "hiczin" not in f["bins"]
        assert "normcc" in f["bins"]  # reads neither
    names = MetaCooler(rich_mapfile).contigs["name"].tolist()
    table = tmp_path / "cov.tsv"
    table.write_text(
        "contig_name\tcov\n" + "".join(f"{n}\t{i + 1}\n" for i, n in enumerate(names))
    )
    annotate(rich_mapfile, contigs_table=table, replace=True, inplace=True)
    with h5py.File(rich_mapfile) as f:
        assert "hiczin" not in f["bins"]
        assert "normcc" in f["bins"]
