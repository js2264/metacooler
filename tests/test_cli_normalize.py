"""CLI tests for ``metacooler normalize`` and the ``--balance`` options."""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
from _cli_output import flat, squash
from _metahic import MetaHiC, make_metahic
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.weights import read_weights


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def community(tmp_path_factory) -> MetaHiC:
    return make_metahic(
        tmp_path_factory.mktemp("community"),
        contigs_per_genome=6,
        n_pairs=12_000,
        seed=21,
    )


@pytest.fixture
def mapfile(tmp_path: Path, community: MetaHiC) -> Path:
    path = tmp_path / "map.metacool"
    shutil.copy2(community.metacool, path)
    return path


def _ok(runner: CliRunner, args: list[str]):
    result = runner.invoke(cli, args)
    assert result.exit_code == 0, result.output
    return result


def _columns(path: Path) -> set[str]:
    with h5py.File(path) as f:
        return set(f["bins"])


def test_help_lists_every_method(runner):
    result = _ok(runner, ["normalize", "--help"])
    for name in ("ice", "metator", "mge", "normcc", "hiczin", "list"):
        assert name in result.output


def test_normalise_is_an_alias(runner, mapfile):
    _ok(runner, ["normalise", "ice", str(mapfile), "-q"])
    assert "weight" in _columns(mapfile)
    assert "normalise" not in _ok(runner, ["--help"]).output


def test_ice_records_its_command(runner, mapfile):
    _ok(runner, ["normalize", "ice", str(mapfile), "--scope", "contig", "--name", "w"])
    with h5py.File(mapfile) as f:
        assert f["bins/w"].attrs["scope"] == "contig"
    commands = MetaCooler(mapfile).metadata["commands"]
    assert "normalize ice" in commands[-1]
    assert "--scope contig" in commands[-1]


def test_existing_column_needs_force(runner, mapfile):
    _ok(runner, ["normalize", "ice", str(mapfile), "-q"])
    result = runner.invoke(cli, ["normalize", "ice", str(mapfile)])
    assert result.exit_code == 1
    assert "already exists" in flat(result.output)
    _ok(runner, ["normalize", "ice", str(mapfile), "--scope", "global", "-f", "-q"])


def test_metator_from_enzymes(runner, mapfile, community):
    _ok(
        runner,
        [
            "normalize",
            "metator",
            str(mapfile),
            "--mode",
            "theoritical_hit",
            "--coverage",
            "cov",
            "--enzyme",
            "DpnII,HinfI",
            "--fasta",
            str(community.fasta),
            "-q",
        ],
    )
    with h5py.File(mapfile) as f:
        attrs = f["bins/metator_theoretical_hit"].attrs
        assert list(attrs["enzymes"]) == ["DpnII", "HinfI"]


@pytest.mark.parametrize(
    ("args", "code", "message"),
    [
        (["--mode", "bogus"], 2, "Unknown MetaTOR normalization"),
        (["--mode", "abundance"], 1, "--coverage"),
        (["--mode", "rs", "--enzyme", "DpnII"], 1, "--fasta"),
        (["--enzyme", ""], 2, "is empty"),
    ],
)
def test_metator_errors(runner, mapfile, args, code, message):
    result = runner.invoke(cli, ["normalize", "metator", str(mapfile), *args])
    assert result.exit_code == code, result.output
    assert message in squash(result.output) or message in flat(result.output)


def test_mge(runner, mapfile):
    _ok(runner, ["normalize", "mge", str(mapfile), "--min-distance", "2kb", "-q"])
    _ok(runner, ["normalize", "mge", str(mapfile), "--scope", "genome", "-q"])
    assert {"mge", "mge_genome"} <= _columns(mapfile)
    result = runner.invoke(
        cli, ["normalize", "mge", str(mapfile), "--min-distance", "far"]
    )
    assert result.exit_code == 2
    assert "--min-distance" in flat(result.output)


def test_normcc_and_hiczin(runner, mapfile, community):
    _ok(runner, ["normalize", "normcc", str(mapfile), "--sites", "RS", "-q"])
    _ok(
        runner,
        [
            "normalize",
            "hiczin",
            str(mapfile),
            "--coverage",
            "cov",
            "--enzyme",
            "DpnII",
            "--fasta",
            str(community.fasta),
            "--min-length",
            "1kb",
            "--max-pairs",
            "40",
            "--seed",
            "2",
            "-q",
        ],
    )
    assert {"normcc", "hiczin"} <= _columns(mapfile)
    with h5py.File(mapfile) as f:
        assert f["bins/hiczin"].attrs["n_pairs_fitted"] == 40
    result = runner.invoke(
        cli, ["normalize", "normcc", str(mapfile), "--min-length", "-5", "--name", "x"]
    )
    assert result.exit_code == 2


def test_list(runner, mapfile):
    result = _ok(runner, ["normalize", "list", str(mapfile)])
    assert "No weight columns" in result.output
    _ok(runner, ["normalize", "ice", str(mapfile), "-q"])
    _ok(runner, ["normalize", "metator", str(mapfile), "-q"])
    n = len(MetaCooler(mapfile).bins())
    with h5py.File(mapfile, "r+") as f:
        f["bins"].create_dataset("foreign", data=np.ones(n))
    result = runner.invoke(
        cli, ["normalize", "list", str(mapfile)], env={"COLUMNS": "200"}
    )
    assert result.exit_code == 0, result.output
    for text in ("weight", "ice", "scope=genome", "metator_empirical_hit", "foreign"):
        assert text in result.output


def test_balance_options_of_the_readers(runner, mapfile, tmp_path):
    _ok(runner, ["normalize", "normcc", str(mapfile), "-q"])
    edges = tmp_path / "edges.tsv"
    _ok(
        runner,
        [
            "network",
            str(mapfile),
            "--level",
            "genome",
            "--balance",
            "normcc",
            "--drop-spurious",
            "-o",
            str(edges),
            "-q",
        ],
    )
    assert edges.read_text()
    pixels = tmp_path / "pixels.tsv"
    _ok(
        runner,
        [
            "dump",
            str(mapfile),
            "-t",
            "pixels",
            "--balance",
            "normcc",
            "-H",
            "-o",
            str(pixels),
        ],
    )
    assert pixels.read_text().splitlines()[0].endswith("balanced")
    result = runner.invoke(
        cli, ["dump", str(mapfile), "-t", "genomes", "--balance", "normcc"]
    )
    assert result.exit_code == 2
    assert "only applies to the pixels table" in flat(result.output)
    image = tmp_path / "map.png"
    _ok(
        runner,
        [
            "show",
            str(mapfile),
            "--balance",
            "normcc",
            "--dpi",
            "40",
            "-o",
            str(image),
            "-q",
        ],
    )
    assert image.stat().st_size > 0
    result = runner.invoke(cli, ["network", str(mapfile), "--balance", "nope"])
    assert result.exit_code == 1
    assert "No weight column 'nope'" in flat(result.output)


def test_weights_match_the_api(runner, mapfile):
    _ok(runner, ["normalize", "metator", str(mapfile), "--mode", "length", "-q"])
    from metacooler.ops.normalize import normalize_metator

    api = normalize_metator(mapfile, mode="length", store=False).weights
    np.testing.assert_allclose(read_weights(MetaCooler(mapfile), "metator_length"), api)
