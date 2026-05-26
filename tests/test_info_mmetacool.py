"""Tests for info/tree/attrs/metadata mmetacool support."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.io.mmetacool import create_mmetacool


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def test_mmetacool(tmp_path: Path, tiny_metacooler: Path) -> Path:
    out = tmp_path / "test.mmetacool"
    create_mmetacool(out, {"sample1": tiny_metacooler, "sample2": tiny_metacooler})
    return out


class TestInfoMmetacool:
    def test_info_mmetacool_without_group_raises(self, runner, test_mmetacool: Path):
        result = runner.invoke(cli, ["info", str(test_mmetacool)])
        assert result.exit_code != 0
        assert (
            "mmetacool" in result.output.lower() or "entries" in result.output.lower()
        )


class TestTreeMmetacool:
    def test_tree_mmetacool(self, runner, test_mmetacool: Path):
        result = runner.invoke(cli, ["tree", str(test_mmetacool)])
        assert result.exit_code == 0
        assert "entries" in result.output.lower()
        assert "sample1" in result.output
        assert "sample2" in result.output

    def test_tree_metacool(self, runner, tiny_metacooler: Path):
        result = runner.invoke(cli, ["tree", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "genomes" in result.output.lower()


class TestAttrsMmetacool:
    def test_attrs_mmetacool_root(self, runner, test_mmetacool: Path):
        result = runner.invoke(cli, ["attrs", str(test_mmetacool)])
        assert result.exit_code == 0
        assert "MMetaCooler" in result.output or "format" in result.output

    def test_attrs_mmetacool_entry(self, runner, test_mmetacool: Path):
        result = runner.invoke(cli, ["attrs", str(test_mmetacool), "entries/sample1"])
        assert result.exit_code == 0


class TestMetadataMmetacool:
    def test_metadata_list_mmetacool_raises(self, runner, test_mmetacool: Path):
        result = runner.invoke(cli, ["metadata", "list", str(test_mmetacool)])
        assert result.exit_code != 0
        assert (
            "mmetacool" in result.output.lower() or "entries" in result.output.lower()
        )
