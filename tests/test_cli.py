"""CLI tests for metacooler."""

from __future__ import annotations

import json
import subprocess

import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli


@pytest.fixture
def runner():
    return CliRunner()


class TestCreate:
    def test_create_help(self, runner):
        result = runner.invoke(cli, ["create", "--help"])
        assert result.exit_code == 0
        assert "cool" in result.output
        assert "pairs" in result.output

    def test_create(self, tiny_metacooler):
        # tiny_metacooler fixture already runs create via create_metacooler
        assert tiny_metacooler.exists()


class TestDump:
    def test_dump_chroms(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["dump", "-t", "chroms", str(tiny_metacooler)])
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) == 4  # 4 contigs, no header

    def test_dump_bins(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["dump", "-t", "bins", str(tiny_metacooler)])
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) == 22  # 22 bins, no header

    def test_dump_pixels(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["dump", "-t", "pixels", str(tiny_metacooler)])
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) > 0

    def test_dump_genomes(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["dump", "-t", "genomes", str(tiny_metacooler)])
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) == 3  # 3 genomes, no header
        assert "__unbinned__" in result.output

    def test_dump_contigs(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["dump", "-t", "contigs", str(tiny_metacooler)])
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) == 4  # 4 contigs, no header

    def test_dump_header(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli, ["dump", "-t", "genomes", "-H", str(tiny_metacooler)]
        )
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) == 4  # 3 genomes + header row

    def test_dump_matches_cooler_for_bins(self, runner, tiny_metacooler):
        mc = runner.invoke(cli, ["dump", "-t", "bins", str(tiny_metacooler)])
        cc = subprocess.run(
            ["cooler", "dump", "-t", "bins", str(tiny_metacooler)],
            capture_output=True,
            text=True,
        )
        assert mc.exit_code == 0
        assert cc.returncode == 0
        assert mc.output == cc.stdout

    def test_dump_matches_cooler_for_pixels(self, runner, tiny_metacooler):
        mc = runner.invoke(cli, ["dump", "-t", "pixels", str(tiny_metacooler)])
        cc = subprocess.run(
            ["cooler", "dump", "-t", "pixels", str(tiny_metacooler)],
            capture_output=True,
            text=True,
        )
        assert mc.exit_code == 0
        assert cc.returncode == 0
        assert mc.output == cc.stdout


class TestInfo:
    def test_info(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["info", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "HDF5::MetaCooler" in result.output
        data = json.loads(result.output)
        assert data["format"] == "HDF5::MetaCooler"
        assert data["format-version"] == "3"  # numpy scalar serialized via default=str
        assert data["metacooler-format-version"] == "1"
        assert data["ngenomes"] == 3
        assert data["ncontigs"] == 4
        assert data["ordering"] == "genome-then-alpha"

    def test_info_table(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["info", "--table", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "Genomes" in result.output


class TestTree:
    def test_tree(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["tree", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "genomes/" in result.output
        assert "contigs/" in result.output


class TestAttrs:
    def test_attrs_root(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["attrs", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "format" in result.output
        assert "metacooler-format" in result.output

    def test_attrs_group(self, runner, tiny_metacooler):
        result = runner.invoke(cli, ["attrs", str(tiny_metacooler), "/genomes"])
        assert result.exit_code == 0


class TestAggregate:
    def test_aggregate_contig(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "agg.cool"
        result = runner.invoke(cli, ["aggregate", str(tiny_metacooler), "-o", str(out)])
        assert result.exit_code == 0
        assert out.exists()

    def test_aggregate_genome(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "agg.cool"
        result = runner.invoke(
            cli,
            ["aggregate", "--level", "genome", str(tiny_metacooler), "-o", str(out)],
        )
        assert result.exit_code == 0
        assert out.exists()


class TestShow:
    def test_show_genome(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        result = runner.invoke(
            cli, ["show", "--level", "genome", "-o", str(out), "-f", str(tiny_metacooler)]
        )
        assert result.exit_code == 0
        assert out.exists()

    def test_show_contig(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        result = runner.invoke(
            cli, ["show", "--level", "contig", "-o", str(out), "-f", str(tiny_metacooler)]
        )
        assert result.exit_code == 0
        assert out.exists()


class TestMerge:
    def test_merge(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        result = runner.invoke(
            cli, ["merge", "-o", str(out), str(tiny_metacooler), str(tiny_metacooler)]
        )
        assert result.exit_code == 0
        assert out.exists()


class TestExtract:
    def test_extract_help(self, runner):
        result = runner.invoke(cli, ["extract", "--help"])
        assert result.exit_code == 0
        assert "--genomes" in result.output

    def test_extract_unknown_genome(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "unknown.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                str(tiny_metacooler),
                "--genomes",
                "NONEXISTENT",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code != 0

    def test_extract_unbinned(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "unbinned.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                str(tiny_metacooler),
                "--genomes",
                "__unbinned__",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
