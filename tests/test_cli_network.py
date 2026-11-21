"""CLI tests for the ``metacooler network`` command."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def a1_rep1_1kb() -> Path:
    p = DATA / "assembly_1_rep1_1kb.metacool"
    assert p.exists(), f"Missing test data: {p}"
    return p


# ---------------------------------------------------------------------------
# TestNetworkCLI — thin fixtures (tiny_metacooler)
# ---------------------------------------------------------------------------


class TestNetworkCLITiny:
    def test_help(self, runner):
        result = runner.invoke(cli, ["network", "--help"])
        assert result.exit_code == 0
        assert "--level" in result.output
        assert "--weight" in result.output
        assert "--output" in result.output

    def test_contig_level_to_file(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "network.tsv"
        result = runner.invoke(
            cli,
            ["network", str(tiny_metacooler), "--level", "contig", "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        # Three tab-separated columns
        for line in lines:
            assert len(line.split("\t")) == 3

    def test_genome_level_to_file(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "network_genome.tsv"
        result = runner.invoke(
            cli,
            ["network", str(tiny_metacooler), "--level", "genome", "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0

    def test_bin_level_to_file(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "network_bin.tsv"
        result = runner.invoke(
            cli,
            ["network", str(tiny_metacooler), "--level", "bin", "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0

    def test_stdout_when_no_output(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli, ["network", str(tiny_metacooler), "--level", "contig"]
        )
        assert result.exit_code == 0, result.output
        # Filter log lines (start with '[') to get TSV rows only
        tsv_lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        assert len(tsv_lines) > 0
        for line in tsv_lines:
            assert len(line.split("\t")) == 3

    def test_bin_level_to_stdout(self, runner, tiny_metacooler, tmp_path):
        """Without -o, --level bin prints what -o writes."""
        out = tmp_path / "bins.tsv"
        result = runner.invoke(
            cli, ["network", str(tiny_metacooler), "--level", "bin", "-q"]
        )
        assert result.exit_code == 0, result.output
        written = runner.invoke(
            cli,
            ["network", str(tiny_metacooler), "--level", "bin", "-q", "-o", str(out)],
        )
        assert written.exit_code == 0, written.output
        assert result.output == out.read_text()

    def test_verbose_flag(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "network_v.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(tiny_metacooler),
                "--level",
                "contig",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_weight_count_explicit(self, runner, tiny_metacooler, tmp_path):
        out = tmp_path / "network_w.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(tiny_metacooler),
                "--level",
                "genome",
                "--weight",
                "count",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_invalid_weight_fails(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli,
            [
                "network",
                str(tiny_metacooler),
                "--level",
                "contig",
                "--weight",
                "nonexistent_column",
            ],
        )
        assert result.exit_code != 0

    def test_contig_nodes_are_contig_names(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli, ["network", str(tiny_metacooler), "--level", "contig"]
        )
        assert result.exit_code == 0, result.output
        lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        node_names = {ln.split("\t")[0] for ln in lines} | {
            ln.split("\t")[1] for ln in lines
        }
        # All node names should be contig names (not integers at contig level)
        assert any("contig" in n for n in node_names)

    def test_genome_nodes_are_genome_names(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli, ["network", str(tiny_metacooler), "--level", "genome"]
        )
        assert result.exit_code == 0, result.output
        lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        node_names = {ln.split("\t")[0] for ln in lines} | {
            ln.split("\t")[1] for ln in lines
        }
        assert any(n.startswith("MAG") or n == "__unbinned__" for n in node_names)

    def test_values_are_positive_integers(self, runner, tiny_metacooler):
        result = runner.invoke(
            cli, ["network", str(tiny_metacooler), "--level", "contig"]
        )
        assert result.exit_code == 0, result.output
        tsv_lines = [
            line
            for line in result.output.strip().split("\n")
            if line and not line.startswith("[")
        ]
        for line in tsv_lines:
            parts = line.split("\t")
            assert len(parts) == 3
            assert int(parts[2]) > 0

    def test_bin_level_equals_dump_pixels(self, runner, tiny_metacooler, tmp_path):
        """network --level bin and dump --table pixels produce identical rows."""
        net_out = tmp_path / "net_bin.tsv"
        dump_out = tmp_path / "dump_pixels.tsv"

        net_result = runner.invoke(
            cli,
            ["network", "--level", "bin", str(tiny_metacooler), "-o", str(net_out)],
        )
        assert net_result.exit_code == 0, net_result.output

        result = runner.invoke(
            cli,
            ["dump", "--table", "pixels", str(tiny_metacooler), "-o", str(dump_out)],
        )
        assert result.exit_code == 0, result.output

        net_rows = sorted(net_out.read_text().strip().split("\n"))
        dump_rows = sorted(dump_out.read_text().strip().split("\n"))
        assert net_rows == dump_rows, (
            "network --level bin and dump --table pixels differ"
        )


# ---------------------------------------------------------------------------
# TestNetworkCLIRealData — real assembly files
# ---------------------------------------------------------------------------


class TestNetworkCLIRealData:
    """Mirrors the three calls from cli_tests.sh network section."""

    def test_network_contig_level(self, runner, a1_rep1_1kb, tmp_path):
        """metacooler network assembly_1_rep1_1kb.metacool --level contig -v -o network.tsv"""
        out = tmp_path / "network.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--level",
                "contig",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        ctgs = []
        cnts = 0
        for line in lines:
            splitted = line.split("\t")
            assert len(splitted) == 3
            cnts += int(splitted[2])
            ctgs.extend(splitted[:2])
        assert cnts == 50000, "Total contact count should be 50000"
        assert len(set(ctgs)) == 44, "Expected 44 unique contigs"

    def test_network_genome_level(self, runner, a1_rep1_1kb, tmp_path):
        """metacooler network assembly_1_rep1_1kb.metacool --level genome -v -o network.tsv"""
        out = tmp_path / "network_genome.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--level",
                "genome",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        node_names = {ln.split("\t")[0] for ln in lines} | {
            ln.split("\t")[1] for ln in lines
        }
        assert any(n.startswith("MAG") for n in node_names)
        gens = []
        cnts = 0
        for line in lines:
            splitted = line.split("\t")
            assert len(splitted) == 3
            cnts += int(splitted[2])
            gens.extend(splitted[:2])
        assert cnts == 50000, "Total contact count should be 50000"
        assert len(set(gens)) == 11, "Expected 11 unique genomes"

    def test_network_genome_level_weight_count(self, runner, a1_rep1_1kb, tmp_path):
        """metacooler network assembly_1_rep1_1kb.metacool --weight count --level genome -v -o network.tsv"""
        out = tmp_path / "network_genome_count.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--weight",
                "count",
                "--level",
                "genome",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        # All weights should be positive
        for line in lines:
            parts = line.split("\t")
            assert int(parts[2]) > 0

    def test_network_contig_node_names_from_assembly(self, runner, a1_rep1_1kb):
        """Node names at contig level come from the contigs table."""
        from metacooler.api import MetaCooler

        mc = MetaCooler(a1_rep1_1kb)
        contig_names = set(mc.contigs["name"].values)

        out_result = runner.invoke(
            cli, ["network", str(a1_rep1_1kb), "--level", "contig"]
        )
        assert out_result.exit_code == 0, out_result.output
        for line in out_result.output.strip().split("\n"):
            if not line or line.startswith("["):
                continue
            parts = line.split("\t")
            assert parts[0] in contig_names, f"Unknown node: {parts[0]}"
            assert parts[1] in contig_names, f"Unknown node: {parts[1]}"

    def test_network_genome_node_names_from_assembly(self, runner, a1_rep1_1kb):
        """Node names at genome level come from the genomes table."""
        from metacooler.api import MetaCooler

        mc = MetaCooler(a1_rep1_1kb)
        genome_names = set(mc.genomes["name"].values)

        out_result = runner.invoke(
            cli, ["network", str(a1_rep1_1kb), "--level", "genome"]
        )
        assert out_result.exit_code == 0, out_result.output
        for line in out_result.output.strip().split("\n"):
            if not line or line.startswith("["):
                continue
            parts = line.split("\t")
            assert parts[0] in genome_names, f"Unknown node: {parts[0]}"
            assert parts[1] in genome_names, f"Unknown node: {parts[1]}"

    def test_network_genome_sum_equals_contig_sum(self, runner, a1_rep1_1kb):
        """Total weight at genome level == total weight at contig level."""
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            genome_out = tmp / "g.tsv"
            contig_out = tmp / "c.tsv"
            runner.invoke(
                cli,
                [
                    "network",
                    str(a1_rep1_1kb),
                    "--level",
                    "genome",
                    "-o",
                    str(genome_out),
                ],
            )
            runner.invoke(
                cli,
                [
                    "network",
                    str(a1_rep1_1kb),
                    "--level",
                    "contig",
                    "-o",
                    str(contig_out),
                ],
            )
            genome_sum = sum(
                int(ln.split("\t")[2])
                for ln in genome_out.read_text().strip().split("\n")
                if ln
            )
            contig_sum = sum(
                int(ln.split("\t")[2])
                for ln in contig_out.read_text().strip().split("\n")
                if ln
            )
            assert genome_sum == contig_sum
