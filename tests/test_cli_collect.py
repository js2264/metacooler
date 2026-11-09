"""CLI tests for ``metacooler collect``."""

from __future__ import annotations

from pathlib import Path

import pytest
from _cli_output import flat as _flat
from _cli_output import squash as _squash
from _cools import BSUB, ECOLI, make_cool
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli

DATA = Path(__file__).parent / "data"


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def cools(tmp_path: Path) -> tuple[Path, Path]:
    return (
        make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1),
        make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2),
    )


def run(runner: CliRunner, *args: object):
    return runner.invoke(cli, [str(a) for a in args])


class TestCollectCLI:
    def test_listed_in_root_help(self, runner):
        result = run(runner, "--help")
        assert result.exit_code == 0
        assert "collect" in result.output

    def test_help(self, runner):
        result = run(runner, "collect", "--help")
        assert result.exit_code == 0, result.output
        for option in (
            "--names",
            "--genomes",
            "--contigs",
            "--filter",
            "--contig-order",
        ):
            assert option in result.output

    def test_basic(self, runner, tmp_path, cools):
        out = tmp_path / "out.metacool"
        result = run(runner, "collect", *cools, "-o", out)
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist() == ["B_subtilis", "E_coli"]
        assert mc.info["ordering"] == "cool"

    def test_records_the_command(self, runner, tmp_path, cools):
        out = tmp_path / "out.metacool"
        result = run(runner, "collect", *cools, "--names", "Ec,Bs", "-o", out)
        assert result.exit_code == 0, result.output
        (command,) = MetaCooler(out).metadata["commands"]
        # CliRunner names the program after the click group ("cli collect").
        assert "collect" in command
        assert "--names Ec --names Bs" in command
        assert str(cools[0]) in command

    def test_names_repeated_or_comma_separated(self, runner, tmp_path, cools):
        for i, names in enumerate(
            (["--names", "Ec,Bs"], ["--names", "Ec", "--names", "Bs"])
        ):
            out = tmp_path / f"out{i}.metacool"
            result = run(runner, "collect", *cools, *names, "-o", out)
            assert result.exit_code == 0, result.output
            assert MetaCooler(out).genomes["name"].tolist() == ["Bs", "Ec"]

    def test_empty_names_rejected(self, runner, tmp_path, cools):
        result = run(
            runner, "collect", *cools, "--names", "", "-o", tmp_path / "o.metacool"
        )
        assert result.exit_code != 0
        assert "--names was given but is empty" in _flat(result.output)

    def test_names_count_mismatch_is_a_clean_error(self, runner, tmp_path, cools):
        result = run(
            runner, "collect", *cools, "--names", "one", "-o", tmp_path / "o.metacool"
        )
        assert result.exit_code == 1
        assert "Traceback" not in result.output
        assert "1 genome name(s) for 2 input(s)" in _flat(result.output)

    def test_requires_an_input(self, runner, tmp_path):
        result = run(runner, "collect", "-o", tmp_path / "o.metacool")
        assert result.exit_code == 2
        assert "Missing argument" in _flat(result.output)

    def test_missing_input(self, runner, tmp_path, cools):
        result = run(
            runner, "collect", cools[0], tmp_path / "nope.cool", "-o", tmp_path / "o"
        )
        assert result.exit_code == 2
        assert "does not exist" in _flat(result.output)

    def test_output_exists_needs_force(self, runner, tmp_path, cools):
        out = tmp_path / "out.metacool"
        assert run(runner, "collect", *cools, "-o", out).exit_code == 0
        result = run(runner, "collect", *cools, "-o", out)
        assert result.exit_code == 1
        assert "already exists" in _flat(result.output)
        result = run(runner, "collect", *cools, "-o", out, "--force")
        assert result.exit_code == 0, result.output

    def test_force_never_deletes_an_input(self, runner, tmp_path, cools):
        before = cools[0].read_bytes()
        result = run(runner, "collect", *cools, "-o", cools[0], "--force")
        assert result.exit_code == 2
        assert "An input is never overwritten" in _flat(result.output)
        assert cools[0].read_bytes() == before

    def test_mcool_root_is_a_clean_error(self, runner, tmp_path, cools):
        result = run(
            runner,
            "collect",
            cools[0],
            DATA / "assembly_1_rep1_1kb.mcool",
            "-o",
            tmp_path / "o.metacool",
        )
        assert result.exit_code == 1
        assert "Traceback" not in result.output
        assert "resolutions/1000" in _squash(result.output)

    def test_mcool_uris(self, runner, tmp_path):
        out = tmp_path / "out.metacool"
        result = run(
            runner,
            "collect",
            f"{DATA / 'assembly_1_rep1_1kb.mcool'}::resolutions/5000",
            f"{DATA / 'assembly_2_rep1_1kb.mcool'}::resolutions/5000",
            "-o",
            out,
        )
        assert result.exit_code == 0, result.output
        assert MetaCooler(out).info["bin-size"] == 5000

    def test_tables_filter_and_metadata(self, runner, tmp_path, cools):
        genomes = tmp_path / "species.tsv"
        genomes.write_text("genome_name\tphylum\nE_coli\tProteobacteria\n")
        contigs = tmp_path / "contigs.tsv"
        contigs.write_text(
            "contig_name\tgenome_name\tcoverage\nchromosome\tE_coli\t42\n"
        )
        metadata = tmp_path / "meta.yaml"
        metadata.write_text("study: demo\n")
        out = tmp_path / "out.metacool"
        result = run(
            runner,
            "collect",
            *cools,
            "--genomes",
            genomes,
            "--contigs",
            contigs,
            "--metadata",
            metadata,
            "--filter",
            "contigs.length >= 300",
            "--assembly-name",
            "two_species",
            "-o",
            out,
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert mc.metadata["study"] == "demo"
        assert mc.info["assembly"] == "two_species"
        assert "E_coli__plasmid1" not in mc.contigs["name"].tolist()
        g = mc.genomes.set_index("name")
        assert g.loc["E_coli", "phylum"] == "Proteobacteria"
        c = mc.contigs.set_index("name")
        assert c.loc["E_coli__chromosome", "coverage"] == 42

    def test_metadata_may_not_set_commands(self, runner, tmp_path, cools):
        metadata = tmp_path / "bad.yaml"
        metadata.write_text("commands:\n  - forged\n")
        result = run(
            runner, "collect", *cools, "--metadata", metadata, "-o", tmp_path / "o"
        )
        assert result.exit_code == 2
        assert "reserved for provenance" in _flat(result.output)

    @pytest.mark.parametrize("strategy", ["auto", "memory", "external"])
    def test_contig_order_and_sort_strategy(self, runner, tmp_path, cools, strategy):
        out = tmp_path / f"{strategy}.metacool"
        result = run(
            runner,
            "collect",
            *cools,
            "--contig-order",
            "genome-then-length",
            "--sort-strategy",
            strategy,
            "--chunksize",
            "50",
            "-o",
            out,
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert mc.info["ordering"] == "genome-then-length"
        assert mc.contigs["name"].tolist()[-2:] == [
            "E_coli__plasmid2",
            "E_coli__plasmid1",
        ]

    def test_no_qc_and_qc_threshold(self, runner, tmp_path, cools):
        out = tmp_path / "noqc.metacool"
        assert run(runner, "collect", *cools, "--no-qc", "-o", out).exit_code == 0
        assert "qc" not in MetaCooler(out).metadata
        out = tmp_path / "qc.metacool"
        result = run(runner, "collect", *cools, "--qc-threshold", "1kb", "-o", out)
        assert result.exit_code == 0, result.output
        assert MetaCooler(out).metadata["qc"]["distance_threshold"] == 1000

    def test_bad_qc_threshold(self, runner, tmp_path, cools):
        result = run(
            runner, "collect", *cools, "--qc-threshold", "far", "-o", tmp_path / "o"
        )
        assert result.exit_code == 2
        assert "--qc-threshold" in _flat(result.output)

    def test_chunksize_must_be_positive(self, runner, tmp_path, cools):
        result = run(
            runner, "collect", *cools, "--chunksize", "0", "-o", tmp_path / "o"
        )
        assert result.exit_code == 2

    def test_quiet(self, runner, tmp_path, cools):
        result = run(runner, "collect", *cools, "-q", "-o", tmp_path / "q.metacool")
        assert result.exit_code == 0
        assert "Created metacool file" not in result.output
