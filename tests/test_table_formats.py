"""--contigs-format / --genomes-format across commands, and create's no-match check."""

from __future__ import annotations

from pathlib import Path

import pytest
from _cools import BSUB, ECOLI, make_cool
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.collect import collect
from metacooler.io.create import from_cool, from_pairs
from metacooler.io.reassign import reassign


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def flat(output: str) -> str:
    """A CLI error message with rich's box drawing and line wrapping undone."""
    return " ".join(output.replace("│", " ").split())


@pytest.fixture
def checkm(tmp_path: Path) -> Path:
    return write(
        tmp_path / "checkm2.tsv",
        "Name\tCompleteness\nMAG_1\t91.5\nMAG_2\t70\n",
    )


@pytest.fixture
def metator(tmp_path: Path) -> Path:
    return write(
        tmp_path / "contig_data_final.txt",
        "ID\tName\tGC_content\tFinal_bin\n"
        "1\tcontig_A\t0.5\tMAG_1\n2\tcontig_B\t0.4\tMAG_1\n3\tcontig_C\t0.6\tMAG_2\n",
    )


class TestCreate:
    def test_create_cool_genomes_format(
        self, tmp_path, tiny_cool, tiny_contigs_table, checkm
    ):
        out = tmp_path / "o.metacool"
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "--genomes",
                str(checkm),
                "--genomes-format",
                "Name",
                "--no-qc",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        genomes = MetaCooler(out).genomes.set_index("name")
        assert genomes.loc["MAG_1", "Completeness"] == 91.5

    def test_create_pairs_by_names(
        self, tmp_path, tiny_pairs, tiny_bins_bed, metator, checkm
    ):
        out = tmp_path / "o.metacool"
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "pairs",
                "--pairs",
                str(tiny_pairs),
                "--bins",
                str(tiny_bins_bed),
                "--contigs",
                str(metator),
                "--contigs-format",
                "Name,Final_bin",
                "--genomes",
                str(checkm),
                "--genomes-format",
                "1",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        contigs = mc.contigs.set_index("name")
        assert contigs.loc["contig_A", "GC_content"] == 0.5
        assert set(mc.genomes["name"]) == {"__unbinned__", "MAG_1", "MAG_2"}

    @pytest.mark.parametrize(
        ("spec", "message"),
        [
            ("0,1", ">= 1"),
            ("Name,Final_bin,noheader", "no header line"),
            ("contig_name", "Unknown --contigs-format 'contig_name'"),
            ("1,2,3", "names 3 columns"),
        ],
    )
    def test_bad_specs_are_refused_while_parsing_options(
        self, tmp_path, tiny_cool, tiny_contigs_table, spec, message
    ):
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "--contigs-format",
                spec,
                "-o",
                str(tmp_path / "o.metacool"),
            ],
        )
        assert result.exit_code == 2, result.output
        assert message in flat(result.output)

    def test_bad_genomes_spec(self, tmp_path, tiny_cool, tiny_contigs_table, checkm):
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "--genomes",
                str(checkm),
                "--genomes-format",
                "1,2",
                "-o",
                str(tmp_path / "o.metacool"),
            ],
        )
        assert result.exit_code == 2
        assert "takes 1: the genome column" in flat(result.output)


class TestNoMatch:
    def test_from_cool_refuses_a_table_naming_no_contig(self, tmp_path, tiny_cool):
        table = write(tmp_path / "c.tsv", "ID\tcontig\tgenome\n1\tcontig_A\tMAG_1\n")
        with pytest.raises(ValueError, match=r"names none of the assembly's 4.*'1,2'"):
            from_cool(
                tiny_cool,
                tmp_path / "o.metacool",
                contig_genome_map=table,
                contig_table_format="1,2",
            )
        assert not (tmp_path / "o.metacool").exists()

    def test_unknown_rows_are_reported(self, tmp_path, tiny_cool, caplog):
        table = write(
            tmp_path / "c.tsv",
            "contig_name\tgenome_name\ncontig_A\tMAG_1\nnot_here\tMAG_1\n",
        )
        with caplog.at_level("WARNING", logger="metacooler"):
            from_cool(
                tiny_cool, tmp_path / "o.metacool", contig_genome_map=table, qc=False
            )
        assert "1 of its 2 contig(s) are not in the assembly" in caplog.text
        assert "'not_here'" in caplog.text

    def test_an_empty_table_unbins_everything(self, tmp_path, tiny_cool, caplog):
        table = write(tmp_path / "c.tsv", "contig_name\tgenome_name\n")
        out = tmp_path / "o.metacool"
        with caplog.at_level("WARNING", logger="metacooler"):
            from_cool(tiny_cool, out, contig_genome_map=table, qc=False)
        assert "has no row" in caplog.text
        assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__"]

    def test_a_dict_naming_no_contig(
        self, tmp_path, tiny_pairs, tiny_bins_bed, tiny_fasta
    ):
        with pytest.raises(ValueError, match="The contig map names none"):
            from_pairs(
                tmp_path / "o.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map={"elsewhere": "MAG_1"},
                fasta=tiny_fasta,
                ps=False,
            )


class TestCollect:
    def test_contig_and_genome_columns_by_name(self, tmp_path):
        ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
        bsub = make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)
        contigs = write(
            tmp_path / "depth.tsv",
            "contigName\tdepth\nE_coli__plasmid1\t310\nB_subtilis__chromosome\t98\n",
        )
        genomes = write(tmp_path / "g.tsv", "x\tName\nhigh\tE_coli\n")
        out = tmp_path / "c.metacool"
        result = CliRunner().invoke(
            cli,
            [
                "collect",
                str(ecoli),
                str(bsub),
                "--contigs",
                str(contigs),
                "--contigs-format",
                "contigName",
                "--genomes",
                str(genomes),
                "--genomes-format",
                "Name",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert mc.contigs.set_index("name").loc["E_coli__plasmid1", "depth"] == 310
        assert mc.genomes.set_index("name").loc["E_coli", "x"] == "high"

    def test_two_column_spec_checks_the_genome_column(self, tmp_path):
        ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
        bsub = make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)
        table = write(
            tmp_path / "t.tsv",
            "contig\tspecies\tgc\nchromosome\tE_coli\t0.51\nchromosome\tB_subtilis\t0.44\n",
        )
        out = tmp_path / "c.metacool"
        collect(
            [ecoli, bsub],
            out,
            contigs_table=table,
            contigs_table_format="contig,species",
        )
        contigs = MetaCooler(out).contigs.set_index("name")
        assert contigs.loc["E_coli__chromosome", "gc"] == 0.51
        assert contigs.loc["B_subtilis__chromosome", "gc"] == 0.44

    def test_headerless_metadata_table(self, tmp_path):
        ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
        table = write(tmp_path / "t.tsv", "plasmid1\t310\n")
        out = tmp_path / "c.metacool"
        collect([ecoli], out, contigs_table=table, contigs_table_format="1,noheader")
        assert (
            MetaCooler(out).contigs.set_index("name").loc["E_coli__plasmid1", "col2"]
            == 310
        )

    def test_missing_contig_column(self, tmp_path):
        ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
        table = write(tmp_path / "t.tsv", "name\tdepth\nplasmid1\t1\n")
        with pytest.raises(ValueError, match="no column named 'contig_name'"):
            collect([ecoli], tmp_path / "c.metacool", contigs_table=table)


def test_reassign_genomes_format(tmp_path, tiny_metacooler, checkm):
    out = tmp_path / "r.metacool"
    reassign(
        tiny_metacooler,
        {"contig_A": "MAG_1"},
        out,
        genomes_table=checkm,
        genomes_table_format="Name",
        qc=False,
    )
    assert (
        MetaCooler(out).genomes.set_index("name").loc["MAG_2", "Completeness"] == 70.0
    )


def test_reassign_cli_genomes_format(
    tmp_path, tiny_metacooler, tiny_contigs_table, checkm
):
    out = tmp_path / "r.metacool"
    result = CliRunner().invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(tiny_contigs_table),
            "--genomes",
            str(checkm),
            "--genomes-format",
            "Name",
            "--no-qc",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Completeness" in MetaCooler(out).genomes.columns


def test_spec_callback_passes_none_through():
    from metacooler.cli._options import _column_spec_callback

    assert _column_spec_callback(("genome",))(None, None, None) is None
