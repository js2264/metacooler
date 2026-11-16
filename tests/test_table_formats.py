"""--contigs-format / --genomes-format across commands, and create's no-match check."""

from __future__ import annotations

from pathlib import Path

import pytest
from _cli_output import flat
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


class TestNoBinningYet:
    """--contigs-format CONTIG,NA: a map with every contig in __unbinned__."""

    @pytest.fixture
    def contig_list(self, tmp_path):
        return write(tmp_path / "contigs.txt", "contig_A\ncontig_B\ncontig_C\n")

    def test_create_cool(self, tmp_path, tiny_cool, contig_list, caplog):
        out = tmp_path / "o.metacool"
        with caplog.at_level("INFO", logger="metacooler"):
            from_cool(
                tiny_cool,
                out,
                contig_genome_map=contig_list,
                contig_table_format="1,NA,noheader",
                qc=False,
            )
        assert "nothing is binned yet" in caplog.text
        genomes = MetaCooler(out).genomes
        assert genomes["name"].tolist() == ["__unbinned__"]
        assert genomes["n_contigs"].tolist() == [4]

    def test_metadata_columns_are_kept(self, tmp_path, tiny_cool):
        table = write(
            tmp_path / "c.tsv", "contig_name\tdepth\ncontig_A\t12\ncontig_D\t3\n"
        )
        out = tmp_path / "o.metacool"
        from_cool(
            tiny_cool,
            out,
            contig_genome_map=table,
            contig_table_format="contig_name,NA",
            qc=False,
        )
        contigs = MetaCooler(out).contigs.set_index("name")
        assert contigs.loc["contig_A", "depth"] == 12
        assert set(MetaCooler(out).genomes["name"]) == {"__unbinned__"}

    def test_an_empty_list_is_refused(self, tmp_path, tiny_cool):
        """No binning is asked for with NA and a list of contigs, not with nothing."""
        table = write(tmp_path / "c.tsv", "contig_name\n")
        with pytest.raises(ValueError, match="There is no row"):
            from_cool(
                tiny_cool,
                tmp_path / "o.metacool",
                contig_genome_map=table,
                contig_table_format="contig_name,NA",
            )

    def test_a_wrong_contig_column_is_still_refused(self, tmp_path, tiny_cool):
        table = write(tmp_path / "c.tsv", "1\tcontig_A\n")
        with pytest.raises(ValueError, match="None of the 1 contig name"):
            from_cool(
                tiny_cool,
                tmp_path / "o.metacool",
                contig_genome_map=table,
                contig_table_format="1,NA,noheader",
            )

    def test_genomes_table_is_refused(self, tmp_path, tiny_cool, contig_list, checkm):
        with pytest.raises(ValueError, match="no genome for it to describe"):
            from_cool(
                tiny_cool,
                tmp_path / "o.metacool",
                contig_genome_map=contig_list,
                contig_table_format="1,NA,noheader",
                genomes_table=checkm,
            )

    def test_drop_unassigned_is_refused(
        self, tmp_path, tiny_pairs, tiny_bins_bed, contig_list
    ):
        with pytest.raises(ValueError, match="would drop every contig"):
            from_pairs(
                tmp_path / "o.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=contig_list,
                contig_table_format="1,NA,noheader",
                drop_unassigned=True,
            )

    def test_create_pairs_cli(self, tmp_path, tiny_pairs, tiny_bins_bed, contig_list):
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
                str(contig_list),
                "--contigs-format",
                "1,NA,noheader",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__"]

    @pytest.mark.parametrize(
        ("option", "spec", "message"),
        [
            ("--contigs-format", "NA,1", "only goes in the genome's place"),
            ("--genomes-format", "NA", "nothing without its genome names"),
        ],
    )
    def test_misplaced_na_is_refused_while_parsing_options(
        self, tmp_path, tiny_cool, tiny_contigs_table, checkm, option, spec, message
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
                "--genomes",
                str(checkm),
                option,
                spec,
                "-o",
                str(tmp_path / "o.metacool"),
            ],
        )
        assert result.exit_code == 2
        assert message in flat(result.output)

    def test_reassign_unbins_the_listed_contigs(self, tmp_path, tiny_metacooler):
        table = write(tmp_path / "c.txt", "contig_A\n")
        out = tmp_path / "r.metacool"
        reassign(
            tiny_metacooler,
            table,
            out,
            contig_table_format="1,NA,noheader",
            layout="keep",
            qc=False,
        )
        mc = MetaCooler(out)
        genome_of = dict(
            zip(
                mc.contigs["name"],
                mc.genomes["name"].to_numpy()[mc.contigs["genome_id"].to_numpy()],
                strict=True,
            )
        )
        assert genome_of["contig_A"] == "__unbinned__"
        assert genome_of["contig_B"] == "MAG_1"

        everything = tmp_path / "all.metacool"
        reassign(
            tiny_metacooler,
            table,
            everything,
            contig_table_format="1,NA,noheader",
            missing="unbinned",
            layout="keep",
            qc=False,
        )
        assert MetaCooler(everything).genomes["name"].tolist() == ["__unbinned__"]

    def test_annotate_reads_contig_na_as_contig(self, tmp_path, tiny_metacooler):
        from metacooler.io.annotate import annotate

        table = write(tmp_path / "c.tsv", "contig_name\tdepth\ncontig_A\t5\n")
        out = tmp_path / "a.metacool"
        annotate(
            tiny_metacooler,
            out,
            contigs_table=table,
            contigs_table_format="contig_name,NA",
        )
        assert MetaCooler(out).contigs.set_index("name").loc["contig_A", "depth"] == 5


class TestNoContigTable:
    """No --contigs at all: the cooler plus the genome layer, nothing binned."""

    def test_create_pairs_cli(self, tmp_path, tiny_pairs, tiny_bins_bed):
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
                "--filter",
                "contigs.length >= 1",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        genomes = MetaCooler(out).genomes
        assert genomes["name"].tolist() == ["__unbinned__"]
        assert genomes["n_contigs"].tolist() == [4]

    def test_create_cool_keeps_the_cool_file(self, tmp_path, tiny_cool, caplog):
        import cooler
        import numpy as np

        out = tmp_path / "o.metacool"
        with caplog.at_level("INFO", logger="metacooler"):
            from_cool(tiny_cool, out, contig_ordering="cool", qc=False)
        assert "No contig table was given: nothing is binned yet" in caplog.text
        src, dst = cooler.Cooler(str(tiny_cool)), cooler.Cooler(str(out))
        assert list(src.chromnames) == list(dst.chromnames)
        cols = ["start", "end"]
        assert src.bins()[:][cols].equals(dst.bins()[:][cols])
        assert np.array_equal(src.pixels()[:].to_numpy(), dst.pixels()[:].to_numpy())
        assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__"]

    def test_from_pairs_without_a_map(self, tmp_path, tiny_pairs, tiny_bins_bed):
        out = tmp_path / "o.metacool"
        from_pairs(out, bins=tiny_bins_bed, pairs=tiny_pairs, ps=False)
        assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__"]

    def test_a_map_of_the_wrong_type(self, tmp_path, tiny_pairs, tiny_bins_bed):
        with pytest.raises(TypeError, match="a dict, a path, or None"):
            from_pairs(
                tmp_path / "o.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=42,  # ty: ignore[invalid-argument-type]
            )

    def test_genomes_and_drop_unassigned_are_refused(self, tmp_path, tiny_cool, checkm):
        with pytest.raises(ValueError, match="no contig table was given, so no"):
            from_cool(tiny_cool, tmp_path / "o.metacool", genomes_table=checkm)
        with pytest.raises(ValueError, match="would drop every contig"):
            from_cool(tiny_cool, tmp_path / "o.metacool", drop_unassigned=True)

    def test_contigs_format_without_contigs_is_refused(self, tmp_path, tiny_cool):
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs-format",
                "metator",
                "-o",
                str(tmp_path / "o.metacool"),
            ],
        )
        assert result.exit_code == 2
        assert "says how to read the --contigs table" in flat(result.output)
