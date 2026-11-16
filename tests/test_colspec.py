"""Tests for the column-spec grammar shared by --contigs-format/--genomes-format."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.io.colspec import (
    ASSIGNMENT_KEYS,
    CONTIGS_FORMAT_ALIASES,
    GENOMES_KEYS,
    alias_listing,
    parse_column_spec,
    resolve_column_spec,
)
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.create import _parse_contig_map
from metacooler.io.genomes import parse_genomes_table


def assignment(spec: str):
    return parse_column_spec(
        spec, keys=ASSIGNMENT_KEYS, aliases=CONTIGS_FORMAT_ALIASES, option="--x"
    )


def metadata(spec: str):
    return parse_column_spec(
        spec,
        keys=ASSIGNMENT_KEYS,
        required=1,
        aliases=CONTIGS_FORMAT_ALIASES,
        option="--x",
    )


def genomes(spec: str):
    return parse_column_spec(spec, keys=GENOMES_KEYS, option="--x")


class TestGrammar:
    def test_positions_names_and_mixtures(self):
        assert assignment("2,5").columns == (2, 5)
        assert assignment("Name,Final_bin").columns == ("Name", "Final_bin")
        assert assignment("1,Final_bin").columns == (1, "Final_bin")
        assert assignment(" 1 , 2 ").header

    def test_noheader(self):
        spec = assignment("1,2,noheader")
        assert spec.columns == (1, 2)
        assert not spec.header
        assert not genomes("3,noheader").header

    def test_aliases_expand_and_show_their_expansion(self):
        spec = assignment("metator")
        assert spec.columns == ("Name", "Final_bin")
        assert spec.text == "metator (= Name,Final_bin)"
        assert assignment("dastool").columns == (1, 2)
        assert not assignment("binette").header
        assert assignment("semibin").columns == ("contig", "bin")
        assert assignment("vamb").columns == ("contigname", "clustername")
        assert assignment("contigs").columns == ("contig_name", "genome_name")

    def test_a_name_needs_a_header(self):
        with pytest.raises(ValueError, match="no header line to find it in"):
            assignment("Name,Final_bin,noheader")

    def test_positions_are_one_based(self):
        with pytest.raises(ValueError, match=">= 1"):
            assignment("0,1")

    def test_empty_column(self):
        with pytest.raises(ValueError, match="empty column"):
            assignment("1,,2")

    def test_one_column_twice(self):
        with pytest.raises(ValueError, match="same column"):
            assignment("contig,contig")

    def test_assignment_needs_two_columns(self):
        with pytest.raises(ValueError, match="Unknown --x 'contig_name'"):
            assignment("contig_name")
        with pytest.raises(ValueError, match="names 3 columns, but this table takes 2"):
            assignment("1,2,3")

    def test_metadata_takes_one_column_or_two(self):
        assert metadata("contigName").columns == ("contigName",)
        assert metadata("1,noheader").columns == (1,)
        assert metadata("contig_name,genome_name").columns == (
            "contig_name",
            "genome_name",
        )
        assert metadata("metator").columns == ("Name", "Final_bin")
        with pytest.raises(ValueError, match="takes 1 to 2"):
            metadata("1,2,3")

    def test_genomes_take_one_column(self):
        assert genomes("user_genome").columns == ("user_genome",)
        with pytest.raises(ValueError, match="takes 1: the genome column"):
            genomes("1,2")

    def test_alias_listing(self):
        assert alias_listing({"a": "1,2"}) == "a (1,2)"


class TestResolution:
    def test_names_are_found_wherever_they_are(self):
        spec = assignment("contig_name,genome_name")
        labels = ["ID", "genome_name", "x", "contig_name"]
        assert resolve_column_spec(
            spec, labels, kind="T", path="t.tsv", keys=ASSIGNMENT_KEYS
        ) == [3, 1]

    def test_missing_name_lists_the_header(self):
        spec = assignment("Name,Final_bin")
        with pytest.raises(ValueError, match="Its header is: Name, bin") as exc:
            resolve_column_spec(
                spec, ["Name", "bin"], kind="T", path="t.tsv", keys=ASSIGNMENT_KEYS
            )
        assert "holds the genome names" in str(exc.value)

    def test_too_narrow_a_table(self):
        with pytest.raises(ValueError, match="needs at least 2 columns"):
            resolve_column_spec(
                assignment("1,2"), ["a"], kind="T", path="t", keys=ASSIGNMENT_KEYS
            )

    def test_position_out_of_range(self):
        with pytest.raises(ValueError, match="requires at least 5 columns"):
            resolve_column_spec(
                assignment("1,5"), ["a", "b"], kind="T", path="t", keys=ASSIGNMENT_KEYS
            )

    def test_name_and_position_on_one_column(self):
        with pytest.raises(ValueError, match="same column"):
            resolve_column_spec(
                assignment("1,a"), ["a", "b"], kind="T", path="t", keys=ASSIGNMENT_KEYS
            )


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


class TestContigTables:
    def test_default_finds_its_columns_by_name(self, tmp_path):
        """An ID column first used to shift a positional read silently."""
        path = write(
            tmp_path / "t.tsv",
            "ID\tcontig_name\tgenome_name\tgc\n1\tc1\tMAG_1\t0.5\n",
        )
        mapping, df = parse_contigs_table(path)
        assert mapping == {"c1": "MAG_1"}
        assert list(df.columns) == ["ID", "contig_name", "genome_id", "gc"]

    def test_metator_alias_is_its_spec(self, tmp_path):
        path = write(tmp_path / "t.tsv", "ID\tName\tFinal_bin\n1\tc1\tMAG_1\n")
        assert (
            parse_contigs_table(path, "metator")[0]
            == parse_contigs_table(path, "Name,Final_bin")[0]
            == {"c1": "MAG_1"}
        )

    def test_contig_genome_keeps_its_extra_columns(self, tmp_path):
        path = write(tmp_path / "t.tsv", "c1\tMAG_1\t1000\nc2\tMAG_2\t800\n")
        mapping, df = parse_contigs_table(path, "contig_genome")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert list(df.columns) == ["contig_name", "genome_id", "col3"]

    @pytest.mark.parametrize("alias", ["dastool", "binette"])
    def test_headerless_contig_to_bin_tables(self, tmp_path, alias):
        path = write(tmp_path / "t.tsv", "c1\tbin.01\nc2\tbin.02\n")
        assert parse_contigs_table(path, alias)[0] == {"c1": "bin.01", "c2": "bin.02"}

    def test_semibin(self, tmp_path):
        path = write(tmp_path / "t.tsv", "contig\tbin\nc1\t0\nc2\t-1\n")
        assert parse_contigs_table(path, "semibin")[0] == {"c1": "0", "c2": "-1"}
        # -1 is SemiBin's "unbinned".
        cg_map, _ = _parse_contig_map(path, "semibin")
        assert cg_map == {"c1": "0", "c2": "__unbinned__"}

    def test_vamb_puts_the_cluster_first(self, tmp_path):
        path = write(
            tmp_path / "t.tsv", "clustername\tcontigname\nvae_1\tc1\nvae_2\tc2\n"
        )
        assert parse_contigs_table(path, "vamb")[0] == {"c1": "vae_1", "c2": "vae_2"}

    def test_wrong_alias_names_the_aliases(self, tmp_path):
        path = write(tmp_path / "t.tsv", "c1\tMAG_1\n")
        with pytest.raises(ValueError, match="not an alias .*dastool"):
            parse_contigs_table(path, "das_tool")


class TestGenomesTables:
    def test_by_name_anywhere(self, tmp_path):
        path = write(
            tmp_path / "g.tsv", "Completeness\tName\tContamination\n90\tMAG_1\t1\n"
        )
        df = parse_genomes_table(path, "Name")
        assert list(df.columns) == ["Completeness", "genome_id", "Contamination"]
        assert df["genome_id"].tolist() == ["MAG_1"]

    def test_by_position(self, tmp_path):
        path = write(tmp_path / "g.tsv", "a\tb\n1\tMAG_1\n")
        assert parse_genomes_table(path, "2")["genome_id"].tolist() == ["MAG_1"]

    def test_without_header(self, tmp_path):
        path = write(tmp_path / "g.tsv", "MAG_1\thigh\nMAG_2\tlow\n")
        df = parse_genomes_table(path, "1,noheader")
        assert list(df.columns) == ["genome_id", "col2"]
        assert df["col2"].tolist() == ["high", "low"]

    def test_default_still_needs_genome_name(self, tmp_path):
        path = write(tmp_path / "g.tsv", "Name\tX\nMAG_1\t1\n")
        with pytest.raises(ValueError, match="no column named 'genome_name'"):
            parse_genomes_table(path)

    def test_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_genomes_table(tmp_path / "nope.tsv")
