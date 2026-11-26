"""Tests for metacooler.io.mag."""

from __future__ import annotations

import pytest

from metacooler.io.contigs import parse_contigs_table


class TestParseContigTable:
    def test_contigs_with_header(self, tmp_path):
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_2\n")
        mapping, df = parse_contigs_table(path, fmt="contigs")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert list(df.columns) == ["contig_name", "genome_id"]

    def test_contigs_extra_columns(self, tmp_path):
        path = tmp_path / "contigs.tsv"
        path.write_text(
            "contig_name\tgenome_name\tgc_content\tcoverage\nc1\tMAG_1\t0.52\t45.3\nc2\tMAG_2\t0.48\t32.1\n"
        )
        mapping, df = parse_contigs_table(path, fmt="contigs")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert list(df.columns) == [
            "contig_name",
            "genome_id",
            "gc_content",
            "coverage",
        ]
        assert df["gc_content"].tolist() == ["0.52", "0.48"]

    def test_contigs_too_few_columns(self, tmp_path):
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\nc1\nc2\n")
        with pytest.raises(ValueError, match="at least 2 columns"):
            parse_contigs_table(path, fmt="contigs")

    def test_contigs_wrong_first_column(self, tmp_path):
        path = tmp_path / "contigs.tsv"
        path.write_text("name\tgenome_name\nc1\tMAG_1\n")
        with pytest.raises(ValueError, match="contig_name"):
            parse_contigs_table(path, fmt="contigs")

    def test_contigs_wrong_second_column(self, tmp_path):
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_id\nc1\tMAG_1\n")
        with pytest.raises(ValueError, match="genome_name"):
            parse_contigs_table(path, fmt="contigs")

    def test_contig_genome(self, tmp_path):
        path = tmp_path / "mag.tsv"
        path.write_text("c1\tMAG_1\nc2\tMAG_2\n")
        mapping, _ = parse_contigs_table(path, fmt="contig_genome")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}

    def test_metator(self, tmp_path):
        path = tmp_path / "mag.tsv"
        path.write_text("Name\tFinal_bin\nc1\tMAG_1\nc2\tMAG_2\n")
        mapping, _ = parse_contigs_table(path, fmt="metator")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_contigs_table("/nonexistent")

    def test_metator_missing_columns(self, tmp_path):
        path = tmp_path / "mag.tsv"
        path.write_text("A\tB\nc1\tMAG_1\n")
        with pytest.raises(ValueError, match=r"'Name'.*Name,Final_cluster"):
            parse_contigs_table(path, fmt="metator")

    def test_unknown_format(self, tmp_path):
        path = tmp_path / "mag.tsv"
        path.write_text("c1\tMAG_1\n")
        with pytest.raises(ValueError, match="Unknown"):
            parse_contigs_table(path, fmt="unknown")  # type: ignore[arg-type]


class TestParseContigsTableColSpec:
    """Tests for the column-index format (e.g. '1,2' or '1,5,noheader')."""

    def test_basic_adjacent_columns_with_header(self, tmp_path):
        """'1,2' parses col 1 as contig_name and col 2 as genome_id (with header)."""
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_2\n")
        mapping, df = parse_contigs_table(path, fmt="1,2")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert "contig_name" in df.columns
        assert "genome_id" in df.columns

    def test_nonadjacent_columns_with_header(self, tmp_path):
        """'1,3' picks columns 1 and 3; column 2 is preserved in DataFrame."""
        path = tmp_path / "contigs.tsv"
        path.write_text("contig\textra\tgenome\nc1\tX\tMAG_1\nc2\tY\tMAG_2\n")
        mapping, df = parse_contigs_table(path, fmt="1,3")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert "contig_name" in df.columns
        assert "genome_id" in df.columns
        assert df["extra"].tolist() == ["X", "Y"]

    def test_header_row_not_included_in_data(self, tmp_path):
        """With header, the column-name row is NOT counted as a data row."""
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_name\nc1\tMAG_1\n")
        mapping, df = parse_contigs_table(path, fmt="1,2")
        assert "contig_name" not in mapping
        assert len(df) == 1

    def test_noheader_includes_first_row_as_data(self, tmp_path):
        """'1,2,noheader' treats the first row as data, not a header."""
        path = tmp_path / "contigs.tsv"
        path.write_text("c1\tMAG_1\nc2\tMAG_2\nc3\tMAG_3\n")
        mapping, df = parse_contigs_table(path, fmt="1,2,noheader")
        assert len(df) == 3
        assert "c1" in mapping
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2", "c3": "MAG_3"}

    def test_nonadjacent_noheader(self, tmp_path):
        """'1,3,noheader' picks columns 1 and 3 from a header-less file."""
        path = tmp_path / "contigs.tsv"
        path.write_text("c1\tX\tMAG_1\nc2\tY\tMAG_2\n")
        mapping, df = parse_contigs_table(path, fmt="1,3,noheader")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert "contig_name" in df.columns
        assert "genome_id" in df.columns

    def test_genome_column_before_contig_column(self, tmp_path):
        """'2,1' treats col 2 as contig_name and col 1 as genome_id."""
        path = tmp_path / "contigs.tsv"
        path.write_text("genome\tcontig\nMAG_1\tc1\nMAG_2\tc2\n")
        mapping, _ = parse_contigs_table(path, fmt="2,1")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}

    def test_zero_index_raises(self, tmp_path):
        """Column index 0 is not valid (must be >= 1)."""
        path = tmp_path / "contigs.tsv"
        path.write_text("c1\tMAG_1\n")
        with pytest.raises(ValueError, match=">= 1"):
            parse_contigs_table(path, fmt="0,1")

    def test_column_out_of_range_raises(self, tmp_path):
        """Requesting column 5 in a 2-column file raises ValueError."""
        path = tmp_path / "contigs.tsv"
        path.write_text("contig\tgenome\nc1\tMAG_1\n")
        with pytest.raises(ValueError, match="5 columns"):
            parse_contigs_table(path, fmt="1,5")

    def test_whitespace_around_indices_is_tolerated(self, tmp_path):
        """' 1 , 2 ' with extra whitespace still parses correctly."""
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_name\nc1\tMAG_1\n")
        mapping, _ = parse_contigs_table(path, fmt=" 1 , 2 ")
        assert mapping == {"c1": "MAG_1"}

    def test_noheader_with_whitespace_flag(self, tmp_path):
        """' 1 , 2 , noheader ' with extra whitespace around noheader flag."""
        path = tmp_path / "contigs.tsv"
        path.write_text("c1\tMAG_1\nc2\tMAG_2\n")
        mapping, _ = parse_contigs_table(path, fmt=" 1 , 2 , noheader ")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}

    def test_full_dataframe_preserved(self, tmp_path):
        """All columns are returned in the DataFrame, not just the two selected."""
        path = tmp_path / "contigs.tsv"
        path.write_text(
            "contig\textra1\textra2\tgenome\nc1\ta\tb\tMAG_1\nc2\tc\td\tMAG_2\n"
        )
        _, df = parse_contigs_table(path, fmt="1,4")
        assert df.shape == (2, 4)
        assert df["extra1"].tolist() == ["a", "c"]
        assert df["extra2"].tolist() == ["b", "d"]


class TestParseContigsTableEdgeCases:
    """Edge cases not covered by the main TestParseContigTable suite."""

    def test_contigs_empty_mapping_when_no_data_rows(self, tmp_path):
        """A contigs file with only a header (no data rows) gives empty mapping."""
        path = tmp_path / "empty.tsv"
        path.write_text("contig_name\tgenome_name\n")
        mapping, df = parse_contigs_table(path, fmt="contigs")
        assert mapping == {}
        assert len(df) == 0

    def test_contig_genome_empty_raises(self, tmp_path):
        """An empty (0-byte) contig table raises a clear, branded error."""
        path = tmp_path / "empty.tsv"
        path.write_text("")
        with pytest.raises(ValueError, match="Contig table is empty"):
            parse_contigs_table(path, fmt="contig_genome")

    def test_metator_preserves_extra_columns(self, tmp_path):
        """metator format keeps columns beyond Name and Final_bin."""
        path = tmp_path / "metator.tsv"
        path.write_text("Name\tFinal_bin\tgc\nc1\tMAG_1\t0.5\nc2\tMAG_2\t0.6\n")
        mapping, df = parse_contigs_table(path, fmt="metator")
        assert mapping == {"c1": "MAG_1", "c2": "MAG_2"}
        assert "gc" in df.columns
        assert df["gc"].tolist() == ["0.5", "0.6"]

    def test_contig_genome_preserves_duplicates(self, tmp_path):
        """contig_genome format: last occurrence wins if duplicates exist."""
        path = tmp_path / "dup.tsv"
        path.write_text("c1\tMAG_1\nc1\tMAG_2\n")
        mapping, _ = parse_contigs_table(path, fmt="contig_genome")
        # dict(zip) keeps last occurrence
        assert mapping["c1"] == "MAG_2"

    def test_contigs_keeps_default_na_false(self, tmp_path):
        """NA-like strings ('NA', 'nan') should NOT be converted to NaN."""
        path = tmp_path / "na.tsv"
        path.write_text("contig_name\tgenome_name\nc_na\tNA\nc_nan\tnan\n")
        mapping, _ = parse_contigs_table(path, fmt="contigs")
        assert mapping["c_na"] == "NA"
        assert mapping["c_nan"] == "nan"
