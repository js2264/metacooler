"""Tests for metacooler.io.dump."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.io.dump import dump_table


class TestDumpTable:
    def test_dump_genomes(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "genomes", out=out)
        result = out.read_text()
        assert "__unbinned__" in result
        assert "MAG_1" in result

    def test_dump_contigs(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "contigs", header=True, out=out)
        result = out.read_text()
        assert "contig_A" in result
        assert "genome_id" in result

    def test_dump_bins(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "bins", header=True, out=out)
        result = out.read_text()
        assert "chrom" in result

    def test_dump_pixels(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "pixels", header=True, out=out)
        result = out.read_text()
        assert "bin1_id" in result

    def test_dump_indexes(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "indexes", header=True, out=out)
        result = out.read_text()
        assert "index" in result

    def test_dump_genomes_filtered(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "genomes", genomes=["MAG_1"], header=True, out=out)
        result = out.read_text()
        # genomes table doesn't have per-genome filtering, just returns all
        assert "name" in result

    def test_dump_no_header(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "genomes", out=out)  # default is now no header
        result = out.read_text()
        assert "name" not in result  # no header

    def test_dump_pixels_range1_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "pixels", range1="contig_A", header=True, out=out)
        result = out.read_text()
        assert "bin1_id" in result

    def test_dump_pixels_range1_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(tiny_metacooler, "pixels", range1="MAG_1", header=True, out=out)
        assert out.exists()
        result = out.read_text()
        assert "bin1_id" in result

    def test_dump_pixels_range1_with_coords(self, tiny_metacooler, tmp_path):
        out_full = tmp_path / "full.tsv"
        out_sub = tmp_path / "sub.tsv"
        dump_table(
            tiny_metacooler, "pixels", range1="contig_A", header=True, out=out_full
        )
        dump_table(
            tiny_metacooler, "pixels", range1="contig_A:0-500", header=True, out=out_sub
        )
        # Sub-range should have fewer or equal rows
        full_lines = out_full.read_text().strip().splitlines()
        sub_lines = out_sub.read_text().strip().splitlines()
        assert len(sub_lines) <= len(full_lines)

    def test_dump_pixels_range1_range2(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        dump_table(
            tiny_metacooler,
            "pixels",
            range1="MAG_1",
            range2="MAG_2",
            header=True,
            out=out,
        )
        assert out.exists()

    def test_dump_pixels_range_unknown_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        with pytest.raises(ValueError, match="not a known genome or contig"):
            dump_table(tiny_metacooler, "pixels", range1="nonexistent", out=out)

    def test_dump_chroms(self, tiny_metacooler, tmp_path):
        out = tmp_path / "chroms.tsv"
        dump_table(tiny_metacooler, "chroms", header=True, out=out)
        result = out.read_text()
        lines = [ln for ln in result.strip().split("\n") if ln]
        # header + 4 contigs
        assert len(lines) == 5
        assert "name" in lines[0]
        assert "length" in lines[0]

    def test_dump_chroms_no_header_count(self, tiny_metacooler, tmp_path):
        out = tmp_path / "chroms.tsv"
        dump_table(tiny_metacooler, "chroms", out=out)
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) == 4

    def test_dump_bins_filtered_by_genome(self, tiny_metacooler, tmp_path):
        out_full = tmp_path / "bins_full.tsv"
        out_filtered = tmp_path / "bins_filtered.tsv"
        dump_table(tiny_metacooler, "bins", out=out_full)
        dump_table(tiny_metacooler, "bins", genomes=["MAG_1"], out=out_filtered)
        full_lines = [ln for ln in out_full.read_text().strip().split("\n") if ln]
        filt_lines = [ln for ln in out_filtered.read_text().strip().split("\n") if ln]
        assert len(filt_lines) < len(full_lines)

    def test_dump_bins_filtered_by_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "bins_contig.tsv"
        dump_table(tiny_metacooler, "bins", contigs=["contig_A"], out=out)
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        # contig_A has 1000/100=10 bins
        assert len(lines) == 10

    def test_dump_contigs_filtered_by_genome(self, tiny_metacooler, tmp_path):
        out_full = tmp_path / "contigs_full.tsv"
        out_filt = tmp_path / "contigs_filt.tsv"
        dump_table(tiny_metacooler, "contigs", out=out_full)
        dump_table(tiny_metacooler, "contigs", genomes=["MAG_1"], out=out_filt)
        full_lines = [ln for ln in out_full.read_text().strip().split("\n") if ln]
        filt_lines = [ln for ln in out_filt.read_text().strip().split("\n") if ln]
        # MAG_1 has 2 contigs (contig_A, contig_B)
        assert len(filt_lines) == 2
        assert len(full_lines) == 4

    def test_dump_pixels_filtered_by_genome(self, tiny_metacooler, tmp_path):
        out_full = tmp_path / "pix_full.tsv"
        out_filt = tmp_path / "pix_filt.tsv"
        dump_table(tiny_metacooler, "pixels", out=out_full)
        dump_table(tiny_metacooler, "pixels", genomes=["MAG_1"], out=out_filt)
        full_count = len([ln for ln in out_full.read_text().strip().split("\n") if ln])
        filt_count = len([ln for ln in out_filt.read_text().strip().split("\n") if ln])
        assert filt_count < full_count

    def test_dump_unknown_genome_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        with pytest.raises(KeyError, match="not found"):
            dump_table(tiny_metacooler, "bins", genomes=["NONEXISTENT"], out=out)

    def test_dump_unknown_contig_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        with pytest.raises(KeyError, match="not found"):
            dump_table(tiny_metacooler, "bins", contigs=["NONEXISTENT"], out=out)

    def test_dump_all_requires_prefix(self, tiny_metacooler, tmp_path):
        with pytest.raises(ValueError, match="--prefix"):
            dump_table(tiny_metacooler, "all", out=tmp_path / "out")

    def test_dump_all_creates_files(self, tiny_metacooler, tmp_path):
        prefix = str(tmp_path / "dump_all")
        dump_table(tiny_metacooler, "all", prefix=prefix)
        for tbl in ("genomes", "contigs", "bins", "pixels"):
            f = Path(f"{prefix}.{tbl}.tsv")
            assert f.exists(), f"Missing: {f}"

    def test_dump_pixels_wide_format(self, tiny_metacooler, tmp_path):
        out = tmp_path / "wide.tsv"
        dump_table(tiny_metacooler, "pixels", pixels_format="wide", out=out)
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        # Wide format: one row per bin row
        assert len(lines) > 0
        # Each row should have as many columns as bins
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        n_bins = len(mc.bins())
        first_row_cols = len(lines[0].split("\t"))
        assert first_row_cols == n_bins

    def test_dump_stdout(self, tiny_metacooler, capsys):
        dump_table(tiny_metacooler, "genomes")
        captured = capsys.readouterr()
        assert len(captured.out.strip()) > 0
        assert "__unbinned__" in captured.out

    def test_dump_contigs_core_columns_first(self, tiny_metacooler, tmp_path):
        """Core columns (name, length, genome_id) appear before extra columns."""
        out = tmp_path / "contigs.tsv"
        dump_table(tiny_metacooler, "contigs", header=True, out=out)
        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "name"
        assert header[1] == "length"
        assert header[2] == "genome_id"

    def test_dump_genomes_core_columns_first(self, tiny_metacooler, tmp_path):
        """Core columns (name, length, n_contigs) appear before extra columns."""
        out = tmp_path / "genomes.tsv"
        dump_table(tiny_metacooler, "genomes", header=True, out=out)
        lines = out.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header[0] == "name"
        assert header[1] == "length"
        assert header[2] == "n_contigs"

    def test_dump_unknown_table_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.tsv"
        with pytest.raises(ValueError, match="Unknown table"):
            dump_table(tiny_metacooler, "bogus_table", out=out)  # ty:ignore

    def test_dump_pixels_wide_invalid_weight_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "wide.tsv"
        with pytest.raises(ValueError, match="Weight column"):
            dump_table(
                tiny_metacooler,
                "pixels",
                pixels_format="wide",
                weight="nonexistent_col",
                out=out,
            )

    def test_dump_pixels_wide_too_many_bins_raises(self, tiny_metacooler, tmp_path):
        from unittest.mock import patch

        out = tmp_path / "wide.tsv"
        # Patch the limit to 0 to force the "too many bins" error path
        with (
            patch("metacooler.io.dump.DENSE_MATRIX_BIN_LIMIT", 0),
            pytest.raises(ValueError, match="Wide format refuses"),
        ):
            dump_table(tiny_metacooler, "pixels", pixels_format="wide", out=out)
