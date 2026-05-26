"""Tests for metacooler.ops.show."""

from __future__ import annotations

import pytest

from metacooler.ops.show import show


class TestShow:
    def test_show_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="genome", out=out)
        assert out.exists()

    def test_show_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="contig", out=out)
        assert out.exists()

    def test_show_bin_requires_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="--genome, --genomes, or --range"):
            show(tiny_metacooler, level="bin", out=out)

    def test_show_bin_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="bin", genome="MAG_1", out=out)
        assert out.exists()

    def test_show_invalid_level(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="Unknown level"):
            show(tiny_metacooler, level="invalid", out=out)  # type: ignore

    def test_show_invalid_transform(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="Unknown transform"):
            show(tiny_metacooler, level="genome", transform="invalid", out=out)  # type: ignore

    def test_show_range1_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="bin", range1="contig_A", out=out)
        assert out.exists()

    def test_show_range1_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="bin", range1="MAG_1", out=out)
        assert out.exists()

    def test_show_range1_with_coords(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="bin", range1="contig_A:0-500", out=out)
        assert out.exists()

    def test_show_range1_range2_inter_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        show(tiny_metacooler, level="bin", range1="MAG_1", range2="MAG_2", out=out)
        assert out.exists()
