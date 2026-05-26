"""Tests for metacooler.ops.range."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.api import MetaCooler
from metacooler.ops.range import parse_metacooler_range


class TestParseMetacoolerRange:
    def test_genome_range_returns_bin_interval(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        i0, i1 = parse_metacooler_range(mc, "MAG_1")
        assert i0 < i1

    def test_genome_range_matches_genome_extent(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        assert parse_metacooler_range(mc, "MAG_1") == mc.genome_extent("MAG_1")

    def test_contig_range_full(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        i0, i1 = parse_metacooler_range(mc, "contig_A")
        # 1000 bp contig with 100 bp bins → 10 bins
        assert i1 - i0 == 10
        assert (i0, i1) == mc.contig_extent("contig_A")

    def test_contig_range_with_coords(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        i0_full, _ = parse_metacooler_range(mc, "contig_A")
        i0, i1 = parse_metacooler_range(mc, "contig_A:0-500")
        # 500 bp with 100 bp bins → 5 bins
        assert i1 - i0 == 5
        assert i0 == i0_full

    def test_unknown_name_raises(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(ValueError, match="not a known genome or contig"):
            parse_metacooler_range(mc, "nonexistent")

    def test_genome_takes_priority_over_contig_if_same_name(
        self, tiny_metacooler: Path
    ):
        """If a name appears in genomes, it is treated as a genome."""
        mc = MetaCooler(tiny_metacooler)
        genome_names = mc.genomes["name"].tolist()
        # Verify genome lookup path is used (result == genome_extent)
        g = genome_names[0]
        result = parse_metacooler_range(mc, g)
        assert result == mc.genome_extent(g)
