"""Tests for metacooler.ops.aggregate."""

from __future__ import annotations

import cooler

from metacooler.api import MetaCooler
from metacooler.ops.aggregate import aggregate


class TestAggregate:
    def test_aggregate_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg_contig.cool"
        aggregate(tiny_metacooler, out, level="contig")
        assert cooler.fileops.is_cooler(str(out))
        clr = cooler.Cooler(str(out))
        assert len(clr.chroms()) == 4

    def test_aggregate_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg_genome.cool"
        aggregate(tiny_metacooler, out, level="genome")
        assert cooler.fileops.is_cooler(str(out))
        clr = cooler.Cooler(str(out))
        assert len(clr.chroms()) == 3

    def test_pixel_sum_preserved(self, tiny_metacooler, tmp_path):
        mc = MetaCooler(tiny_metacooler)
        original_sum = mc.pixels()["count"].sum()

        out = tmp_path / "agg.cool"
        aggregate(tiny_metacooler, out, level="contig")
        clr = cooler.Cooler(str(out))
        agg_sum = clr.pixels()[:]["count"].sum()
        assert agg_sum == original_sum
