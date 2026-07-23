"""Tests for binnify operation."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.api import MetaCooler
from metacooler.ops.binnify import binnify


class TestBinnify:
    def test_binnify_doubles_binsize(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        assert out.exists()
        mc2 = MetaCooler(out)
        assert mc2.shape[0] > 0
        # Should have fewer bins
        assert mc2.shape[0] < mc.shape[0]

    def test_binnify_preserves_genomes(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.genomes) == len(mc.genomes)
        assert set(mc2.genomes["name"]) == set(mc.genomes["name"])

    def test_binnify_preserves_contigs(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.contigs) == len(mc.contigs)

    def test_binnify_smaller_binsize_raises(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "fail.metacool"
        with pytest.raises(ValueError, match="must be strictly larger"):
            binnify(tiny_metacooler, out, binsize=original_binsize // 2)

    def test_binnify_same_binsize_raises(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "fail.metacool"
        with pytest.raises(ValueError, match="must be strictly larger"):
            binnify(tiny_metacooler, out, binsize=original_binsize)

    def test_binnify_has_pixels(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.pixels()) > 0

    def test_binnify_genome_level_ligates_to_one_contig_per_genome(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        mc2 = MetaCooler(out)
        # Output "contigs" are the ligated genomes: one per genome, named after it.
        assert len(mc2.contigs) == len(mc.genomes)
        assert set(mc2.contigs["name"]) == set(mc.genomes["name"])
        assert set(mc2.genomes["n_contigs"].unique()) == {1}
        # Each ligated contig length equals its genome's total contig length.
        for _, row in mc2.contigs.iterrows():
            gid = int(mc.genomes.index[mc.genomes["name"] == row["name"]][0])
            expected = int(
                mc.contigs.loc[mc.contigs["genome_id"] == gid, "length"].sum()
            )
            assert int(row["length"]) == expected

    def test_binnify_genome_level_preserves_total_counts(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        import cooler

        mc = MetaCooler(tiny_metacooler)
        src_sum = int(mc._clr.pixels()[:]["count"].sum())

        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        c = cooler.Cooler(str(out))
        assert int(c.pixels()[:]["count"].sum()) == src_sum

    def test_binnify_genome_level_is_valid_metacool(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        import cooler

        from metacooler.schema import validate

        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        validate(out)  # raises if the schema is violated
        assert cooler.fileops.is_cooler(str(out))  # still openable by cooler
