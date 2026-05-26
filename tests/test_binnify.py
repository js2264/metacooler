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
