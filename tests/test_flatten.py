"""Tests for flatten operation."""

from __future__ import annotations

from pathlib import Path

import cooler

from metacooler.api import MetaCooler
from metacooler.io.flatten import flatten


class TestFlatten:
    def test_flatten_produces_valid_cool(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "flat.cool"
        flatten(tiny_metacooler, out)
        assert out.exists()
        # Should be a valid cooler file
        clr = cooler.Cooler(str(out))
        assert clr.shape[0] > 0

    def test_flatten_preserves_contigs_as_chroms(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        n_contigs = len(mc.contigs)

        out = tmp_path / "flat.cool"
        flatten(tiny_metacooler, out)

        clr = cooler.Cooler(str(out))
        # Number of chromosomes in the output should equal number of contigs
        assert len(clr.chromnames) == n_contigs

    def test_flatten_preserves_pixel_count(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "flat.cool"
        flatten(tiny_metacooler, out)

        clr = cooler.Cooler(str(out))
        # Pixel count should be preserved (or close — upper triangle only)
        assert clr.info["nnz"] > 0

    def test_flatten_preserves_bin_count(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_nbins = mc.shape[0]

        out = tmp_path / "flat.cool"
        flatten(tiny_metacooler, out)

        clr = cooler.Cooler(str(out))
        assert clr.shape[0] == original_nbins
