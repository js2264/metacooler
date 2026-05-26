"""Tests for zoomify operation."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.io.mmetacool import list_entries
from metacooler.ops.zoomify import zoomify
from metacooler.schema import is_mmetacooler


class TestZoomify:
    def test_zoomify_creates_mmetacool(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, levels=3)
        assert out.exists()
        assert is_mmetacooler(out)

    def test_zoomify_has_multiple_entries(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, levels=3)
        entries = list_entries(out)
        # Should have at least 2 entries (base + 1 coarser)
        assert len(entries) >= 2

    def test_zoomify_base_resolution_present(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        base_binsize = mc.binsize

        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, levels=3)

        entries = list_entries(out)
        # Base resolution should be present as an entry name
        assert str(base_binsize) in entries

    def test_zoomify_single_level(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, levels=1)
        entries = list_entries(out)
        # With levels=1, should have only the base resolution
        assert len(entries) == 1

    def test_zoomify_variable_bin_raises(
        self, tmp_path: Path, tiny_variable_bin_metacooler: Path
    ):
        """Variable-bin metacooler should raise ValueError."""
        out = tmp_path / "zoomed.mmetacool"
        with pytest.raises(ValueError, match="variable-bin"):
            zoomify(tiny_variable_bin_metacooler, out)
