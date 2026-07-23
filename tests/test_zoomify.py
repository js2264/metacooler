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

    def test_zoomify_explicit_resolutions(self, tmp_path: Path, tiny_metacooler: Path):
        # base binsize is 100bp; explicit cooler-style list.
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, resolutions="100,200,400")
        entries = set(list_entries(out))
        assert {"100", "200", "400"} <= entries

    def test_zoomify_nice_progression(self, tmp_path: Path, tiny_metacooler: Path):
        # '100N' is the 1-2-5 nice progression: 100,200,500,1000.
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, resolutions="100N")
        entries = set(list_entries(out))
        assert {"100", "200", "500"} <= entries

    def test_zoomify_binary_progression(self, tmp_path: Path, tiny_metacooler: Path):
        # '100B' is the binary x2 progression: 100,200,400,800.
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, resolutions="100B")
        entries = set(list_entries(out))
        assert {"100", "200", "400"} <= entries
        assert "500" not in entries  # 500 is a 'nice' step, not binary

    def test_zoomify_ignores_finer_than_base(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        # Resolutions below the base (100) are dropped; base is always kept.
        out = tmp_path / "zoomed.mmetacool"
        zoomify(tiny_metacooler, out, resolutions="50,200")
        entries = set(list_entries(out))
        assert "50" not in entries
        assert {"100", "200"} <= entries

    def test_zoomify_invalid_resolution_raises(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        out = tmp_path / "zoomed.mmetacool"
        with pytest.raises(ValueError, match="Invalid resolution token"):
            zoomify(tiny_metacooler, out, resolutions="oops")

    def test_zoomify_entry_uri_rejected(self, tmp_path: Path, tiny_mmetacool: Path):
        out = tmp_path / "zoomed.mmetacool"
        with pytest.raises(ValueError, match="entry URIs"):
            zoomify(f"{tiny_mmetacool}::entries/entry1", out)
