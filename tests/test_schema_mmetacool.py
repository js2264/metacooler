"""Tests for mmetacooler schema functions."""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest

from metacooler.schema import (
    MMETACOOLER_FORMAT,
    MMETACOOLER_FORMAT_VERSION,
    MetaCoolerFormatError,
    is_mmetacooler,
    list_metacoolers,
    validate_mmetacooler,
)


@pytest.fixture
def tiny_mmetacool(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """Create a 2-entry .mmetacool from the tiny_metacooler fixture."""
    out = tmp_path / "tiny.mmetacool"
    with h5py.File(out, "w") as dst:
        dst.attrs["format"] = MMETACOOLER_FORMAT
        dst.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
        grp = dst.create_group("entries")
        with h5py.File(tiny_metacooler, "r") as src:
            src.copy(src["/"], grp, name="entry1")
            src.copy(src["/"], grp, name="entry2")
    return out


class TestIsMmetacooler:
    def test_valid_mmetacool(self, tiny_mmetacool: Path):
        assert is_mmetacooler(tiny_mmetacool) is True

    def test_metacool_is_not_mmetacool(self, tiny_metacooler: Path):
        assert is_mmetacooler(tiny_metacooler) is False

    def test_nonexistent_file(self):
        assert is_mmetacooler("/nonexistent/file.mmetacool") is False

    def test_empty_hdf5(self, tmp_path: Path):
        p = tmp_path / "empty.h5"
        with h5py.File(p, "w"):
            pass
        assert is_mmetacooler(p) is False

    def test_wrong_format(self, tmp_path: Path):
        p = tmp_path / "wrong.h5"
        with h5py.File(p, "w") as f:
            f.attrs["format"] = "WrongFormat"
            f.create_group("entries")
        assert is_mmetacooler(p) is False


class TestValidateMmetacooler:
    def test_valid(self, tiny_mmetacool: Path):
        validate_mmetacooler(tiny_mmetacool)  # should not raise

    def test_missing_file(self):
        with pytest.raises(MetaCoolerFormatError, match="File not found"):
            validate_mmetacooler("/nonexistent.mmetacool")

    def test_wrong_format(self, tmp_path: Path):
        p = tmp_path / "wrong.h5"
        with h5py.File(p, "w") as f:
            f.attrs["format"] = "WrongFormat"
            f.create_group("entries")
        with pytest.raises(MetaCoolerFormatError, match="Expected format"):
            validate_mmetacooler(p)

    def test_missing_format_version(self, tmp_path: Path):
        p = tmp_path / "nover.h5"
        with h5py.File(p, "w") as f:
            f.attrs["format"] = MMETACOOLER_FORMAT
            f.create_group("entries")
        with pytest.raises(MetaCoolerFormatError, match="format-version"):
            validate_mmetacooler(p)

    def test_missing_entries(self, tmp_path: Path):
        p = tmp_path / "noentries.h5"
        with h5py.File(p, "w") as f:
            f.attrs["format"] = MMETACOOLER_FORMAT
            f.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
        with pytest.raises(MetaCoolerFormatError, match="/entries"):
            validate_mmetacooler(p)

    def test_entry_missing_table(self, tmp_path: Path, tiny_metacooler: Path):
        p = tmp_path / "bad_entry.mmetacool"
        with h5py.File(p, "w") as dst:
            dst.attrs["format"] = MMETACOOLER_FORMAT
            dst.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
            grp = dst.create_group("entries")
            # Create a partial entry (missing pixels)
            entry = grp.create_group("bad")
            entry.create_group("bins")
            entry.create_group("contigs")
            entry.create_group("genomes")
            entry.create_group("indexes")
            # missing pixels
        with pytest.raises(MetaCoolerFormatError, match="missing required group"):
            validate_mmetacooler(p)


class TestListMetacoolers:
    def test_list_entries(self, tiny_mmetacool: Path):
        entries = list_metacoolers(tiny_mmetacool)
        assert "/entries/entry1" in entries
        assert "/entries/entry2" in entries
        assert len(entries) == 2

    def test_empty_mmetacool(self, tmp_path: Path):
        p = tmp_path / "empty.mmetacool"
        with h5py.File(p, "w") as f:
            f.attrs["format"] = MMETACOOLER_FORMAT
            f.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
            f.create_group("entries")
        assert list_metacoolers(p) == []
