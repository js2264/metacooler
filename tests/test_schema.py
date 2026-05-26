"""Tests for metacooler.schema."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from metacooler.schema import (
    MMETACOOLER_FORMAT,
    MMETACOOLER_FORMAT_VERSION,
    MetaCoolerFormatError,
    is_metacooler,
    is_mmetacooler,
    list_metacoolers,
    read_table,
    validate,
)


class TestIsMetacooler:
    def test_false_for_missing_file(self):
        assert is_metacooler("/nonexistent/file.metacool") is False

    def test_true_for_metacooler(self, tiny_metacooler):
        assert is_metacooler(tiny_metacooler) is True

    def test_false_for_mmetacool(self, tmp_path, tiny_metacooler):
        out = tmp_path / "test.mmetacool"
        with h5py.File(out, "w") as dst:
            dst.attrs["format"] = MMETACOOLER_FORMAT
            dst.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
            dst.create_group("entries")
        assert is_metacooler(out) is False

    def test_false_for_plain_hdf5(self, tmp_path):
        p = tmp_path / "plain.h5"
        with h5py.File(p, "w") as f:
            f.create_dataset("data", data=np.array([1, 2, 3]))
        assert is_metacooler(p) is False


class TestIsMmetacooler:
    def test_false_for_metacooler(self, tiny_metacooler):
        assert is_mmetacooler(tiny_metacooler) is False

    def test_false_for_missing_file(self):
        assert is_mmetacooler("/nonexistent.mmetacool") is False


class TestValidate:
    def test_passes_for_valid(self, tiny_metacooler):
        validate(tiny_metacooler)

    def test_fails_for_missing_file(self):
        with pytest.raises(MetaCoolerFormatError, match="File not found"):
            validate("/nonexistent")


class TestReadTable:
    def test_read_genomes_table(self, tiny_metacooler):
        df = read_table(tiny_metacooler, "genomes")
        assert "name" in df.columns
        assert "length" in df.columns
        assert "n_contigs" in df.columns
        assert len(df) == 3  # MAG_1, MAG_2, __unbinned__

    def test_read_contigs_table(self, tiny_metacooler):
        df = read_table(tiny_metacooler, "contigs")
        assert "name" in df.columns
        assert "length" in df.columns
        assert "genome_id" in df.columns
        assert len(df) == 4

    def test_byte_strings_decoded(self, tiny_metacooler):
        df = read_table(tiny_metacooler, "genomes")
        # All name values must be str, not bytes
        for name in df["name"]:
            assert isinstance(name, str)

    def test_read_nonexistent_group_raises(self, tiny_metacooler):
        with pytest.raises(KeyError):
            read_table(tiny_metacooler, "nonexistent_group")


class TestListMetacoolers:
    def test_list_on_mmetacool(self, tmp_path, tiny_metacooler):
        out = tmp_path / "test.mmetacool"
        with h5py.File(out, "w") as dst:
            dst.attrs["format"] = MMETACOOLER_FORMAT
            dst.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
            grp = dst.create_group("entries")
            with h5py.File(tiny_metacooler, "r") as src:
                src.copy(src["/"], grp, name="entry1")
                src.copy(src["/"], grp, name="entry2")
        entries = list_metacoolers(out)
        assert "/entries/entry1" in entries
        assert "/entries/entry2" in entries
        assert len(entries) == 2

    def test_empty_mmetacool_returns_empty_list(self, tmp_path):
        out = tmp_path / "empty.mmetacool"
        with h5py.File(out, "w") as f:
            f.create_group("entries")
        entries = list_metacoolers(out)
        assert entries == []
