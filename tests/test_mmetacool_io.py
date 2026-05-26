"""Tests for mmetacool I/O operations."""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest

from metacooler.io.mmetacool import (
    add_entry,
    copy_entry,
    copy_entry_cross_file,
    create_mmetacool,
    list_entries,
    move_entry,
    remove_entry,
)
from metacooler.schema import MMETACOOLER_FORMAT, MMETACOOLER_FORMAT_VERSION


class TestCreateMmetacool:
    def test_create_from_single_entry(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "single.mmetacool"
        create_mmetacool(out, {"sample1": tiny_metacooler})
        with h5py.File(out, "r") as f:
            assert f.attrs["format"] == MMETACOOLER_FORMAT
            assert f.attrs["format-version"] == MMETACOOLER_FORMAT_VERSION
            assert "entries" in f
            assert "sample1" in f["entries"]
            assert "bins" in f["entries/sample1"]

    def test_create_from_multiple_entries(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "multi.mmetacool"
        create_mmetacool(out, {"entry1": tiny_metacooler, "entry2": tiny_metacooler})
        entries = list_entries(out)
        assert sorted(entries) == ["entry1", "entry2"]

    def test_create_missing_source(self, tmp_path: Path):
        out = tmp_path / "fail.mmetacool"
        with pytest.raises(FileNotFoundError):
            create_mmetacool(out, {"bad": Path("/nonexistent.metacool")})


class TestAddEntry:
    def test_add_new_entry(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "add.mmetacool"
        create_mmetacool(out, {"first": tiny_metacooler})
        add_entry(out, "second", tiny_metacooler)
        entries = list_entries(out)
        assert sorted(entries) == ["first", "second"]

    def test_add_duplicate_raises(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "dup.mmetacool"
        create_mmetacool(out, {"first": tiny_metacooler})
        with pytest.raises(KeyError, match="already exists"):
            add_entry(out, "first", tiny_metacooler)

    def test_add_overwrite(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "overwrite.mmetacool"
        create_mmetacool(out, {"first": tiny_metacooler})
        add_entry(out, "first", tiny_metacooler, overwrite=True)
        entries = list_entries(out)
        assert entries == ["first"]

    def test_add_missing_source(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "addfail.mmetacool"
        create_mmetacool(out, {"first": tiny_metacooler})
        with pytest.raises(FileNotFoundError):
            add_entry(out, "bad", Path("/nonexistent.metacool"))


class TestRemoveEntry:
    def test_remove_existing(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "rm.mmetacool"
        create_mmetacool(out, {"entry1": tiny_metacooler, "entry2": tiny_metacooler})
        remove_entry(out, "entry1")
        entries = list_entries(out)
        assert entries == ["entry2"]

    def test_remove_nonexistent(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "rmbad.mmetacool"
        create_mmetacool(out, {"entry1": tiny_metacooler})
        with pytest.raises(KeyError, match="not found"):
            remove_entry(out, "nonexistent")


class TestCopyEntry:
    def test_copy_within_file(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "cp.mmetacool"
        create_mmetacool(out, {"original": tiny_metacooler})
        copy_entry(out, "original", "copy")
        entries = list_entries(out)
        assert sorted(entries) == ["copy", "original"]

    def test_copy_same_name_raises(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "cpsame.mmetacool"
        create_mmetacool(out, {"entry": tiny_metacooler})
        with pytest.raises(ValueError, match="same"):
            copy_entry(out, "entry", "entry")

    def test_copy_overwrite(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "cpover.mmetacool"
        create_mmetacool(out, {"src": tiny_metacooler, "dst": tiny_metacooler})
        copy_entry(out, "src", "dst", overwrite=True)
        entries = list_entries(out)
        assert sorted(entries) == ["dst", "src"]

    def test_copy_dst_exists_no_overwrite(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "cpfail.mmetacool"
        create_mmetacool(out, {"src": tiny_metacooler, "dst": tiny_metacooler})
        with pytest.raises(KeyError, match="already exists"):
            copy_entry(out, "src", "dst")

    def test_copy_nonexistent_src(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "cpbad.mmetacool"
        create_mmetacool(out, {"entry": tiny_metacooler})
        with pytest.raises(KeyError, match="not found"):
            copy_entry(out, "nonexistent", "new")


class TestMoveEntry:
    def test_move_within_file(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "mv.mmetacool"
        create_mmetacool(out, {"original": tiny_metacooler})
        move_entry(out, "original", "renamed")
        entries = list_entries(out)
        assert entries == ["renamed"]

    def test_move_same_name_raises(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "mvsame.mmetacool"
        create_mmetacool(out, {"entry": tiny_metacooler})
        with pytest.raises(ValueError, match="same"):
            move_entry(out, "entry", "entry")


class TestListEntries:
    def test_list_multiple(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "list.mmetacool"
        create_mmetacool(out, {"a": tiny_metacooler, "b": tiny_metacooler})
        entries = list_entries(out)
        assert sorted(entries) == ["a", "b"]

    def test_list_empty(self, tmp_path: Path):
        out = tmp_path / "empty.mmetacool"
        with h5py.File(out, "w") as f:
            f.attrs["format"] = MMETACOOLER_FORMAT
            f.attrs["format-version"] = MMETACOOLER_FORMAT_VERSION
            f.create_group("entries")
        assert list_entries(out) == []


class TestCrossFileCopy:
    def test_cross_file_copy(self, tmp_path: Path, tiny_metacooler: Path):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"entry1": tiny_metacooler})
        create_mmetacool(dst_file, {"existing": tiny_metacooler})
        copy_entry_cross_file(src_file, "entry1", dst_file, "copied")
        entries = list_entries(dst_file)
        assert sorted(entries) == ["copied", "existing"]

    def test_cross_file_copy_nonexistent_src(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"entry1": tiny_metacooler})
        create_mmetacool(dst_file, {})
        with pytest.raises(KeyError, match="not found"):
            copy_entry_cross_file(src_file, "nonexistent", dst_file, "new")

    def test_cross_file_copy_with_overwrite(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"entry1": tiny_metacooler})
        create_mmetacool(dst_file, {"existing": tiny_metacooler})
        # Copy once
        copy_entry_cross_file(src_file, "entry1", dst_file, "existing", overwrite=True)
        entries = list_entries(dst_file)
        assert "existing" in entries

    def test_cross_file_copy_dst_exists_no_overwrite_raises(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"entry1": tiny_metacooler})
        create_mmetacool(dst_file, {"existing": tiny_metacooler})
        with pytest.raises((KeyError, ValueError)):
            copy_entry_cross_file(
                src_file, "entry1", dst_file, "existing", overwrite=False
            )


class TestMoveEntry_2:
    def test_move_entry_basic(self, tmp_path: Path, tiny_metacooler: Path):
        mmetacool = tmp_path / "move.mmetacool"
        create_mmetacool(mmetacool, {"old_name": tiny_metacooler})
        move_entry(mmetacool, "old_name", "new_name")
        entries = list_entries(mmetacool)
        assert "new_name" in entries
        assert "old_name" not in entries

    def test_move_nonexistent_raises(self, tmp_path: Path, tiny_metacooler: Path):
        mmetacool = tmp_path / "move.mmetacool"
        create_mmetacool(mmetacool, {"entry1": tiny_metacooler})
        with pytest.raises((KeyError, ValueError)):
            move_entry(mmetacool, "ghost", "new_name")

    def test_move_to_existing_with_overwrite(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mmetacool = tmp_path / "move.mmetacool"
        create_mmetacool(
            mmetacool, {"entry1": tiny_metacooler, "entry2": tiny_metacooler}
        )
        move_entry(mmetacool, "entry1", "entry2", overwrite=True)
        entries = list_entries(mmetacool)
        assert "entry2" in entries
        assert "entry1" not in entries

    def test_move_preserves_data_integrity(self, tmp_path: Path, tiny_metacooler: Path):
        mmetacool = tmp_path / "integrity.mmetacool"
        create_mmetacool(mmetacool, {"orig": tiny_metacooler})
        move_entry(mmetacool, "orig", "renamed")
        with h5py.File(mmetacool, "r") as f:
            assert "bins" in f["entries/renamed"]


class TestListEntries_2:
    def test_list_empty_mmetacool(self, tmp_path: Path):
        mmetacool = tmp_path / "empty.mmetacool"
        create_mmetacool(mmetacool, {})
        assert list_entries(mmetacool) == []

    def test_list_entries_sorted(self, tmp_path: Path, tiny_metacooler: Path):
        mmetacool = tmp_path / "sorted.mmetacool"
        create_mmetacool(
            mmetacool,
            {"zz": tiny_metacooler, "aa": tiny_metacooler, "mm": tiny_metacooler},
        )
        entries = list_entries(mmetacool)
        assert entries == sorted(entries)

    def test_list_entries_uses_groups_key(self, tmp_path: Path, tiny_metacooler: Path):
        """Verify list_entries reads from entries/groups subgroup (not flat)."""
        mmetacool = tmp_path / "grp.mmetacool"
        create_mmetacool(mmetacool, {"s1": tiny_metacooler, "s2": tiny_metacooler})
        assert set(list_entries(mmetacool)) == {"s1", "s2"}
