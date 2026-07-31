"""Tests for mmetacool management CLI commands (ls, cp, mv, rm)."""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.io.mmetacool import create_mmetacool, list_entries
from metacooler.schema import COOLER_FORMAT


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def two_entry_mmetacool(tmp_path: Path, tiny_metacooler: Path) -> Path:
    out = tmp_path / "test.mmetacool"
    create_mmetacool(out, {"entry1": tiny_metacooler, "entry2": tiny_metacooler})
    return out


@pytest.fixture
def tiny_mcool(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """A tiny .mcool file with one resolution."""
    out = tmp_path / "tiny.mcool"
    with h5py.File(out, "w") as dst:
        dst.attrs["format"] = COOLER_FORMAT
        res = dst.create_group("resolutions")
        with h5py.File(tiny_metacooler, "r") as src:
            src.copy(src["/"], res, name="1000")
    return out


class TestLsCommand:
    def test_ls_mmetacool(self, runner, two_entry_mmetacool: Path):
        result = runner.invoke(cli, ["ls", str(two_entry_mmetacool)])
        assert result.exit_code == 0
        assert "entry1" in result.output
        assert "entry2" in result.output

    def test_ls_metacool(self, runner, tiny_metacooler: Path):
        result = runner.invoke(cli, ["ls", str(tiny_metacooler)])
        assert result.exit_code == 0
        assert "::" in result.output

    def test_ls_mcool(self, runner, tiny_mcool: Path):
        result = runner.invoke(cli, ["ls", str(tiny_mcool)])
        assert result.exit_code == 0
        assert "resolutions/1000" in result.output


class TestCpCommand:
    def test_cp_within_file(self, runner, two_entry_mmetacool: Path):
        src = f"{two_entry_mmetacool}::entries/entry1"
        dst = f"{two_entry_mmetacool}::entries/entry3"
        result = runner.invoke(cli, ["cp", src, dst])
        assert result.exit_code == 0
        entries = list_entries(two_entry_mmetacool)
        assert "entry3" in entries

    def test_cp_cross_file(self, runner, tmp_path: Path, tiny_metacooler: Path):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"original": tiny_metacooler})
        create_mmetacool(dst_file, {"existing": tiny_metacooler})

        src = f"{src_file}::entries/original"
        dst = f"{dst_file}::entries/copied"
        result = runner.invoke(cli, ["cp", src, dst])
        assert result.exit_code == 0
        entries = list_entries(dst_file)
        assert "copied" in entries


class TestMvCommand:
    def test_mv_within_file(self, runner, two_entry_mmetacool: Path):
        src = f"{two_entry_mmetacool}::entries/entry1"
        dst = f"{two_entry_mmetacool}::entries/renamed"
        result = runner.invoke(cli, ["mv", src, dst])
        assert result.exit_code == 0
        entries = list_entries(two_entry_mmetacool)
        assert "renamed" in entries
        assert "entry1" not in entries

    def test_mv_cross_file_fails(self, runner, tmp_path: Path, tiny_metacooler: Path):
        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        create_mmetacool(src_file, {"entry": tiny_metacooler})
        create_mmetacool(dst_file, {})

        src = f"{src_file}::entries/entry"
        dst = f"{dst_file}::entries/entry"
        result = runner.invoke(cli, ["mv", src, dst])
        assert result.exit_code != 0
        assert (
            "same file" in result.output.lower()
            or "same file" in str(result.exception).lower()
        )


class TestRmCommand:
    def test_rm_with_force(self, runner, two_entry_mmetacool: Path):
        target = f"{two_entry_mmetacool}::entries/entry1"
        result = runner.invoke(cli, ["rm", "-f", target])
        assert result.exit_code == 0
        entries = list_entries(two_entry_mmetacool)
        assert "entry1" not in entries
        assert "entry2" in entries

    def test_rm_nonexistent_entry(self, runner, two_entry_mmetacool: Path):
        """A missing group is rejected by MetaCoolURI at parse time (exit 2).

        This is a Click UsageError raised while converting TARGET, not an
        exit-1 error from the command body.
        """
        target = f"{two_entry_mmetacool}::entries/nonexistent"
        result = runner.invoke(cli, ["rm", "-f", target])
        assert result.exit_code == 2, result.output
        flat = _flat(result.output)
        assert "Invalid value for 'TARGET'" in flat
        assert "Group '/entries/nonexistent' not found" in flat
        # nothing was removed
        assert list_entries(two_entry_mmetacool) == ["entry1", "entry2"]

    def test_rm_with_confirmation(self, runner, two_entry_mmetacool: Path):
        target = f"{two_entry_mmetacool}::entries/entry1"
        result = runner.invoke(cli, ["rm", target], input="y\n")
        assert result.exit_code == 0
        entries = list_entries(two_entry_mmetacool)
        assert "entry1" not in entries

    def test_rm_declined_exits_nonzero(self, runner, two_entry_mmetacool: Path):
        """Declining the confirmation must fail, not silently report success."""
        target = f"{two_entry_mmetacool}::entries/entry1"
        result = runner.invoke(cli, ["rm", target], input="n\n")
        assert result.exit_code != 0
        assert "Aborted" in result.output
        assert "entry1" in list_entries(two_entry_mmetacool)

    def test_rm_bare_mmetacool_rejected(self, runner, two_entry_mmetacool: Path):
        """A collection root is not a valid rm target (UsageError, exit 2)."""
        result = runner.invoke(cli, ["rm", "-f", str(two_entry_mmetacool)])
        assert result.exit_code == 2, result.output
        assert "is an mmetacool collection" in _flat(result.output)
        assert list_entries(two_entry_mmetacool) == ["entry1", "entry2"]
