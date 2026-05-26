"""Tests for stack CLI command."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.io.mmetacool import list_entries


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def two_metacools(tmp_path: Path, tiny_metacooler: Path) -> tuple[Path, Path]:
    """Two copies of tiny_metacooler for stacking tests."""
    mc1 = tmp_path / "sample1.metacool"
    mc2 = tmp_path / "sample2.metacool"
    shutil.copy(tiny_metacooler, mc1)
    shutil.copy(tiny_metacooler, mc2)
    return mc1, mc2


@pytest.fixture
def three_metacools(tmp_path: Path, tiny_metacooler: Path) -> tuple[Path, Path, Path]:
    """Three copies of tiny_metacooler for multi-file stacking tests."""
    mc1 = tmp_path / "sample1.metacool"
    mc2 = tmp_path / "sample2.metacool"
    mc3 = tmp_path / "sample3.metacool"
    shutil.copy(tiny_metacooler, mc1)
    shutil.copy(tiny_metacooler, mc2)
    shutil.copy(tiny_metacooler, mc3)
    return mc1, mc2, mc3


class TestStackCommand:
    def test_stack_auto_named(self, runner, two_metacools, tmp_path):
        """Stack with auto-generated names from filenames."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1), str(mc2)])
        assert result.exit_code == 0
        assert "Created" in result.output
        entries = list_entries(output)
        assert "sample1" in entries
        assert "sample2" in entries

    def test_stack_custom_names(self, runner, two_metacools, tmp_path):
        """Stack with custom entry names."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli,
            [
                "stack",
                "-o",
                str(output),
                "--names",
                "MAG_A,MAG_B",
                str(mc1),
                str(mc2),
            ],
        )
        assert result.exit_code == 0
        entries = list_entries(output)
        assert "MAG_A" in entries
        assert "MAG_B" in entries

    def test_stack_custom_names_with_spaces(self, runner, two_metacools, tmp_path):
        """Stack with custom entry names that have spaces (trimmed)."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli,
            [
                "stack",
                "-o",
                str(output),
                "--names",
                "  genome_1  ,  genome_2  ",
                str(mc1),
                str(mc2),
            ],
        )
        assert result.exit_code == 0
        entries = list_entries(output)
        assert "genome_1" in entries
        assert "genome_2" in entries

    def test_stack_names_count_mismatch(self, runner, two_metacools, tmp_path):
        """Stack with mismatched --names count."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli,
            [
                "stack",
                "-o",
                str(output),
                "--names",
                "MAG_A",  # Only 1 name for 2 files
                str(mc1),
                str(mc2),
            ],
        )
        assert result.exit_code != 0
        assert "must match" in result.output.lower()

    def test_stack_requires_two_inputs(self, runner, two_metacools, tmp_path):
        """Stack requires at least 2 inputs."""
        mc1, _ = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1)])
        assert result.exit_code != 0
        assert "least 2" in result.output.lower()

    def test_stack_no_overwrite_by_default(self, runner, two_metacools, tmp_path):
        """Stack fails if output exists without --force."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        # Create output file first
        output.touch()
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1), str(mc2)])
        assert result.exit_code != 0
        assert "already exists" in result.output.lower()

    def test_stack_force_overwrite(self, runner, two_metacools, tmp_path):
        """Stack overwrites with --force."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        # Create output file first
        output.touch()
        result = runner.invoke(
            cli, ["stack", "-o", str(output), "-f", str(mc1), str(mc2)]
        )
        assert result.exit_code == 0
        assert "Created" in result.output

    def test_stack_name_collision(self, runner, tmp_path, tiny_metacooler):
        """Stack detects name collisions."""
        # Create two files with same stem name in different directories
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()

        mc1 = dir1 / "sample.metacool"
        mc2 = dir2 / "sample.metacool"
        shutil.copy(tiny_metacooler, mc1)
        shutil.copy(tiny_metacooler, mc2)

        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1), str(mc2)])
        assert result.exit_code != 0
        assert "collision" in result.output.lower()

    def test_stack_verbose_logging(self, runner, two_metacools, tmp_path):
        """Stack with verbose flag produces debug output."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli, ["stack", "-v", "-o", str(output), str(mc1), str(mc2)]
        )
        assert result.exit_code == 0
        # Verbose should produce some output
        assert "Created" in result.output

    def test_stack_nonexistent_input(self, runner, tmp_path):
        """Stack fails if input doesn't exist."""
        mc1 = tmp_path / "missing1.metacool"
        mc2 = tmp_path / "missing2.metacool"
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1), str(mc2)])
        assert result.exit_code != 0
        assert "does not exist" in result.output.lower()

    def test_stack_three_files(self, runner, three_metacools, tmp_path):
        """Stack multiple files (more than 2)."""
        mc1, mc2, mc3 = three_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli, ["stack", "-o", str(output), str(mc1), str(mc2), str(mc3)]
        )
        assert result.exit_code == 0
        entries = list_entries(output)
        assert len(entries) == 3
        assert "sample1" in entries
        assert "sample2" in entries
        assert "sample3" in entries

    def test_stack_entries_are_readable(self, runner, two_metacools, tmp_path):
        """Verify stacked entries can be read back via ls command."""
        mc1, mc2 = two_metacools
        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(cli, ["stack", "-o", str(output), str(mc1), str(mc2)])
        assert result.exit_code == 0

        # Now try to list entries via ls command
        result_ls = runner.invoke(cli, ["ls", str(output)])
        assert result_ls.exit_code == 0
        assert "sample1" in result_ls.output
        assert "sample2" in result_ls.output

    def test_stack_output_path_required(self, runner, two_metacools):
        """Stack requires --output argument."""
        mc1, mc2 = two_metacools
        result = runner.invoke(cli, ["stack", str(mc1), str(mc2)])
        assert result.exit_code != 0
        assert "output" in result.output.lower() or "required" in result.output.lower()

    def test_stack_many_files(self, runner, tmp_path, tiny_metacooler):
        """Stack with many files (5+)."""
        files = []
        for i in range(5):
            mc = tmp_path / f"sample{i}.metacool"
            shutil.copy(tiny_metacooler, mc)
            files.append(mc)

        output = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli, ["stack", "-o", str(output)] + [str(f) for f in files]
        )
        assert result.exit_code == 0
        entries = list_entries(output)
        assert len(entries) == 5
        for i in range(5):
            assert f"sample{i}" in entries
