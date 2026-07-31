"""Tests for the integrity report (``metacooler.schema.check`` / ``metacooler check``)."""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.schema import check

DATA_DIR = Path(__file__).parent / "data"
SOURCE = DATA_DIR / "assembly_1_1kb.metacool"
COLLECTION = DATA_DIR / "assembly_1_1kb.mmetacool"
PLAIN_COOL = DATA_DIR / "assembly_1_rep1_1kb.cool"

_runner = CliRunner()


def _by_name(results) -> dict[str, bool]:
    """Map each check name to whether it passed."""
    return {r.name: r.ok for r in results}


def _failed(results) -> list[str]:
    return [r.name for r in results if not r.ok]


def _detail(results, name: str) -> str:
    return next(r.detail for r in results if r.name == name)


@pytest.fixture
def scratch_copy(tmp_path) -> Path:
    """A writable copy of the shipped metacool, for corruption tests."""
    dest = tmp_path / "copy.metacool"
    shutil.copy(SOURCE, dest)
    return dest


class TestCheckSound:
    def test_a_good_file_passes_every_check(self):
        results = check(SOURCE)
        assert _failed(results) == []
        assert set(_by_name(results)) == {
            "readable",
            "format",
            "structure",
            "dtypes",
            "cooler-compat",
        }

    def test_deep_adds_a_full_pixel_scan(self):
        results = check(SOURCE, deep=True)
        assert _failed(results) == []
        assert "pixels" in _by_name(results)

    def test_a_collection_is_checked_entry_by_entry(self):
        results = check(COLLECTION)
        assert _failed(results) == []
        names = _by_name(results)
        assert "HDF5::MMetaCooler" in next(
            r.detail for r in results if r.name == "format"
        )
        # One dtypes/cooler-compat pair per entry, each tagged with its name.
        per_entry = [n for n in names if n.startswith("dtypes [")]
        assert len(per_entry) == 2

    def test_an_entry_uri_checks_only_that_entry(self):
        entry = f"{COLLECTION}::entries/assembly_1_rep1_1kb"
        results = check(entry)
        assert _failed(results) == []
        assert _by_name(results)["format"] is True
        assert not any(r.name.startswith("dtypes [") for r in results)


class TestCheckReportsProblems:
    def test_missing_file(self, tmp_path):
        results = check(tmp_path / "ghost.metacool")
        assert _failed(results) == ["readable"]
        assert "not found" in results[0].detail.lower()

    def test_not_hdf5(self, tmp_path):
        path = tmp_path / "text.metacool"
        path.write_text("not hdf5\n")
        results = check(path)
        assert _failed(results) == ["readable"]

    def test_missing_group_in_a_uri(self):
        results = check(f"{COLLECTION}::entries/nope")
        assert _failed(results) == ["readable"]
        assert "not found" in results[0].detail

    def test_a_plain_cool_is_not_a_metacooler(self):
        results = check(PLAIN_COOL)
        assert _failed(results) == ["format"]
        # The report says what it *is*, so the fix is obvious.
        assert "plain cooler" in results[-1].detail
        # Nothing downstream is attempted once the format is wrong.
        assert len(results) == 2

    def test_a_bad_dtype_is_caught_where_validate_does_not_look(self, scratch_copy):
        """``validate`` casts on read, so a float64 count column slips past it."""
        with h5py.File(scratch_copy, "r+") as f:
            values = f["genomes/n_bins"][:]
            del f["genomes/n_bins"]
            f["genomes"].create_dataset("n_bins", data=values.astype(np.float64))
        results = check(scratch_copy)
        assert _failed(results) == ["dtypes"]
        assert "/genomes/n_bins is float64" in _detail(results, "dtypes")

    def test_a_wrong_declared_sum_is_caught_only_by_deep(self, scratch_copy):
        with h5py.File(scratch_copy, "r+") as f:
            f.attrs["sum"] = int(f.attrs["sum"]) + 7
        assert _failed(check(scratch_copy)) == []
        assert _failed(check(scratch_copy, deep=True)) == ["pixels"]
        assert "attrs['sum']" in _detail(check(scratch_copy, deep=True), "pixels")

    def test_scrambled_pixels_are_caught_by_deep(self, scratch_copy):
        """The sampled windows of ``validate`` can miss a local reordering."""
        with h5py.File(scratch_copy, "r+") as f:
            b2 = f["pixels/bin2_id"]
            # Swap two neighbours that share a bin1_id, breaking the (bin1, bin2)
            # order without moving any row out of its bin1_offset slot.
            b1 = f["pixels/bin1_id"][:]
            same = np.flatnonzero(b1[:-1] == b1[1:])
            i = int(same[len(same) // 2])
            b2[i], b2[i + 1] = b2[i + 1], b2[i]
        results = check(scratch_copy, deep=True)
        assert "pixels" in _failed(results)
        assert "not sorted" in _detail(results, "pixels")

    def test_a_broken_structure_is_reported_not_raised(self, scratch_copy):
        with h5py.File(scratch_copy, "r+") as f:
            n_contigs = f["contigs/name"].shape[0]
            f.attrs["nchroms"] = n_contigs + 1
        results = check(scratch_copy)
        assert "structure" in _failed(results)
        assert "nchroms" in _detail(results, "structure")
        # The later checks still ran, so one problem does not hide the others.
        assert "cooler-compat" in _by_name(results)


class TestCheckCLI:
    def test_help(self):
        r = _runner.invoke(cli, ["check", "--help"])
        assert r.exit_code == 0
        assert "--deep" in r.output

    def test_sound_file_exits_zero(self):
        r = _runner.invoke(cli, ["check", str(SOURCE)])
        assert r.exit_code == 0, r.output
        assert "OK (5 checks passed)" in r.output

    def test_a_collection_root_is_accepted(self):
        r = _runner.invoke(cli, ["check", str(COLLECTION)])
        assert r.exit_code == 0, r.output

    def test_failure_exits_non_zero_and_names_the_checks(self, scratch_copy):
        with h5py.File(scratch_copy, "r+") as f:
            del f["indexes/genome_offset"]
        r = _runner.invoke(cli, ["check", str(scratch_copy)])
        assert r.exit_code == 1, r.output
        assert "structure" in r.output

    def test_plain_cool_exits_non_zero(self):
        r = _runner.invoke(cli, ["check", str(PLAIN_COOL)])
        assert r.exit_code == 1, r.output

    def test_top_level_export(self):
        from metacooler import CheckResult as _CheckResult
        from metacooler import check as _check
        from metacooler.schema import CheckResult

        assert _check is check
        assert _CheckResult is CheckResult
