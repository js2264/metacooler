"""The metadata keys metacooler writes itself cannot be set or removed by hand."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import pytest
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli._provenance import RESERVED_METADATA_KEYS
from metacooler.cli.main import cli


def flat(output: str) -> str:
    """A CLI error message with rich's box drawing and line wrapping undone."""
    return " ".join(output.replace("│", " ").split())


@pytest.fixture
def mc(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """A copy of the tiny map with a provenance record, as the CLI writes one."""
    path = tmp_path / "edit.metacool"
    path.write_bytes(tiny_metacooler.read_bytes())
    with h5py.File(path, "a") as f:
        meta = json.loads(f.attrs.get("metadata", "{}"))
        meta["commands"] = ["metacooler create pairs --pairs p.pairs -o edit.metacool"]
        f.attrs["metadata"] = json.dumps(meta)
    return path


@pytest.mark.parametrize("key", sorted(RESERVED_METADATA_KEYS))
def test_add_refuses_a_reserved_key(mc, key):
    before = MetaCooler(mc).metadata
    result = CliRunner().invoke(
        cli, ["metadata", "add", str(mc), "-k", key, "--value", "[]"]
    )
    assert result.exit_code == 2
    assert f"Cannot set {key!r}" in flat(result.output)
    assert MetaCooler(mc).metadata == before


def test_remove_refuses_a_reserved_key(mc):
    before = MetaCooler(mc).metadata
    assert before["commands"]
    result = CliRunner().invoke(cli, ["metadata", "remove", str(mc), "-k", "commands"])
    assert result.exit_code == 2
    assert "Cannot remove 'commands'" in flat(result.output)
    assert MetaCooler(mc).metadata == before


def test_reserved_keys_stay_readable(mc):
    result = CliRunner().invoke(cli, ["metadata", "get", str(mc), "-k", "commands"])
    assert result.exit_code == 0, result.output
    assert "create pairs" in result.output


def test_other_keys_are_still_editable(mc):
    runner = CliRunner()
    add = ["metadata", "add", str(mc), "-k", "sample", "--value", "gut_T0"]
    assert runner.invoke(cli, add).exit_code == 0
    assert MetaCooler(mc).metadata["sample"] == "gut_T0"
    remove = ["metadata", "remove", str(mc), "-k", "sample"]
    assert runner.invoke(cli, remove).exit_code == 0
    assert "sample" not in MetaCooler(mc).metadata


@pytest.mark.parametrize("key", ["reassign", "annotate"])
def test_load_refuses_the_history_keys(mc, tmp_path, key):
    meta = tmp_path / "meta.yaml"
    meta.write_text(f"{key}: overwritten\n")
    result = CliRunner().invoke(cli, ["metadata", "load", str(meta), str(mc)])
    assert result.exit_code != 0
    assert f"must not contain a {key!r} key" in flat(result.output)
