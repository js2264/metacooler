"""The log of a command that writes a map or a file (CONVENTIONS.md § 11).

``<sample>.metacooler_<YYYYMMDDHHMMSS>.log`` beside what the command writes,
with hicstuff's header, a line per record from INFO up (from DEBUG with -v,
whatever -q says), and a last line saying what the command wrote and how long
it took; a command that only reads logs to stderr only.
"""

from __future__ import annotations

import logging
import re
import shutil
import time
from pathlib import Path

import pytest
from click.testing import CliRunner

from metacooler._version import __version__
from metacooler.cli import logging_utils
from metacooler.cli.logging_utils import (
    duration,
    log_file_name,
    run_log,
    setup_logging,
)
from metacooler.cli.main import cli

DATA = Path(__file__).parent / "data"
MAP = DATA / "assembly_1_1kb.metacool"
#: A line of the log after its header.
RECORD = re.compile(r"^\d{4}-\d{2}-\d{2},\d{2}:\d{2}:\d{2} :: [A-Z]+ :: ")


def logs(directory: Path) -> list[Path]:
    """The logs metacooler wrote in *directory*."""
    return sorted(directory.glob("*.metacooler_*.log"))


def one_log(directory: Path, sample: str) -> list[str]:
    """The lines of the one log of *sample* in *directory*."""
    found = logs(directory)
    assert len(found) == 1, found
    assert re.fullmatch(
        rf"{re.escape(sample)}\.metacooler_\d{{14}}\.log", found[0].name
    )
    return found[0].read_text().splitlines()


@pytest.fixture
def a_map(tmp_path: Path) -> Path:
    """The test data's map, copied where its command's log may go."""
    path = tmp_path / "in" / "gut T0.metacool"
    path.parent.mkdir()
    shutil.copy(MAP, path)
    return path


class TestTheLogOfACommand:
    def test_beside_the_output_with_hicstuffs_header(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "out dir" / "flat.cool"
        out.parent.mkdir()
        result = CliRunner().invoke(cli, ["flatten", str(a_map), "-o", str(out)])
        assert result.exit_code == 0, result.output
        lines = one_log(out.parent, "flat")
        assert lines[0] == f"## metacooler: v{__version__} log file"
        assert re.fullmatch(r"## date: \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", lines[1])
        # The command line as given, quoted where it needs to be.
        assert lines[2] == f"## cmd: metacooler flatten '{a_map}' -o '{out}'"
        header = lines[: lines.index("---")]
        assert f"## metacool_path: {a_map}" in header
        assert f"## output: {out}" in header
        assert "## force: False" in header
        assert "## quiet: False" in header
        records = lines[lines.index("---") + 1 :]
        assert records and all(RECORD.match(line) for line in records)
        assert re.search(
            rf":: INFO :: Cool file written after 0h 0m \d+s: {re.escape(str(out))}$",
            records[-1],
        )
        # Nothing of the log on stderr: it shows what it always showed.
        assert "## metacooler" not in result.stderr
        assert "Cool file written after" not in result.stderr

    def test_a_command_that_only_reads_writes_no_log(self, a_map: Path):
        for args in (["info", str(a_map)], ["dump", "-t", "genomes", str(a_map)]):
            result = CliRunner().invoke(cli, args)
            assert result.exit_code == 0, result.output
        assert logs(a_map.parent) == []

    def test_a_table_to_a_file_is_logged(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "genomes.tsv.gz"
        result = CliRunner().invoke(
            cli, ["dump", "-t", "genomes", str(a_map), "-o", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert "Table written after" in one_log(tmp_path, "genomes")[-1]

    def test_quiet_quietens_stderr_never_the_file(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "flat.cool"
        result = CliRunner().invoke(cli, ["flatten", str(a_map), "-o", str(out), "-q"])
        assert result.exit_code == 0, result.output
        assert "[INFO]" not in result.stderr
        lines = one_log(tmp_path, "flat")
        assert "## quiet: True" in lines
        assert any(" :: INFO :: Flattening " in line for line in lines)

    def test_an_edit_in_place_logs_beside_the_map(self, a_map: Path):
        result = CliRunner().invoke(
            cli, ["metadata", "add", str(a_map), "-k", "enzyme", "--value", "DpnII"]
        )
        assert result.exit_code == 0, result.output
        lines = one_log(a_map.parent, "gut T0")
        assert "## key: enzyme" in lines
        assert "Metadata written after" in lines[-1]

    def test_an_entry_of_a_collection_names_its_log(self, tiny_mmetacool: Path):
        entry = f"{tiny_mmetacool}::entries/entry1"
        result = CliRunner().invoke(
            cli, ["metadata", "add", entry, "-k", "enzyme", "--value", "DpnII"]
        )
        assert result.exit_code == 0, result.output
        assert "Metadata written after" in one_log(tiny_mmetacool.parent, "entry1")[-1]

    def test_several_values_are_joined(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "two.cool"
        result = CliRunner().invoke(
            cli,
            [
                "extract",
                str(a_map),
                "--genomes",
                "MAG_01",
                "--genomes",
                "MAG_02",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "## genomes: MAG_01,MAG_02" in one_log(tmp_path, "two")

    def test_a_refused_output_leaves_no_log(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "flat.cool"
        out.write_text("there already")
        result = CliRunner().invoke(cli, ["flatten", str(a_map), "-o", str(out)])
        assert result.exit_code != 0
        assert logs(tmp_path) == []

    def test_a_missing_directory_fails_as_it_did(self, a_map: Path, tmp_path: Path):
        out = tmp_path / "missing" / "flat.cool"
        result = CliRunner().invoke(cli, ["flatten", str(a_map), "-o", str(out)])
        assert result.exit_code != 0
        assert not out.parent.exists()


class TestRunLog:
    """:func:`run_log` itself, outside a command line."""

    def test_its_levels(self, tmp_path: Path):
        logger = logging.getLogger("metacooler")
        setup_logging(verbose=False)
        with run_log(tmp_path / "a.tsv", "Table written") as path:
            logger.debug("hidden")
            logger.info("shown")
        text = Path(str(path)).read_text()
        assert ":: INFO :: shown" in text and "hidden" not in text
        setup_logging(verbose=True)
        with run_log(tmp_path / "b.tsv", "Table written") as path:
            logger.debug("detail")
        assert ":: DEBUG :: detail" in Path(str(path)).read_text()
        setup_logging(quiet=True)
        assert logger.level == logging.INFO

    def test_a_failure_ends_the_log(self, tmp_path: Path):
        setup_logging()
        with (
            pytest.raises(ValueError, match="boom"),
            run_log(tmp_path / "a.cool", "Cool file written") as path,
        ):
            raise ValueError("boom")
        last = Path(str(path)).read_text().splitlines()[-1]
        assert re.search(r":: ERROR :: Failed after 0h 0m \d+s: boom$", last)
        # The handler is gone: the next record goes nowhere near that file.
        logging.getLogger("metacooler").info("a record once the run is over")
        assert "once the run is over" not in Path(str(path)).read_text()

    def test_an_exit_with_status_0_is_no_failure(self, tmp_path: Path):
        setup_logging()
        with (
            pytest.raises(SystemExit),
            run_log(tmp_path / "a.cool", "Cool file written") as path,
        ):
            raise SystemExit(0)
        assert "Cool file written after" in Path(str(path)).read_text()

    def test_nothing_for_stdout(self):
        with run_log(None, "Table written") as path:
            assert path is None

    def test_nothing_where_the_directory_is_missing(self, tmp_path: Path):
        with run_log(tmp_path / "missing" / "a.cool", "Cool file written") as path:
            assert path is None
        assert not (tmp_path / "missing").exists()
        with run_log(tmp_path / "made" / "a", "FASTA written", mkdir=True) as path:
            assert path is not None and path.parent == tmp_path / "made"

    def test_a_log_that_cannot_be_opened_is_said(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog
    ):
        def refuse(*args, **kwargs):
            raise PermissionError("no writing here")

        monkeypatch.setattr(logging_utils.logging, "FileHandler", refuse)
        with (
            caplog.at_level(logging.WARNING, logger="metacooler"),
            run_log(tmp_path / "a.cool", "Cool file written") as path,
        ):
            assert path is None
        assert "No log file for this run: no writing here" in caplog.text

    def test_its_name(self, tmp_path: Path):
        started = time.mktime((2026, 11, 26, 10, 15, 30, 0, 0, -1))
        assert log_file_name(tmp_path / "gut.metacool", started) == (
            tmp_path / "gut.metacooler_20261126101530.log"
        )
        assert log_file_name(f"{tmp_path}/c.mmetacool::entries/r1", started) == (
            tmp_path / "r1.metacooler_20261126101530.log"
        )
        assert log_file_name(
            tmp_path / "report", started, directory=tmp_path / "report", sample="gut"
        ) == (tmp_path / "report" / "gut.metacooler_20261126101530.log")
        # Without a sample's name, the tool's.
        assert log_file_name(tmp_path / "x", started, sample="").name == (
            "metacooler_20261126101530.log"
        )

    def test_durations(self):
        assert duration(0) == "0h 0m 0s"
        assert duration(1051.4) == "0h 17m 31s"
        assert duration(3600 * 2 + 61) == "2h 1m 1s"
