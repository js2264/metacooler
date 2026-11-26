"""metacooler report: the report of a map (CONVENTIONS.md § 12).

The test data's report.txt is tests/data/assembly_1_1kb.report.txt, and is
compared with what the command writes: a change to it is a change to the
report's format, which raises ops.report.FORMAT_VERSION when a key, a column,
a unit or a MultiQC column changes or goes. To write it again, after a
change made on purpose:

    uv run metacooler report tests/data/assembly_1_1kb.metacool -o tests/data -f
    rm tests/data/assembly_1_1kb.report.html tests/data/assembly_1_1kb.report_mqc.yaml
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

import metacooler
from metacooler.cli.main import cli
from metacooler.ops import report as ops_report
from metacooler.ops.report import (
    FORMAT_VERSION,
    GENOME_COLUMNS,
    MQC_COLUMNS,
    PS_COLUMNS,
    _length_steps,
    _number,
    _pct,
    collect,
    n50,
)

DATA = Path(__file__).parent / "data"
MAP = DATA / "assembly_1_1kb.metacool"
REFERENCE = DATA / "assembly_1_1kb.report.txt"
#: The keys of the section report whose values change from run to run.
_RUN_KEYS = re.compile(r"^(map|written|metacooler_version)( +).*$", re.MULTILINE)


def _stable(text: str) -> str:
    """report.txt without what changes from run to run."""
    return _RUN_KEYS.sub(r"\1\2-", text)


@pytest.fixture
def a_map(tmp_path: Path) -> Path:
    path = tmp_path / "maps" / "assembly_1_1kb.metacool"
    path.parent.mkdir()
    shutil.copy(MAP, path)
    return path


def run(*args: str):
    result = CliRunner().invoke(cli, ["report", *args])
    return result


class TestTheCommand:
    def test_three_files_and_nothing_else(self, a_map: Path, tmp_path: Path):
        outdir = tmp_path / "report"
        result = run(str(a_map), "-o", str(outdir))
        assert result.exit_code == 0, result.output
        # No log file: metacooler logs to the terminal only.
        assert sorted(path.name for path in outdir.iterdir()) == [
            "assembly_1_1kb.report.html",
            "assembly_1_1kb.report.txt",
            "assembly_1_1kb.report_mqc.yaml",
        ]
        # The map is read only.
        assert sorted(path.name for path in a_map.parent.iterdir()) == [a_map.name]

    def test_report_txt_is_the_reference(self, a_map: Path, tmp_path: Path):
        result = run(str(a_map), "-o", str(tmp_path))
        assert result.exit_code == 0, result.output
        written = (tmp_path / "assembly_1_1kb.report.txt").read_text()
        assert _stable(written) == _stable(REFERENCE.read_text())
        assert f"format_version      {FORMAT_VERSION}\n" in written

    def test_read_report_reads_back_what_was_written(self, a_map: Path, tmp_path: Path):
        txt, _, _ = metacooler.write_report(a_map, tmp_path)
        sections, tables = metacooler.read_report(txt)
        assert list(sections) == [
            "report",
            "run",
            "commands",
            "genomes",
            "contigs",
            "bins",
            "contacts",
            "qc",
        ]
        assert sections["report"]["sample"] == "assembly_1_1kb"
        assert sections["genomes"]["genomes"] == "10"
        assert sections["contacts"]["contacts"] == "100000"
        assert sections["qc"]["stored"] == "false"
        assert sections["qc"]["inter_genome_pct"] == "28.61"
        assert list(tables["genomes"].columns) == list(GENOME_COLUMNS)
        assert len(tables["genomes"]) == 11
        assert list(tables["ps"].columns) == list(PS_COLUMNS)

    def test_its_multiqc_row(self, a_map: Path, tmp_path: Path):
        _, _, mqc = metacooler.write_report(a_map, tmp_path)
        content = yaml.safe_load(mqc.read_text())
        assert content["id"] == "metacooler"
        assert content["plot_type"] == "generalstats"
        assert [next(iter(column)) for column in content["pconfig"]] == [
            key for key, _, _ in MQC_COLUMNS
        ]
        row = content["data"]["assembly_1_1kb"]
        assert row == {
            "genomes": 10,
            "contigs": 44,
            "contacts": 100000,
            "inter_genome_pct": pytest.approx(28.608),
        }

    def test_its_page(self, a_map: Path, tmp_path: Path):
        _, html, _ = metacooler.write_report(a_map, tmp_path)
        page = html.read_text()
        for section in (
            "overview",
            "genomes",
            "contigs",
            "contacts",
            "ps",
            "run",
            "files",
        ):
            assert f'id="{section}"' in page
        found = re.search(
            r'<script id="cx-report-figures" type="application/json">(.*?)</script>',
            page,
            re.S,
        )
        assert found is not None
        figures = json.loads(found.group(1))
        assert sorted(figures) == ["plot-contigs", "plot-genomes", "plot-ps", "plot-qc"]
        assert 'href="assembly_1_1kb.report.txt"' in page
        assert 'href="assembly_1_1kb.report_mqc.yaml"' in page
        assert "metacooler merge" in page  # the map's provenance, in Run

    def test_an_earlier_report_is_refused_without_force(
        self, a_map: Path, tmp_path: Path
    ):
        assert run(str(a_map), "-o", str(tmp_path)).exit_code == 0
        written = sorted(path.name for path in tmp_path.iterdir())
        result = run(str(a_map), "-o", str(tmp_path))
        assert result.exit_code != 0
        assert "exists" in result.output
        assert sorted(path.name for path in tmp_path.iterdir()) == written
        result = run(str(a_map), "-o", str(tmp_path), "-f")
        assert result.exit_code == 0, result.output

    def test_an_entry_is_reported_under_its_name(
        self, tiny_mmetacool: Path, tmp_path: Path
    ):
        outdir = tmp_path / "reports"
        result = run(f"{tiny_mmetacool}::entries/entry2", "-o", str(outdir))
        assert result.exit_code == 0, result.output
        assert (outdir / "entry2.report.txt").is_file()
        assert (outdir / "entry2.report_mqc.yaml").is_file()

    def test_a_collection_is_refused(self, tiny_mmetacool: Path, tmp_path: Path):
        result = run(str(tiny_mmetacool), "-o", str(tmp_path / "r"))
        assert result.exit_code == 2
        assert "entries/entry1" in result.output
        assert not (tmp_path / "r").exists()


class TestWhatItHolds:
    def test_a_stored_qc_block_is_used(self, tiny_metacooler: Path):
        before = collect(tiny_metacooler).report.sections["qc"]
        result = CliRunner().invoke(cli, ["qc", str(tiny_metacooler), "--store"])
        assert result.exit_code == 0, result.output
        after = collect(tiny_metacooler).report.sections["qc"]
        assert after["stored"] is True and after["source"] == "pixels"
        # Measured between bin midpoints either way: the same breakdown.
        for key in ("intra_contig_short", "inter_genome", "trans_contig_ratio"):
            assert after[key] == before[key]

    def test_the_contacts_without_a_sum_attribute(self, a_map: Path):
        expected = collect(a_map).report.sections["contacts"]["contacts"]
        with h5py.File(a_map, "r+") as f:
            del f.attrs["sum"]
        assert collect(a_map, chunksize=1000).report.sections["contacts"][
            "contacts"
        ] == (expected)

    def test_a_map_without_a_ps(
        self, a_map: Path, monkeypatch: pytest.MonkeyPatch, caplog
    ):
        import metacooler.ops.ps as ops_ps

        def empty(*args, **kwargs):
            raise ValueError("no contig selected")

        monkeypatch.setattr(ops_ps, "contact_probability", empty)
        collected = collect(a_map)
        assert collected.report.tables["ps"].empty
        assert "No P(s) for this map: no contig selected" in caplog.text
        assert collected.figures["plot-ps"]["data"][0]["x"] == []


class TestItsNumbers:
    def test_n50(self):
        assert n50(np.array([2, 3, 4, 5, 6])) == (5, 2)
        assert n50(np.array([10])) == (10, 1)
        assert n50(np.array([], dtype=np.int64)) == (0, 0)

    @pytest.mark.parametrize(
        ("value", "number"),
        [
            (None, None),
            (float("nan"), None),
            (np.int64(3), 3),
            (2.5, 2.5),
            ("1000", 1000),
            (b"12", 12),
            ("fixed", None),
            (object(), None),
        ],
    )
    def test_number(self, value, number):
        assert _number(value) == number

    def test_pct(self):
        assert _pct(1, 4) == 25
        assert _pct(None, 4) is None
        assert _pct(1, None) is None
        assert _pct(1, 0) is None

    def test_length_steps(self):
        assert _length_steps(np.array([], dtype=np.int64)) == ([], [])
        edges, steps = _length_steps(np.array([1000, 1000]))
        assert edges[0] < 1000 < edges[-1] and sum(steps[:-1]) == 2
        edges, steps = _length_steps(np.array([1000, 10_000, 100_000]))
        assert len(edges) == len(steps) == 25 and sum(steps[:-1]) == 3

    def test_the_tool(self):
        tool = ops_report.tool()
        assert tool.id == "metacooler" and tool.format_version == FORMAT_VERSION
        assert Path(str(tool.template)).is_file()

    def test_the_genomes_table_follows_the_map(self, tiny_metacooler: Path):
        table = collect(tiny_metacooler).report.tables["genomes"]
        genomes = metacooler.MetaCooler(tiny_metacooler).genomes
        assert table["name"].tolist() == genomes["name"].astype(str).tolist()
        pd.testing.assert_series_equal(
            table["length_bp"], genomes["length"], check_names=False, check_dtype=False
        )
