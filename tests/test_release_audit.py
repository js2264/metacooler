"""Regression tests from the pre-release audit.

Every test in this module pins down a defect found while auditing the package
for its first public release. As in ``test_regressions.py``, each docstring
states what went wrong before the fix, and the tests assert on observable
behaviour (file contents, exit codes, raised messages) rather than on how the
fix is implemented.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.io._utils import check_output_path
from metacooler.schema import validate


def _run(*args: str):
    return CliRunner().invoke(cli, [str(a) for a in args])


class TestOutputNeverDestroysInput:
    """``-o <input> -f`` deleted the input, then failed with "File not found".

    The ``--force`` guard unlinked an existing output before the command ran,
    and nothing compared the output against the inputs, so naming the input as
    the output destroyed it.
    """

    @pytest.mark.parametrize(
        "command",
        [
            ["binnify", "{src}", "--binsize", "200"],
            ["merge", "{src}", "{src}"],
            ["zoomify", "{src}", "--levels", "2"],
            ["aggregate", "{src}", "--level", "genome"],
            ["flatten", "{src}"],
            ["extract", "{src}", "--genomes", "MAG_1"],
            ["dump", "-t", "contigs", "{src}"],
            ["network", "{src}"],
            ["qc", "{src}", "--table"],
        ],
        ids=lambda c: c[0],
    )
    def test_input_survives_being_named_as_output(
        self, tmp_path: Path, tiny_metacooler: Path, command: list[str]
    ):
        src = tmp_path / "input.metacool"
        shutil.copy(tiny_metacooler, src)
        before = src.read_bytes()
        args = [a.format(src=src) for a in command] + ["-o", str(src), "-f"]
        result = _run(*args)
        assert result.exit_code != 0
        assert "same file as the input" in _flat(result.output)
        assert src.read_bytes() == before

    def test_entry_uri_input_is_compared_by_file(
        self, tmp_path: Path, tiny_mmetacool: Path
    ):
        """An ``::entries/x`` URI names the same file as its bare path."""
        before = tiny_mmetacool.read_bytes()
        result = _run(
            "binnify",
            f"{tiny_mmetacool}::entries/entry1",
            "--binsize",
            "200",
            "-o",
            tiny_mmetacool,
            "-f",
        )
        assert result.exit_code != 0
        assert "same file as the input" in _flat(result.output)
        assert tiny_mmetacool.read_bytes() == before

    def test_api_guard_compares_inputs(self, tmp_path: Path):
        f = tmp_path / "x.metacool"
        f.write_bytes(b"x")
        with pytest.raises(ValueError, match="same file as the input"):
            check_output_path(f, force=True, inputs=[str(f)])
        assert f.exists()

    def test_api_guard_ignores_stdin_and_missing_inputs(self, tmp_path: Path):
        f = tmp_path / "x.metacool"
        f.write_bytes(b"x")
        check_output_path(f, force=True, inputs=["-", None, tmp_path / "nope"])
        assert not f.exists()


class TestForceKeepsPreviousOutputUntilSuccess:
    """``create --force`` deleted the previous map before doing any work.

    ``create`` stages its output under ``<name>.partial`` and renames it into
    place only once validation passes, precisely so a failed run cannot destroy
    a finished file — but the CLI's ``--force`` guard unlinked the old output
    (and its ``.ps.tsv``) up front, so a run that then failed left nothing.
    """

    def test_failed_recreate_keeps_old_map_and_ps(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ):
        out = tmp_path / "out.metacool"
        base = [
            "create",
            "pairs",
            "--pairs",
            tiny_pairs,
            "--bins",
            tiny_bins_bed,
            "-o",
            out,
        ]
        assert _run(*base, "--contigs", tiny_contigs_table).exit_code == 0
        ps = tmp_path / "out.ps.tsv"
        before_map, before_ps = out.read_bytes(), ps.read_bytes()

        bad = tmp_path / "bad.tsv"
        bad.write_text("contig_name\tgenome_name\ncontig_A\tMAG_é\n")
        result = _run(*base, "--contigs", bad, "-f")
        assert result.exit_code != 0
        assert out.read_bytes() == before_map
        assert ps.read_bytes() == before_ps
        validate(out)

    def test_successful_recreate_still_overwrites(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ):
        out = tmp_path / "out.metacool"
        args = ["create", "pairs", "--pairs", tiny_pairs, "--bins", tiny_bins_bed]
        args += ["--contigs", tiny_contigs_table, "-o", out]
        assert _run(*args).exit_code == 0
        result = _run(*args, "-f", "--assembly-name", "second")
        assert result.exit_code == 0, result.output
        import metacooler

        assert metacooler.MetaCooler(out).info["assembly"] == "second"
