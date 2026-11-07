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

import cooler
import h5py
import numpy as np
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io._utils import check_output_path, cooler_writer_module
from metacooler.io.merge import merge
from metacooler.ops.binnify import binnify
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
        assert MetaCooler(out).info["assembly"] == "second"


def _scale_counts_to_float(path: Path, factor: float) -> None:
    """Rewrite ``pixels/count`` as float64 times *factor* (a rescaled map)."""
    with h5py.File(path, "a") as f:
        counts = f["pixels"]["count"][:].astype(np.float64) * factor
        del f["pixels"]["count"]
        f["pixels"].create_dataset("count", data=counts)
        f.attrs["sum"] = float(counts.sum())


def _pixel_sum(uri: str | Path) -> float:
    return float(cooler.Cooler(str(uri)).pixels()[:]["count"].sum())


class TestMerge:
    def test_float_counts_are_not_truncated(self, tmp_path: Path, tiny_metacooler):
        """``merge`` forced ``count`` to int64, flooring every float count."""
        a, b = tmp_path / "a.metacool", tmp_path / "b.metacool"
        for p in (a, b):
            shutil.copy(tiny_metacooler, p)
            _scale_counts_to_float(p, 0.3)
        out = tmp_path / "m.metacool"
        merge([a, b], out)
        assert _pixel_sum(out) == pytest.approx(_pixel_sum(a) + _pixel_sum(b))

    def test_integer_counts_widen_to_int64(self, tmp_path: Path, tiny_metacooler):
        out = tmp_path / "m.metacool"
        merge([tiny_metacooler, tiny_metacooler], out)
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"].dtype == np.int64

    def test_cli_does_not_restore_the_stale_qc_block(
        self, tmp_path: Path, tiny_metacooler
    ):
        """The CLI re-merged every input's metadata, putting ``qc`` back.

        ``io.merge`` strips the block (it describes one input, not the sum), but
        the CLI then copied every input key into the output: a merged file of
        twice the contacts reported the first input's ``n_contacts``.
        """
        assert "qc" in MetaCooler(tiny_metacooler).metadata
        out = tmp_path / "m.metacool"
        result = _run("merge", tiny_metacooler, tiny_metacooler, "-o", out)
        assert result.exit_code == 0, result.output
        meta = MetaCooler(out).metadata
        assert "qc" not in meta
        assert "merge" in meta["commands"][-1]

    def test_bin_size_survives_single_bin_contigs(
        self, tmp_path: Path, tiny_metacooler
    ):
        """With every contig one bin, cooler inferred ``bin-size=None``.

        The merged file then claimed variable bins, so ``zoomify`` refused it
        and ``binnify`` applied the variable-bin resolution rule.
        """
        coarse = tmp_path / "coarse.metacool"
        binnify(tiny_metacooler, coarse, binsize=2000)
        assert MetaCooler(coarse).info["bin-size"] == 2000
        out = tmp_path / "m.metacool"
        merge([coarse, coarse], out)
        info = MetaCooler(out).info
        assert info["bin-size"] == 2000
        assert info["bin-type"] == "fixed"

    def test_entry_uris_are_accepted(self, tmp_path: Path, tiny_mmetacool: Path):
        """The CLI accepts entry URIs, but ``merge`` opened them as plain paths."""
        e1 = f"{tiny_mmetacool}::entries/entry1"
        e2 = f"{tiny_mmetacool}::entries/entry2"
        out = tmp_path / "m.metacool"
        result = _run("merge", e1, e2, "-o", out)
        assert result.exit_code == 0, result.output
        assert _pixel_sum(out) == _pixel_sum(e1) + _pixel_sum(e2)
        validate(out)

    def test_goes_through_the_writer_shims(
        self, tmp_path: Path, tiny_metacooler, monkeypatch
    ):
        """``merge_coolers`` ran ``cooler.create`` without the shims.

        Every other writer suppresses cooler's per-contig enum build and
        ``get_binsize`` loop; merge alone paid both, taking 65 s at 150 k
        contigs where every other command took under 10 s.
        """
        module = cooler_writer_module()
        real = module.write_bins
        seen: list[bool] = []

        def spy(*args, chrom_as_enum=True, **kwargs):
            seen.append(chrom_as_enum)
            return real(*args, chrom_as_enum=chrom_as_enum, **kwargs)

        def forbidden(bins):
            raise AssertionError("cooler's per-contig get_binsize was called")

        monkeypatch.setattr(module, "write_bins", spy)
        monkeypatch.setattr(module, "get_binsize", forbidden)
        merge([tiny_metacooler, tiny_metacooler], tmp_path / "m.metacool")
        assert seen and not any(seen)

    def test_api_refuses_an_input_as_output(self, tmp_path: Path, tiny_metacooler):
        a = tmp_path / "a.metacool"
        shutil.copy(tiny_metacooler, a)
        before = a.read_bytes()
        with pytest.raises(ValueError, match="same file as the input"):
            merge([a, tiny_metacooler], a)
        assert a.read_bytes() == before
