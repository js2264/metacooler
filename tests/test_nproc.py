"""``--nproc``: parallel work must not change a single output byte that matters."""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import metacooler.io.create as create_mod
from metacooler.io.create import _map_in_order, _PairsDropStats, from_pairs
from metacooler.ops.zoomify import _parents, zoomify


def _tables(path: Path, root: str = "/") -> dict:
    out = {}
    with h5py.File(path, "r") as f:
        grp = f[root]

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                out[name] = obj[()].tolist()

        grp.visititems(visit)
        meta = json.loads(grp.attrs["metadata"])
        meta.pop("commands", None)
        out["metadata"] = meta
    return out


class TestMapInOrder:
    def test_results_come_back_in_submission_order(self):
        def slow_then_fast(x):
            time.sleep(0.05 if x % 2 == 0 else 0)
            return x * 10

        assert list(_map_in_order(slow_then_fast, iter(range(9)), nproc=3)) == [
            x * 10 for x in range(9)
        ]

    def test_serial(self):
        assert list(_map_in_order(lambda x: x + 1, [1, 2], nproc=1)) == [2, 3]

    def test_worker_error_propagates(self):
        def boom(x):
            raise RuntimeError(f"bad {x}")

        with pytest.raises(RuntimeError, match="bad 0"):
            list(_map_in_order(boom, iter(range(4)), nproc=2))


def test_drop_stats_merge_keeps_the_first_names_in_order():
    total = _PairsDropStats()
    for names in (["a", "b"], ["c"], [f"x{i}" for i in range(20)]):
        batch = _PairsDropStats(total=5, unknown_contig=len(names))
        batch.note_unknown_names(names)
        total.merge(batch)
    assert total.total == 15
    assert total.unknown_contig == 23
    assert list(total.unknown_names) == ["a", "b", "c", *[f"x{i}" for i in range(7)]]


@pytest.fixture
def small_batches(monkeypatch):
    """Re-cut the pairs stream into 7-row batches, so several are in flight."""
    real = create_mod._stream_pairs

    def rechunked(*args, **kwargs):
        for chunk in real(*args, **kwargs):
            for lo in range(0, len(chunk), 7):
                yield chunk.iloc[lo : lo + 7]

    monkeypatch.setattr(create_mod, "_stream_pairs", rechunked)


@pytest.mark.usefixtures("small_batches")
@pytest.mark.parametrize("sort_in_memory", [True, False])
@pytest.mark.parametrize("one_bin_per_contig", [False, True])
def test_create_pairs_nproc_is_invisible(
    tmp_path: Path,
    tiny_bins_bed: Path,
    tiny_contigs_table: Path,
    tiny_fasta: Path,
    tiny_pairs: Path,
    caplog,
    sort_in_memory: bool,
    one_bin_per_contig: bool,
):
    # Unknown contigs and out-of-range positions scattered across batches.
    lines = tiny_pairs.read_text().splitlines(keepends=True)
    for i in range(3, len(lines), 9):
        lines[i] = f"bad_{i}\tcontig_Z{i}\t5\tcontig_A\t10\t+\t+\n"
    for i in range(5, len(lines), 13):
        lines[i] = f"oob_{i}\tcontig_D\t9999\tcontig_A\t10\t+\t+\n"
    pairs = tmp_path / "mixed.pairs"
    pairs.write_text("".join(lines))

    messages: dict[int, list[str]] = {}
    tables: dict[int, dict] = {}
    ps: dict[int, pd.DataFrame] = {}
    for nproc in (1, 3):
        out = tmp_path / f"n{nproc}.metacool"
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            from_pairs(
                uri=out,
                bins=None if one_bin_per_contig else tiny_bins_bed,
                pairs=pairs,
                contig_genome_map=tiny_contigs_table,
                contig_table_format="contigs",
                fasta=tiny_fasta,
                sort_in_memory=sort_in_memory,
                nproc=nproc,
            )
        messages[nproc] = [r.getMessage() for r in caplog.records]
        tables[nproc] = _tables(out)
        ps[nproc] = pd.read_csv(out.with_suffix(".ps.tsv"), sep="\t")
    assert tables[1] == tables[3]
    pd.testing.assert_frame_equal(ps[1], ps[3])
    discarded = [m for m in messages[1] if m.startswith("Discarded")]
    assert discarded
    assert discarded == [m for m in messages[3] if m.startswith("Discarded")]


def test_create_pairs_rejects_zero_nproc(tmp_path: Path):
    with pytest.raises(ValueError, match="nproc must be at least 1"):
        from_pairs(uri=tmp_path / "x.metacool", nproc=0)


class TestZoomifyNproc:
    def test_parents(self):
        assert _parents([1000, 2000, 5000, 10000, 20000]) == {
            2000: 1000,
            5000: 1000,
            10000: 5000,
            20000: 10000,
        }

    def test_parallel_ladder_matches_serial(self, tiny_metacooler: Path, tmp_path):
        ladder = "100,200,500,1000"
        serial = tmp_path / "serial.mmetacool"
        zoomify(tiny_metacooler, serial, resolutions=ladder)
        scratch = tmp_path / "scratch"
        parallel = tmp_path / "parallel.mmetacool"
        zoomify(
            tiny_metacooler, parallel, resolutions=ladder, nproc=2, temp_dir=scratch
        )
        with h5py.File(serial, "r") as a, h5py.File(parallel, "r") as b:
            assert list(a["entries"]) == list(b["entries"])
            names = list(a["entries"])
        for name in names:
            assert _tables(serial, f"entries/{name}") == _tables(
                parallel, f"entries/{name}"
            )
        assert list(scratch.iterdir()) == []  # every level file was cleaned up

    def test_rejects_zero_nproc(self, tiny_metacooler: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="nproc must be at least 1"):
            zoomify(tiny_metacooler, tmp_path / "z.mmetacool", nproc=0)

    # The workers only raise, so the fork-while-threaded hazard is moot here.
    @pytest.mark.filterwarnings("ignore:This process .* is multi-threaded")
    def test_worker_failure_removes_the_output(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch
    ):
        import multiprocessing

        import metacooler.ops.zoomify as zoomify_mod

        # fork, so that the patched level builder is the one the workers run.
        real_context = multiprocessing.get_context
        monkeypatch.setattr(
            multiprocessing, "get_context", lambda method=None: real_context("fork")
        )
        monkeypatch.setattr(zoomify_mod, "_build_level", _failing_build_level)
        out = tmp_path / "z.mmetacool"
        with pytest.raises(RuntimeError, match=r"level \d+ failed"):
            zoomify(tiny_metacooler, out, resolutions="100,200,500", nproc=2)
        assert not out.exists()


def _failing_build_level(parent, output, binsize, source_uri, source_metadata):
    raise RuntimeError(f"level {binsize} failed")


def test_cli_rejects_zero_nproc(tiny_metacooler: Path, tmp_path: Path):
    from metacooler.cli.main import cli

    result = CliRunner().invoke(
        cli,
        [
            "zoomify",
            str(tiny_metacooler),
            "-o",
            str(tmp_path / "z.mmetacool"),
            "--nproc",
            "0",
        ],
    )
    assert result.exit_code != 0
    assert "--nproc" in result.output


def test_build_level_in_process(tiny_metacooler: Path, tmp_path: Path):
    """The worker entry point, run here so that coverage sees it."""
    from metacooler.ops.binnify import binnify
    from metacooler.ops.zoomify import _build_level

    level = tmp_path / "level.metacool"
    _build_level(str(tiny_metacooler), str(level), 200, str(tiny_metacooler), {})
    direct = tmp_path / "direct.metacool"
    binnify(tiny_metacooler, direct, binsize=200)
    assert _tables(level)["pixels/count"] == _tables(direct)["pixels/count"]


def test_positions_to_int64():
    from metacooler.io.create import _positions_to_int64

    assert _positions_to_int64(np.array([1.0, 2.0]), "pos1").dtype == np.int64
    with pytest.raises(ValueError, match="non-integer positions"):
        _positions_to_int64(np.array([1.5]), "pos1")
    with pytest.raises(ValueError, match="non-numeric dtype"):
        _positions_to_int64(np.array(["a"]), "pos1")
