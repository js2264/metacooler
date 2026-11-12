"""Row-block re-binning (``ops/binnify``) and zoomify's level-from-level ladder."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import metacooler.ops.binnify as binnify_mod
from metacooler.api import MetaCooler
from metacooler.ops.binnify import _row_blocks, binnify
from metacooler.ops.zoomify import zoomify
from metacooler.schema import validate


def _pixels(uri) -> pd.DataFrame:
    return MetaCooler(uri).pixels()[["bin1_id", "bin2_id", "count"]]


class TestRowBlocks:
    def test_blocks_cover_whole_output_rows(self):
        old_to_new = np.array([0, 0, 1, 1, 2, 2, 3])
        bin1_offset = np.array([0, 3, 6, 9, 12, 15, 18, 21])
        blocks = _row_blocks(old_to_new, bin1_offset, chunksize=7)
        assert blocks[0][0] == 0 and blocks[-1][1] == len(old_to_new)
        for (_, hi), (lo, _) in zip(blocks[:-1], blocks[1:], strict=True):
            assert hi == lo
            assert old_to_new[lo] != old_to_new[lo - 1]  # cut between rows

    def test_dropped_bin_never_splits_a_row(self):
        old_to_new = np.array([0, 0, -1, 0, 1, 1])
        bin1_offset = np.arange(0, 70, 10)
        blocks = _row_blocks(old_to_new, bin1_offset, chunksize=10)
        assert (0, 4) in blocks  # output row 0 (with its dropped bin) stays whole

    def test_one_huge_row_is_one_block(self):
        old_to_new = np.array([0, 0, 0, 1])
        bin1_offset = np.array([0, 100, 200, 300, 301])
        assert _row_blocks(old_to_new, bin1_offset, chunksize=10) == [(0, 3), (3, 4)]

    def test_empty(self):
        assert _row_blocks(np.array([], dtype=np.int64), np.array([0]), 10) == []


class TestOrderedRebinning:
    @pytest.mark.parametrize("chunksize", [1, 7, 10_000])
    def test_row_blocks_match_unordered_path(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch, chunksize: int
    ):
        ordered = tmp_path / "ordered.metacool"
        binnify(tiny_metacooler, ordered, binsize=300, chunksize=chunksize)

        monkeypatch.setattr(binnify_mod, "_is_monotonic", lambda values: False)
        unordered = tmp_path / "unordered.metacool"
        binnify(tiny_metacooler, unordered, binsize=300, chunksize=chunksize)

        pd.testing.assert_frame_equal(_pixels(ordered), _pixels(unordered))
        validate(ordered)

    def test_square_source_takes_the_unordered_path(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch
    ):
        reference = tmp_path / "reference.metacool"
        binnify(tiny_metacooler, reference, binsize=300)
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["storage-mode"] = "square"
        called = []
        real = binnify_mod._remap_pixel_chunks
        monkeypatch.setattr(
            binnify_mod,
            "_remap_pixel_chunks",
            lambda *a, **k: called.append(1) or real(*a, **k),
        )
        out = tmp_path / "square.metacool"
        binnify(tiny_metacooler, out, binsize=300)
        assert called
        pd.testing.assert_frame_equal(_pixels(out), _pixels(reference))

    def test_zero_pixel_source(self, tmp_path: Path, tiny_bins_bed, tiny_fasta):
        from metacooler.io.create import from_pairs

        pairs = tmp_path / "none.pairs"
        pairs.write_text(
            "#columns: readID chrom1 pos1 chrom2 pos2\nr1\tnope\t1\tnope\t2\n"
        )
        table = tmp_path / "t.tsv"
        table.write_text("contig_name\tgenome_name\ncontig_A\tMAG_1\n")
        src = tmp_path / "empty.metacool"
        from_pairs(
            uri=src,
            bins=tiny_bins_bed,
            pairs=pairs,
            contig_genome_map=table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
        )
        out = tmp_path / "coarse.metacool"
        binnify(src, out, binsize=300)
        assert MetaCooler(out).info["nnz"] == 0
        validate(out)


class TestZoomifyLadder:
    def test_each_level_equals_rebinning_the_base(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        out = tmp_path / "z.mmetacool"
        zoomify(tiny_metacooler, out, resolutions="100,200,400,800")
        for res in (200, 400, 800):
            direct = tmp_path / f"direct_{res}.metacool"
            binnify(tiny_metacooler, direct, binsize=res)
            entry = f"{out}::entries/{res}"
            pd.testing.assert_frame_equal(_pixels(entry), _pixels(direct))
            pd.testing.assert_frame_equal(
                MetaCooler(entry).bins(), MetaCooler(direct).bins()
            )
            with h5py.File(out, "r") as f:
                meta = json.loads(f[f"entries/{res}"].attrs["metadata"])
            # Provenance names the base, as when every level came from it.
            assert meta["binnify"] == [
                {"source": str(tiny_metacooler), "binsize": res, "level": "contig"}
            ]

    def test_levels_come_from_the_nearest_finer_level(
        self, tiny_metacooler: Path, tmp_path: Path, monkeypatch
    ):
        parents: list[tuple[int | None, int]] = []
        real = binnify_mod._binnify

        def spy(mc, output, binsize, **kwargs):
            parents.append((mc.binsize, binsize))
            return real(mc, output, binsize, **kwargs)

        monkeypatch.setattr(binnify_mod, "_binnify", spy)
        zoomify(
            tiny_metacooler, tmp_path / "z.mmetacool", resolutions="100,200,500,1000"
        )
        assert parents == [(100, 200), (100, 500), (500, 1000)]


def test_dropped_source_bins_lose_their_pixels(tiny_metacooler: Path):
    """A ``-1`` in the bin map drops every pixel touching that source bin."""
    from metacooler.ops.binnify import _coarsen_row_blocks

    mc = MetaCooler(tiny_metacooler)
    n = mc.shape[0]
    old_to_new = np.arange(n, dtype=np.int64) // 2
    old_to_new[[1, 4]] = -1
    got = pd.concat(list(_coarsen_row_blocks(mc, old_to_new, chunksize=5)))

    px = mc.pixels()
    b1 = old_to_new[px["bin1_id"].to_numpy()]
    b2 = old_to_new[px["bin2_id"].to_numpy()]
    keep = (b1 >= 0) & (b2 >= 0)
    expected = (
        pd.DataFrame(
            {"bin1_id": b1[keep], "bin2_id": b2[keep], "count": px["count"][keep]}
        )
        .groupby(["bin1_id", "bin2_id"], as_index=False)["count"]
        .sum()
    )
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), expected, check_dtype=False
    )
