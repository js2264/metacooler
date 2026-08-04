"""Tests for binnify operation."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.api import MetaCooler
from metacooler.ops.binnify import binnify


class TestBinnify:
    def test_binnify_doubles_binsize(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        assert out.exists()
        mc2 = MetaCooler(out)
        assert mc2.shape[0] > 0
        # Should have fewer bins
        assert mc2.shape[0] < mc.shape[0]

    def test_binnify_preserves_genomes(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.genomes) == len(mc.genomes)
        assert set(mc2.genomes["name"]) == set(mc.genomes["name"])

    def test_binnify_preserves_contigs(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.contigs) == len(mc.contigs)

    def test_binnify_smaller_binsize_raises(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "fail.metacool"
        with pytest.raises(ValueError, match="must be strictly larger"):
            binnify(tiny_metacooler, out, binsize=original_binsize // 2)

    def test_binnify_same_binsize_raises(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "fail.metacool"
        with pytest.raises(ValueError, match="must be strictly larger"):
            binnify(tiny_metacooler, out, binsize=original_binsize)

    def test_binnify_has_pixels(self, tmp_path: Path, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        original_binsize = mc.binsize
        assert original_binsize is not None

        out = tmp_path / "rebinned.metacool"
        binnify(tiny_metacooler, out, binsize=original_binsize * 2)

        mc2 = MetaCooler(out)
        assert len(mc2.pixels()) > 0

    def test_binnify_genome_level_ligates_to_one_contig_per_genome(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        mc = MetaCooler(tiny_metacooler)
        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        mc2 = MetaCooler(out)
        # Output "contigs" are the ligated genomes: one per genome, named after it.
        assert len(mc2.contigs) == len(mc.genomes)
        assert set(mc2.contigs["name"]) == set(mc.genomes["name"])
        assert set(mc2.genomes["n_contigs"].unique()) == {1}
        # Each ligated contig length equals its genome's total contig length.
        for _, row in mc2.contigs.iterrows():
            gid = int(mc.genomes.index[mc.genomes["name"] == row["name"]][0])
            expected = int(
                mc.contigs.loc[mc.contigs["genome_id"] == gid, "length"].sum()
            )
            assert int(row["length"]) == expected

    def test_binnify_genome_level_preserves_total_counts(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        import cooler

        mc = MetaCooler(tiny_metacooler)
        src_sum = int(mc._clr.pixels()[:]["count"].sum())

        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        c = cooler.Cooler(str(out))
        assert int(c.pixels()[:]["count"].sum()) == src_sum

    def test_binnify_genome_level_is_valid_metacool(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        import cooler

        from metacooler.schema import validate

        out = tmp_path / "genome.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")

        validate(out)  # raises if the schema is violated
        assert cooler.fileops.is_cooler(str(out))  # still openable by cooler


class TestMapOldBinsPerContig:
    """The vectorised bin mapper must agree with the per-contig loop it replaced.

    That loop rebuilt an ``old_chrom == cid`` mask per contig, making the whole
    function O(n_contigs x n_bins) -- unrunnable past a few tens of thousands of
    contigs. These cases pin the behaviour the rewrite has to preserve.
    """

    @staticmethod
    def _reference(old_bins, bin_chrom_codes, n_contigs, new_bins):
        """The original per-contig implementation, kept as the oracle."""
        import numpy as np

        old_midpoints = (old_bins["start"].values + old_bins["end"].values) / 2
        lookup = {}
        for cid, grp in new_bins.groupby("chrom", sort=False):
            lookup[int(cid)] = (
                grp["start"].values.astype(np.int64),
                grp.index.values.astype(np.int64),
            )
        out = np.full(len(old_bins), -1, dtype=np.int64)
        for cid in range(n_contigs):
            mask = bin_chrom_codes == cid
            if not mask.any():
                continue
            starts, ids = lookup.get(
                cid, (np.array([], np.int64), np.array([], np.int64))
            )
            if len(starts) == 0:
                continue
            idx = np.searchsorted(starts, old_midpoints[mask], side="right") - 1
            valid = idx >= 0
            pos = np.where(mask)[0]
            out[pos[valid]] = ids[idx[valid]]
        return out

    @pytest.mark.parametrize(
        ("lengths", "label"),
        [
            ([10_000] * 20, "uniform contigs"),
            ([1, 999, 12_345, 50_000, 7], "ragged contigs"),
            ([500, 900, 1_500], "contigs shorter than one target bin"),
            ([0, 4_000, 0, 9_000], "zero-length contigs interleaved"),
            ([123_456], "single contig"),
            ([0, 0, 0], "every contig zero-length"),
        ],
    )
    @pytest.mark.parametrize(
        ("src", "tgt"), [(1_000, 5_000), (500, 10_000), (1_000, 1_000)]
    )
    def test_matches_per_contig_loop(self, lengths, label, src, tgt):
        import numpy as np

        from metacooler.ops.binnify import _map_old_bins_per_contig, _tile_uniform

        arr = np.asarray(lengths, dtype=np.int64)
        old_bins = _tile_uniform(arr, src)
        codes = old_bins["chrom"].values.astype(np.int64)
        new_bins = _tile_uniform(arr, tgt)
        expected = self._reference(old_bins, codes, len(arr), new_bins)
        got = _map_old_bins_per_contig(old_bins, codes, len(arr), new_bins)
        assert np.array_equal(got, expected), f"{label} {src}->{tgt}"

    def test_midpoint_past_stale_contig_length_clamps_to_last_bin(self):
        """A bin reaching past /contigs/length must not leak into the next contig."""
        import numpy as np
        import pandas as pd

        from metacooler.ops.binnify import _map_old_bins_per_contig, _tile_uniform

        old_bins = pd.DataFrame(
            {"chrom": [0, 0, 1], "start": [0, 900, 0], "end": [100, 5_000, 10]}
        )
        codes = np.array([0, 0, 1], dtype=np.int64)
        new_bins = _tile_uniform(np.array([200, 50]), 100)
        got = _map_old_bins_per_contig(old_bins, codes, 2, new_bins)
        # contig 0 owns new bins 0-1, contig 1 owns new bin 2.
        assert got.tolist() == [0, 1, 2]
