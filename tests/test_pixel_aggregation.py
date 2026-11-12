"""Key-based pixel aggregation (``io/_utils``) and dictionary-encoded pair names."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metacooler.io._utils import (
    aggregate_pixels,
    pixel_keys,
    pixels_from_keys,
    reduce_sorted_keys,
)
from metacooler.io.create import _contig_codes, _decode_pairs_chunk, _PairsDropStats


def _random_frames(seed: int, n_frames: int, rows: int, dtype=np.int32):
    rng = np.random.default_rng(seed)
    frames = []
    for _ in range(n_frames):
        b1 = rng.integers(0, 50, rows)
        b2 = b1 + rng.integers(0, 10, rows)
        frames.append(
            pd.DataFrame(
                {
                    "bin1_id": b1,
                    "bin2_id": b2,
                    "count": rng.integers(1, 5, rows).astype(dtype),
                }
            )
        )
    return frames


def _groupby_reference(frames):
    return (
        pd.concat(frames, ignore_index=True)
        .groupby(["bin1_id", "bin2_id"], sort=True, as_index=False)["count"]
        .sum()
    )


class TestPixelKeys:
    def test_round_trip_and_order(self):
        b1 = np.array([0, 3, 3, 2**31 - 1])
        b2 = np.array([5, 1, 7, 2**31 - 1])
        keys = pixel_keys(b1, b2)
        back = pixels_from_keys(keys, np.ones(4, dtype=np.int32))
        np.testing.assert_array_equal(back["bin1_id"], b1)
        np.testing.assert_array_equal(back["bin2_id"], b2)
        assert list(np.argsort(keys)) == [0, 1, 2, 3]

    @pytest.mark.parametrize("bad", [-1, 2**31])
    def test_out_of_range_refused(self, bad):
        with pytest.raises(OverflowError, match="2\\*\\*31"):
            pixel_keys(np.array([0, bad]), np.array([1, 1]))

    def test_empty(self):
        assert pixel_keys(np.array([]), np.array([])).dtype == np.int64
        keys, sums = reduce_sorted_keys(np.empty(0, np.int64), None)
        assert len(keys) == 0 and sums.dtype == np.int32
        keys, sums = reduce_sorted_keys(np.empty(0, np.int64), np.empty(0, np.int64))
        assert sums.dtype == np.int64

    def test_unit_counts_are_run_lengths(self):
        keys, sums = reduce_sorted_keys(np.array([1, 1, 2, 5, 5, 5]), None)
        np.testing.assert_array_equal(keys, [1, 2, 5])
        np.testing.assert_array_equal(sums, [2, 1, 3])
        assert sums.dtype == np.int32


class TestAggregatePixels:
    @pytest.mark.parametrize("dtype", [np.int32, np.int64, np.float64])
    @pytest.mark.parametrize("presorted", [False, True])
    def test_matches_groupby(self, dtype, presorted):
        frames = _random_frames(1, 4, 300, dtype)
        if presorted:
            frames = [aggregate_pixels([f]) for f in frames]
        got = aggregate_pixels(frames)
        expected = _groupby_reference(frames)
        pd.testing.assert_frame_equal(got, expected, check_dtype=False)
        assert got["count"].dtype == dtype

    def test_single_sorted_frame_passes_through(self):
        (frame,) = _random_frames(2, 1, 200)
        once = aggregate_pixels([frame])
        pd.testing.assert_frame_equal(aggregate_pixels([once]), once)

    def test_unique_frame_keeps_counts(self):
        frame = pd.DataFrame({"bin1_id": [2, 0], "bin2_id": [3, 1], "count": [7, 9]})
        got = aggregate_pixels([frame])
        assert got["count"].tolist() == [9, 7]

    def test_consume_empties_the_list(self):
        frames = _random_frames(3, 3, 50)
        expected = _groupby_reference(frames)
        got = aggregate_pixels(frames, consume=True)
        assert frames == []
        pd.testing.assert_frame_equal(got, expected, check_dtype=False)


class TestDictionaryEncodedNames:
    NAMES = {"a": 0, "b": 1}

    def test_categorical_and_plain_agree(self):
        raw = pd.Series(["a", "zz", "b", None, "a"])
        plain = _contig_codes(raw, self.NAMES)
        cat = _contig_codes(raw.astype("category"), self.NAMES)
        np.testing.assert_array_equal(plain, [0, -1, 1, -1, 0])
        np.testing.assert_array_equal(cat, plain)
        assert cat.dtype == np.int32

    def test_drop_stats_exact(self):
        chunk = pd.DataFrame(
            {
                "chrom1": pd.Series(["a", "zz", "b", "a", "a"], dtype="category"),
                "pos1": [1, 2, 3, 4, 500],
                "chrom2": pd.Series(["b", "a", "yy", "a", "b"], dtype="category"),
                "pos2": [5, 6, 7, 8, 9],
            }
        )
        stats = _PairsDropStats()
        out = _decode_pairs_chunk(
            chunk, self.NAMES, stats, contig_length=np.array([100, 100])
        )
        assert out is not None
        assert stats.total == 5
        assert stats.unknown_contig == 2
        assert stats.out_of_bounds == 1
        assert list(stats.unknown_names) == ["zz", "yy"]
        np.testing.assert_array_equal(out[0], [0, 0])


def test_from_cool_sums_duplicate_source_pixels(tmp_path):
    """``create cool`` only permutes bins, but a cooler written without cooler's
    duplicate check can still repeat a pixel; it is summed, never written twice.
    """
    import cooler

    from metacooler.api import MetaCooler
    from metacooler.io.create import from_cool
    from metacooler.schema import validate

    bins = pd.DataFrame(
        {"chrom": ["a", "a", "b"], "start": [0, 100, 0], "end": [100, 200, 100]}
    )
    pixels = pd.DataFrame(
        {"bin1_id": [0, 0, 0, 1], "bin2_id": [1, 1, 2, 2], "count": [2, 3, 1, 4]}
    )
    src = tmp_path / "dup.cool"
    cooler.create_cooler(str(src), bins, pixels, ordered=True, dupcheck=False)
    assert cooler.Cooler(str(src)).info["nnz"] == 4  # the repeat is stored

    table = tmp_path / "t.tsv"
    table.write_text("contig_name\tgenome_name\na\tG1\nb\tG2\n")
    out = tmp_path / "out.metacool"
    from_cool(src, out, table, contig_ordering="cool", qc=False)
    validate(out)
    got = MetaCooler(out).pixels()
    assert got[["bin1_id", "bin2_id", "count"]].values.tolist() == [
        [0, 1, 5],
        [0, 2, 1],
        [1, 2, 4],
    ]


def test_sort_kind_follows_run_length():
    from metacooler.io._utils import _sort_kind

    runs = np.concatenate([np.arange(i, i + 50) for i in range(0, 5000, 7)])
    assert _sort_kind(runs) == "stable"
    scrambled = np.random.default_rng(0).permutation(5000)
    assert _sort_kind(scrambled) == "quicksort"
    assert _sort_kind(scrambled, runs_known_sorted=True) == "stable"


class TestRemapBins:
    def test_remap_bins_by_names_or_categories(self):
        from metacooler.io.create import _remap_bins

        bins = pd.DataFrame(
            {
                "chrom": ["a", "a", "b", "c", "c"],
                "start": [0, 10, 0, 0, 10],
                "end": [10, 20, 5, 10, 15],
                "weight": [0.1, 0.2, 0.3, 0.4, 0.5],
            }
        )
        order = ["c", "a"]  # b is dropped
        for frame in (bins, bins.astype({"chrom": "category"})):
            remapped, old_to_new = _remap_bins(frame, order, keep_extra_columns=True)
            assert remapped["chrom"].tolist() == [0, 0, 1, 1]
            assert remapped["start"].tolist() == [0, 10, 0, 10]
            assert remapped["weight"].tolist() == [0.4, 0.5, 0.1, 0.2]
            assert old_to_new.tolist() == [2, 3, -1, 0, 1]
            assert remapped["chrom"].dtype == np.int32
            bare, _ = _remap_bins(frame, order)
            assert list(bare.columns) == ["chrom", "start", "end"]

    def test_from_cool_reads_a_grouped_uri(
        self, tmp_path, tiny_cool, tiny_contigs_table
    ):
        import cooler

        from metacooler.io.create import from_cool

        clr = cooler.Cooler(str(tiny_cool))
        grouped = f"{tmp_path / 'multi.mcool'}::/resolutions/100"
        cooler.create_cooler(grouped, clr.bins()[:], clr.pixels()[:], ordered=True)
        a, b = tmp_path / "plain.metacool", tmp_path / "grouped.metacool"
        from_cool(tiny_cool, a, tiny_contigs_table, qc=False)
        from_cool(grouped, b, tiny_contigs_table, qc=False)
        from metacooler.api import MetaCooler

        pd.testing.assert_frame_equal(MetaCooler(a).pixels(), MetaCooler(b).pixels())
        pd.testing.assert_frame_equal(MetaCooler(a).bins(), MetaCooler(b).bins())
