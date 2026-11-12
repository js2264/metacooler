"""Genome/contig retrieval through the offset indexes (``api.py``).

Every selector is checked against a brute-force answer computed from the full
bin and pixel tables, under a genome-grouped ordering and under one that
scatters a genome, and with the many-range and gap-merging read paths forced.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import metacooler.api as api
from metacooler.api import MetaCooler, _in_ranges, _merge_ranges, _ranges_to_ids
from metacooler.io._utils import chrom_codes
from metacooler.io.create import from_pairs
from metacooler.schema import MetaCoolerFormatError, validate


@pytest.fixture
def scattered_metacooler(
    tmp_path: Path, tiny_fasta: Path, tiny_bins_bed: Path, tiny_pairs: Path
) -> Path:
    """MAG_1 owns contig_A and contig_C, which ``fasta`` order separates."""
    table = tmp_path / "scattered.tsv"
    table.write_text(
        "contig_name\tgenome_name\ncontig_A\tMAG_1\ncontig_B\tMAG_2\ncontig_C\tMAG_1\n"
    )
    out = tmp_path / "scattered.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=table,
        contig_table_format="contigs",
        fasta=tiny_fasta,
        contig_ordering="fasta",
        chunksize=1000,
    )
    return out


def _brute_bin_ids(mc: MetaCooler, genome: str) -> np.ndarray:
    gid = int(np.flatnonzero(mc.genomes["name"] == genome)[0])
    cids = mc.contigs.index[mc.contigs["genome_id"] == gid]
    return np.flatnonzero(np.isin(chrom_codes(mc.bins()), cids))


def _brute_pair(mc: MetaCooler, g1: str, g2: str) -> pd.DataFrame:
    a, b = _brute_bin_ids(mc, g1), _brute_bin_ids(mc, g2)
    px = mc.pixels()
    b1, b2 = px["bin1_id"], px["bin2_id"]
    sel = (b1.isin(a) & b2.isin(b)) | (b1.isin(b) & b2.isin(a))
    return px[sel].reset_index(drop=True)


@pytest.fixture(params=["grouped", "scattered"])
def any_metacooler(request, tiny_metacooler: Path, scattered_metacooler: Path):
    return tiny_metacooler if request.param == "grouped" else scattered_metacooler


@pytest.fixture(params=["default", "forced"])
def read_paths(request, monkeypatch):
    """``forced``: bulk reads for every range set, merged pixel blocks, tiny chunks."""
    if request.param == "forced":
        monkeypatch.setattr(api, "_MAX_RANGE_READS", 0)
        monkeypatch.setattr(api, "_PIXEL_READ_GAP", 10**9)
    return request.param


class TestRangeHelpers:
    def test_merge_ranges(self):
        lo = np.array([10, 0, 5, 30, 20, 40])
        hi = np.array([12, 5, 8, 30, 25, 41])
        np.testing.assert_array_equal(
            _merge_ranges(lo, hi), [[0, 8], [10, 12], [20, 25], [40, 41]]
        )
        np.testing.assert_array_equal(
            _merge_ranges(lo, hi, gap=2), [[0, 12], [20, 25], [40, 41]]
        )
        assert _merge_ranges(np.array([3]), np.array([3])).shape == (0, 2)

    def test_in_ranges(self):
        ranges = np.array([[2, 4], [7, 9]])
        values = np.arange(11)
        np.testing.assert_array_equal(
            np.flatnonzero(_in_ranges(values, ranges)), [2, 3, 7, 8]
        )
        np.testing.assert_array_equal(
            np.flatnonzero(_in_ranges(values, ranges[:1])), [2, 3]
        )
        assert not _in_ranges(values, np.empty((0, 2), dtype=np.int64)).any()

    def test_ranges_to_ids(self):
        np.testing.assert_array_equal(
            _ranges_to_ids(np.array([[2, 4], [7, 9], [9, 10]])), [2, 3, 7, 8, 9]
        )
        assert _ranges_to_ids(np.empty((0, 2), dtype=np.int64)).size == 0


@pytest.mark.usefixtures("read_paths")
class TestSelectorsMatchBruteForce:
    def test_genome_selectors(self, any_metacooler: Path):
        mc = MetaCooler(any_metacooler)
        for genome in mc.genomes["name"]:
            expected = _brute_bin_ids(mc, genome)
            np.testing.assert_array_equal(mc.genome_bin_ids(genome), expected)
            pd.testing.assert_frame_equal(
                mc.bins_for_genome(genome), mc.bins().iloc[expected]
            )

    def test_genome_pairs(self, any_metacooler: Path):
        mc = MetaCooler(any_metacooler)
        names = mc.genomes["name"].tolist()
        for g1 in names:
            for g2 in names:
                pd.testing.assert_frame_equal(
                    mc.pixels_for_genome_pair(g1, g2), _brute_pair(mc, g1, g2)
                )

    def test_contig_selectors(self, any_metacooler: Path):
        mc = MetaCooler(any_metacooler)
        codes = chrom_codes(mc.bins())
        for cid, contig in enumerate(mc.chromnames):
            ids = np.flatnonzero(codes == cid)
            assert mc.contig_extent(contig) == (int(ids[0]), int(ids[-1]) + 1)
            pd.testing.assert_frame_equal(
                mc.bins_for_contig(contig), mc.bins().iloc[ids]
            )

    def test_chunked_reads(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        rows = np.array([[0, mc.shape[0]]])
        chunks = list(mc._iter_pixel_rows(rows, chunksize=7))
        assert max(len(c) for c in chunks) <= 7
        pd.testing.assert_frame_equal(pd.concat(chunks), mc.pixels())


class TestExtents:
    def test_scattered_genome_extent_raises(self, scattered_metacooler: Path):
        mc = MetaCooler(scattered_metacooler)
        assert len(mc.genome_bin_ranges("MAG_1")) == 2
        with pytest.raises(ValueError, match="not contiguous"):
            mc.genome_extent("MAG_1")
        # A genome that happens to be contiguous still has an extent.
        lo, hi = mc.genome_extent("MAG_2")
        np.testing.assert_array_equal(np.arange(lo, hi), _brute_bin_ids(mc, "MAG_2"))

    def test_unknown_names(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="Genome 'nope'"):
            mc.genome_bin_ids("nope")
        with pytest.raises(KeyError, match="Contig 'nope'"):
            mc.contig_extent("nope")

    def test_empty_pair_keeps_columns(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        empty = mc._empty_pixels()
        assert list(empty.columns) == list(mc.pixels().columns)
        assert len(empty) == 0

    def test_no_rows_reads_nothing(self, tiny_metacooler: Path):
        mc = MetaCooler(tiny_metacooler)
        assert list(mc._iter_pixel_rows(np.array([[3, 3]]))) == []
        bins = mc._bins_in_ranges(np.empty((0, 2), dtype=np.int64))
        assert len(bins) == 0
        assert list(bins.columns) == list(mc.bins().columns)


class TestStructuralOpen:
    """An open checks structure only; ``validate()`` still checks contents."""

    def _corrupt_bins(self, path: Path) -> None:
        with h5py.File(path, "r+") as f:
            start = f["bins"]["start"][:]
            start[1] = -5
            f["bins"]["start"][:] = start

    def test_open_skips_content_checks(self, tiny_metacooler: Path):
        self._corrupt_bins(tiny_metacooler)
        MetaCooler(tiny_metacooler)  # structural checks pass
        validate(tiny_metacooler, full=False)
        with pytest.raises(MetaCoolerFormatError, match="non-negative"):
            validate(tiny_metacooler)

    def test_open_checks_offset_endpoints(self, tiny_metacooler: Path):
        with h5py.File(tiny_metacooler, "r+") as f:
            off = f["indexes"]["bin1_offset"][:]
            off[-1] += 1
            f["indexes"]["bin1_offset"][:] = off
        with pytest.raises(MetaCoolerFormatError, match="bin1_offset endpoints"):
            MetaCooler(tiny_metacooler)

    def test_open_checks_offset_presence(self, tiny_metacooler: Path):
        with h5py.File(tiny_metacooler, "r+") as f:
            del f["indexes"]["genome_offset"]
        with pytest.raises(MetaCoolerFormatError, match="genome_offset is missing"):
            MetaCooler(tiny_metacooler)

    def test_open_checks_offset_length(self, tiny_metacooler: Path):
        with h5py.File(tiny_metacooler, "r+") as f:
            del f["indexes"]["genome_offset"]
            f["indexes"].create_dataset("genome_offset", data=np.array([0]))
        with pytest.raises(MetaCoolerFormatError, match="genome_offset length"):
            MetaCooler(tiny_metacooler)


class TestDumpAndExtractRowBlocks:
    def test_dump_range_matches_full_scan(self, any_metacooler: Path, tmp_path: Path):
        from metacooler.io.dump import dump_table

        mc = MetaCooler(any_metacooler)
        px = mc.pixels()
        for contig in mc.chromnames:
            out = tmp_path / f"{contig}.tsv"
            dump_table(
                any_metacooler,
                "pixels",
                range1=contig,
                range2="contig_A",
                header=True,
                out=out,
            )
            i0, i1 = mc.contig_extent(contig)
            j0, j1 = mc.contig_extent("contig_A")
            b1, b2 = px["bin1_id"], px["bin2_id"]
            sel = (b1.between(i0, i1 - 1) & b2.between(j0, j1 - 1)) | (
                b1.between(j0, j1 - 1) & b2.between(i0, i1 - 1)
            )
            got = pd.read_csv(out, sep="\t")
            pd.testing.assert_frame_equal(
                got, px[sel].reset_index(drop=True), check_dtype=False
            )

    def test_dump_selection_without_range(
        self, scattered_metacooler: Path, tmp_path: Path
    ):
        from metacooler.io.dump import dump_table

        mc = MetaCooler(scattered_metacooler)
        out = tmp_path / "mag1.tsv"
        dump_table(
            scattered_metacooler, "pixels", genomes=["MAG_1"], header=True, out=out
        )
        pd.testing.assert_frame_equal(
            pd.read_csv(out, sep="\t"),
            _brute_pair(mc, "MAG_1", "MAG_1"),
            check_dtype=False,
        )

    @pytest.mark.parametrize("join", [False, True])
    def test_dump_empty_selection_keeps_header(
        self, tiny_metacooler: Path, tmp_path: Path, join: bool, monkeypatch
    ):
        from metacooler.io.dump import dump_table

        # Nothing matches: no chunk is ever produced.
        monkeypatch.setattr(MetaCooler, "_iter_pixel_rows", lambda *a, **k: iter(()))
        out = tmp_path / "empty.tsv"
        dump_table(
            tiny_metacooler,
            "pixels",
            range1="contig_D:1-1",
            range2="contig_A:1-1",
            header=True,
            join=join,
            out=out,
        )
        text = out.read_text().splitlines()
        assert len(text) == 1
        assert ("contig1" in text[0]) == join

    def test_extract_square_source(self, tiny_metacooler: Path, tmp_path: Path):
        """A ``square`` source takes the folding, merging write path."""
        import cooler

        from metacooler.io.extract import extract

        plain = tmp_path / "plain.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=plain)
        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs["storage-mode"] = "square"
        square = tmp_path / "square.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=square)
        pd.testing.assert_frame_equal(
            pd.DataFrame(cooler.Cooler(str(plain)).pixels()[:]),
            pd.DataFrame(cooler.Cooler(str(square)).pixels()[:]),
        )
