"""
Formal correctness tests for metacooler.io.create.from_pairs and from_cool.

Tests are structured in three layers:
1. Unit tests for individual internal helpers (_bin_pairs_chunk, _build_bins, etc.)
2. Round-trip invariant tests (total count conservation, pixel symmetry, index integrity)
3. Regression tests for known bugs (drop_unassigned crash, fragment underflow, fixed-bin overflow)

Run with:
    pytest tests/test_create_correctness.py -v
"""

from __future__ import annotations

from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from metacooler.io.create import (
    _bin_pairs_chunk,
    from_cool,
    from_pairs,
)
from metacooler.schema import is_metacooler, validate

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


def _write_pairs(path: Path, rows: list[tuple], chromsizes: dict[str, int]) -> None:
    """Write a minimal 4DN pairs file (uncompressed)."""
    with open(path, "w") as f:
        f.write("## pairs format v1.0\n")
        for name, size in chromsizes.items():
            f.write(f"#chromsize: {name} {size}\n")
        f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
        for i, (c1, p1, c2, p2) in enumerate(rows):
            f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")


def _write_contigs_table(path: Path, rows: list[tuple[str, str]]) -> None:
    """Write a contigs TSV (contig_name, genome_name)."""
    with open(path, "w") as f:
        f.write("contig_name\tgenome_name\n")
        for contig, genome in rows:
            f.write(f"{contig}\t{genome}\n")


def _write_cool(
    path: Path,
    chromsizes: dict[str, int],
    binsize: int,
    pairs: list[tuple[str, int, str, int]],
) -> None:
    """Create a minimal cooler file directly from a list of contacts."""
    chroms = list(chromsizes.keys())
    bins_rows = []
    for chrom in chroms:
        length = chromsizes[chrom]
        for start in range(0, length, binsize):
            bins_rows.append(
                {"chrom": chrom, "start": start, "end": min(start + binsize, length)}
            )
    bins_df = pd.DataFrame(bins_rows)

    # Build bin lookup
    chrom_offset = {}
    cumulative = 0
    for chrom in chroms:
        chrom_offset[chrom] = cumulative
        cumulative += sum(1 for b in bins_rows if b["chrom"] == chrom)

    pixel_rows = []
    for c1, p1, c2, p2 in pairs:
        b1 = chrom_offset[c1] + p1 // binsize
        b2 = chrom_offset[c2] + p2 // binsize
        if b1 > b2:
            b1, b2 = b2, b1
        pixel_rows.append({"bin1_id": b1, "bin2_id": b2, "count": 1})

    pixels_df = pd.DataFrame(pixel_rows)
    pixels_df = pixels_df.groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()

    cooler.create_cooler(
        str(path), bins=bins_df, pixels=pixels_df, ordered=True, symmetric_upper=True
    )


# ---------------------------------------------------------------------------
# Layer 1: Unit tests for internal helpers
# ---------------------------------------------------------------------------


class TestBinPairsChunk:
    """Tests for _bin_pairs_chunk: the innermost binning function."""

    def _make_chunk(self, rows: list[tuple[str, int, str, int]]) -> pd.DataFrame:
        return pd.DataFrame(rows, columns=["chrom1", "pos1", "chrom2", "pos2"])

    def test_fixed_bin_basic_count(self):
        """Each pair produces exactly one count unit."""
        # Contig A: length 1000, binsize 100 -> 10 bins (ids 0..9)
        # Contig B: length 600,  binsize 100 ->  6 bins (ids 10..15)
        contig_to_id = {"A": 0, "B": 1}
        chrom_offset = np.array([0, 10, 16], dtype=np.int64)
        binsize = 100

        # 3 contacts, all in different (bin1, bin2) pairs -> 3 pixels, total count=3
        chunk = self._make_chunk(
            [
                ("A", 50, "A", 150),  # bin 0 vs bin 1 (intra A)
                ("A", 250, "B", 50),  # bin 2 vs bin 10 (inter A-B)
                ("B", 100, "B", 300),  # bin 11 vs bin 13 (intra B)
            ]
        )
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert result["count"].sum() == 3
        assert len(result) == 3

    def test_fixed_bin_dedup_within_chunk(self):
        """Multiple pairs mapping to the same bin pair are summed."""
        contig_to_id = {"A": 0}
        chrom_offset = np.array([0, 10], dtype=np.int64)
        binsize = 100

        # 5 contacts all mapping to the same (bin 0, bin 1) pair
        chunk = self._make_chunk(
            [
                ("A", 10, "A", 110),
                ("A", 20, "A", 120),
                ("A", 30, "A", 130),
                ("A", 40, "A", 140),
                ("A", 50, "A", 150),
            ]
        )
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert len(result) == 1  # all collapsed to one pixel
        assert result["count"].sum() == 5

    def test_fixed_bin_upper_triangle_enforced(self):
        """bin1_id <= bin2_id is always enforced regardless of input order."""
        contig_to_id = {"A": 0, "B": 1}
        chrom_offset = np.array([0, 10, 16], dtype=np.int64)
        binsize = 100

        # Inter pair stored as (B, A) — should be swapped to (A, B)
        chunk = self._make_chunk([("B", 50, "A", 50)])
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert (result["bin1_id"] <= result["bin2_id"]).all()

    def test_fixed_bin_self_pair_stays(self):
        """Pairs where both reads fall in the same bin -> single diagonal pixel."""
        contig_to_id = {"A": 0}
        chrom_offset = np.array([0, 10], dtype=np.int64)
        binsize = 100
        chunk = self._make_chunk([("A", 50, "A", 75)])  # both in bin 0
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert result.iloc[0]["bin1_id"] == result.iloc[0]["bin2_id"] == 0

    def test_unknown_contigs_dropped(self):
        """Pairs referencing contigs not in contig_to_id are silently dropped."""
        contig_to_id = {"A": 0}
        chrom_offset = np.array([0, 10], dtype=np.int64)
        binsize = 100

        chunk = self._make_chunk(
            [
                ("A", 50, "A", 150),  # valid
                ("UNKNOWN", 50, "A", 50),  # invalid chrom1 -> dropped
                ("A", 50, "UNKNOWN", 50),  # invalid chrom2 -> dropped
            ]
        )
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert result["count"].sum() == 1  # only the valid pair

    def test_all_unknown_returns_none(self):
        """If all pairs are on unknown contigs, returns None (not empty DF)."""
        contig_to_id = {"A": 0}
        chrom_offset = np.array([0, 10], dtype=np.int64)
        binsize = 100
        chunk = self._make_chunk([("UNKNOWN", 50, "UNKNOWN", 150)])
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is None

    # --- Bug regression: fixed-size overflow ---
    def test_fixed_bin_position_at_contig_boundary(self):
        """A read at exactly the contig length (off-by-one edge case) stays in last bin."""
        # Contig A: length=100, binsize=100 -> 1 bin (id 0)
        contig_to_id = {"A": 0}
        chrom_offset = np.array([0, 1], dtype=np.int64)
        binsize = 100

        # pos=99 (0-based after conversion from 1-based 100) -> bin 0 (correct)
        chunk = self._make_chunk(
            [("A", 100, "A", 100)]
        )  # 1-based pos=100 -> 0-based 99
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        assert result.iloc[0]["bin1_id"] == 0
        assert result.iloc[0]["bin2_id"] == 0

    def test_fixed_bin_overflow_bug(self):
        """
        REGRESSION: A 1-based position beyond the contig length must be clamped
        to the contig's last bin, not overflow into the next contig's bin IDs.

        Contig A has length=200, binsize=100, so 2 bins (ids 0, 1).
        chrom_offset = [0, 2, ...].
        A read at 1-based pos=201 -> 0-based 200 -> 0 + 200//100 = 2 unclamped,
        which is contig B's first bin. It must be clamped to bin 1 (A's last).
        """
        contig_to_id = {"A": 0, "B": 1}
        chrom_offset = np.array([0, 2, 4], dtype=np.int64)  # A: bins 0-1, B: bins 2-3
        binsize = 100

        chunk = self._make_chunk([("A", 1, "A", 201)])
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, binsize, None)
        assert result is not None
        # Must be clamped to bin 1 (last bin of A), NOT overflow into B (bin 2).
        assert result.iloc[0]["bin2_id"] == 1

    # --- Bug regression: fragment-level underflow ---
    def test_fragment_bin_position_before_first_fragment(self):
        """
        BUG: With restriction-enzyme fragment bins that do not start at 0,
        a read before the first fragment underflows into the previous contig's bins.

        Contig A fragments: [50, 200], contig B fragments: [0, 150]
        A read at position 30 (0-based) on contig A -> underflows to last bin of... nothing,
        or the last bin of the previous contig if A is not first.

        This test uses a two-contig setup where A is second (global bins 2, 3)
        and position 30 on A underflows to global bin 1 (B's last bin).
        """
        # Simulate: B=[0,150], A=[50,200] as fragment bins
        # Global bin_starts: [0, 150, 50, 200] (B then A in contig order)
        # chrom_offset: B starts at 0, A starts at 2
        contig_to_id = {"B": 0, "A": 1}
        # B: bins at [0, 150], A: bins at [50, 200]
        bin_starts = np.array([0, 150, 50, 200], dtype=np.int64)
        chrom_offset = np.array([0, 2, 4], dtype=np.int64)

        # Read on A at pos 30 (0-based), which is before the first fragment (start=50)
        chunk = self._make_chunk([("A", 31, "A", 31)])  # 1-based -> 0-based = 30
        result = _bin_pairs_chunk(chunk, contig_to_id, chrom_offset, None, bin_starts)
        assert result is not None
        b1 = result.iloc[0]["bin1_id"]
        # Must be clamped to A's first bin (global bin 2), NOT underflow into
        # contig B's last bin (global bin 1).
        assert b1 >= 2  # must be in A's bin range [2, 3]


# ---------------------------------------------------------------------------
# Layer 2: Round-trip invariant tests (from_pairs)
# ---------------------------------------------------------------------------


class TestFromPairsInvariants:
    """
    Formal invariants that must hold for every valid from_pairs output.

    These do NOT depend on specific counts — they test structural correctness
    regardless of input size.
    """

    @pytest.fixture
    def simple_setup(self, tmp_path):
        """
        Deterministic 3-contig, 2-genome setup with known contacts.

        Contacts (1-based positions, fixed binsize=100):
            A(100) - A(200) : intra-A, bins (1,2)          count 2
            A(50)  - B(50)  : inter A-B, bins (0, 10)      count 3
            B(150) - C(50)  : inter B-C, bins (11, 20)     count 1
            C(100) - C(200) : intra-C, bins (20, 21)       count 1
        Total pairs = 2+3+1+1 = 7
        """
        chromsizes = {"A": 500, "B": 400, "C": 300}
        pairs_rows = (
            [("A", 100, "A", 200)] * 2
            + [("A", 50, "B", 50)] * 3
            + [("B", 150, "C", 50)]
            + [("C", 100, "C", 200)]
        )
        contigs = [("A", "MAG1"), ("B", "MAG1"), ("C", "MAG2")]
        n_pairs = len(pairs_rows)

        pairs_path = tmp_path / "test.pairs"
        _write_pairs(pairs_path, pairs_rows, chromsizes)
        contigs_path = tmp_path / "contigs.tsv"
        _write_contigs_table(contigs_path, contigs)
        out_path = tmp_path / "out.metacool"

        from_pairs(
            uri=out_path,
            pairs=pairs_path,
            bins=None,
            binsize=100,
            contig_genome_map=contigs_path,
            contig_table_format="contigs",
            assembly_name="test",
            contig_ordering="genome-then-alpha",
            chunksize=2,  # tiny chunksize to exercise cross-chunk dedup
        )
        return out_path, n_pairs

    def test_output_is_valid_metacooler(self, simple_setup):
        out_path, _ = simple_setup
        assert is_metacooler(out_path)
        validate(out_path)  # must not raise

    def test_pixel_sum_equals_pair_count(self, simple_setup):
        """
        Invariant: sum(pixels["count"]) == number of input pairs.

        This is the most fundamental correctness check: every pair must
        be counted exactly once, no pair dropped, no pair double-counted.
        """
        out_path, n_pairs = simple_setup
        with h5py.File(out_path, "r") as f:
            total = int(f["pixels"]["count"][:].sum())
        assert total == n_pairs, (
            f"Expected {n_pairs} total contacts, got {total}. "
            "Either pairs are being dropped or double-counted."
        )

    def test_pixels_upper_triangle(self, simple_setup):
        """Invariant: all pixels satisfy bin1_id <= bin2_id."""
        out_path, _ = simple_setup
        with h5py.File(out_path, "r") as f:
            b1 = f["pixels"]["bin1_id"][:]
            b2 = f["pixels"]["bin2_id"][:]
        assert (b1 <= b2).all(), "Upper-triangle violation: found bin1_id > bin2_id"

    def test_pixels_sorted_by_bin1_bin2(self, simple_setup):
        """Invariant: pixels are sorted by (bin1_id, bin2_id) for index validity."""
        out_path, _ = simple_setup
        with h5py.File(out_path, "r") as f:
            b1 = f["pixels"]["bin1_id"][:]
            b2 = f["pixels"]["bin2_id"][:]
        pairs_arr = np.stack([b1, b2], axis=1)
        sorted_pairs = pairs_arr[np.lexsort((b2, b1))]
        np.testing.assert_array_equal(
            pairs_arr,
            sorted_pairs,
            err_msg="Pixels are not sorted by (bin1_id, bin2_id)",
        )

    def test_pixels_no_duplicate_bin_pairs(self, simple_setup):
        """Invariant: no two pixels share the same (bin1_id, bin2_id)."""
        out_path, _ = simple_setup
        with h5py.File(out_path, "r") as f:
            b1 = f["pixels"]["bin1_id"][:]
            b2 = f["pixels"]["bin2_id"][:]
        pairs_arr = list(zip(b1.tolist(), b2.tolist(), strict=False))
        assert len(pairs_arr) == len(set(pairs_arr)), (
            "Duplicate (bin1_id, bin2_id) found"
        )

    def test_chrom_offset_consistent_with_bins(self, simple_setup):
        """Invariant: indexes/chrom_offset matches the actual bin/chrom distribution."""
        out_path, _ = simple_setup
        with h5py.File(out_path, "r") as f:
            chrom_offset = f["indexes"]["chrom_offset"][:]
            chrom_codes = f["bins"]["chrom"][:]
        n_contigs = len(chrom_offset) - 1
        expected_offset = np.zeros(n_contigs + 1, dtype=np.int64)
        counts = np.bincount(chrom_codes.astype(np.int64), minlength=n_contigs)
        expected_offset[1:] = np.cumsum(counts)
        np.testing.assert_array_equal(
            chrom_offset,
            expected_offset,
            err_msg="chrom_offset does not match actual bin distribution",
        )

    def test_bin1_offset_consistent(self, simple_setup):
        """Invariant: indexes/bin1_offset matches the actual pixel/bin1_id distribution."""
        out_path, _ = simple_setup
        with h5py.File(out_path, "r") as f:
            bin1_offset = f["indexes"]["bin1_offset"][:]
            b1 = f["pixels"]["bin1_id"][:]
            n_bins = len(f["bins"]["chrom"])
        expected = np.zeros(n_bins + 1, dtype=np.int64)
        counts = np.bincount(b1.astype(np.int64), minlength=n_bins)
        expected[1:] = np.cumsum(counts)
        np.testing.assert_array_equal(
            bin1_offset,
            expected,
            err_msg="bin1_offset does not match actual pixel distribution",
        )

    def test_genomes_table_correct(self, simple_setup):
        """Invariant: genomes table entries are consistent with contigs."""
        out_path, _ = simple_setup
        from metacooler.api import MetaCooler

        mc = MetaCooler(out_path)
        genomes = mc.genomes
        contigs = mc.contigs

        for _, genome_row in genomes.iterrows():
            g_name = genome_row["name"]
            g_contigs = contigs[contigs["genome_id"] == genome_row.name]
            assert len(g_contigs) == genome_row["n_contigs"], (
                f"Genome {g_name}: n_contigs={genome_row['n_contigs']} "
                f"but found {len(g_contigs)} rows in contigs table"
            )
            assert int(g_contigs["length"].sum()) == genome_row["length"], (
                f"Genome {g_name}: length={genome_row['length']} "
                f"but sum of contig lengths={g_contigs['length'].sum()}"
            )

    def test_cross_chunk_dedup(self, tmp_path):
        """
        The same (bin1, bin2) pair spanning multiple chunks must be merged,
        not double-counted. Uses chunksize=1 to force every pair into its own chunk.
        """
        chromsizes = {"A": 1000}
        # 5 pairs all mapping to (bin 0, bin 1) in 100-bp bins
        pairs_rows = [("A", 50, "A", 150)] * 5
        pairs_path = tmp_path / "pairs.pairs"
        _write_pairs(pairs_path, pairs_rows, chromsizes)
        contigs_path = tmp_path / "contigs.tsv"
        _write_contigs_table(contigs_path, [("A", "G1")])
        out = tmp_path / "out.metacool"

        from_pairs(
            uri=out,
            pairs=pairs_path,
            binsize=100,
            contig_genome_map=contigs_path,
            contig_table_format="contigs",
            chunksize=1,  # one pair per chunk -> cross-chunk dedup required
        )
        with h5py.File(out, "r") as f:
            total = int(f["pixels"]["count"][:].sum())
            n_pixels = len(f["pixels"]["count"])
        assert total == 5, f"Expected count=5, got {total}"
        assert n_pixels == 1, (
            f"Expected 1 pixel, got {n_pixels} (cross-chunk dedup failed)"
        )

    def test_fragment_bins_count_invariant(self, tmp_path):
        """
        Fragment-level (BED) bins: total count must equal input pair count.

        Fragments starting at position 0 on all contigs, so no underflow possible.
        """
        chromsizes = {"A": 500, "B": 300}
        # BED bins: irregular sizes, all starting at 0
        bed_lines = [
            "A\t0\t200\n",
            "A\t200\t500\n",
            "B\t0\t150\n",
            "B\t150\t300\n",
        ]
        bed_path = tmp_path / "frags.bed"
        bed_path.write_text("".join(bed_lines))

        pairs_rows = [
            ("A", 100, "A", 300),  # frag 0 vs frag 1 of A
            ("A", 50, "B", 100),  # frag 0 of A vs frag 0 of B
            ("B", 200, "B", 250),  # frag 1 vs frag 1 of B (diagonal)
        ]
        n_pairs = len(pairs_rows)
        pairs_path = tmp_path / "frags.pairs"
        _write_pairs(pairs_path, pairs_rows, chromsizes)
        contigs_path = tmp_path / "contigs.tsv"
        _write_contigs_table(contigs_path, [("A", "G1"), ("B", "G2")])
        out = tmp_path / "out.metacool"

        from_pairs(
            uri=out,
            pairs=pairs_path,
            bins=bed_path,
            contig_genome_map=contigs_path,
            contig_table_format="contigs",
        )
        with h5py.File(out, "r") as f:
            total = int(f["pixels"]["count"][:].sum())
        assert total == n_pairs

    def test_intra_vs_inter_counts(self, tmp_path):
        """
        Intra-genome and inter-genome contacts must both be preserved exactly.

        Known setup: 2 genomes (G1=A, G2=B), 3 intra-A + 2 inter + 4 intra-B.
        After binning, the per-(genome pair) sums must match.
        """
        chromsizes = {"A": 1000, "B": 1000}
        pairs_rows = (
            [("A", 50, "A", 550)] * 3  # intra-G1
            + [("A", 50, "B", 50)] * 2  # inter G1-G2
            + [("B", 50, "B", 550)] * 4  # intra-G2
        )
        pairs_path = tmp_path / "p.pairs"
        _write_pairs(pairs_path, pairs_rows, chromsizes)
        contigs_path = tmp_path / "c.tsv"
        _write_contigs_table(contigs_path, [("A", "G1"), ("B", "G2")])
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            pairs=pairs_path,
            binsize=500,
            contig_genome_map=contigs_path,
            contig_table_format="contigs",
        )

        from metacooler.api import MetaCooler

        mc = MetaCooler(out)
        bins = mc.bins()
        contigs = mc.contigs
        pixels = mc.pixels()

        bin_to_contig = bins["chrom"].cat.codes.values
        contig_to_genome = contigs["genome_id"].values
        bin_to_genome = contig_to_genome[bin_to_contig]

        g1 = pixels["bin1_id"].map(lambda i: bin_to_genome[i])
        g2 = pixels["bin2_id"].map(lambda i: bin_to_genome[i])
        total = pixels["count"].sum()
        intra_g1 = pixels.loc[(g1 == 0) & (g2 == 0), "count"].sum()
        intra_g2 = pixels.loc[(g1 == 1) & (g2 == 1), "count"].sum()
        inter = pixels.loc[g1 != g2, "count"].sum()

        assert total == 9
        assert intra_g1 == 3
        assert intra_g2 == 4
        assert inter == 2


# ---------------------------------------------------------------------------
# Layer 3: Round-trip invariant tests (from_cool)
# ---------------------------------------------------------------------------


class TestFromCoolInvariants:
    """
    Invariants for from_cool: pixel counts must be perfectly preserved
    after the bin-remapping step.
    """

    @pytest.fixture
    def simple_cool(self, tmp_path):
        """
        A cool file with known contact counts.
        Chromosomes: chrA (1000 bp), chrB (800 bp), binsize=200.
        Contacts: 5 intra-A, 3 inter A-B, 2 intra-B.
        """
        chromsizes = {"chrA": 1000, "chrB": 800}
        pairs = (
            [("chrA", 100, "chrA", 600)] * 5
            + [("chrA", 100, "chrB", 100)] * 3
            + [("chrB", 100, "chrB", 600)] * 2
        )
        cool_path = tmp_path / "source.cool"
        _write_cool(cool_path, chromsizes, binsize=200, pairs=pairs)
        contigs_path = tmp_path / "contigs.tsv"
        _write_contigs_table(contigs_path, [("chrA", "G1"), ("chrB", "G2")])
        return cool_path, contigs_path, 10  # 10 total contacts

    def test_pixel_sum_preserved(self, simple_cool, tmp_path):
        """
        Invariant: sum(pixels["count"]) must equal the source cooler's sum.
        This verifies no contacts are lost or duplicated during remapping.
        """
        cool_path, contigs_path, n_contacts = simple_cool
        out = tmp_path / "out.metacool"
        from_cool(cool_uri=cool_path, output=out, contig_genome_map=contigs_path)

        src = cooler.Cooler(str(cool_path))
        src_sum = src.pixels()[:]["count"].sum()

        with h5py.File(out, "r") as f:
            mc_sum = int(f["pixels"]["count"][:].sum())
        assert mc_sum == src_sum == n_contacts

    def test_nnz_preserved(self, simple_cool, tmp_path):
        """Invariant: number of non-zero pixels must equal the source cooler's nnz."""
        cool_path, contigs_path, _ = simple_cool
        out = tmp_path / "out.metacool"
        from_cool(cool_uri=cool_path, output=out, contig_genome_map=contigs_path)

        src = cooler.Cooler(str(cool_path))
        with h5py.File(out, "r") as f:
            mc_nnz = len(f["pixels"]["count"])
        assert mc_nnz == src.info["nnz"]

    def test_per_chrom_pair_counts_preserved(self, simple_cool, tmp_path):
        """
        Invariant: contacts per (chrom1, chrom2) pair must be identical
        in the source cooler and the output metacooler (after flattening back).

        This is the strongest form of count preservation: not just totals,
        but per-contig-pair subtotals.
        """
        from metacooler.io.flatten import flatten

        cool_path, contigs_path, _ = simple_cool
        out = tmp_path / "out.metacool"
        flat = tmp_path / "flat.cool"
        from_cool(
            cool_uri=cool_path,
            output=out,
            contig_genome_map=contigs_path,
            contig_ordering="cool",
        )  # preserve order for easy comparison
        flatten(out, flat)

        def _pair_counts(clr: cooler.Cooler) -> dict:
            bins = clr.bins()[:][["chrom"]]
            pixels = clr.pixels()[:]
            c1 = bins["chrom"].iloc[pixels["bin1_id"].values].values
            c2 = bins["chrom"].iloc[pixels["bin2_id"].values].values
            df = pd.DataFrame({"c1": c1, "c2": c2, "count": pixels["count"].values})
            # Normalize order so (A,B)==(B,A)
            df["key"] = [
                (min(a, b), max(a, b)) for a, b in zip(df["c1"], df["c2"], strict=False)
            ]
            return df.groupby("key")["count"].sum().to_dict()

        src_counts = _pair_counts(cooler.Cooler(str(cool_path)))
        mc_counts = _pair_counts(cooler.Cooler(str(flat)))
        assert src_counts == mc_counts, (
            f"Per-(chrom,chrom) counts differ after round-trip.\n"
            f"Source: {src_counts}\nMeta: {mc_counts}"
        )

    def test_reordering_preserves_counts(self, simple_cool, tmp_path):
        """
        When contig_ordering changes the order of contigs relative to the source
        cooler, pixel sums must still be preserved (remapping must be correct).
        """
        cool_path, contigs_path, n_contacts = simple_cool
        for ordering in ("genome-then-alpha", "genome-then-length", "alpha"):
            out = tmp_path / f"out_{ordering}.metacool"
            from_cool(
                cool_uri=cool_path,
                output=out,
                contig_genome_map=contigs_path,
                contig_ordering=ordering,
            )
            with h5py.File(out, "r") as f:
                mc_sum = int(f["pixels"]["count"][:].sum())
            assert mc_sum == n_contacts, (
                f"ordering={ordering!r}: expected {n_contacts} contacts, got {mc_sum}"
            )

    def test_upper_triangle_preserved_after_reorder(self, simple_cool, tmp_path):
        """After remapping, all pixels must still satisfy bin1_id <= bin2_id."""
        cool_path, contigs_path, _ = simple_cool
        out = tmp_path / "out.metacool"
        from_cool(
            cool_uri=cool_path,
            output=out,
            contig_genome_map=contigs_path,
            contig_ordering="alpha",
        )  # alpha may reverse contig order
        with h5py.File(out, "r") as f:
            b1 = f["pixels"]["bin1_id"][:]
            b2 = f["pixels"]["bin2_id"][:]
        assert (b1 <= b2).all(), (
            "Upper-triangle violated after contig reordering in from_cool"
        )

    def test_is_valid_metacooler(self, simple_cool, tmp_path):
        """Output must pass metacooler schema validation."""
        cool_path, contigs_path, _ = simple_cool
        out = tmp_path / "out.metacool"
        from_cool(cool_uri=cool_path, output=out, contig_genome_map=contigs_path)
        assert is_metacooler(out)
        validate(out)  # must not raise


# ---------------------------------------------------------------------------
# Layer 3b: Bug regression tests for from_cool
# ---------------------------------------------------------------------------


class TestFromCoolBugRegressions:
    def test_drop_unassigned_crashes(self, tmp_path):
        """
        REGRESSION: from_cool(..., drop_unassigned=True) must succeed when the
        cool file has bins for contigs that are being dropped.

        Dropped contigs' bins map to a -1 sentinel and their pixels are filtered
        out (rather than raising ValueError / IndexError).
        """
        chromsizes = {"chrA": 1000, "chrB": 800, "chrC": 600}
        pairs = [("chrA", 100, "chrA", 600)] * 3 + [("chrB", 100, "chrC", 100)] * 2
        cool_path = tmp_path / "src.cool"
        _write_cool(cool_path, chromsizes, binsize=200, pairs=pairs)

        # Only assign chrA and chrB; chrC is absent from the table -> will be dropped
        contigs_path = tmp_path / "contigs.tsv"
        _write_contigs_table(contigs_path, [("chrA", "G1"), ("chrB", "G2")])
        out = tmp_path / "out.metacool"

        from_cool(
            cool_uri=cool_path,
            output=out,
            contig_genome_map=contigs_path,
            drop_unassigned=True,
        )

        # validate the result and check that the dropped contig is gone
        assert is_metacooler(out)
        from metacooler.api import MetaCooler

        mc = MetaCooler(out)
        assert "chrC" not in mc.contigs["name"].values
        # chrA-chrA contacts should be preserved; chrB-chrC contacts should be dropped
        # since chrC was dropped -> only the chrA intra-contacts remain
        with h5py.File(out, "r") as f:
            total = int(f["pixels"]["count"][:].sum())
        assert total == 3  # only the intra-chrA contacts

    def test_nproc_parallel_count_preserved(self, tmp_path):
        """
        from_cool(..., nproc=N) must produce the same pixel sum as nproc=1.
        The parallel path uses ThreadPoolExecutor — a race condition or ordering
        bug would cause count discrepancies.
        """
        chromsizes = {"chrA": 2000, "chrB": 1500}
        # 50 contacts spread across both contigs
        rng = np.random.default_rng(42)
        pairs = []
        for _ in range(30):
            pairs.append(
                ("chrA", int(rng.integers(1, 2000)), "chrA", int(rng.integers(1, 2000)))
            )
        for _ in range(20):
            pairs.append(
                ("chrA", int(rng.integers(1, 2000)), "chrB", int(rng.integers(1, 1500)))
            )
        cool_path = tmp_path / "src.cool"
        _write_cool(cool_path, chromsizes, binsize=200, pairs=pairs)
        contigs_path = tmp_path / "c.tsv"
        _write_contigs_table(contigs_path, [("chrA", "G1"), ("chrB", "G2")])

        out1 = tmp_path / "out_nproc1.metacool"
        out4 = tmp_path / "out_nproc4.metacool"
        from_cool(
            cool_uri=cool_path, output=out1, contig_genome_map=contigs_path, nproc=1
        )
        from_cool(
            cool_uri=cool_path, output=out4, contig_genome_map=contigs_path, nproc=4
        )

        with h5py.File(out1, "r") as f:
            sum1 = int(f["pixels"]["count"][:].sum())
        with h5py.File(out4, "r") as f:
            sum4 = int(f["pixels"]["count"][:].sum())
        assert sum1 == sum4, (
            f"nproc=1 sum={sum1} != nproc=4 sum={sum4}: parallel path has a count bug"
        )


# ---------------------------------------------------------------------------
# Efficiency notes (not executable, but documented for review)
# ---------------------------------------------------------------------------
"""
EFFICIENCY ISSUES IN from_pairs AND from_cool (not bugs, but scalability limits):

1. pixel_stream_sorted() in from_pairs collects ALL pixel chunks into RAM before
   sorting (line: `frames = list(pixel_stream())`). For a metagenome with 10^8
   pixels this requires ~2.4 GB just for the pixel table (3 int64 cols × 8 bytes).

   Better approach: pass ordered=False to cooler.create_cooler and let cooler
   do its own disk-based external mergesort. cooler handles this correctly and
   only buffers `mergebuf` rows at a time (default 20M). This makes from_pairs
   O(chunksize) in memory regardless of file size.

   Proposed fix:
       cooler.create_cooler(
           str(uri), bins=bins_df,
           pixels=pixel_stream(),   # not pixel_stream_sorted()
           ordered=False,           # let cooler sort
           symmetric_upper=True,
           ...
       )

2. Same issue in from_cool: _remap_pixel_stream_sorted() collects all remapped
   chunks in RAM. The source cooler is already sorted, and remapping preserves
   sortedness IFF the new contig order is consistent (which it is when
   ordering='cool'). For other orderings, passing ordered=False is correct and
   avoids the RAM spike.

3. _stream_pairs block_size heuristic (line: `block_size = max(chunksize*120, 8MB)`):
   120 bytes/row is a rough estimate. Pairs rows with long contig names (e.g.
   ">100-char metagenomic scaffold names) will be underestimated, causing pyarrow
   to return more rows per batch than `chunksize`. This doesn't break correctness
   but makes the chunksize parameter misleading.

4. from_cool with nproc > 1 submits ALL chunk futures eagerly before consuming any:
       futures = [pool.submit(...) for lo, hi in ranges]  # all submitted at once
       for f in futures: yield f.result()
   For large files this fills the thread queue and may OOM if completed futures
   are not garbage-collected fast enough. Better: use a sliding window:
       for lo, hi in ranges:
           futures.append(pool.submit(_remap_chunk, lo, hi))
           if len(futures) >= 2 * nproc:
               yield futures.pop(0).result()
   Or use concurrent.futures.as_completed() with a bounded semaphore.
"""
