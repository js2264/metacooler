"""Tests for metacooler.io.extract."""

from __future__ import annotations

import cooler
import pytest

from metacooler.io.extract import extract


class TestExtract:
    def test_extract_single_genome(self, tiny_metacooler, tmp_path):
        out = tmp_path / "mag1.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=out)
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A", "contig_B"}

    def test_extract_multiple_genomes(self, tiny_metacooler, tmp_path):
        out = tmp_path / "multi.cool"
        extract(tiny_metacooler, genomes=["MAG_1", "MAG_2"], output=out)
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A", "contig_B", "contig_C"}

    def test_extract_unbinned(self, tiny_metacooler, tmp_path):
        out = tmp_path / "unbinned.cool"
        extract(tiny_metacooler, genomes=["__unbinned__"], output=out)
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_D"}

    def test_extract_contig(self, tiny_metacooler, tmp_path):
        out = tmp_path / "contig_a.cool"
        extract(tiny_metacooler, contigs=["contig_A"], output=out)
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A"}

    def test_extract_genomes_and_contigs_combined(self, tiny_metacooler, tmp_path):
        out = tmp_path / "combined.cool"
        extract(
            tiny_metacooler,
            genomes=["MAG_1"],
            contigs=["contig_C"],
            output=out,
        )
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A", "contig_B", "contig_C"}

    def test_extract_split(self, tiny_metacooler, tmp_path):
        out = tmp_path / "split"
        extract(
            tiny_metacooler,
            genomes=["MAG_1", "MAG_2"],
            output=out,
            split=True,
        )
        mag1 = tmp_path / "split_MAG_1.cool"
        mag2 = tmp_path / "split_MAG_2.cool"
        assert mag1.exists()
        assert mag2.exists()
        c1 = cooler.Cooler(str(mag1))
        c2 = cooler.Cooler(str(mag2))
        assert set(c1.chromnames) == {"contig_A", "contig_B"}
        assert set(c2.chromnames) == {"contig_C"}

    def test_extract_split_with_contigs_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "split"
        # --split is genome-only: combining it with --contigs is an error.
        with pytest.raises(ValueError, match="cannot be combined with --split"):
            extract(
                tiny_metacooler,
                genomes=["MAG_1"],
                contigs=["contig_C"],
                output=out,
                split=True,
            )

    def test_extract_split_contigs_only_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "split"
        # --split with only --contigs is rejected (split is genome-only).
        with pytest.raises(ValueError, match="cannot be combined with --split"):
            extract(tiny_metacooler, contigs=["contig_C"], output=out, split=True)

    def test_extract_unknown_genome_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "nope.cool"
        with pytest.raises(KeyError, match="NONEXISTENT"):
            extract(tiny_metacooler, genomes=["NONEXISTENT"], output=out)

    def test_extract_unknown_contig_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "nope.cool"
        with pytest.raises(KeyError, match="NONEXISTENT"):
            extract(tiny_metacooler, contigs=["NONEXISTENT"], output=out)

    def test_extract_no_genomes_or_contigs_raises(self, tiny_metacooler, tmp_path):
        out = tmp_path / "nope.cool"
        with pytest.raises(ValueError, match="At least one of"):
            extract(tiny_metacooler, output=out)

    def test_extract_is_valid_cooler(self, tiny_metacooler, tmp_path):
        out = tmp_path / "mag1.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=out)
        assert cooler.fileops.is_cooler(str(out))

    def test_extracted_bins_within_selection(self, tiny_metacooler, tmp_path):
        out = tmp_path / "mag1.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=out)
        c = cooler.Cooler(str(out))
        bins = c.bins()[:]
        assert set(bins["chrom"].unique()) == {"contig_A", "contig_B"}

    def test_extracted_pixels_upper_triangle(self, tiny_metacooler, tmp_path):
        out = tmp_path / "mag1.cool"
        extract(tiny_metacooler, genomes=["MAG_1"], output=out)
        c = cooler.Cooler(str(out))
        pixels = c.pixels()[:]
        assert (pixels["bin1_id"] <= pixels["bin2_id"]).all()

    def test_extract_all_genomes(self, tiny_metacooler, tmp_path):
        out = tmp_path / "all.cool"
        extract(
            tiny_metacooler,
            genomes=["MAG_1", "MAG_2", "__unbinned__"],
            output=out,
        )
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A", "contig_B", "contig_C", "contig_D"}

    def test_extract_contig_with_no_pixels(
        self, tmp_path, tiny_fasta, tiny_contigs_table, tiny_bins_bed
    ):
        """Extracting a contig with zero pixels should produce a valid empty cool file."""
        pairs_path = tmp_path / "only_A.pairs"
        with open(pairs_path, "w") as f:
            f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
            for i in range(10):
                f.write(f"read_{i}\tcontig_A\t{10 + i}\tcontig_A\t{20 + i}\t+\t+\n")

        metacool_path = tmp_path / "sparse.metacool"
        from metacooler.io.create import from_pairs

        from_pairs(
            uri=metacool_path,
            bins=tiny_bins_bed,
            pairs=pairs_path,
            contig_genome_map=tiny_contigs_table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
            assembly_name="sparse_test",
            chunksize=1000,
        )

        out = tmp_path / "contig_b.cool"
        extract(metacool_path, contigs=["contig_B"], output=out)
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert set(c.chromnames) == {"contig_B"}
        assert len(c.pixels()[:]) == 0
