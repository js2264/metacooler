"""Tests for metacooler.io.merge."""

from __future__ import annotations

import pytest

from metacooler.io.merge import merge


class TestMerge:
    def test_merge_identical(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        merge([tiny_metacooler, tiny_metacooler], out)
        assert out.exists()

    def test_merge_single_file_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        with pytest.raises(ValueError, match="At least two"):
            merge([tiny_metacooler], out)

    def test_merge_different_files_fails(self, tiny_metacooler, tmp_path):
        # Create a different metacool file
        import tempfile
        from pathlib import Path

        from metacooler.io.create import from_pairs

        with tempfile.TemporaryDirectory() as td:
            t = Path(td)
            fa = t / "tiny.fa"
            fa.write_text(">contig_X\n" + "ACGT" * 10 + "\n")
            mag = t / "mag.tsv"
            mag.write_text("contig_name\tgenome_name\ncontig_X\tMAG_X\n")
            bed = t / "bins.bed"
            bed.write_text("contig_X\t0\t40\n")
            pairs = t / "pairs.pairs"
            pairs.write_text("#c\nread_0\tcontig_X\t1\tcontig_X\t1\t+\t+\n")
            other = t / "other.metacooler"
            from_pairs(
                other,
                bins=bed,
                pairs=pairs,
                contig_genome_map=mag,
                contig_table_format="contigs",
                fasta=fa,
                chunksize=1000,
            )
            out = tmp_path / "merged.metacool"
            with pytest.raises(ValueError, match="Contigs table mismatch"):
                merge([tiny_metacooler, other], out)

    def test_merge_bins_chrom_dtype_is_not_enum(self, tiny_metacooler, tmp_path):
        """bins/chrom must be plain int32 (no HDF5 enum), not an HDF5 enum dtype."""
        import h5py

        out = tmp_path / "merged.metacool"
        merge([tiny_metacooler, tiny_metacooler], out)
        with h5py.File(out, "r") as f:
            chrom_dtype = f["bins"]["chrom"].dtype
        assert chrom_dtype == "int32", (
            f"Expected int32 but got {chrom_dtype} — enum dtype not stripped"
        )

    def test_merge_pixel_sums_double(self, tiny_metacooler, tmp_path):
        """Merging identical files should double all pixel counts."""
        import h5py

        out = tmp_path / "merged.metacool"
        merge([tiny_metacooler, tiny_metacooler], out)
        with h5py.File(tiny_metacooler, "r") as f:
            orig_sum = int(f["pixels"]["count"][:].sum())
        with h5py.File(out, "r") as f:
            merged_sum = int(f["pixels"]["count"][:].sum())
        assert merged_sum == 2 * orig_sum, (
            f"Expected {2 * orig_sum} total counts but got {merged_sum}"
        )
