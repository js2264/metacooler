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
