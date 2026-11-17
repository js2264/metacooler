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


class TestMetadataColumns:
    """Only identity columns decide whether two maps merge (annotate's knock-on)."""

    def _copy(self, src, dst):
        dst.write_bytes(src.read_bytes())
        return dst

    def test_an_annotated_replicate_still_merges(
        self, tiny_metacooler, tmp_path, caplog
    ):
        from metacooler.api import MetaCooler
        from metacooler.io.annotate import annotate

        rep = self._copy(tiny_metacooler, tmp_path / "rep.metacool")
        table = tmp_path / "g.tsv"
        table.write_text("genome_name\tcompleteness\nMAG_1\t0.9\nMAG_2\t0.7\n")
        annotate(rep, inplace=True, genomes_table=table)
        out = tmp_path / "merged.metacool"
        with caplog.at_level("WARNING", logger="metacooler"):
            merge([tiny_metacooler, rep], out)
        assert "genomes.completeness (only here)" in caplog.text
        # The first input's metadata is the one kept.
        assert "completeness" not in MetaCooler(out).genomes.columns

        out2 = tmp_path / "merged2.metacool"
        caplog.clear()
        with caplog.at_level("WARNING", logger="metacooler"):
            merge([rep, tiny_metacooler], out2)
        assert "genomes.completeness (absent here)" in caplog.text
        assert "completeness" in MetaCooler(out2).genomes.columns

    def test_differing_values_are_reported(self, tiny_metacooler, tmp_path, caplog):
        from metacooler.io.annotate import annotate

        a = self._copy(tiny_metacooler, tmp_path / "a.metacool")
        b = self._copy(tiny_metacooler, tmp_path / "b.metacool")
        for path, value in ((a, "0.9"), (b, "0.8")):
            table = tmp_path / f"{path.stem}.tsv"
            table.write_text(f"genome_name\tcompleteness\nMAG_1\t{value}\n")
            annotate(path, inplace=True, genomes_table=table)
        with caplog.at_level("WARNING", logger="metacooler"):
            merge([a, b], tmp_path / "m.metacool")
        assert "genomes.completeness (other values)" in caplog.text

    def test_a_user_qc_column_of_contigs_is_compared(
        self, tiny_metacooler, tmp_path, caplog
    ):
        """Only the two QC columns of /contigs are skipped; qc_pass is metadata."""
        from metacooler.io.annotate import annotate

        a = self._copy(tiny_metacooler, tmp_path / "a.metacool")
        b = self._copy(tiny_metacooler, tmp_path / "b.metacool")
        for path, value in ((a, "1"), (b, "0")):
            table = tmp_path / f"{path.stem}.tsv"
            table.write_text(f"contig_name\tqc_pass\ncontig_A\t{value}\n")
            annotate(path, inplace=True, contigs_table=table)
        with caplog.at_level("WARNING", logger="metacooler"):
            merge([a, b], tmp_path / "m.metacool")
        assert "contigs.qc_pass (other values)" in caplog.text
        assert "qc_intra_contig_short" not in caplog.text

    def test_identity_hash_ignores_integer_width(self, tiny_metacooler, tmp_path):
        import h5py
        import numpy as np

        from metacooler.io.merge import _hash_group

        wide = self._copy(tiny_metacooler, tmp_path / "wide.metacool")
        with h5py.File(wide, "a") as f:
            lengths = f["contigs/length"][:]
            del f["contigs/length"]
            f["contigs"].create_dataset("length", data=lengths.astype(np.int64))
        with h5py.File(tiny_metacooler, "r") as f0, h5py.File(wide, "r") as f1:
            assert f0["contigs/length"].dtype != f1["contigs/length"].dtype
            assert _hash_group(f0, "contigs") == _hash_group(f1, "contigs")

    def test_a_different_assignment_is_still_refused(self, tiny_metacooler, tmp_path):
        from metacooler.io.reassign import reassign

        moved = tmp_path / "moved.metacool"
        reassign(
            tiny_metacooler,
            {"contig_A": "MAG_2"},
            moved,
            layout="keep",
            qc=False,
        )
        with pytest.raises(ValueError, match="Contigs table mismatch"):
            merge([tiny_metacooler, moved], tmp_path / "m.metacool")
