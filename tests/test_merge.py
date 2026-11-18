"""Tests for metacooler.io.merge."""

from __future__ import annotations

import shutil

import pytest

from metacooler.io.merge import merge


def _replica(path, tmp_path):
    """A byte-identical second replicate: merge refuses the *same* file twice."""
    copy = tmp_path / f"replica_{path.name}"
    shutil.copy(path, copy)
    return copy


def _cooler_info_decodes_bytes() -> bool:
    """Whether cooler's ``info()`` decodes an attribute stored as bytes."""
    import io

    import h5py
    import numpy as np
    from cooler.api import info

    with h5py.File(io.BytesIO(), "w") as f:
        f.attrs["storage-mode"] = np.bytes_(b"square")
        return info(f)["storage-mode"] == "square"


class TestMerge:
    def test_merge_identical(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        merge([tiny_metacooler, _replica(tiny_metacooler, tmp_path)], out)
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
        merge([tiny_metacooler, _replica(tiny_metacooler, tmp_path)], out)
        with h5py.File(out, "r") as f:
            chrom_dtype = f["bins"]["chrom"].dtype
        assert chrom_dtype == "int32", (
            f"Expected int32 but got {chrom_dtype} — enum dtype not stripped"
        )

    def test_merge_pixel_sums_double(self, tiny_metacooler, tmp_path):
        """Merging identical files should double all pixel counts."""
        import h5py

        out = tmp_path / "merged.metacool"
        merge([tiny_metacooler, _replica(tiny_metacooler, tmp_path)], out)
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


class TestMergeWithoutPixels:
    """Maps that hold no pixel merge into a map that holds none.

    ``cooler.merge_coolers`` concatenates the inputs' pixel chunks, and it
    failed with "No objects to concatenate" when no input had any: libraries
    with no contact left, or the per-library maps of an empty run.
    """

    @staticmethod
    def _map(tmp_path, name, pairs_body=""):
        from metacooler.io.create import from_pairs

        directory = tmp_path / name
        directory.mkdir()
        bed = directory / "bins.bed"
        bed.write_text("c1\t0\t1000\nc2\t0\t500\nc3\t0\t300\n")
        table = directory / "contigs.tsv"
        table.write_text("contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_1\nc3\tMAG_2\n")
        pairs = directory / "lib.pairs"
        pairs.write_text(
            "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n"
            + pairs_body
        )
        path = tmp_path / f"{name}.metacool"
        from_pairs(path, bins=bed, pairs=pairs, contig_genome_map=table, ps=False)
        return path

    def test_maps_without_pixels_merge(self, tmp_path):
        import cooler

        from metacooler.api import MetaCooler
        from metacooler.schema import read_table, validate

        a = self._map(tmp_path, "a")
        b = self._map(tmp_path, "b")
        out = tmp_path / "merged.metacool"
        merge([a, b], out)

        validate(out, full=True)
        clr = cooler.Cooler(str(out))
        assert clr.info["nnz"] == 0
        assert clr.info["sum"] == 0
        assert clr.pixels()[:].empty
        source = cooler.Cooler(str(a))
        assert clr.bins()[:].equals(source.bins()[:])
        for table in ("contigs", "genomes"):
            assert read_table(out, table).equals(read_table(a, table))
        assert MetaCooler(out).metadata["qc"]["n_contacts"] == 0

    @staticmethod
    def _square(path):
        """*path*, declared square, as a map that stores both triangles is."""
        import h5py

        with h5py.File(path, "r+") as f:
            f.attrs["storage-mode"] = "square"
        return path

    @pytest.mark.parametrize("mode", ["symmetric-upper", "square"])
    def test_a_merge_without_pixels_is_laid_out_as_any_merge(self, tmp_path, mode):
        """The same groups, columns, dtypes and attribute values as a merge of
        maps that have pixels, over the same contigs, but for what counts the
        pixels and when the file was made."""
        import h5py

        def made(name, body=""):
            path = self._map(tmp_path, name, body)
            return self._square(path) if mode == "square" else path

        empty = tmp_path / "empty.metacool"
        merge([made("a"), made("b")], empty)
        pair = "r1\tc1\t10\tc2\t20\nr2\tc3\t5\tc3\t250\n"
        full = tmp_path / "full.metacool"
        merge([made("c", pair), made("d", pair)], full)

        def layout(path):
            with h5py.File(path, "r") as f:
                datasets = {}
                f.visititems(
                    lambda name, node: (
                        datasets.__setitem__(name, node.dtype)
                        if isinstance(node, h5py.Dataset)
                        else None
                    )
                )
                attrs = {
                    key: f.attrs[key]
                    for key in f.attrs
                    if key not in ("creation-date", "nnz", "sum", "metadata")
                }
                return datasets, attrs

        datasets, attrs = layout(empty)
        assert attrs["storage-mode"] == mode
        assert (datasets, attrs) == layout(full)

    @pytest.mark.parametrize("with_pixels", [False, True])
    def test_maps_of_two_storage_modes_are_refused(self, tmp_path, with_pixels):
        """With or without pixels, before anything is written: cooler refused
        them only with pixels, and the empty merge took the first one's mode."""
        pair = "r1\tc1\t10\tc2\t20\n" if with_pixels else ""
        a = self._map(tmp_path, "a", pair)
        b = self._square(self._map(tmp_path, "b", pair))
        out = tmp_path / "merged.metacool"
        with pytest.raises(ValueError, match="storage-mode mismatch"):
            merge([a, b], out)
        assert not out.exists()

    @pytest.mark.parametrize("with_pixels", [False, True])
    def test_attributes_stored_as_bytes_are_read_as_text(
        self, tmp_path, with_pixels, request
    ):
        """Another writer may store ``storage-mode`` and ``contig-ordering`` as
        bytes, which ``validate`` accepts: compared raw, they refused a merge
        with a map that stores them as text, and made the empty merge square.

        With pixels, cooler's own ``merge_coolers`` compares the storage modes
        too, as its ``info()`` reads them: cooler 0.10's reads them raw, and
        refuses to merge the two maps (devel does the same)."""
        import cooler
        import h5py
        import numpy as np

        if with_pixels and not _cooler_info_decodes_bytes():
            request.applymarker(
                pytest.mark.xfail(
                    strict=True,
                    raises=ValueError,
                    reason="this cooler's merge_coolers compares the raw "
                    "storage-mode: a map storing it as bytes and one storing "
                    "it as text are not merged",
                )
            )

        pair = "r1\tc1\t10\tc2\t20\n" if with_pixels else ""
        a = self._map(tmp_path, "a", pair)
        b = self._map(tmp_path, "b", pair)
        for path in (a, b) if not with_pixels else (b,):
            with h5py.File(path, "r+") as f:
                for key in ("storage-mode", "contig-ordering"):
                    f.attrs[key] = np.bytes_(str(f.attrs[key]).encode())
        out = tmp_path / "merged.metacool"
        merge([a, b], out)
        assert cooler.Cooler(str(out)).info["storage-mode"] == "symmetric-upper"

    def test_the_cli_merges_maps_without_pixels(self, tmp_path):
        from click.testing import CliRunner

        from metacooler.api import MetaCooler
        from metacooler.cli.main import cli
        from metacooler.schema import validate

        a = self._map(tmp_path, "a")
        b = self._map(tmp_path, "b")
        out = tmp_path / "merged.metacool"
        result = CliRunner().invoke(cli, ["merge", str(a), str(b), "-o", str(out)])
        assert result.exit_code == 0, result.output
        validate(out, full=True)
        assert MetaCooler(out).metadata["commands"]
