"""Tests for metacooler.schema validation paths."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from metacooler.schema import (
    COOLER_FORMAT,
    METACOOLER_FORMAT,
    MetaCoolerFormatError,
    validate,
)


class TestValidate:
    def test_missing_file(self):
        with pytest.raises(MetaCoolerFormatError, match="File not found"):
            validate("/nonexistent")

    def test_wrong_format_attr(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = "wrong"
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            f.create_group("genomes")
            f.create_group("contigs")
            f.create_group("bins")
            f.create_group("pixels")
            f.create_group("indexes")
        with pytest.raises(MetaCoolerFormatError, match="Expected format"):
            validate(path)

    def test_missing_metacooler_format(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.create_group("genomes")
            f.create_group("contigs")
            f.create_group("bins")
            f.create_group("pixels")
            f.create_group("indexes")
        with pytest.raises(MetaCoolerFormatError, match="metacooler-format"):
            validate(path)

    def test_missing_group(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            f.create_group("genomes")
            f.create_group("contigs")
            f.create_group("bins")
            f.create_group("pixels")
            # missing indexes
        with pytest.raises(MetaCoolerFormatError, match="Missing required group"):
            validate(path)

    def test_missing_column(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            # missing length and n_contigs
            f.create_group("contigs")
            f.create_group("bins")
            f.create_group("pixels")
            f.create_group("indexes")
        with pytest.raises(MetaCoolerFormatError, match="Missing column"):
            validate(path)

    def test_bad_chrom_offset(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1]))
            c.create_dataset("genome_id", data=np.array([0]))
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([0]))
            p.create_dataset("bin2_id", data=np.array([0]))
            p.create_dataset("count", data=np.array([1]))
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 1, 2]))  # wrong length
            idx.create_dataset("bin1_offset", data=np.array([0, 1]))
        with pytest.raises(MetaCoolerFormatError, match="chrom_offset length"):
            validate(path)

    def test_bad_genome_id_range(self, tmp_path):
        path = tmp_path / "bad.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1]))
            c.create_dataset("genome_id", data=np.array([5]))  # out of range
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([0]))
            p.create_dataset("bin2_id", data=np.array([0]))
            p.create_dataset("count", data=np.array([1]))
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 1]))
            idx.create_dataset("bin1_offset", data=np.array([0, 1]))
            idx.create_dataset("genome_offset", data=np.array([0, 1]))
        with pytest.raises(
            MetaCoolerFormatError, match="genome_id values out of range"
        ):
            validate(path)

    def test_empty_pixels_warns(self, tmp_path):
        path = tmp_path / "empty_pix.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1000]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1000]))
            c.create_dataset("genome_id", data=np.array([0]))
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1000]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([], dtype=np.int64))
            p.create_dataset("bin2_id", data=np.array([], dtype=np.int64))
            p.create_dataset("count", data=np.array([], dtype=np.int32))
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 0]))
            idx.create_dataset("bin1_offset", data=np.array([0, 0]))
            idx.create_dataset("genome_offset", data=np.array([0, 0]))
        validate(path)

    def test_bin1_offset_wrong_length_raises(self, tmp_path):
        path = tmp_path / "bad_bin1off.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1000]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1000]))
            c.create_dataset("genome_id", data=np.array([0]))
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1000]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([0]))
            p.create_dataset("bin2_id", data=np.array([0]))
            p.create_dataset("count", data=np.array([1]))
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 1]))
            # bin1_offset should have length n_bins+1=2 but has 3
            idx.create_dataset("bin1_offset", data=np.array([0, 1, 2]))
            idx.create_dataset("genome_offset", data=np.array([0, 1]))
        with pytest.raises(MetaCoolerFormatError, match="bin1_offset length"):
            validate(path)

    def test_genome_offset_wrong_length_raises(self, tmp_path):
        """genome_offset is required and must be n_genomes+1."""
        path = tmp_path / "bad_goff.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1000]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1000]))
            c.create_dataset("genome_id", data=np.array([0]))
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1000]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([0]))
            p.create_dataset("bin2_id", data=np.array([0]))
            p.create_dataset("count", data=np.array([1]))
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 1]))
            idx.create_dataset("bin1_offset", data=np.array([0, 1]))
            # genome_offset should have length n_genomes+1=2 but has 3
            idx.create_dataset("genome_offset", data=np.array([0, 1, 2]))
        with pytest.raises(MetaCoolerFormatError, match="genome_offset length"):
            validate(path)

    def test_missing_pixels_columns(self, tmp_path):
        path = tmp_path / "no_count.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1000]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1000]))
            c.create_dataset("genome_id", data=np.array([0]))
            b = f.create_group("bins")
            b.create_dataset("chrom", data=np.array([0]))
            b.create_dataset("start", data=np.array([0]))
            b.create_dataset("end", data=np.array([1000]))
            p = f.create_group("pixels")
            p.create_dataset("bin1_id", data=np.array([0]))
            p.create_dataset("bin2_id", data=np.array([0]))
            # missing count
            idx = f.create_group("indexes")
            idx.create_dataset("chrom_offset", data=np.array([0, 1]))
            idx.create_dataset("bin1_offset", data=np.array([0, 1]))
            idx.create_dataset("genome_offset", data=np.array([0, 1]))
        with pytest.raises(MetaCoolerFormatError, match="Missing column in /pixels"):
            validate(path)

    def test_missing_contigs_column(self, tmp_path):
        path = tmp_path / "no_genome_id.metacool"
        with h5py.File(path, "w") as f:
            f.attrs["format"] = COOLER_FORMAT
            f.attrs["metacooler-format"] = METACOOLER_FORMAT
            g = f.create_group("genomes")
            g.create_dataset("name", data=np.array([b"g1"]))
            g.create_dataset("length", data=np.array([1000]))
            g.create_dataset("n_contigs", data=np.array([1]))
            g.create_dataset("n_bins", data=np.array([1]))
            c = f.create_group("contigs")
            c.create_dataset("name", data=np.array([b"c1"]))
            c.create_dataset("length", data=np.array([1000]))
            # missing genome_id — validation should fail here after genomes pass
            f.create_group("bins")
            f.create_group("pixels")
            f.create_group("indexes")
        with pytest.raises(MetaCoolerFormatError, match="Missing column in /contigs"):
            validate(path)
