"""Tests for metacooler.api.MetaCooler."""

from __future__ import annotations

import cooler
import numpy as np
import pytest

from metacooler.api import MetaCooler
from metacooler.schema import is_metacooler


class TestMetaCooler:
    def test_init(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        assert mc is not None

    def test_cooler_compat_is_cooler(self, tiny_metacooler):
        assert cooler.fileops.is_cooler(str(tiny_metacooler)) is True

    def test_is_metacooler(self, tiny_metacooler):
        assert is_metacooler(tiny_metacooler) is True

    def test_bins(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        bins = mc.bins()
        assert len(bins) > 0
        assert {"chrom", "start", "end"}.issubset(bins.columns)

    def test_pixels(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        pix = mc.pixels()
        assert len(pix) > 0
        assert {"bin1_id", "bin2_id", "count"}.issubset(pix.columns)

    def test_chroms(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        chroms = mc.chroms()
        assert len(chroms) == 4  # 4 contigs

    def test_genomes(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        genomes = mc.genomes
        assert len(genomes) == 3  # MAG_1, MAG_2, __unbinned__
        assert "__unbinned__" in genomes["name"].values

    def test_contigs(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        contigs = mc.contigs
        assert len(contigs) == 4
        assert "genome_id" in contigs.columns

    def test_info(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        info = mc.info
        assert info["ngenomes"] == 3
        assert info["ncontigs"] == 4
        assert info["format"] == "HDF5::MetaCooler"
        assert info["format-version"] == 3  # cooler format version
        assert info["metacooler-format-version"] == 1
        assert info["ordering"] == "genome-then-alpha"

    def test_bins_for_genome(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        bins = mc.bins_for_genome("MAG_1")
        assert len(bins) > 0

    def test_bins_for_contig(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        bins = mc.bins_for_contig("contig_A")
        assert len(bins) == 10  # 1000 bp / 100 bp bins

    def test_genome_extent(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        start, end = mc.genome_extent("MAG_1")
        assert end > start

    def test_contig_extent(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        start, end = mc.contig_extent("contig_A")
        assert end - start == 10  # 10 bins

    def test_aggregate_matrix_contig(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        mat = mc.aggregate_matrix(level="contig", sparse=False)
        assert mat.shape == (4, 4)
        assert np.allclose(mat, mat.T)  # symmetric  # type: ignore

    def test_aggregate_matrix_genome(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        mat = mc.aggregate_matrix(level="genome", sparse=False)
        assert mat.shape == (3, 3)
        assert np.allclose(mat, mat.T)  # type: ignore

    def test_to_cooler_contig(self, tiny_metacooler, tmp_path):
        mc = MetaCooler(tiny_metacooler)
        out = tmp_path / "out.cool"
        mc.to_cooler(out, level="contig")
        assert cooler.fileops.is_cooler(str(out))
        clr = cooler.Cooler(str(out))
        assert len(clr.chroms()) == 4

    def test_to_cooler_genome(self, tiny_metacooler, tmp_path):
        mc = MetaCooler(tiny_metacooler)
        out = tmp_path / "out.cool"
        mc.to_cooler(out, level="genome")
        assert cooler.fileops.is_cooler(str(out))
        clr = cooler.Cooler(str(out))
        assert len(clr.chroms()) == 3


class TestMetaCoolerErrors:
    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            MetaCooler(tmp_path / "nonexistent.metacool")

    def test_group_path_on_non_mmetacool_raises(self, tiny_metacooler):
        uri = f"{tiny_metacooler}::some/group"
        with pytest.raises(ValueError, match="not an mmetacool file"):
            MetaCooler(uri)

    def test_bins_for_genome_unknown(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="not found"):
            mc.bins_for_genome("NONEXISTENT")

    def test_bins_for_contig_unknown(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="not found"):
            mc.bins_for_contig("NONEXISTENT")

    def test_genome_extent_unknown(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="not found"):
            mc.genome_extent("NONEXISTENT")

    def test_contig_extent_unknown(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="not found"):
            mc.contig_extent("NONEXISTENT")

    def test_aggregate_matrix_invalid_level(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(ValueError, match="level must be"):
            mc.aggregate_matrix(level="invalid")  # type: ignore

    def test_pixels_for_genome_pair_unknown_genome(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with pytest.raises(KeyError, match="not found"):
            mc.pixels_for_genome_pair("MAG_1", "NONEXISTENT")


class TestMetaCoolerProperties:
    def test_filename_property(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        assert mc.filename == str(tiny_metacooler)

    def test_store_property(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        from pathlib import Path

        assert isinstance(mc.store, Path)
        assert mc.store == tiny_metacooler

    def test_root_property(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        assert mc.root == "/"

    def test_uri_property(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        assert "::" in mc.uri
        assert str(tiny_metacooler) in mc.uri

    def test_binsize_property(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        bs = mc.binsize
        # Should be None (variable bins) or 100 (fixed 100 bp bins)
        assert bs is None or bs == 100

    def test_chromnames(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        names = mc.chromnames
        assert isinstance(names, list)
        assert len(names) == 4

    def test_chromsizes(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        import pandas as pd

        sizes = mc.chromsizes
        assert isinstance(sizes, pd.Series)
        assert len(sizes) == 4

    def test_shape(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        s = mc.shape
        assert isinstance(s, tuple)
        assert len(s) == 2
        assert s[0] == s[1]  # square matrix

    def test_extent(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        ext = mc.extent("contig_A")
        assert isinstance(ext, tuple)
        assert len(ext) == 2
        assert ext[1] > ext[0]

    def test_offset(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        off = mc.offset("contig_A")
        assert isinstance(off, (int, __import__("numpy").integer))
        assert off >= 0

    def test_matrix(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        mat = mc.matrix(balance=False)
        assert mat is not None

    def test_info_has_required_keys(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        info = mc.info
        for key in (
            "format",
            "format-version",
            "ngenomes",
            "ncontigs",
            "nbins",
            "nnz",
            "ordering",
        ):
            assert key in info, f"Missing key in info: {key}"

    def test_info_nbins_correct(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        info = mc.info
        bins = mc.bins()
        assert info["nbins"] == len(bins)

    def test_info_ncontigs_correct(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        info = mc.info
        contigs = mc.contigs
        assert info["ncontigs"] == len(contigs)

    def test_info_ngenomes_correct(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        info = mc.info
        genomes = mc.genomes
        assert info["ngenomes"] == len(genomes)


class TestMetaCoolerPixelsForGenomePair:
    def test_same_genome(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        pix = mc.pixels_for_genome_pair("MAG_1", "MAG_1")
        # Intra-genome pixels should exist
        assert len(pix) >= 0  # may be 0 if no contacts, but should not raise

    def test_cross_genome(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        pix = mc.pixels_for_genome_pair("MAG_1", "MAG_2")
        assert isinstance(pix, __import__("pandas").DataFrame)
        assert {"bin1_id", "bin2_id", "count"}.issubset(pix.columns)

    def test_cross_genome_symmetric(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        pix_12 = mc.pixels_for_genome_pair("MAG_1", "MAG_2")
        pix_21 = mc.pixels_for_genome_pair("MAG_2", "MAG_1")
        assert len(pix_12) == len(pix_21)


class TestMetaCoolerAggregateMatrix:
    def test_sparse_returns_sparse_matrix(self, tiny_metacooler):
        from scipy import sparse as sp

        mc = MetaCooler(tiny_metacooler)
        mat = mc.aggregate_matrix(level="contig", sparse=True)
        assert sp.issparse(mat)

    def test_dense_symmetric(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        mat = mc.aggregate_matrix(level="genome", sparse=False)
        assert np.allclose(mat, mat.T)  # type: ignore

    def test_genome_shape_matches_ngenomes(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        n = len(mc.genomes)
        mat = mc.aggregate_matrix(level="genome", sparse=False)
        assert mat.shape == (n, n)

    def test_contig_diagonal_nonnegative(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        mat = mc.aggregate_matrix(level="contig", sparse=False)
        assert (np.diag(mat) >= 0).all()  # type: ignore


class TestMetaCoolerOpenContextManager:
    def test_open_returns_h5py_file(self, tiny_metacooler):
        import h5py

        mc = MetaCooler(tiny_metacooler)
        with mc.open("r") as f:
            assert isinstance(f, h5py.File)
            assert "genomes" in f
            assert "contigs" in f

    def test_open_read_mode(self, tiny_metacooler):
        mc = MetaCooler(tiny_metacooler)
        with mc.open("r") as f:
            names = f["genomes"]["name"][:]
            assert len(names) == 3


class TestMetaCoolerMmetacoolURI:
    def test_mmetacool_uri_opens_entry(self, tmp_path, tiny_metacooler):
        from metacooler.io.mmetacool import create_mmetacool

        mmetacool = tmp_path / "test.mmetacool"
        create_mmetacool(mmetacool, {"entry1": tiny_metacooler})
        uri = f"{mmetacool}::entries/entry1"
        mc = MetaCooler(uri)
        assert mc is not None
        assert len(mc.genomes) == 3

    def test_mmetacool_uri_missing_group_raises(self, tmp_path, tiny_metacooler):
        from metacooler.io.mmetacool import create_mmetacool

        mmetacool = tmp_path / "test.mmetacool"
        create_mmetacool(mmetacool, {"entry1": tiny_metacooler})
        uri = f"{mmetacool}::entries/nonexistent"
        with pytest.raises(KeyError):
            MetaCooler(uri)

    def test_mmetacool_entry_info(self, tmp_path, tiny_metacooler):
        from metacooler.io.mmetacool import create_mmetacool

        mmetacool = tmp_path / "test.mmetacool"
        create_mmetacool(mmetacool, {"myentry": tiny_metacooler})
        uri = f"{mmetacool}::entries/myentry"
        mc = MetaCooler(uri)
        info = mc.info
        assert info["format"] == "HDF5::MetaCooler"
        assert info["ncontigs"] == 4
