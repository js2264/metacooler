"""Edge-case tests for metacooler.io.create."""

from __future__ import annotations

import bz2
import gzip
from pathlib import Path

import h5py
import pandas as pd
import pytest

from metacooler.io.create import (
    _build_bins,
    _build_contig_order,
    _build_contigs_table,
    _build_genomes_table,
    _drop_nonpositive_contigs,
    _parse_pairs_header,
    _read_contig_lengths,
    _stream_pairs,
    _tile_bins,
    _total_time_str,
    _validate_bin_coords,
    from_cool,
    from_pairs,
)
from metacooler.io.genomes import parse_genomes_table


class TestBuildContigOrder:
    def test_fasta_order(self):
        lengths = {"c": 100, "a": 200, "b": 150}
        cg_map = {"a": "MAG_1", "b": "MAG_1", "c": "MAG_2"}
        order = _build_contig_order(lengths, cg_map, ordering="fasta")
        assert order == ["c", "a", "b"]

    def test_alpha_order(self):
        lengths = {"c": 100, "a": 200, "b": 150}
        cg_map = {"a": "MAG_1", "b": "MAG_1", "c": "MAG_2"}
        order = _build_contig_order(lengths, cg_map, ordering="alpha")
        assert order == ["a", "b", "c"]

    def test_unknown_ordering_raises(self):
        with pytest.raises(ValueError, match="Unknown ordering"):
            _build_contig_order({"a": 100}, {"a": "MAG_1"}, ordering="bogus")  # ty:ignore[invalid-argument-type]


class TestBuildGenomesTable:
    def test_unbinned_absent_when_all_contigs_assigned(self):
        """No phantom __unbinned__ row when every contig has an explicit genome."""
        contigs = ["a", "b"]
        lengths = {"a": 100, "b": 200}
        cg_map = {"a": "MAG_1", "b": "MAG_1"}
        df = _build_genomes_table(contigs, lengths, cg_map)
        assert "__unbinned__" not in df["name"].values
        assert list(df.columns[:3]) == ["name", "length", "n_contigs"]

    def test_unbinned_present_when_contig_unassigned(self):
        """__unbinned__ row is injected when at least one contig has no assignment."""
        contigs = ["a", "b", "c"]
        lengths = {"a": 100, "b": 200, "c": 50}
        cg_map = {"a": "MAG_1", "b": "MAG_1"}  # "c" is unassigned
        df = _build_genomes_table(contigs, lengths, cg_map)
        assert "__unbinned__" in df["name"].values

    def test_extra_columns_from_genomes_file(self):
        contigs = ["a", "b", "c"]
        lengths = {"a": 100, "b": 200, "c": 50}
        cg_map = {"a": "MAG_1", "b": "MAG_1", "c": "MAG_2"}
        extra = pd.DataFrame(
            {
                "genome_id": ["MAG_1", "MAG_2"],
                "col_1": ["high", "medium"],
            }
        )
        df = _build_genomes_table(contigs, lengths, cg_map, extra_genomes_df=extra)
        assert "col_1" in df.columns
        mag1_row = df[df["name"] == "MAG_1"].iloc[0]
        assert mag1_row["col_1"] == "high"

    def test_extra_columns_no_match(self):
        contigs = ["a"]
        lengths = {"a": 100}
        cg_map = {"a": "MAG_1"}
        extra = pd.DataFrame(
            {
                "genome_id": ["MAG_X"],
                "col_1": ["value"],
            }
        )
        df = _build_genomes_table(contigs, lengths, cg_map, extra_genomes_df=extra)
        # MAG_1 has no match → col_1 is NaN
        assert "col_1" in df.columns

    def test_extra_columns_numeric_conversion(self):
        """Test that numeric string columns are converted to float."""
        contigs = ["a", "b"]
        lengths = {"a": 100, "b": 200}
        cg_map = {"a": "MAG_1", "b": "MAG_1"}
        extra = pd.DataFrame(
            {
                "genome_id": ["MAG_1"],
                "completeness": ["95.5"],
            }
        )
        df = _build_genomes_table(contigs, lengths, cg_map, extra_genomes_df=extra)
        assert "completeness" in df.columns
        mag1_row = df[df["name"] == "MAG_1"].iloc[0]
        assert mag1_row["completeness"] == pytest.approx(95.5)


class TestParseGenomesTable:
    def test_single_column(self, tmp_path):
        path = tmp_path / "genomes.tsv"
        path.write_text("genome_name\nMAG_1\nMAG_2\n")
        df = parse_genomes_table(path)
        assert list(df.columns) == ["genome_id"]
        assert df["genome_id"].tolist() == ["MAG_1", "MAG_2"]

    def test_multiple_columns(self, tmp_path):
        path = tmp_path / "genomes.tsv"
        path.write_text(
            "genome_name\tquality\tphylum\nMAG_1\thigh\tBacteroidetes\nMAG_2\tmedium\tFirmicutes\n"
        )
        df = parse_genomes_table(path)
        assert list(df.columns) == ["genome_id", "quality", "phylum"]
        assert df["genome_id"].tolist() == ["MAG_1", "MAG_2"]
        assert df["phylum"].tolist() == ["Bacteroidetes", "Firmicutes"]

    def test_wrong_first_column(self, tmp_path):
        path = tmp_path / "genomes.tsv"
        path.write_text("genome_id\nMAG_1\n")
        with pytest.raises(ValueError, match="genome_name"):
            parse_genomes_table(path)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_genomes_table("/nonexistent.tsv")

    def test_genomes_table_in_create(
        self, tmp_path, tiny_fasta, tiny_bins_bed, tiny_pairs, tiny_contigs_table
    ):
        genomes_path = tmp_path / "genomes.tsv"
        genomes_path.write_text(
            "genome_name\tquality\nMAG_1\thigh_quality\nMAG_2\tmedium_quality\n"
        )
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
            genomes_table=genomes_path,
            assembly_name="test",
            chunksize=1000,
        )
        assert out.exists()
        from metacooler.api import MetaCooler

        mc = MetaCooler(out)
        genomes_df = mc.genomes
        assert "quality" in genomes_df.columns
        assert "high_quality" in genomes_df["quality"].values


class TestReadContigLengths:
    def test_from_bins(self):
        import pandas as pd

        bins = pd.DataFrame(
            {"chrom": ["a", "a", "b"], "start": [0, 50, 0], "end": [50, 100, 80]}
        )
        lengths = _read_contig_lengths(bins=bins)
        assert lengths == {"a": 100, "b": 80}

    def test_no_source_raises(self):
        with pytest.raises(ValueError, match="Contig lengths required"):
            _read_contig_lengths()

    def test_from_gzip_fasta(self, tmp_path):
        """Test _read_contig_lengths with a gzipped FASTA."""
        fa_path = tmp_path / "test.fa.gz"
        with gzip.open(fa_path, "wt") as fh:
            fh.write(">contig_a\n" + "A" * 1000 + "\n")
            fh.write(">contig_b\n" + "T" * 500 + "\n")
        lengths = _read_contig_lengths(fasta=fa_path)
        assert lengths == {"contig_a": 1000, "contig_b": 500}


class TestParsePairsHeader:
    def test_parse_chromsizes(self, tmp_path):
        path = tmp_path / "pairs.pairs"
        path.write_text(
            "#chromsize: chr1 1000\n#chromsize: chr2 500\nread_0\tchr1\t1\tchr2\t1\n"
        )
        result = _parse_pairs_header(path)
        assert result == {"chr1": 1000, "chr2": 500}

    def test_parse_empty(self, tmp_path):
        path = tmp_path / "pairs.pairs"
        path.write_text("read_0\tchr1\t1\tchr2\t1\n")
        result = _parse_pairs_header(path)
        assert result is None

    def test_parse_gzipped(self, tmp_path):
        path = tmp_path / "pairs.pairs.gz"
        with gzip.open(path, "wt") as f:
            f.write("#chromsize: chr1 1000\nread_0\tchr1\t1\tchr1\t1\n")
        result = _parse_pairs_header(path)
        assert result == {"chr1": 1000}


class TestStreamPairs:
    def test_stream_plain(self, tmp_path):
        path = tmp_path / "pairs.pairs"
        path.write_text("#c\nread_0\tchr1\t10\tchr2\t20\t+\t+\n")
        chunks = list(_stream_pairs(path, chunksize=1))
        assert len(chunks) == 1
        assert chunks[0]["chrom1"].iloc[0] == "chr1"
        assert chunks[0]["pos1"].iloc[0] == 10

    def test_stream_gzipped(self, tmp_path):
        path = tmp_path / "pairs.pairs.gz"
        with gzip.open(path, "wt") as f:
            f.write("#c\nread_0\tchr1\t10\tchr2\t20\t+\t+\n")
        chunks = list(_stream_pairs(path, chunksize=1))
        assert len(chunks) == 1
        assert chunks[0]["chrom1"].iloc[0] == "chr1"

    def test_stream_bz2(self, tmp_path):
        """Test _stream_pairs with a bz2-compressed pairs file."""
        path = tmp_path / "pairs.pairs.bz2"
        with bz2.open(path, "wt") as fh:
            fh.write("#c\nread_0\tchr1\t10\tchr2\t20\t+\t+\n")
        chunks = list(_stream_pairs(path, chunksize=1))
        assert len(chunks) == 1
        assert chunks[0]["chrom1"].iloc[0] == "chr1"


class TestTileBins:
    def test_tile(self):
        order = ["a", "b"]
        lengths = {"a": 250, "b": 100}
        name_to_id = {"a": 0, "b": 1}
        df = _tile_bins(order, lengths, name_to_id, binsize=100)
        assert len(df) == 4  # a: 3 bins, b: 1 bin
        assert df["chrom"].tolist() == [0, 0, 0, 1]
        assert df["end"].tolist() == [100, 200, 250, 100]


class TestBuildContigsTableExtra:
    """Tests for _build_contigs_table with extra columns."""

    def test_extra_columns_numeric_conversion(self):
        """Test that numeric string columns in the contigs table are converted."""
        genomes_df = pd.DataFrame(
            {"name": ["MAG_1", "__unbinned__"], "length": [300, 0], "n_contigs": [2, 0]}
        )
        contig_df = pd.DataFrame(
            {
                "contig_name": ["a", "b"],
                "genome_id": ["MAG_1", "MAG_1"],
                "gc_content": ["0.52", "0.48"],
            }
        )
        df = _build_contigs_table(
            contig_order=["a", "b"],
            contig_lengths={"a": 200, "b": 100},
            contig_genome_map={"a": "MAG_1", "b": "MAG_1"},
            genomes=genomes_df,
            contig_df=contig_df,
        )
        assert "gc_content" in df.columns
        # Numeric string "0.52" should be converted to float
        assert df.loc[df["name"] == "a", "gc_content"].iloc[0] == pytest.approx(0.52)

    def test_extra_columns_non_numeric_string(self):
        """Test that non-numeric extra columns remain as strings."""
        genomes_df = pd.DataFrame(
            {"name": ["MAG_1", "__unbinned__"], "length": [300, 0], "n_contigs": [2, 0]}
        )
        contig_df = pd.DataFrame(
            {
                "contig_name": ["a"],
                "genome_id": ["MAG_1"],
                "taxonomy": ["Bacteria"],
            }
        )
        df = _build_contigs_table(
            contig_order=["a"],
            contig_lengths={"a": 200},
            contig_genome_map={"a": "MAG_1"},
            genomes=genomes_df,
            contig_df=contig_df,
        )
        assert "taxonomy" in df.columns
        assert df.loc[df["name"] == "a", "taxonomy"].iloc[0] == "Bacteria"


class TestFromPairsValidation:
    """Tests for argument validation and error paths in from_pairs."""

    def test_invalid_contig_genome_map_type(self, tmp_path, tiny_bins_bed, tiny_pairs):
        with pytest.raises(TypeError, match="contig_genome_map must be"):
            from_pairs(
                uri=tmp_path / "out.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=42,  # type: ignore
            )

    def test_contig_genome_map_as_dict(self, tmp_path, tiny_bins_bed, tiny_pairs):
        """from_pairs accepts a dict for contig_genome_map."""
        cg_map = {
            "contig_A": "MAG_1",
            "contig_B": "MAG_1",
            "contig_C": "MAG_2",
        }
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=cg_map,
        )
        assert out.exists()

    def test_bins_as_dataframe(self, tmp_path, tiny_pairs, tiny_contigs_table):
        """Test that bins can be passed as a pre-built DataFrame."""
        bins_df = pd.DataFrame(
            {
                "chrom": ["contig_A", "contig_A", "contig_B", "contig_C", "contig_D"],
                "start": [0, 500, 0, 0, 0],
                "end": [500, 1000, 600, 400, 200],
            }
        )
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=bins_df,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
        )
        assert out.exists()

    def test_invalid_bins_type(self, tmp_path, tiny_pairs, tiny_contigs_table):
        with pytest.raises(TypeError, match="bins must be"):
            from_pairs(
                uri=tmp_path / "out.metacool",
                bins=42,  # type: ignore
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
            )

    def test_binsize_required_when_bins_none(
        self, tmp_path, tiny_fasta, tiny_pairs, tiny_contigs_table
    ):
        """When bins=None and binsize=None, from_pairs should raise."""
        with pytest.raises(ValueError, match="binsize required"):
            from_pairs(
                uri=tmp_path / "out.metacool",
                bins=None,
                fasta=tiny_fasta,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                binsize=None,
            )

    def test_contigs_extra_columns_via_from_pairs(
        self, tmp_path, tiny_bins_bed, tiny_pairs
    ):
        """Contigs table with extra numeric columns should be stored correctly."""
        ctg_table = tmp_path / "ctg_extras.tsv"
        ctg_table.write_text(
            "contig_name\tgenome_name\tgc_content\n"
            "contig_A\tMAG_1\t0.52\n"
            "contig_B\tMAG_1\t0.48\n"
            "contig_C\tMAG_2\t0.55\n"
        )
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=ctg_table,
        )
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "gc_content" in f["contigs"]


class TestFromPairsEmptyChunks:
    """Tests for the empty-chunk fast-path in pixel_stream."""

    def test_empty_chunks_binsize_path(self, tmp_path, tiny_fasta, tiny_contigs_table):
        """Chunks with unknown contigs are skipped in the binsize path."""
        pairs_path = tmp_path / "mixed.pairs"
        # Two reads: one on unknown contig (→ empty chunk), one valid
        pairs_path.write_text(
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
            "read_0\tunknown_ctg\t100\tunknown_ctg\t200\t+\t+\n"
            "read_1\tcontig_A\t100\tcontig_A\t200\t+\t+\n"
        )
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=None,
            binsize=1000,
            fasta=tiny_fasta,
            pairs=pairs_path,
            contig_genome_map=tiny_contigs_table,
            chunksize=1,  # one read per chunk so we get at least one empty chunk
        )
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"][:].sum() >= 1

    def test_empty_chunks_fragment_path(
        self, tmp_path, tiny_bins_bed, tiny_contigs_table
    ):
        """Chunks with unknown contigs are skipped in the fragment-level path."""
        pairs_path = tmp_path / "mixed.pairs"
        pairs_path.write_text(
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
            "read_0\tunknown_ctg\t100\tunknown_ctg\t200\t+\t+\n"
            "read_1\tcontig_A\t100\tcontig_A\t200\t+\t+\n"
        )
        out = tmp_path / "out.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=pairs_path,
            contig_genome_map=tiny_contigs_table,
            chunksize=1,  # one read per chunk so we get at least one empty chunk
        )
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"][:].sum() >= 1


_COOL_PATH = Path("tests/data/assembly_1_rep1_1kb.cool")
_CONTIGS_TABLE = Path("tests/data/contigs_table_1.tsv")


class TestFromCoolExtras:
    """Tests for optional from_cool parameters."""

    def test_from_cool_with_genomes_table(self, tmp_path):
        """from_cool with genomes_table stores extra genome metadata."""
        genomes_table = tmp_path / "genomes.tsv"
        genomes_table.write_text("genome_name\tquality\nMAG_01\thigh\nMAG_02\tmedium\n")
        out = tmp_path / "out.metacool"
        from_cool(
            cool_uri=_COOL_PATH,
            output=out,
            contig_genome_map=_CONTIGS_TABLE,
            genomes_table=genomes_table,
        )
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "quality" in f["genomes"]

    def test_from_cool_with_assembly_name(self, tmp_path):
        """from_cool stores assembly_name in file attrs."""
        out = tmp_path / "out.metacool"
        from_cool(
            cool_uri=_COOL_PATH,
            output=out,
            contig_genome_map=_CONTIGS_TABLE,
            assembly_name="my_assembly",
        )
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert f.attrs["genome-assembly"] == "my_assembly"

    def test_from_cool_drop_unassigned_succeeds(self, tmp_path):
        """drop_unassigned=True should create a file with only assigned contigs."""
        out = tmp_path / "out.metacool"
        from_cool(
            cool_uri=_COOL_PATH,
            output=out,
            contig_genome_map=_CONTIGS_TABLE,
            drop_unassigned=True,
        )
        assert out.exists()
        import cooler

        clr = cooler.Cooler(str(out))
        # All chroms in the output must come from the contigs table
        import pandas as pd

        contigs_table = pd.read_csv(_CONTIGS_TABLE, sep="\t")
        assigned_contigs = set(contigs_table.iloc[:, 0].astype(str))
        output_chroms = set(clr.chroms()[:]["name"].astype(str))
        assert output_chroms <= assigned_contigs


class TestValidateBinCoords:
    def test_negative_coords_raise(self):
        df = pd.DataFrame({"chrom": ["c"], "start": [-1], "end": [100]})
        with pytest.raises(ValueError, match="non-negative"):
            _validate_bin_coords(df)

    def test_inverted_coords_raise(self):
        df = pd.DataFrame({"chrom": ["c"], "start": [100], "end": [50]})
        with pytest.raises(ValueError, match="end < start"):
            _validate_bin_coords(df)

    def test_zero_width_warns_but_ok(self, caplog):
        df = pd.DataFrame({"chrom": ["c"], "start": [100], "end": [100]})
        with caplog.at_level("WARNING"):
            _validate_bin_coords(df)  # must not raise
        assert any("zero-width" in r.message for r in caplog.records)

    def test_valid_coords_pass(self):
        df = pd.DataFrame({"chrom": ["c", "c"], "start": [0, 100], "end": [100, 200]})
        _validate_bin_coords(df)  # no raise


class TestBuildBinsErrors:
    def test_unknown_contig_raises(self):
        bins = pd.DataFrame({"chrom": ["ghost"], "start": [0], "end": [100]})
        with pytest.raises(ValueError, match="unknown contig"):
            _build_bins(bins, {"c": 0})


class TestParsePairsHeaderErrors:
    def test_malformed_chromsize_raises(self, tmp_path: Path):
        p = tmp_path / "bad.pairs"
        p.write_text(
            "#chromsize: contigA notanumber\nread1\tcontigA\t1\tcontigA\t2\t+\t+\n"
        )
        with pytest.raises(ValueError, match="Malformed #chromsize"):
            _parse_pairs_header(p)

    def test_no_chromsize_returns_none(self, tmp_path: Path):
        p = tmp_path / "plain.pairs"
        p.write_text("read1\tcontigA\t1\tcontigA\t2\t+\t+\n")
        assert _parse_pairs_header(p) is None


class TestTotalTimeStr:
    def test_seconds(self):
        import time

        assert _total_time_str(time.perf_counter()).endswith("s")

    def test_minutes(self):
        import time

        s = _total_time_str(time.perf_counter() - 90)
        assert "m:" in s and "h:" not in s

    def test_hours(self):
        import time

        assert "h:" in _total_time_str(time.perf_counter() - 3700)


class TestCreateCleanupOnFailure:
    def test_from_pairs_removes_partial_file_on_validate_error(
        self,
        tmp_path: Path,
        tiny_fasta: Path,
        tiny_contigs_table: Path,
        tiny_bins_bed: Path,
        tiny_pairs: Path,
        monkeypatch,
    ):
        # A validation failure mid-creation must not leave a half-written file.
        import metacooler.io.create as create_mod

        def _boom(*args, **kwargs):
            raise ValueError("simulated validation failure")

        monkeypatch.setattr(create_mod, "validate", _boom)
        out = tmp_path / "should_be_removed.metacool"
        with pytest.raises(ValueError, match="simulated validation failure"):
            from_pairs(
                uri=out,
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                contig_table_format="contigs",
                fasta=tiny_fasta,
                chunksize=1000,
            )
        assert not out.exists()


class TestDropNonpositiveContigs:
    def test_drops_zero_and_negative_length(self, caplog):
        lengths = {"good": 100, "empty": 0, "neg": -5}
        with caplog.at_level("WARNING"):
            out_lengths, _ = _drop_nonpositive_contigs(lengths)
        assert out_lengths == {"good": 100}
        assert any("non-positive length" in r.message for r in caplog.records)

    def test_filters_bins_in_lockstep(self):
        lengths = {"good": 100, "empty": 0}
        bins = pd.DataFrame(
            {
                "chrom": ["good", "good", "empty"],
                "start": [0, 50, 0],
                "end": [50, 100, 0],
            }
        )
        out_lengths, out_bins = _drop_nonpositive_contigs(lengths, bins)
        assert out_lengths == {"good": 100}
        assert out_bins is not None
        assert set(out_bins["chrom"]) == {"good"}

    def test_no_change_when_all_positive(self):
        lengths = {"a": 100, "b": 200}
        out_lengths, out_bins = _drop_nonpositive_contigs(lengths, None)
        assert out_lengths == lengths
        assert out_bins is None


class TestFromPairsGzip:
    def test_gzipped_pairs_roundtrip(
        self,
        tmp_path: Path,
        tiny_fasta: Path,
        tiny_contigs_table: Path,
        tiny_bins_bed: Path,
        tiny_pairs: Path,
    ):
        # A .pairs.gz input exercises the gzip reader and the size-based
        # sort_in_memory auto-detection (ISIZE trailer branch).
        gz = tmp_path / "pairs.pairs.gz"
        with open(tiny_pairs, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())

        out = tmp_path / "from_gz.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=gz,
            contig_genome_map=tiny_contigs_table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
            chunksize=1000,
        )
        assert out.exists()
        from metacooler.api import MetaCooler

        assert MetaCooler(out).shape[0] > 0


class TestBufferPixelsFolding:
    """``_buffer_pixels`` must decide on deduplicated pixels, not raw rows.

    The sort-strategy decision used to be driven by the sum of per-chunk row
    counts.  Bin pairs recur across chunks and Hi-C pixel counts saturate on a
    fixed bin grid, so that proxy overestimates the real pixel table badly on a
    deep library — tripping the budget and paying for an external sort (one
    ``cooler.create`` call per chunk, seconds each on a metagenome) that the
    data never needed.  Folding the buffer as it fills makes the budget measure
    the thing it is supposed to bound.
    """

    @staticmethod
    def _chunks(n_chunks: int, rows: int, n_bins: int, seed: int = 0):
        """Build aggregated pixel chunks drawn from a bin grid of *n_bins*."""
        import numpy as np

        rng = np.random.default_rng(seed)
        out = []
        for _ in range(n_chunks):
            b1 = rng.integers(0, n_bins, rows)
            b2 = rng.integers(0, n_bins, rows)
            frame = pd.DataFrame(
                {
                    "bin1_id": np.minimum(b1, b2),
                    "bin2_id": np.maximum(b1, b2),
                    "count": np.ones(rows, dtype="int64"),
                }
            )
            out.append(
                frame.groupby(["bin1_id", "bin2_id"], sort=False, as_index=False)[
                    "count"
                ].sum()
            )
        return out

    def test_saturating_stream_stays_in_memory(self):
        """Raw rows over budget but distinct pixels under it must not spill."""
        from metacooler.io.create import _buffer_pixels

        chunks = self._chunks(n_chunks=40, rows=5000, n_bins=200)
        raw_rows = sum(len(c) for c in chunks)
        distinct = len(
            pd.concat(chunks)
            .groupby(["bin1_id", "bin2_id"], as_index=False)["count"]
            .sum()
        )
        budget = 50_000
        # The premise of the test: the two measures straddle the budget.
        assert distinct < budget < raw_rows

        frames, rest = _buffer_pixels(iter(chunks), max_rows=budget)
        assert rest is None, "a saturating stream must not fall back to external sort"
        assert len(frames) == 1, "the buffer should be folded to a single frame"

    def test_folding_preserves_counts_exactly(self):
        """Folding must be arithmetically identical to one global aggregate."""
        from metacooler.io.create import _buffer_pixels

        chunks = self._chunks(n_chunks=40, rows=5000, n_bins=200, seed=1)
        expected = (
            pd.concat(chunks)
            .groupby(["bin1_id", "bin2_id"], as_index=False)["count"]
            .sum()
            .sort_values(["bin1_id", "bin2_id"])
            .reset_index(drop=True)
        )
        frames, _ = _buffer_pixels(iter(chunks), max_rows=50_000)
        got = (
            frames[0]
            .sort_values(["bin1_id", "bin2_id"])
            .reset_index(drop=True)
            .astype(expected.dtypes.to_dict())
        )
        pd.testing.assert_frame_equal(got, expected)
        assert got["count"].sum() == sum(int(c["count"].sum()) for c in chunks)

    def test_non_saturating_stream_still_spills(self):
        """The external-sort safety net must survive the optimisation.

        Folding narrows *when* the fallback fires; it must not remove it, or a
        genuinely huge pixel table would be held entirely in RAM.
        """
        from metacooler.io.create import _buffer_pixels

        chunks = self._chunks(n_chunks=40, rows=5000, n_bins=10_000_000)
        frames, rest = _buffer_pixels(iter(chunks), max_rows=50_000)
        assert rest is not None
        assert len(frames) == 1, "the partial accumulator is handed back as one frame"

    def test_empty_stream(self):
        from metacooler.io.create import _buffer_pixels

        assert _buffer_pixels(iter([]), max_rows=1000) == ([], None)

    def test_single_oversized_chunk_spills(self):
        """One chunk larger than the whole budget must still be detected."""
        from metacooler.io.create import _buffer_pixels

        chunks = self._chunks(n_chunks=1, rows=20_000, n_bins=10_000_000)
        frames, rest = _buffer_pixels(iter(chunks), max_rows=100)
        assert rest is not None
        assert len(frames) == 1
