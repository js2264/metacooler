"""Comprehensive CLI integration tests using real data files.

Covers all major CLI commands from cli_tests.sh:
- create (pairs, cool)
- info / tree / attrs
- dump (tables, range queries)
- stack / merge
- metadata (list, add, remove, get, load)
- ls / cp / mv / rm
- flatten (round-trip)
- extract (by genome, by contig, with --split)
- binnify
- zoomify
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import cooler
import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.io.mmetacool import list_entries
from metacooler.schema import (
    COOLER_FORMAT,
    METACOOLER_FORMAT,
    METACOOLER_FORMAT_VERSION,
)

# ---------------------------------------------------------------------------
# Constants – expected values derived from running cli_tests.sh
# ---------------------------------------------------------------------------

DATA = Path(__file__).parent / "data"

# assembly_1_rep1_1kb.metacool properties
A1_REP1_1KB_NGENOMES = 11
A1_REP1_1KB_NCONTIGS = 44
A1_REP1_1KB_NBINS = 661
A1_REP1_1KB_NNZ = 21633
A1_REP1_1KB_SUM = 50000
A1_REP1_1KB_BINSIZE = 1000

# assembly_1_rep2_1kb merged properties
MERGED_NNZ = 36249
MERGED_SUM = 100000

# ---------------------------------------------------------------------------
# Exact expected table content for assembly_1_rep1_1kb.metacool
# (verified by running `metacooler dump -t <table> assembly_1_rep1_1kb.metacool`)
# Each tuple represents one tab-separated row.
# ---------------------------------------------------------------------------

# All 11 genome rows (genome_name, length, n_contigs, n_bins, [phylum, …])
# __unbinned__ has empty metadata columns (tabs present but values empty).
EXPECTED_GENOMES_ROWS: list[tuple] = [
    ("__unbinned__", "42000", "5", "42", "", "", "", ""),
    ("MAG_01", "60000", "3", "60", "Bacteroidetes", "0.63", "0.1", "high"),
    ("MAG_02", "124000", "4", "124", "Firmicutes", "0.92", "0.14", "low"),
    ("MAG_03", "80000", "5", "80", "Proteobacteria", "0.73", "0.03", "low"),
    ("MAG_04", "24000", "3", "24", "Firmicutes", "0.99", "0.24", "medium"),
    ("MAG_05", "41000", "4", "41", "Proteobacteria", "0.89", "0.27", "high"),
    ("MAG_06", "83000", "5", "83", "Proteobacteria", "0.56", "0.49", "high"),
    ("MAG_07", "30000", "3", "30", "Proteobacteria", "0.89", "0.37", "high"),
    ("MAG_08", "49000", "4", "49", "Firmicutes", "0.9", "0.47", "high"),
    ("MAG_09", "103000", "5", "103", "Proteobacteria", "0.95", "0.04", "low"),
    ("MAG_10", "25000", "3", "25", "Proteobacteria", "0.99", "0.42", "medium"),
]

# First 20 contig rows (name, length, genome, genome_id, gc_content, coverage).
# Contigs are laid out grouped by genome: __unbinned__ (genome_id 0) first, then
# MAG_01 (1), MAG_02 (2), … — longest contig first within each genome.  The
# genome column holds the MAG name; genome_id keeps the internal integer index.
EXPECTED_CONTIGS_ROWS_20: list[tuple] = [
    ("contig_1_010", "12000", "__unbinned__", "0", "0.5", "91.35"),
    ("contig_1_035", "9000", "__unbinned__", "0", "0.54", "63.57"),
    ("contig_1_032", "8000", "__unbinned__", "0", "0.41", "87.76"),
    ("contig_1_015", "7000", "__unbinned__", "0", "0.47", "72.0"),
    ("contig_1_003", "6000", "__unbinned__", "0", "0.52", "8.28"),
    ("contig_1_023", "25000", "MAG_01", "1", "0.42", "66.93"),
    ("contig_1_036", "22000", "MAG_01", "1", "0.56", "37.79"),
    ("contig_1_012", "13000", "MAG_01", "1", "0.49", "45.37"),
    ("contig_1_005", "50000", "MAG_02", "2", "0.56", "63.69"),
    ("contig_1_034", "50000", "MAG_02", "2", "0.41", "28.51"),
    ("contig_1_001", "12000", "MAG_02", "2", "0.59", "75.43"),
    ("contig_1_020", "12000", "MAG_02", "2", "0.57", "35.01"),
    ("contig_1_029", "27000", "MAG_03", "3", "0.54", "11.43"),
    ("contig_1_019", "23000", "MAG_03", "3", "0.51", "10.62"),
    ("contig_1_027", "10000", "MAG_03", "3", "0.44", "74.74"),
    ("contig_1_030", "10000", "MAG_03", "3", "0.58", "67.79"),
    ("contig_1_033", "10000", "MAG_03", "3", "0.58", "59.19"),
    ("contig_1_000", "18000", "MAG_04", "4", "0.51", "99.4"),
    ("contig_1_002", "3000", "MAG_04", "4", "0.45", "8.81"),
    ("contig_1_037", "3000", "MAG_04", "4", "0.6", "93.22"),
]

# First 20 chrom rows (name, length) – same order as contigs
EXPECTED_CHROMS_ROWS_20: list[tuple] = [
    ("contig_1_010", "12000"),
    ("contig_1_035", "9000"),
    ("contig_1_032", "8000"),
    ("contig_1_015", "7000"),
    ("contig_1_003", "6000"),
    ("contig_1_023", "25000"),
    ("contig_1_036", "22000"),
    ("contig_1_012", "13000"),
    ("contig_1_005", "50000"),
    ("contig_1_034", "50000"),
    ("contig_1_001", "12000"),
    ("contig_1_020", "12000"),
    ("contig_1_029", "27000"),
    ("contig_1_019", "23000"),
    ("contig_1_027", "10000"),
    ("contig_1_030", "10000"),
    ("contig_1_033", "10000"),
    ("contig_1_000", "18000"),
    ("contig_1_002", "3000"),
    ("contig_1_037", "3000"),
]

# First 20 bin rows (chrom_name, start, end) – 1 kb bins on the first contigs
EXPECTED_BINS_ROWS_20: list[tuple] = [
    ("contig_1_010", "0", "1000"),
    ("contig_1_010", "1000", "2000"),
    ("contig_1_010", "2000", "3000"),
    ("contig_1_010", "3000", "4000"),
    ("contig_1_010", "4000", "5000"),
    ("contig_1_010", "5000", "6000"),
    ("contig_1_010", "6000", "7000"),
    ("contig_1_010", "7000", "8000"),
    ("contig_1_010", "8000", "9000"),
    ("contig_1_010", "9000", "10000"),
    ("contig_1_010", "10000", "11000"),
    ("contig_1_010", "11000", "12000"),
    ("contig_1_035", "0", "1000"),
    ("contig_1_035", "1000", "2000"),
    ("contig_1_035", "2000", "3000"),
    ("contig_1_035", "3000", "4000"),
    ("contig_1_035", "4000", "5000"),
    ("contig_1_035", "5000", "6000"),
    ("contig_1_035", "6000", "7000"),
    ("contig_1_035", "7000", "8000"),
]

# First 20 pixel rows (bin1_id, bin2_id, count)
EXPECTED_PIXELS_ROWS_20: list[tuple] = [
    ("0", "0", "1"),
    ("0", "1", "1"),
    ("0", "2", "1"),
    ("0", "3", "1"),
    ("0", "4", "2"),
    ("0", "5", "2"),
    ("0", "6", "4"),
    ("0", "7", "3"),
    ("0", "8", "2"),
    ("0", "9", "2"),
    ("0", "11", "1"),
    ("0", "12", "4"),
    ("0", "13", "3"),
    ("0", "14", "5"),
    ("0", "15", "3"),
    ("0", "16", "4"),
    ("0", "17", "3"),
    ("0", "18", "2"),
    ("0", "19", "3"),
    ("0", "20", "5"),
]

# First 20 pixels for range query contig_1_023:5000-15000.
# contig_1_023 occupies bins 42–66; bins 47–56 = positions 5000–15000.
EXPECTED_PIXELS_RANGE_ROWS_20: list[tuple] = [
    ("47", "47", "1"),
    ("47", "48", "1"),
    ("47", "50", "1"),
    ("47", "52", "2"),
    ("47", "54", "2"),
    ("47", "55", "1"),
    ("47", "56", "2"),
    ("48", "48", "1"),
    ("48", "49", "2"),
    ("48", "50", "1"),
    ("48", "51", "5"),
    ("48", "52", "2"),
    ("48", "53", "2"),
    ("48", "55", "2"),
    ("48", "56", "5"),
    ("49", "49", "1"),
    ("49", "50", "2"),
    ("49", "51", "1"),
    ("49", "52", "2"),
    ("49", "53", "1"),
]


def _parse_rows(output: str) -> list[tuple[str, ...]]:
    """Split TSV dump output into a list of field tuples (one per line)."""
    return [tuple(line.split("\t")) for line in output.strip().split("\n") if line]


# ---------------------------------------------------------------------------
# Module-scoped fixtures using real static test data
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def a1_rep1_1kb() -> Path:
    """Read-only path to assembly_1_rep1_1kb.metacool."""
    p = DATA / "assembly_1_rep1_1kb.metacool"
    assert p.exists(), f"Missing test data: {p}"
    return p


@pytest.fixture(scope="module")
def a1_rep2_1kb() -> Path:
    p = DATA / "assembly_1_rep2_1kb.metacool"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def a1_rep1_DpnII() -> Path:
    p = DATA / "assembly_1_rep1_DpnII.metacool"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def a1_rep2_DpnII() -> Path:
    p = DATA / "assembly_1_rep2_DpnII.metacool"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def a1_1kb_mmetacool() -> Path:
    p = DATA / "assembly_1_1kb.mmetacool"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def a1_rep1_cool() -> Path:
    p = DATA / "assembly_1_rep1_1kb.cool"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def contigs_table_1() -> Path:
    p = DATA / "contigs_table_1.tsv"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def mag_table_1() -> Path:
    p = DATA / "mag_table_1.tsv"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def dpnii_bed_1() -> Path:
    p = DATA / "assembly_1_DpnII-digested.bed"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def pairs_1_rep1() -> Path:
    p = DATA / "contacts_1_rep1.pairs.gz"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def pairs_1_rep2() -> Path:
    p = DATA / "contacts_1_rep2.pairs.gz"
    assert p.exists()
    return p


@pytest.fixture(scope="module")
def metadata_yaml() -> Path:
    p = DATA / "metadata.yaml"
    assert p.exists()
    return p


# ---------------------------------------------------------------------------
# TestCreatePairs
# ---------------------------------------------------------------------------


class TestCreatePairs:
    def test_create_pairs_fixed_bins(
        self,
        runner,
        tmp_path,
        pairs_1_rep1,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "created.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--pairs",
                str(pairs_1_rep1),
                "--bins",
                "1000",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_create_pairs_produces_correct_counts(
        self,
        runner,
        tmp_path,
        pairs_1_rep1,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "created.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--pairs",
                str(pairs_1_rep1),
                "--bins",
                "1000",
                "--output",
                str(out),
            ],
        )
        info = runner.invoke(cli, ["info", str(out)])
        data = json.loads(info.output)
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["ncontigs"] == A1_REP1_1KB_NCONTIGS
        assert data["nbins"] == A1_REP1_1KB_NBINS
        assert data["nnz"] == A1_REP1_1KB_NNZ
        assert data["sum"] == A1_REP1_1KB_SUM
        assert data["bin-size"] == A1_REP1_1KB_BINSIZE

    def test_create_pairs_bed_bins(
        self,
        runner,
        tmp_path,
        pairs_1_rep1,
        contigs_table_1,
        mag_table_1,
        dpnii_bed_1,
    ):
        out = tmp_path / "created_dpnii.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--pairs",
                str(pairs_1_rep1),
                "--bins",
                str(dpnii_bed_1),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        info = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert info["bin-type"] == "variable"

    def test_create_pairs_format_is_metacooler(
        self,
        runner,
        tmp_path,
        pairs_1_rep1,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "created.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--pairs",
                str(pairs_1_rep1),
                "--bins",
                "1000",
                "--output",
                str(out),
            ],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["format"] == METACOOLER_FORMAT
        assert data["format-version"] == METACOOLER_FORMAT_VERSION

    def test_create_pairs_no_overwrite_without_force(
        self,
        runner,
        tmp_path,
        pairs_1_rep1,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "created.metacool"
        args = [
            "create",
            "pairs",
            "--contigs",
            str(contigs_table_1),
            "--genomes",
            str(mag_table_1),
            "--pairs",
            str(pairs_1_rep1),
            "--bins",
            "1000",
            "--output",
            str(out),
        ]
        runner.invoke(cli, args)
        result2 = runner.invoke(cli, args)
        assert result2.exit_code != 0


# ---------------------------------------------------------------------------
# TestCreateCool
# ---------------------------------------------------------------------------


class TestCreateCool:
    def test_create_from_cool(
        self,
        runner,
        tmp_path,
        a1_rep1_cool,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "from_cool.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_create_from_cool_correct_counts(
        self,
        runner,
        tmp_path,
        a1_rep1_cool,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "from_cool.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(out),
            ],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["ncontigs"] == A1_REP1_1KB_NCONTIGS
        assert data["nbins"] == A1_REP1_1KB_NBINS
        assert data["nnz"] == A1_REP1_1KB_NNZ
        assert data["sum"] == A1_REP1_1KB_SUM

    def test_create_from_cool_ordering_is_genome_then_length(
        self,
        runner,
        tmp_path,
        a1_rep1_cool,
        contigs_table_1,
        mag_table_1,
    ):
        out = tmp_path / "from_cool.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(out),
            ],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["ordering"] == "genome-then-length"


# ---------------------------------------------------------------------------
# TestInfo
# ---------------------------------------------------------------------------


class TestInfo:
    def test_info_metacool(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["info", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data["format"] == METACOOLER_FORMAT
        assert data["format-version"] == METACOOLER_FORMAT_VERSION
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["ncontigs"] == A1_REP1_1KB_NCONTIGS
        assert data["nbins"] == A1_REP1_1KB_NBINS
        assert data["nnz"] == A1_REP1_1KB_NNZ
        assert data["sum"] == A1_REP1_1KB_SUM
        assert data["bin-size"] == A1_REP1_1KB_BINSIZE
        assert data["storage-mode"] == "symmetric-upper"

    def test_info_table_flag(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["info", "--table", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        assert "Genomes" in result.output

    def test_info_via_mmetacool_uri(self, runner, a1_1kb_mmetacool):
        uri = f"{a1_1kb_mmetacool}::entries/assembly_1_rep1_1kb"
        result = runner.invoke(cli, ["info", uri])
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data["format"] == METACOOLER_FORMAT
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["nbins"] == A1_REP1_1KB_NBINS

    def test_info_mmetacool_second_entry(self, runner, a1_1kb_mmetacool):
        uri = f"{a1_1kb_mmetacool}::entries/assembly_1_rep2_1kb"
        result = runner.invoke(cli, ["info", uri])
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data["format"] == METACOOLER_FORMAT
        # rep2 has different nnz than rep1
        assert data["nbins"] == A1_REP1_1KB_NBINS
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES


# ---------------------------------------------------------------------------
# TestTree
# ---------------------------------------------------------------------------


class TestTree:
    def test_tree_metacool_sections(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["tree", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        assert "genomes/" in result.output
        assert "contigs/" in result.output
        assert "bins/" in result.output
        assert "pixels/" in result.output
        assert "indexes/" in result.output

    def test_tree_metacool_genome_columns(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["tree", str(a1_rep1_1kb)])
        assert "name" in result.output
        assert "length" in result.output
        assert "n_contigs" in result.output

    def test_tree_metacool_pixel_columns(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["tree", str(a1_rep1_1kb)])
        assert "bin1_id" in result.output
        assert "bin2_id" in result.output
        assert "count" in result.output

    def test_tree_cool_file(self, runner, a1_rep1_cool):
        result = runner.invoke(cli, ["tree", str(a1_rep1_cool)])
        assert result.exit_code == 0, result.output
        assert "chroms/" in result.output
        assert "bins/" in result.output
        assert "pixels/" in result.output


# ---------------------------------------------------------------------------
# TestAttrs
# ---------------------------------------------------------------------------


class TestAttrs:
    def test_attrs_metacool_root(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["attrs", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        assert "format" in result.output
        assert "metacooler-format" in result.output
        assert COOLER_FORMAT in result.output

    def test_attrs_metacool_has_metadata_key(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["attrs", str(a1_rep1_1kb)])
        assert "metadata" in result.output
        assert "commands" in result.output

    def test_attrs_mmetacool_format(self, runner, a1_1kb_mmetacool):
        result = runner.invoke(cli, ["attrs", str(a1_1kb_mmetacool)])
        assert result.exit_code == 0, result.output
        assert "HDF5::MMetaCooler" in result.output
        assert "format-version" in result.output

    def test_attrs_mmetacool_has_entries(self, runner, a1_1kb_mmetacool):
        result = runner.invoke(cli, ["attrs", str(a1_1kb_mmetacool)])
        assert "assembly_1_rep1_1kb" in result.output
        assert "assembly_1_rep2_1kb" in result.output


# ---------------------------------------------------------------------------
# TestDump
# ---------------------------------------------------------------------------


class TestDump:
    def test_dump_genomes_count(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "genomes", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) == A1_REP1_1KB_NGENOMES

    def test_dump_genomes_includes_unbinned(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "genomes", str(a1_rep1_1kb)])
        assert "__unbinned__" in result.output

    def test_dump_genomes_header(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "genomes", "-H", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        # header + ngenomes rows
        assert len(lines) == A1_REP1_1KB_NGENOMES + 1
        assert "name" in lines[0]

    def test_dump_contigs_count(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "contigs", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) == A1_REP1_1KB_NCONTIGS

    def test_dump_chroms_matches_contigs(self, runner, a1_rep1_1kb):
        chroms = runner.invoke(cli, ["dump", "-t", "chroms", str(a1_rep1_1kb)])
        contigs = runner.invoke(cli, ["dump", "-t", "contigs", str(a1_rep1_1kb)])
        assert chroms.exit_code == 0
        assert len(chroms.output.strip().split("\n")) == A1_REP1_1KB_NCONTIGS
        # chroms and contigs tables have same number of rows
        assert len(chroms.output.strip().split("\n")) == len(
            contigs.output.strip().split("\n")
        )

    def test_dump_bins_count(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "bins", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) == A1_REP1_1KB_NBINS

    def test_dump_pixels_count(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "pixels", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) == A1_REP1_1KB_NNZ

    def test_dump_pixels_columns(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "pixels", "-H", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        header = result.output.strip().split("\n")[0]
        assert "bin1_id" in header
        assert "bin2_id" in header
        assert "count" in header

    def test_dump_pixels_range(self, runner, a1_rep1_1kb):
        result = runner.invoke(
            cli,
            ["dump", "-t", "pixels", "-r", "contig_1_023:5000-15000", str(a1_rep1_1kb)],
        )
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) > 0
        # range-restricted should return fewer pixels than full dump
        assert len(lines) < A1_REP1_1KB_NNZ

    def test_dump_pixels_range2(self, runner, a1_rep1_1kb):
        result = runner.invoke(
            cli,
            [
                "dump",
                "-t",
                "pixels",
                "-r",
                "contig_1_023:5000-15000",
                "-r2",
                "contig_1_036:0-10000",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output

    def test_dump_indexes(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["dump", "-t", "indexes", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = [ln for ln in result.output.strip().split("\n") if ln]
        assert len(lines) > 0

    def test_dump_bins_matches_cooler(self, runner, a1_rep1_1kb):
        """metacooler dump bins should match cooler dump bins."""
        import subprocess

        mc = runner.invoke(cli, ["dump", "-t", "bins", str(a1_rep1_1kb)])
        cc = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "bins", str(a1_rep1_1kb)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        assert mc.exit_code == 0
        assert cc.returncode == 0
        assert mc.output == cc.stdout

    def test_dump_pixels_matches_cooler(self, runner, a1_rep1_1kb):
        """metacooler dump pixels should match cooler dump pixels."""
        import subprocess

        mc = runner.invoke(cli, ["dump", "-t", "pixels", str(a1_rep1_1kb)])
        cc = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "pixels", str(a1_rep1_1kb)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        assert mc.exit_code == 0
        assert cc.returncode == 0
        assert mc.output == cc.stdout

    # --- Exact row-content checks -------------------------------------------

    def test_dump_genomes_exact_rows(self, runner, a1_rep1_1kb):
        """All genome rows match the expected order and field values exactly."""
        result = runner.invoke(cli, ["dump", "-t", "genomes", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        assert len(rows) == len(EXPECTED_GENOMES_ROWS)
        for i, (actual, expected) in enumerate(
            zip(rows, EXPECTED_GENOMES_ROWS, strict=False)
        ):
            assert actual == expected, (
                f"Genome row {i} mismatch:\n  got:      {actual}\n  expected: {expected}"
            )

    def test_dump_genomes_header_row(self, runner, a1_rep1_1kb):
        """Header row fields are correct."""
        result = runner.invoke(cli, ["dump", "-t", "genomes", "-H", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        lines = result.output.strip().split("\n")
        header = tuple(lines[0].split("\t"))
        assert header[0] == "name"
        assert header[1] == "length"
        assert header[2] == "n_contigs"
        assert header[3] == "n_bins"

    def test_dump_contigs_exact_first_20_rows(self, runner, a1_rep1_1kb):
        """First 20 contig rows match expected order and values exactly."""
        result = runner.invoke(cli, ["dump", "-t", "contigs", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        for i, expected in enumerate(EXPECTED_CONTIGS_ROWS_20):
            assert rows[i] == expected, (
                f"Contig row {i} mismatch:\n  got:      {rows[i]}\n  expected: {expected}"
            )

    def test_dump_contigs_header_row(self, runner, a1_rep1_1kb):
        """Contigs header contains correct column names in correct order."""
        result = runner.invoke(cli, ["dump", "-t", "contigs", "-H", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        header = tuple(result.output.strip().split("\n")[0].split("\t"))
        assert header[0] == "name"
        assert header[1] == "length"
        assert header[2] == "genome"
        assert header[3] == "genome_id"
        assert header[4] == "gc_content"
        assert header[5] == "coverage"

    def test_dump_chroms_exact_first_20_rows(self, runner, a1_rep1_1kb):
        """First 20 chroms rows match contigs ordering (name, length)."""
        result = runner.invoke(cli, ["dump", "-t", "chroms", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        for i, expected in enumerate(EXPECTED_CHROMS_ROWS_20):
            assert rows[i] == expected, (
                f"Chrom row {i} mismatch:\n  got:      {rows[i]}\n  expected: {expected}"
            )

    def test_dump_bins_exact_first_20_rows(self, runner, a1_rep1_1kb):
        """First 20 bin rows are 1 kb tiles on contig_0_2 in the correct order."""
        result = runner.invoke(cli, ["dump", "-t", "bins", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        for i, expected in enumerate(EXPECTED_BINS_ROWS_20):
            assert rows[i] == expected, (
                f"Bin row {i} mismatch:\n  got:      {rows[i]}\n  expected: {expected}"
            )

    def test_dump_bins_first_contig_is_correct(self, runner, a1_rep1_1kb):
        """All bins for the first contig (contig_1_023, 25 kb) appear first."""
        result = runner.invoke(cli, ["dump", "-t", "bins", str(a1_rep1_1kb)])
        rows = _parse_rows(result.output)
        contig_1_023_bins = [r for r in rows if r[0] == "contig_1_023"]
        assert len(contig_1_023_bins) == 25  # 25 × 1 kb bins
        # Starts must be 0, 1000, …, 24000 in order
        starts = [int(r[1]) for r in contig_1_023_bins]
        assert starts == list(range(0, 25000, 1000))

    def test_dump_pixels_exact_first_20_rows(self, runner, a1_rep1_1kb):
        """First 20 pixel rows match expected bin IDs and counts exactly."""
        result = runner.invoke(cli, ["dump", "-t", "pixels", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        for i, expected in enumerate(EXPECTED_PIXELS_ROWS_20):
            assert rows[i] == expected, (
                f"Pixel row {i} mismatch:\n  got:      {rows[i]}\n  expected: {expected}"
            )

    def test_dump_pixels_range_exact_first_20_rows(self, runner, a1_rep1_1kb):
        """Range query contig_1_023:5000-15000 returns the correct pixels in order."""
        result = runner.invoke(
            cli,
            ["dump", "-t", "pixels", "-r", "contig_1_023:5000-15000", str(a1_rep1_1kb)],
        )
        assert result.exit_code == 0, result.output
        rows = _parse_rows(result.output)
        assert len(rows) >= len(EXPECTED_PIXELS_RANGE_ROWS_20)
        for i, expected in enumerate(EXPECTED_PIXELS_RANGE_ROWS_20):
            assert rows[i] == expected, (
                f"Range pixel row {i} mismatch:\n  got:      {rows[i]}\n  expected: {expected}"
            )

    def test_dump_pixels_range_bin_ids_within_range(self, runner, a1_rep1_1kb):
        """All pixels from range query have bin1_id in [5, 14] (contig_1_023:5000-15000)."""
        result = runner.invoke(
            cli,
            ["dump", "-t", "pixels", "-r", "contig_1_023:5000-15000", str(a1_rep1_1kb)],
        )
        rows = _parse_rows(result.output)
        # Compute the expected bin range from the file's own extent so the test
        # is independent of where contig_1_023 lands in the contig ordering.
        from metacooler.api import MetaCooler

        lo, hi = MetaCooler(a1_rep1_1kb).extent("contig_1_023:5000-15000")
        for r in rows:
            bin1, bin2 = int(r[0]), int(r[1])
            assert lo <= bin1 < hi or lo <= bin2 < hi


# ---------------------------------------------------------------------------
# TestStack
# ---------------------------------------------------------------------------


class TestStack:
    def test_stack_two_metacools(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "stacked.mmetacool"
        result = runner.invoke(
            cli,
            ["stack", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert "Created" in result.output
        assert "2 entries" in result.output

    def test_stack_entry_names_from_filenames(
        self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb
    ):
        out = tmp_path / "stacked.mmetacool"
        runner.invoke(
            cli,
            ["stack", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        entries = list_entries(out)
        assert "assembly_1_rep1_1kb" in entries
        assert "assembly_1_rep2_1kb" in entries

    def test_stack_entries_are_accessible(
        self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb
    ):
        out = tmp_path / "stacked.mmetacool"
        runner.invoke(
            cli,
            ["stack", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        ls_result = runner.invoke(cli, ["ls", str(out)])
        assert "assembly_1_rep1_1kb" in ls_result.output
        assert "assembly_1_rep2_1kb" in ls_result.output

    def test_stack_multi_file(
        self,
        runner,
        tmp_path,
        a1_rep1_1kb,
        a1_rep2_1kb,
        a1_rep1_DpnII,
        a1_rep2_DpnII,
    ):
        out = tmp_path / "big_stack.mmetacool"
        result = runner.invoke(
            cli,
            [
                "stack",
                "--output",
                str(out),
                str(a1_rep1_1kb),
                str(a1_rep2_1kb),
                str(a1_rep1_DpnII),
                str(a1_rep2_DpnII),
            ],
        )
        assert result.exit_code == 0, result.output
        entries = list_entries(out)
        assert len(entries) == 4

    def test_stack_entry_info_via_uri(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "stacked.mmetacool"
        runner.invoke(
            cli,
            ["stack", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        uri = f"{out}::entries/assembly_1_rep1_1kb"
        info = runner.invoke(cli, ["info", uri])
        assert info.exit_code == 0, info.output
        data = json.loads(info.output)
        assert data["nbins"] == A1_REP1_1KB_NBINS
        assert data["nnz"] == A1_REP1_1KB_NNZ

    def test_stack_entry_names(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "stacked_names.mmetacool"
        runner.invoke(
            cli,
            [
                "stack",
                "--output",
                str(out),
                str(a1_rep1_1kb),
                str(a1_rep2_1kb),
                "--names",
                "rep1,rep2",
            ],
        )
        ls = runner.invoke(cli, ["ls", str(out)])
        assert ls.exit_code == 0, ls.output
        assert "rep1" in ls.output
        assert "rep2" in ls.output

        out = tmp_path / "stacked_names.mmetacool"
        stacked = runner.invoke(
            cli,
            [
                "stack",
                "--output",
                str(out),
                str(a1_rep1_1kb),
                str(a1_rep2_1kb),
                "--names",
                "rep1,rep2,repX",
            ],
        )
        assert stacked.exit_code == 1, stacked.output


# ---------------------------------------------------------------------------
# TestMerge
# ---------------------------------------------------------------------------


class TestMerge:
    def test_merge_two_replicates(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "merged.metacool"
        result = runner.invoke(
            cli,
            ["merge", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_merge_sums_contacts(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "merged.metacool"
        runner.invoke(
            cli,
            ["merge", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["sum"] == MERGED_SUM

    def test_merge_preserves_bins_and_contigs(
        self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb
    ):
        out = tmp_path / "merged.metacool"
        runner.invoke(
            cli,
            ["merge", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["ncontigs"] == A1_REP1_1KB_NCONTIGS
        assert data["nbins"] == A1_REP1_1KB_NBINS

    def test_merge_nnz_equals_union(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        out = tmp_path / "merged.metacool"
        runner.invoke(
            cli,
            ["merge", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["nnz"] == MERGED_NNZ

    def test_merge_is_valid_metacool(self, runner, tmp_path, a1_rep1_1kb, a1_rep2_1kb):
        from metacooler.api import MetaCooler

        out = tmp_path / "merged.metacool"
        runner.invoke(
            cli,
            ["merge", "--output", str(out), str(a1_rep1_1kb), str(a1_rep2_1kb)],
        )
        mc = MetaCooler(str(out))
        assert mc.info["ngenomes"] == A1_REP1_1KB_NGENOMES


# ---------------------------------------------------------------------------
# TestMetadata
# ---------------------------------------------------------------------------


class TestMetadata:
    @pytest.fixture
    def editable_metacool(self, tmp_path, a1_rep1_1kb) -> Path:
        """A writable copy of the metacool file."""
        out = tmp_path / "editable.metacool"
        shutil.copy(a1_rep1_1kb, out)
        return out

    def test_metadata_list(self, runner, editable_metacool):
        result = runner.invoke(cli, ["metadata", "list", str(editable_metacool)])
        assert result.exit_code == 0, result.output
        # output is JSON
        data = json.loads(result.output)
        assert "commands" in data

    def test_metadata_add(self, runner, editable_metacool):
        result = runner.invoke(
            cli,
            [
                "metadata",
                "add",
                "--key",
                "testkey",
                "--value",
                "testvalue",
                str(editable_metacool),
            ],
        )
        assert result.exit_code == 0, result.output
        listed = json.loads(
            runner.invoke(cli, ["metadata", "list", str(editable_metacool)]).output
        )
        assert listed["testkey"] == "testvalue"

    def test_metadata_get(self, runner, editable_metacool):
        runner.invoke(
            cli,
            [
                "metadata",
                "add",
                "--key",
                "mykey",
                "--value",
                "myvalue",
                str(editable_metacool),
            ],
        )
        result = runner.invoke(
            cli,
            ["metadata", "get", "--key", "mykey", str(editable_metacool)],
        )
        assert result.exit_code == 0, result.output
        assert json.loads(result.output) == "myvalue"

    def test_metadata_remove(self, runner, editable_metacool):
        runner.invoke(
            cli,
            [
                "metadata",
                "add",
                "--key",
                "removeme",
                "--value",
                "x",
                str(editable_metacool),
            ],
        )
        result = runner.invoke(
            cli,
            [
                "metadata",
                "remove",
                "--key",
                "removeme",
                str(editable_metacool),
            ],
        )
        assert result.exit_code == 0, result.output
        listed = json.loads(
            runner.invoke(cli, ["metadata", "list", str(editable_metacool)]).output
        )
        assert "removeme" not in listed

    def test_metadata_load_yaml(self, runner, editable_metacool, metadata_yaml):
        result = runner.invoke(
            cli, ["metadata", "load", str(metadata_yaml), str(editable_metacool)]
        )
        assert result.exit_code == 0, result.output
        listed = json.loads(
            runner.invoke(cli, ["metadata", "list", str(editable_metacool)]).output
        )
        assert "sample" in listed
        assert listed["sample"]["name"] == "gut_microbiome_T0"
        assert "sequencing" in listed
        assert "assembly" in listed

    def test_metadata_load_yaml_preserves_existing(
        self, runner, editable_metacool, metadata_yaml
    ):
        runner.invoke(
            cli,
            [
                "metadata",
                "add",
                "--key",
                "keep",
                "--value",
                "yes",
                str(editable_metacool),
            ],
        )
        runner.invoke(
            cli, ["metadata", "load", str(metadata_yaml), str(editable_metacool)]
        )
        listed = json.loads(
            runner.invoke(cli, ["metadata", "list", str(editable_metacool)]).output
        )
        assert "keep" in listed
        assert "sample" in listed

    def test_metadata_get_nonexistent_key_fails(self, runner, editable_metacool):
        result = runner.invoke(
            cli,
            [
                "metadata",
                "get",
                "--key",
                "nonexistent_key_xyz",
                str(editable_metacool),
            ],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# TestLs
# ---------------------------------------------------------------------------


class TestLs:
    def test_ls_metacool(self, runner, a1_rep1_1kb):
        result = runner.invoke(cli, ["ls", str(a1_rep1_1kb)])
        assert result.exit_code == 0, result.output
        assert "::" in result.output

    def test_ls_mmetacool(self, runner, a1_1kb_mmetacool):
        result = runner.invoke(cli, ["ls", str(a1_1kb_mmetacool)])
        assert result.exit_code == 0, result.output
        assert "assembly_1_rep1_1kb" in result.output
        assert "assembly_1_rep2_1kb" in result.output

    def test_ls_mmetacool_each_entry_is_uri(self, runner, a1_1kb_mmetacool):
        result = runner.invoke(cli, ["ls", str(a1_1kb_mmetacool)])
        for line in result.output.strip().split("\n"):
            assert "::" in line

    def test_ls_cool_file(self, runner, a1_rep1_cool):
        result = runner.invoke(cli, ["ls", str(a1_rep1_cool)])
        assert result.exit_code == 0, result.output
        assert "::" in result.output

    def test_ls_mcool_file(self, runner):
        mcool = DATA / "assembly_1_rep1_1kb.mcool"
        result = runner.invoke(cli, ["ls", str(mcool)])
        assert result.exit_code == 0, result.output
        assert "resolutions" in result.output


# ---------------------------------------------------------------------------
# TestCpMvRm
# ---------------------------------------------------------------------------


class TestCpMvRm:
    @pytest.fixture
    def editable_mmetacool(self, tmp_path, a1_1kb_mmetacool) -> Path:
        out = tmp_path / "editable.mmetacool"
        shutil.copy(a1_1kb_mmetacool, out)
        return out

    def test_cp_within_file(self, runner, editable_mmetacool):
        src = f"{editable_mmetacool}::entries/assembly_1_rep1_1kb"
        dst = f"{editable_mmetacool}::entries/copied"
        result = runner.invoke(cli, ["cp", src, dst])
        assert result.exit_code == 0, result.output
        entries = list_entries(editable_mmetacool)
        assert "copied" in entries
        # original still exists
        assert "assembly_1_rep1_1kb" in entries

    def test_cp_cross_file(self, runner, tmp_path, a1_1kb_mmetacool, a1_rep1_1kb):
        from metacooler.io.mmetacool import create_mmetacool

        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        shutil.copy(a1_1kb_mmetacool, src_file)
        create_mmetacool(dst_file, {})

        src = f"{src_file}::entries/assembly_1_rep1_1kb"
        dst = f"{dst_file}::entries/copied_entry"
        result = runner.invoke(cli, ["cp", src, dst])
        assert result.exit_code == 0, result.output
        entries = list_entries(dst_file)
        assert "copied_entry" in entries

    def test_mv_renames_within_file(self, runner, editable_mmetacool):
        src = f"{editable_mmetacool}::entries/assembly_1_rep1_1kb"
        dst = f"{editable_mmetacool}::entries/renamed"
        result = runner.invoke(cli, ["mv", src, dst])
        assert result.exit_code == 0, result.output
        entries = list_entries(editable_mmetacool)
        assert "renamed" in entries
        assert "assembly_1_rep1_1kb" not in entries

    def test_mv_cross_file_fails(self, runner, tmp_path, a1_1kb_mmetacool):
        from metacooler.io.mmetacool import create_mmetacool

        src_file = tmp_path / "src.mmetacool"
        dst_file = tmp_path / "dst.mmetacool"
        shutil.copy(a1_1kb_mmetacool, src_file)
        create_mmetacool(dst_file, {})

        src = f"{src_file}::entries/assembly_1_rep1_1kb"
        dst = f"{dst_file}::entries/moved"
        result = runner.invoke(cli, ["mv", src, dst])
        assert result.exit_code != 0

    def test_rm_removes_entry(self, runner, editable_mmetacool):
        src = f"{editable_mmetacool}::entries/assembly_1_rep1_1kb"
        result = runner.invoke(cli, ["rm", "-f", src])
        assert result.exit_code == 0, result.output
        entries = list_entries(editable_mmetacool)
        assert "assembly_1_rep1_1kb" not in entries
        # other entry still present
        assert "assembly_1_rep2_1kb" in entries

    def test_rm_without_force_prompts(self, runner, editable_mmetacool):
        """Without -f flag rm should either prompt or fail non-interactively."""
        src = f"{editable_mmetacool}::entries/assembly_1_rep1_1kb"
        # non-interactive runner will get a no-input response
        runner.invoke(cli, ["rm", src], input="n\n")
        # should not have removed the entry
        entries = list_entries(editable_mmetacool)
        assert "assembly_1_rep1_1kb" in entries


# ---------------------------------------------------------------------------
# TestFlatten
# ---------------------------------------------------------------------------


class TestFlatten:
    @pytest.fixture
    def round_trip_cool(
        self, tmp_path, runner, a1_rep1_cool, contigs_table_1, mag_table_1
    ) -> Path:
        """Build metacool from cool, flatten back to cool; return the flat path."""
        metacool = tmp_path / "from_cool.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(metacool),
                "--contig-order",
                "cool",
            ],
        )
        flat = tmp_path / "flat.cool"
        result = runner.invoke(cli, ["flatten", str(metacool), "-o", str(flat)])
        assert result.exit_code == 0, result.output
        return flat

    def test_flatten_produces_cool(self, round_trip_cool):
        assert round_trip_cool.exists()

    def test_flatten_round_trip_bins(
        self, runner, tmp_path, a1_rep1_cool, contigs_table_1, mag_table_1
    ):
        """Bins table survives metacool round trip unchanged."""
        import subprocess

        metacool = tmp_path / "from_cool.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(metacool),
                "--contig-order",
                "cool",
            ],
        )
        flat = tmp_path / "flat.cool"
        runner.invoke(cli, ["flatten", str(metacool), "-o", str(flat)])

        orig = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "bins", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        round_trip = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "bins", str(flat)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        assert orig == round_trip

    def test_flatten_round_trip_pixels(
        self, runner, tmp_path, a1_rep1_cool, contigs_table_1, mag_table_1
    ):
        """Pixel table survives metacool round trip unchanged."""
        import subprocess

        metacool = tmp_path / "from_cool.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--contigs",
                str(contigs_table_1),
                "--genomes",
                str(mag_table_1),
                "--cool",
                str(a1_rep1_cool),
                "--output",
                str(metacool),
                "--contig-order",
                "cool",
            ],
        )
        flat = tmp_path / "flat.cool"
        runner.invoke(cli, ["flatten", str(metacool), "-o", str(flat)])

        orig = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "pixels", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        round_trip = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "pixels", str(flat)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        assert orig == round_trip

    def test_flatten_cool_counts_match(self, round_trip_cool):
        c = cooler.Cooler(str(round_trip_cool))
        assert c.info["nbins"] == A1_REP1_1KB_NBINS
        assert c.info["nnz"] == A1_REP1_1KB_NNZ
        assert c.info["sum"] == A1_REP1_1KB_SUM

    # --- Exact row-content checks after round trip ---------------------------

    def test_round_trip_chroms_exact_first_20_rows(self, round_trip_cool, a1_rep1_cool):
        """First 20 chrom rows of flattened cool match original cool exactly."""
        import subprocess

        orig = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "chroms", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        rt = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "chroms", str(round_trip_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        assert orig.returncode == 0 and rt.returncode == 0
        orig_rows = _parse_rows(orig.stdout)[:20]
        rt_rows = _parse_rows(rt.stdout)[:20]
        for i, (expected, got) in enumerate(zip(orig_rows, rt_rows, strict=False)):
            assert got == expected, (
                f"Round-trip chrom row {i} mismatch:\n  got:      {got}\n  expected: {expected}"
            )

    def test_round_trip_bins_exact_first_20_rows(self, round_trip_cool, a1_rep1_cool):
        """First 20 bin rows of flattened cool match original cool exactly."""
        import subprocess

        orig = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "bins", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        rt = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "bins", str(round_trip_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        assert orig.returncode == 0 and rt.returncode == 0
        orig_rows = _parse_rows(orig.stdout)[:20]
        rt_rows = _parse_rows(rt.stdout)[:20]
        for i, (expected, got) in enumerate(zip(orig_rows, rt_rows, strict=False)):
            assert got == expected, (
                f"Round-trip bin row {i} mismatch:\n  got:      {got}\n  expected: {expected}"
            )

    def test_round_trip_pixels_exact_first_20_rows(self, round_trip_cool, a1_rep1_cool):
        """First 20 pixel rows of flattened cool match original cool exactly."""
        import subprocess

        orig = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "pixels", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        rt = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "pixels", str(round_trip_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        )
        assert orig.returncode == 0 and rt.returncode == 0
        orig_rows = _parse_rows(orig.stdout)[:20]
        rt_rows = _parse_rows(rt.stdout)[:20]
        for i, (expected, got) in enumerate(zip(orig_rows, rt_rows, strict=False)):
            assert got == expected, (
                f"Round-trip pixel row {i} mismatch:\n  got:      {got}\n  expected: {expected}"
            )

    def test_round_trip_bins_contig_order_matches_original(
        self, round_trip_cool, a1_rep1_cool
    ):
        """Contig order in flattened cool matches original cool exactly."""
        import subprocess

        orig_chroms = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "chroms", str(a1_rep1_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        rt_chroms = subprocess.run(
            ["uv", "run", "cooler", "dump", "-t", "chroms", str(round_trip_cool)],
            capture_output=True,
            text=True,
            cwd=str(DATA),
        ).stdout
        orig_names = [r.split("\t")[0] for r in orig_chroms.strip().split("\n") if r]
        rt_names = [r.split("\t")[0] for r in rt_chroms.strip().split("\n") if r]
        assert orig_names == rt_names

    def test_round_trip_pixel_sum_conserved(self, round_trip_cool):
        """Total contact count is unchanged after round trip."""
        c = cooler.Cooler(str(round_trip_cool))
        assert c.info["sum"] == A1_REP1_1KB_SUM


# ---------------------------------------------------------------------------
# TestExtract
# ---------------------------------------------------------------------------


class TestExtract:
    def test_extract_single_genome(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "MAG_01.cool"
        result = runner.invoke(
            cli,
            ["extract", "--output", str(out), "--genomes", "MAG_01", str(a1_rep1_1kb)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_extract_single_genome_correct_counts(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "MAG_01.cool"
        runner.invoke(
            cli,
            ["extract", "--output", str(out), "--genomes", "MAG_01", str(a1_rep1_1kb)],
        )
        c = cooler.Cooler(str(out))
        assert c.info["nbins"] == 60
        assert c.info["nnz"] == 1503

    def test_extract_multiple_genomes(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "MAG_01-02.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--genomes",
                "MAG_01,MAG_02",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        c = cooler.Cooler(str(out))
        # MAG_01 has 3 contigs, MAG_02 has 4 contigs
        assert c.info["nchroms"] == 7

    def test_extract_genome_split(self, runner, tmp_path, a1_rep1_1kb):
        out_prefix = tmp_path / "split"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out_prefix),
                "--split",
                "--genomes",
                "MAG_01,MAG_02",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "split_MAG_01.cool").exists()
        assert (tmp_path / "split_MAG_02.cool").exists()

    def test_extract_genome_split_correct_counts(self, runner, tmp_path, a1_rep1_1kb):
        out_prefix = tmp_path / "split"
        runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out_prefix),
                "--split",
                "--genomes",
                "MAG_01,MAG_02",
                str(a1_rep1_1kb),
            ],
        )
        mag01 = cooler.Cooler(str(tmp_path / "split_MAG_01.cool"))
        mag02 = cooler.Cooler(str(tmp_path / "split_MAG_02.cool"))
        assert mag01.info["nbins"] == 60
        assert mag01.info["nnz"] == 1503
        assert mag02.info["nbins"] == 124

    def test_extract_single_contig(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "contig_1_023.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--contigs",
                "contig_1_023",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert c.info["nchroms"] == 1
        assert c.info["nbins"] == 25
        assert c.info["nnz"] == 243

    def test_extract_multiple_contigs(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "two_contigs.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--contigs",
                "contig_1_023,contig_1_036",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert c.info["nchroms"] == 2

    def test_extract_contig_split(self, runner, tmp_path, a1_rep1_1kb):
        """--split with --contigs (no --genomes) writes to a single output file."""
        out = tmp_path / "contigs.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--split",
                "--contigs",
                "contig_1_023,contig_1_036",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        # split without --genomes falls back to single file
        assert out.exists()
        c = cooler.Cooler(str(out))
        assert c.info["nchroms"] == 2

    def test_extract_unknown_genome_fails(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "bad.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--genomes",
                "NONEXISTENT_GENOME",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code != 0

    def test_extract_unbinned_genome(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "unbinned.cool"
        result = runner.invoke(
            cli,
            [
                "extract",
                "--output",
                str(out),
                "--genomes",
                "__unbinned__",
                str(a1_rep1_1kb),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


# ---------------------------------------------------------------------------
# TestBinnify
# ---------------------------------------------------------------------------


class TestBinnify:
    def test_binnify_fixed_to_larger(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "rebinned_5kb.metacool"
        result = runner.invoke(
            cli,
            ["binnify", "--binsize", "5000", "--output", str(out), str(a1_rep1_1kb)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_binnify_reduces_bins(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "rebinned_5kb.metacool"
        runner.invoke(
            cli,
            ["binnify", "--binsize", "5000", "--output", str(out), str(a1_rep1_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["nbins"] < A1_REP1_1KB_NBINS
        assert data["nbins"] == 147

    def test_binnify_preserves_sum(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "rebinned_5kb.metacool"
        runner.invoke(
            cli,
            ["binnify", "--binsize", "5000", "--output", str(out), str(a1_rep1_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["sum"] == A1_REP1_1KB_SUM

    def test_binnify_preserves_genome_and_contig_count(
        self, runner, tmp_path, a1_rep1_1kb
    ):
        out = tmp_path / "rebinned_5kb.metacool"
        runner.invoke(
            cli,
            ["binnify", "--binsize", "5000", "--output", str(out), str(a1_rep1_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["ngenomes"] == A1_REP1_1KB_NGENOMES
        assert data["ncontigs"] == A1_REP1_1KB_NCONTIGS

    def test_binnify_updates_binsize_attr(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "rebinned_5kb.metacool"
        runner.invoke(
            cli,
            ["binnify", "--binsize", "5000", "--output", str(out), str(a1_rep1_1kb)],
        )
        data = json.loads(runner.invoke(cli, ["info", str(out)]).output)
        assert data["bin-size"] == 5000

    def test_binnify_from_variable_bins(self, runner, tmp_path, a1_rep1_DpnII):
        out = tmp_path / "rebinned_5kb.metacool"
        result = runner.invoke(
            cli,
            [
                "binnify",
                "--binsize",
                "5000",
                "--output",
                str(out),
                str(a1_rep1_DpnII),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_binnify_no_overwrite_without_force(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "rebinned_5kb.metacool"
        args = [
            "binnify",
            "--binsize",
            "5000",
            "--output",
            str(out),
            str(a1_rep1_1kb),
        ]
        runner.invoke(cli, args)
        result2 = runner.invoke(cli, args)
        assert result2.exit_code != 0


# ---------------------------------------------------------------------------
# TestZoomify
# ---------------------------------------------------------------------------


class TestZoomify:
    def test_zoomify_creates_mmetacool(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "zoomed.mmetacool"
        result = runner.invoke(
            cli, ["zoomify", "--levels", "3", "--output", str(out), str(a1_rep1_1kb)]
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_zoomify_correct_number_of_levels(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "zoomed.mmetacool"
        runner.invoke(
            cli, ["zoomify", "--levels", "3", "--output", str(out), str(a1_rep1_1kb)]
        )
        entries = list_entries(out)
        assert len(entries) == 3

    def test_zoomify_base_resolution_included(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "zoomed.mmetacool"
        runner.invoke(
            cli, ["zoomify", "--levels", "3", "--output", str(out), str(a1_rep1_1kb)]
        )
        entries = list_entries(out)
        # base resolution should be present
        assert "1000" in entries

    def test_zoomify_coarser_levels_have_fewer_bins(
        self, runner, tmp_path, a1_rep1_1kb
    ):
        out = tmp_path / "zoomed.mmetacool"
        runner.invoke(
            cli, ["zoomify", "--levels", "3", "--output", str(out), str(a1_rep1_1kb)]
        )
        entries = list_entries(out)
        resolutions = sorted(int(e) for e in entries)
        assert resolutions[0] < resolutions[-1]

        uri_fine = f"{out}::entries/{resolutions[0]}"
        uri_coarse = f"{out}::entries/{resolutions[-1]}"
        fine = json.loads(runner.invoke(cli, ["info", uri_fine]).output)
        coarse = json.loads(runner.invoke(cli, ["info", uri_coarse]).output)
        assert fine["nbins"] > coarse["nbins"]

    def test_zoomify_preserves_sum_at_all_levels(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "zoomed.mmetacool"
        runner.invoke(
            cli, ["zoomify", "--levels", "3", "--output", str(out), str(a1_rep1_1kb)]
        )
        for entry in list_entries(out):
            uri = f"{out}::entries/{entry}"
            data = json.loads(runner.invoke(cli, ["info", uri]).output)
            assert data["sum"] == A1_REP1_1KB_SUM

    def test_zoomify_no_overwrite_without_force(self, runner, tmp_path, a1_rep1_1kb):
        out = tmp_path / "zoomed.mmetacool"
        args = [
            "zoomify",
            "--levels",
            "3",
            "--output",
            str(out),
            str(a1_rep1_1kb),
        ]
        runner.invoke(cli, args)
        result2 = runner.invoke(cli, args)
        assert result2.exit_code != 0


# ---------------------------------------------------------------------------
# TestNetwork — mirrors the cli_tests.sh network section
# ---------------------------------------------------------------------------


class TestNetwork:
    """Tests matching the three ``metacooler network`` calls in cli_tests.sh."""

    def test_network_contig_level(self, runner, tmp_path, a1_rep1_1kb):
        """metacooler network assembly_1_rep1_1kb.metacool --level contig -v -o network.tsv"""
        out = tmp_path / "network.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--level",
                "contig",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        for line in lines:
            assert len(line.split("\t")) == 3

    def test_network_genome_level(self, runner, tmp_path, a1_rep1_1kb):
        """metacooler network assembly_1_rep1_1kb.metacool --level genome -v -o network.tsv"""
        out = tmp_path / "network_genome.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--level",
                "genome",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0

    def test_network_genome_level_weight_count(self, runner, tmp_path, a1_rep1_1kb):
        """metacooler network assembly_1_rep1_1kb.metacool --weight count --level genome -v -o network.tsv"""
        out = tmp_path / "network_genome_count.tsv"
        result = runner.invoke(
            cli,
            [
                "network",
                str(a1_rep1_1kb),
                "--weight",
                "count",
                "--level",
                "genome",
                "-v",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        lines = [ln for ln in out.read_text().strip().split("\n") if ln]
        assert len(lines) > 0
        # Values must be positive integers
        for line in lines:
            parts = line.split("\t")
            assert int(parts[2]) > 0

    def test_network_contig_level_no_output(self, runner, a1_rep1_1kb):
        """Without -o the output goes to stdout."""
        result = runner.invoke(
            cli,
            ["network", str(a1_rep1_1kb), "--level", "contig"],
        )
        assert result.exit_code == 0, result.output
        # Filter out log lines (start with '[') to get only TSV rows
        tsv_lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        assert len(tsv_lines) > 0
        for line in tsv_lines:
            assert len(line.split("\t")) == 3

    def test_network_contig_node_count_leq_ncontigs(self, runner, a1_rep1_1kb):
        """Number of distinct nodes ≤ number of contigs."""
        result = runner.invoke(cli, ["network", str(a1_rep1_1kb), "--level", "contig"])
        assert result.exit_code == 0, result.output
        lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        nodes = {ln.split("\t")[0] for ln in lines} | {
            ln.split("\t")[1] for ln in lines
        }
        assert len(nodes) <= A1_REP1_1KB_NCONTIGS

    def test_network_genome_node_count_leq_ngenomes(self, runner, a1_rep1_1kb):
        """Number of distinct nodes ≤ number of genomes."""
        result = runner.invoke(cli, ["network", str(a1_rep1_1kb), "--level", "genome"])
        assert result.exit_code == 0, result.output
        lines = [
            ln
            for ln in result.output.strip().split("\n")
            if ln and not ln.startswith("[")
        ]
        nodes = {ln.split("\t")[0] for ln in lines} | {
            ln.split("\t")[1] for ln in lines
        }
        assert len(nodes) <= A1_REP1_1KB_NGENOMES
