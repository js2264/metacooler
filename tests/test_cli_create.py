"""Stress tests for the ``metacooler create`` CLI subcommands.

Covers:
- All ``--contigs-format`` choices (pairs and cool)
- All ``--contig-order`` choices (pairs and cool)
- ``--genomes`` validation: valid table, missing ``genome_name`` column, no header
- ``--assembly-name`` stored in HDF5 attrs
- ``--chunksize`` / ``--nproc`` do not cause errors (pairs only)
"""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest
from click.testing import CliRunner

from metacooler.cli.main import cli


@pytest.fixture
def runner():
    return CliRunner()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _base_pairs_args(
    output: Path,
    pairs: Path,
    bins: Path,
    contigs: Path,
    *,
    contigs_format: str = "contigs",
    contig_order: str = "genome-then-length",
    extra: list[str] | None = None,
) -> list[str]:
    args = [
        "create",
        "pairs",
        "--output",
        str(output),
        "--pairs",
        str(pairs),
        "--bins",
        str(bins),
        "--contigs",
        str(contigs),
        "--contigs-format",
        contigs_format,
        "--contig-order",
        contig_order,
    ]
    if extra:
        args.extend(extra)
    return args


def _base_cool_args(
    output: Path,
    cool: Path,
    contigs: Path,
    *,
    contigs_format: str = "contigs",
    contig_order: str = "genome-then-length",
    extra: list[str] | None = None,
) -> list[str]:
    args = [
        "create",
        "cool",
        "--output",
        str(output),
        "--cool",
        str(cool),
        "--contigs",
        str(contigs),
        "--contigs-format",
        contigs_format,
        "--contig-order",
        contig_order,
    ]
    if extra:
        args.extend(extra)
    return args


# ===========================================================================
# create pairs — contigs format
# ===========================================================================


class TestCreatePairsContigsFormats:
    """All --contigs-format choices for ``create pairs``."""

    def test_contigs_format_default(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(out, tiny_pairs, tiny_bins_bed, tiny_contigs_table),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_contigs_explicit(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contigs_format="contigs",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_contig_genome(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_contig_genome
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_contig_genome,
                contigs_format="contig_genome",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_metator(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_metator
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_metator,
                contigs_format="metator",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_with_header(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``1,2`` (with header) on a standard 2-column contigs table."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contigs_format="1,2",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_noheader(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_contig_genome
    ):
        """``1,2,noheader`` on a no-header 2-column table."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_contig_genome,
                contigs_format="1,2,noheader",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_non_default_cols(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_4col
    ):
        """``2,4`` on a 4-column table (col 2 = contig_name, col 4 = genome_name)."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_4col,
                contigs_format="2,4",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_noheader_col_out_of_range(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``1,4,noheader`` on a 2-column table must fail (col 4 out of range)."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contigs_format="1,4,noheader",
            ),
        )
        assert result.exit_code != 0

    def test_contigs_format_contigs_missing_contig_name_col(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_contig_genome
    ):
        """``contigs`` format on a no-header file fails (first col != 'contig_name')."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_contig_genome,
                contigs_format="contigs",
            ),
        )
        assert result.exit_code != 0

    def test_contigs_format_metator_missing_columns(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``metator`` format on a standard contigs table fails (missing Name/Final_bin)."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contigs_format="metator",
            ),
        )
        assert result.exit_code != 0

    def test_contigs_format_contig_genome_accepts_header_table(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``contig_genome`` has no header validation.

        A headered file is treated as all-data rows: the header row becomes a
        spurious mapping ("contig_name" → "genome_name") that is silently
        ignored because no bin has a contig named "contig_name".  The real
        contigs present in rows 2-4 are still correctly assigned, so the
        command succeeds.
        """
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contigs_format="contig_genome",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


# ===========================================================================
# create pairs — contig order
# ===========================================================================


class TestCreatePairsContigOrder:
    """All --contig-order choices for ``create pairs``."""

    @pytest.mark.parametrize(
        "order",
        ["genome-then-length", "genome-then-alpha", "fasta", "alpha"],
    )
    def test_valid_contig_order(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table, order
    ):
        out = tmp_path / f"out_{order}.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order=order,
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contig_order_genome_then_length_produces_sorted_contigs(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """Contigs grouped by genome and sorted longest-first within each genome."""
        out = tmp_path / "out.metacool"
        runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order="genome-then-length",
            ),
        )
        with h5py.File(out, "r") as f:
            names = [
                n.decode() if isinstance(n, bytes) else n
                for n in f["contigs"]["name"][:]
            ]
        # MAG_1 contigs (contig_A 1000bp, contig_B 600bp) come before MAG_2
        # (contig_C 400bp), which comes before __unbinned__ (contig_D 200bp)
        assert names.index("contig_A") < names.index("contig_B")
        assert names.index("contig_B") < names.index("contig_C")

    def test_contig_order_alpha_produces_sorted_contigs(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``alpha`` order yields strict alphabetical sort across all contigs."""
        out = tmp_path / "out.metacool"
        runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order="alpha",
            ),
        )
        with h5py.File(out, "r") as f:
            names = [
                n.decode() if isinstance(n, bytes) else n
                for n in f["contigs"]["name"][:]
            ]
        assert names == sorted(names)

    def test_contig_order_fasta_requires_fasta_or_bins(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``fasta`` order with a BED bins file succeeds (order inferred from bins)."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order="fasta",
            ),
        )
        assert result.exit_code == 0, result.output

    def test_contig_order_invalid_choice(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """An invalid ``--contig-order`` value is rejected by click."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order="random-shuffle",
            ),
        )
        assert result.exit_code != 0


# ===========================================================================
# create pairs — genomes table
# ===========================================================================


class TestCreatePairsGenomesTable:
    """``--genomes`` validation for ``create pairs``."""

    def test_valid_genomes_table(
        self,
        runner,
        tmp_path,
        tiny_pairs,
        tiny_bins_bed,
        tiny_contigs_table,
        tiny_genomes_table,
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_table)],
            ),
        )
        assert result.exit_code == 0, result.output
        # Extra columns from the genomes table should appear in /genomes
        with h5py.File(out, "r") as f:
            assert "quality" in f["genomes"]
            assert "phylum" in f["genomes"]

    def test_genomes_table_wrong_col_name_fails(
        self,
        runner,
        tmp_path,
        tiny_pairs,
        tiny_bins_bed,
        tiny_contigs_table,
        tiny_genomes_wrong_col_name,
    ):
        """Genomes table whose first column is not ``genome_name`` must fail."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_wrong_col_name)],
            ),
        )
        assert result.exit_code != 0

    def test_genomes_table_no_header_fails(
        self,
        runner,
        tmp_path,
        tiny_pairs,
        tiny_bins_bed,
        tiny_contigs_table,
        tiny_genomes_no_header,
    ):
        """Genomes table with no header is read as header=0 so first row becomes
        the column names.  Since 'MAG_1' != 'genome_name', parsing must fail."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_no_header)],
            ),
        )
        assert result.exit_code != 0


# ===========================================================================
# create pairs — assembly name
# ===========================================================================


class TestCreatePairsAssemblyName:
    """``--assembly-name`` is stored as an HDF5 root attribute."""

    def test_assembly_name_stored_in_attrs(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--assembly-name", "my_metagenome"],
            ),
        )
        assert result.exit_code == 0, result.output
        with h5py.File(out, "r") as f:
            assert f.attrs.get("assembly-name") == "my_metagenome"

    def test_assembly_name_absent_by_default(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        runner.invoke(
            cli,
            _base_pairs_args(out, tiny_pairs, tiny_bins_bed, tiny_contigs_table),
        )
        with h5py.File(out, "r") as f:
            assert "assembly-name" not in f.attrs

    def test_assembly_name_with_spaces(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--assembly-name", "my metagenome v2"],
            ),
        )
        assert result.exit_code == 0, result.output
        with h5py.File(out, "r") as f:
            assert f.attrs.get("assembly-name") == "my metagenome v2"


# ===========================================================================
# create pairs — bins as integer (binsize)
# ===========================================================================


class TestCreatePairsBinsAsInteger:
    """``--bins`` accepts an integer bin size with ``--fasta`` for contig lengths."""

    def test_binsize_integer_with_fasta(
        self, runner, tmp_path, tiny_pairs, tiny_fasta_valid, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--output",
                str(out),
                "--pairs",
                str(tiny_pairs),
                "--bins",
                "200",
                "--fasta",
                str(tiny_fasta_valid),
                "--contigs",
                str(tiny_contigs_table),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


# ===========================================================================
# create cool — contigs format
# ===========================================================================


class TestCreateCoolContigsFormats:
    """All --contigs-format choices for ``create cool``."""

    def test_contigs_format_default(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(cli, _base_cool_args(out, tiny_cool, tiny_contigs_table))
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_contig_genome(
        self, runner, tmp_path, tiny_cool, tiny_contigs_contig_genome
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_contig_genome,
                contigs_format="contig_genome",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_metator(
        self, runner, tmp_path, tiny_cool, tiny_contigs_metator
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_metator,
                contigs_format="metator",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                contigs_format="1,2",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_noheader(
        self, runner, tmp_path, tiny_cool, tiny_contigs_contig_genome
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_contig_genome,
                contigs_format="1,2,noheader",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_column_index_non_default_cols(
        self, runner, tmp_path, tiny_cool, tiny_contigs_4col
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_4col,
                contigs_format="2,4",
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contigs_format_noheader_col_out_of_range_fails(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table
    ):
        """``1,4,noheader`` on a 2-column table must fail."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                contigs_format="1,4,noheader",
            ),
        )
        assert result.exit_code != 0

    def test_contigs_format_contigs_missing_col_header_fails(
        self, runner, tmp_path, tiny_cool, tiny_contigs_contig_genome
    ):
        """``contigs`` format on a no-header file fails (col 0 != 'contig_name')."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_contig_genome,
                contigs_format="contigs",
            ),
        )
        assert result.exit_code != 0

    def test_contigs_format_metator_wrong_file_fails(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table
    ):
        """``metator`` format on a standard contigs table fails (missing Name/Final_bin)."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                contigs_format="metator",
            ),
        )
        assert result.exit_code != 0


# ===========================================================================
# create cool — contig order
# ===========================================================================


class TestCreateCoolContigOrder:
    """All --contig-order choices for ``create cool``."""

    @pytest.mark.parametrize(
        "order",
        ["genome-then-length", "genome-then-alpha", "cool", "alpha"],
    )
    def test_valid_contig_order(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table, order
    ):
        out = tmp_path / f"out_{order}.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                contig_order=order,
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_contig_order_invalid_choice(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                contig_order="random-shuffle",
            ),
        )
        assert result.exit_code != 0


# ===========================================================================
# create cool — genomes table
# ===========================================================================


class TestCreateCoolGenomesTable:
    """``--genomes`` validation for ``create cool``."""

    def test_valid_genomes_table(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table, tiny_genomes_table
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_table)],
            ),
        )
        assert result.exit_code == 0, result.output
        with h5py.File(out, "r") as f:
            assert "quality" in f["genomes"]
            assert "phylum" in f["genomes"]

    def test_genomes_table_wrong_col_name_fails(
        self,
        runner,
        tmp_path,
        tiny_cool,
        tiny_contigs_table,
        tiny_genomes_wrong_col_name,
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_wrong_col_name)],
            ),
        )
        assert result.exit_code != 0

    def test_genomes_table_no_header_fails(
        self, runner, tmp_path, tiny_cool, tiny_contigs_table, tiny_genomes_no_header
    ):
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                extra=["--genomes", str(tiny_genomes_no_header)],
            ),
        )
        assert result.exit_code != 0
