"""Stress tests for the ``metacooler create`` CLI subcommands.

Covers:
- All ``--contigs-format`` choices (pairs and cool)
- All ``--contig-order`` choices (pairs and cool)
- ``--genomes`` validation: valid table, missing ``genome_name`` column, no header
- ``--assembly-name`` stored in HDF5 attrs
- ``--chunksize`` / ``--nproc`` do not cause errors (pairs only)
- ``create pairs`` without ``--bins``: one bin per contig, as ``cooler cload
  pairs`` bins over a BED of the contigs
- Bins + pixel correctness: ``create pairs --contig-order fasta`` matches
  ``cooler cload pairs``, for fixed-width and digest (variable-width) bins,
  across every chunking and sort path
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import ClassVar, cast

import h5py
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner
from numpy.testing import assert_array_equal

from metacooler.api import MetaCooler
from metacooler.cli.main import cli

# cooler's CLI, next to the interpreter (see tests/test_cooler_compat.py).
_COOLER = shutil.which("cooler", path=str(Path(sys.executable).parent)) or str(
    Path(sys.executable).parent / "cooler"
)


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

    def test_contigs_format_metator_reads_final_cluster(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_metator
    ):
        """metaTOR 2.0's table (Final_cluster) gives the map that 1.x's
        (Final_bin) gives."""
        renamed = tmp_path / "contig_data_final.txt"
        renamed.write_text(
            tiny_contigs_metator.read_text().replace("Final_bin", "Final_cluster", 1)
        )
        maps = {}
        for name, table in (("1x", tiny_contigs_metator), ("2.0", renamed)):
            maps[name] = tmp_path / f"{name}.metacool"
            result = runner.invoke(
                cli,
                _base_pairs_args(
                    maps[name],
                    tiny_pairs,
                    tiny_bins_bed,
                    table,
                    contigs_format="metator",
                ),
            )
            assert result.exit_code == 0, result.output
        old, new = (MetaCooler(maps[name]) for name in ("1x", "2.0"))
        pd.testing.assert_frame_equal(old.genomes, new.genomes)
        pd.testing.assert_frame_equal(old.contigs, new.contigs)

    def test_contigs_format_help_names_the_metator_fallback(self):
        """The help's alias list says what each alias reads, and the column
        the metator alias reads in a table of MetaTOR 1.x."""
        from metacooler.cli._options import CONTIGS_FORMAT_HELP

        assert (
            "[yellow]metator[/yellow] (Name,Final_cluster; Final_bin if the table "
            "has no Final_cluster)"
        ) in CONTIGS_FORMAT_HELP
        assert "[yellow]dastool[/yellow] (1,2,noheader)," in CONTIGS_FORMAT_HELP

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
        """``metator`` format on a standard contigs table fails: it has no
        ``Name``, nor ``Final_cluster`` or ``Final_bin``."""
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
        spurious mapping ("contig_name" → "genome_name") that is ignored, with
        a warning, because no bin has a contig named "contig_name".  The real
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
        assert "were ignored: 'contig_name'" in " ".join(result.output.split())


# ===========================================================================
# create pairs — contig order
# ===========================================================================


class TestCreatePairsContigOrder:
    """All --contig-order choices for ``create pairs``."""

    @pytest.mark.parametrize(
        "order",
        ["genome-then-length", "genome-then-alpha", "alpha"],
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

    def test_contig_order_fasta_with_fasta(
        self,
        runner,
        tmp_path,
        tiny_pairs,
        tiny_bins_bed,
        tiny_contigs_table,
        tiny_fasta,
    ):
        """``fasta`` order works when a FASTA is supplied to define the order."""
        out = tmp_path / "out_fasta.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                contig_order="fasta",
                extra=["--fasta", str(tiny_fasta)],
            ),
        )
        assert result.exit_code == 0, result.output
        with h5py.File(out, "r") as f:
            names = [
                n.decode() if isinstance(n, bytes) else n
                for n in f["contigs"]["name"][:]
            ]
        # FASTA order, NOT alphabetical-by-accident (they coincide here, so
        # assert the recorded ordering too).
        assert names == ["contig_A", "contig_B", "contig_C", "contig_D"]
        with h5py.File(out, "r") as f:
            assert f.attrs["contig-ordering"] == "fasta"

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

    def test_contig_order_fasta_without_reference_is_refused(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """``fasta`` order is refused when nothing defines the FASTA order.

        Regression test for review item #5: with lengths coming from a BED and
        no FASTA and no ``#chromsize:`` pairs header, there is no FASTA order to
        preserve.  The old code silently fell back to *alphabetical* while still
        recording ``contig-ordering="fasta"``, so two files with different
        layouts both claimed the same ordering and ``merge`` would scramble them.
        """
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
        assert result.exit_code != 0
        assert "fasta" in result.output.lower()
        assert not out.exists()

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
            assert f.attrs.get("genome-assembly") == "my_metagenome"

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
            assert f.attrs.get("genome-assembly") == "my metagenome v2"


# ===========================================================================
# create pairs — bins as integer (binsize)
# ===========================================================================


class TestCreatePairsBinsAsInteger:
    """``--bins`` accepts an integer bin size with ``--fasta`` for contig lengths."""

    def test_binsize_integer_with_fasta(
        self, runner, tmp_path, tiny_pairs, tiny_fasta, tiny_contigs_table
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
                str(tiny_fasta),
                "--contigs",
                str(tiny_contigs_table),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


class TestCreatePairsOneBinPerContig:
    """Without ``--bins``, ``create pairs`` makes one bin per contig."""

    DATA = Path("tests/data")

    @pytest.mark.parametrize("nproc", ["1", "2"])
    def test_matches_cooler_cload_of_contigs(self, runner, tmp_path, nproc):
        """Bins and pixels equal ``cooler cload pairs`` over one bin per contig.

        The contig lengths come from the pairs header: no ``--fasta``.
        """
        sizes = pd.read_csv(
            self.DATA / "assembly_1.chrom.sizes",
            sep="\t",
            header=None,
            names=["chrom", "end"],
        )
        bed = tmp_path / "contigs.bed"
        sizes.assign(start=0)[["chrom", "start", "end"]].to_csv(
            bed, sep="\t", header=False, index=False
        )
        pairs = self.DATA / "contacts_1_rep1.pairs.gz"
        reference = tmp_path / "contigs.cool"
        # cooler's own CLI: a third-party tool, so a subprocess.
        cload = subprocess.run(
            [
                _COOLER,
                "cload",
                "pairs",
                "-c1",
                "2",
                "-p1",
                "3",
                "-c2",
                "4",
                "-p2",
                "5",
                str(bed),
                str(pairs),
                str(reference),
            ],
            capture_output=True,
            text=True,
        )
        assert cload.returncode == 0, cload.stderr

        out = tmp_path / "contigs.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--output",
                str(out),
                "--pairs",
                str(pairs),
                "--contig-order",
                "fasta",
                "--no-ps",
                "--nproc",
                nproc,
            ],
        )
        assert result.exit_code == 0, result.output

        import cooler

        ref, got = cooler.Cooler(str(reference)), cooler.Cooler(str(out))
        cols = ["chrom", "start", "end"]
        bins_ref = ref.bins()[:][cols].astype({"chrom": str})
        bins_got = got.bins()[:][cols].astype({"chrom": str})
        pd.testing.assert_frame_equal(bins_ref, bins_got)
        assert len(bins_got) == len(sizes)
        pix_ref = cast(pd.DataFrame, ref.pixels()[:])
        pix_got = cast(pd.DataFrame, got.pixels()[:])
        pd.testing.assert_frame_equal(pix_ref, pix_got)
        assert pix_got["count"].sum() > 0

    def test_with_contigs_and_fasta(
        self, runner, tmp_path, tiny_pairs, tiny_fasta, tiny_contigs_table
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
                "--fasta",
                str(tiny_fasta),
                "--contigs",
                str(tiny_contigs_table),
                "--no-ps",
            ],
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(str(out))
        assert mc.binsize is None
        assert len(mc.bins()) == len(mc.contigs) == 4
        assert mc.genomes["n_bins"].tolist() == mc.genomes["n_contigs"].tolist()

    def test_needs_contig_lengths(self, runner, tmp_path, tiny_pairs):
        """No ``#chromsize`` lines and no ``--fasta``: nothing to tile."""
        out = tmp_path / "out.metacool"
        result = runner.invoke(
            cli,
            ["create", "pairs", "--output", str(out), "--pairs", str(tiny_pairs)],
        )
        assert result.exit_code != 0
        assert "Contig lengths required" in flat(result.output)
        assert not out.exists()


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
        """``metator`` format on a standard contigs table fails: it has no
        ``Name``, nor ``Final_cluster`` or ``Final_bin``."""
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


class TestCreatePairsPixelCorrectness:
    """``create pairs --contig-order fasta`` must reproduce ``cooler cload pairs``.

    Both the **bins** table and the **pixel** table must match the reference
    cooler exactly.  Comparing pixels alone is not enough: bin ids are only
    meaningful relative to a bin table, so two files could agree pixel-for-pixel
    while binning the genome differently.

    The reference coolers in ``tests/data/`` are built by ``raw-to-cool.sh``::

        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \\
            assembly_1.chrom.sizes:1000 contacts_1_rep1.pairs.gz \\
            assembly_1_rep1_1kb.cool
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \\
            assembly_1_DpnII-digested.bed contacts_1_rep1.pairs.gz \\
            assembly_1_rep1_DpnII.cool

    ``--contig-order fasta`` preserves the chromsize order from the pairs
    header, which matches both ``assembly_1.chrom.sizes`` and the digest BED, so
    bin numbering is identical and the tables are directly comparable.
    """

    DATA = Path("tests/data")

    # (label, --bins value, reference cooler) — fixed-width and variable-width.
    # The variable-width case is the one that exercises a caller-supplied bin
    # table rather than metacooler's own tiling.
    BIN_SPECS: ClassVar = [
        ("fixed-1kb", "1000", "assembly_1_rep1_1kb.cool"),
        (
            "variable-dpnii",
            str(Path("tests/data") / "assembly_1_DpnII-digested.bed"),
            "assembly_1_rep1_DpnII.cool",
        ),
    ]

    # (label, read_block, chunksize, --sort-strategy, pixel budget, expect_spill)
    #
    # None of these knobs may perturb the output by a single count. The
    # multi-chunk rows matter most — they are the only end-to-end check that
    # ``_buffer_pixels``' folding, and cooler's external sort, agree with a
    # plain single-pass ``cload`` to the row.
    #
    # Forcing many chunks needs BOTH levers: pyarrow's read block has an 8 MiB
    # floor and this 50 k-contact fixture fits inside one block, so lowering
    # ``--chunksize`` alone still yields exactly one chunk (the batches are
    # coalesced straight back) and the multi-chunk cases would pass vacuously.
    # Capping ``_PAIRS_READ_BLOCK_MAX_BYTES`` to 64 KiB gives ~39 real chunks.
    # ``expect_spill`` is None where ``--sort-strategy`` is explicit, since
    # ``_buffer_pixels`` is then bypassed entirely.
    STRATEGIES: ClassVar = [
        ("default-single-pass", None, None, "auto", None, False),
        ("many-chunks-in-memory", 1 << 16, 1000, "memory", None, None),
        ("many-chunks-external-sort", 1 << 16, 1000, "external", None, None),
        ("many-chunks-folded-in-memory", 1 << 16, 1000, "auto", 50_000, False),
        ("many-chunks-folded-then-spilled", 1 << 16, 1000, "auto", 500, True),
    ]

    @pytest.fixture
    def runner(self):
        return CliRunner()

    @pytest.mark.parametrize(
        ("bin_label", "bins_arg", "reference"),
        BIN_SPECS,
        ids=[b[0] for b in BIN_SPECS],
    )
    @pytest.mark.parametrize(
        ("label", "read_block", "chunksize", "strategy", "budget", "expect_spill"),
        STRATEGIES,
        ids=[s[0] for s in STRATEGIES],
    )
    def test_matches_cooler_cload(
        self,
        runner,
        tmp_path,
        monkeypatch,
        bin_label,
        bins_arg,
        reference,
        label,
        read_block,
        chunksize,
        strategy,
        budget,
        expect_spill,
    ):
        """Bins and pixels must equal the ``cload`` reference on every path."""
        import cooler
        import pandas as pd

        from metacooler.io import create as create_mod

        if budget is not None:
            monkeypatch.setattr(create_mod, "_SORT_IN_MEMORY_PIXEL_THRESHOLD", budget)

        # Spy on the code paths so a future change that quietly collapses these
        # cases into "one chunk, no folding" fails loudly instead of passing.
        n_chunks = 0
        spilled: bool | None = None
        real_stream = create_mod._stream_binned_pairs
        real_buffer = create_mod._buffer_pixels

        def counting_stream(*args, **kwargs):
            nonlocal n_chunks
            for frame in real_stream(*args, **kwargs):
                n_chunks += 1
                yield frame

        def spying_buffer(pixel_iter, max_rows, fold_rows=None):
            nonlocal spilled
            frames, rest = real_buffer(pixel_iter, max_rows, fold_rows)
            spilled = rest is not None
            return frames, rest

        monkeypatch.setattr(create_mod, "_stream_binned_pairs", counting_stream)
        monkeypatch.setattr(create_mod, "_buffer_pixels", spying_buffer)
        if read_block is not None:
            monkeypatch.setattr(create_mod, "_PAIRS_READ_BLOCK_MAX_BYTES", read_block)

        out = tmp_path / f"{bin_label}_{label}.metacool"
        args = [
            "create",
            "pairs",
            "--output",
            str(out),
            "--pairs",
            str(self.DATA / "contacts_1_rep1.pairs.gz"),
            "--bins",
            bins_arg,
            "--contigs",
            str(self.DATA / "contigs_table_1.tsv"),
            "--genomes",
            str(self.DATA / "mag_table_1.tsv"),
            "--contig-order",
            "fasta",
            "--sort-strategy",
            strategy,
        ]
        if chunksize is not None:
            args += ["--chunksize", str(chunksize)]

        result = runner.invoke(cli, args)
        assert result.exit_code == 0, result.output
        assert out.exists()

        # The case really did take the path its label claims.
        if read_block is None:
            assert n_chunks == 1, "baseline case should be a single chunk"
        else:
            assert n_chunks > 10, (
                f"expected many chunks, got {n_chunks} — the read-block cap and "
                "--chunksize no longer force chunking, so this case is vacuous"
            )
        if expect_spill is not None:
            assert spilled is expect_spill

        ref = cooler.Cooler(str(self.DATA / reference))
        got = cooler.Cooler(str(out))

        # --- bins: same tiling, same order, same widths ---
        cols = ["chrom", "start", "end"]
        bins_ref = ref.bins()[:][cols].reset_index(drop=True).astype({"chrom": str})
        bins_got = got.bins()[:][cols].reset_index(drop=True).astype({"chrom": str})
        pd.testing.assert_frame_equal(bins_ref, bins_got)

        # --- pixels ---
        pix_ref = cast(pd.DataFrame, ref.pixels()[:]).reset_index(drop=True)
        pix_got = cast(pd.DataFrame, got.pixels()[:]).reset_index(drop=True)
        pd.testing.assert_frame_equal(pix_ref, pix_got, check_like=False)

        # Guard the guard: an empty reference would satisfy everything above.
        assert len(bins_got) > 0
        assert len(pix_got) > 0
        assert pix_got["count"].sum() == pix_ref["count"].sum()
        if bin_label == "variable-dpnii":
            widths = bins_got["end"] - bins_got["start"]
            assert widths.nunique() > 1, "digest bins should not be uniform width"


# ===========================================================================
# create pairs & cool — sort_in_memory option
# ===========================================================================


class TestCreatePairsSortInMemory:
    """Test --sort-strategy options for create pairs."""

    def test_sort_in_memory_flag(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """--sort-strategy memory forces in-memory sorting."""
        out = tmp_path / "out_inmem.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--sort-strategy", "memory"],
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f

    def test_sort_external_flag(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """--sort-strategy external forces cooler external sort."""
        out = tmp_path / "out_external.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(
                out,
                tiny_pairs,
                tiny_bins_bed,
                tiny_contigs_table,
                extra=["--sort-strategy", "external"],
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f

    def test_sort_default_auto(
        self, runner, tmp_path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """Default (no flag) uses auto-detect."""
        out = tmp_path / "out_auto.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(out, tiny_pairs, tiny_bins_bed, tiny_contigs_table),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f

    def test_sort_auto_spills_to_external(
        self,
        monkeypatch,
        runner,
        tmp_path,
        tiny_pairs,
        tiny_bins_bed,
        tiny_contigs_table,
    ):
        """auto falls back to the external sort once the row budget is exceeded.

        The spilled path must produce exactly the same pixels as the in-memory
        one: the frames buffered before the threshold was crossed are replayed
        ahead of the rest of the stream, so none may be dropped or duplicated.
        """
        from metacooler.io import create as create_mod

        ref = tmp_path / "out_auto_ref.metacool"
        assert (
            runner.invoke(
                cli,
                _base_pairs_args(
                    ref,
                    tiny_pairs,
                    tiny_bins_bed,
                    tiny_contigs_table,
                    extra=["--sort-strategy", "memory"],
                ),
            ).exit_code
            == 0
        )

        # A budget of 0 rows forces the very first chunk over the threshold.
        monkeypatch.setattr(create_mod, "_SORT_IN_MEMORY_PIXEL_THRESHOLD", 0)
        out = tmp_path / "out_auto_spill.metacool"
        result = runner.invoke(
            cli,
            _base_pairs_args(out, tiny_pairs, tiny_bins_bed, tiny_contigs_table),
        )
        assert result.exit_code == 0, result.output

        with h5py.File(ref, "r") as a, h5py.File(out, "r") as b:
            assert a.attrs["nnz"] == b.attrs["nnz"]
            assert a.attrs["sum"] == b.attrs["sum"]
            for col in ("bin1_id", "bin2_id", "count"):
                assert_array_equal(a["pixels"][col][:], b["pixels"][col][:])


class TestCreateCoolSortInMemory:
    """Test --sort-strategy options for create cool."""

    def test_sort_in_memory_flag(self, runner, tmp_path, tiny_cool, tiny_contigs_table):
        """--sort-strategy memory forces in-memory sorting for create cool."""
        out = tmp_path / "out_inmem.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out, tiny_cool, tiny_contigs_table, extra=["--sort-strategy", "memory"]
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f

    def test_sort_external_flag(self, runner, tmp_path, tiny_cool, tiny_contigs_table):
        """--sort-strategy external forces cooler external sort for create cool."""
        out = tmp_path / "out_external.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(
                out,
                tiny_cool,
                tiny_contigs_table,
                extra=["--sort-strategy", "external"],
            ),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f

    def test_sort_default_auto(self, runner, tmp_path, tiny_cool, tiny_contigs_table):
        """Default (no flag) uses auto-detect for create cool."""
        out = tmp_path / "out_auto.metacool"
        result = runner.invoke(
            cli,
            _base_cool_args(out, tiny_cool, tiny_contigs_table),
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "pixels" in f
