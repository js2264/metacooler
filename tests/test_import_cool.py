"""Tests for cool to metacool import and round-trip conversion."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import cooler
import pytest
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli


class TestImportCool:
    """Tests for cool -> metacool import."""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_import_cool_help(self, runner):
        """Test that the cool subcommand help works."""
        result = runner.invoke(cli, ["create", "cool", "--help"])
        assert result.exit_code == 0
        assert "--contigs" in result.output
        assert "--output" in result.output

    def test_import_cool_basic(self, tmp_path: Path, runner):
        """Test basic cool to metacool conversion."""
        cool_path = Path("tests/data/assembly_1_rep1_1kb.cool")
        contigs_table = Path("tests/data/contigs_table_1.tsv")
        imported_metacool_path = tmp_path / "imported.metacool"

        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(cool_path),
                "--contigs",
                str(contigs_table),
                "-o",
                str(imported_metacool_path),
            ],
        )

        assert result.exit_code == 0, result.output
        assert imported_metacool_path.exists()

        mc = MetaCooler(imported_metacool_path)
        assert len(mc.genomes) > 0
        assert len(mc.contigs) > 0
        assert len(mc.bins()) > 0

    def test_roundtrip_cool_metacool_cool(self, tmp_path: Path, runner):
        """Test round-trip: cool -> metacool -> cool preserves data."""
        from metacooler.io.flatten import flatten

        cool1 = Path("tests/data/assembly_1_rep1_1kb.cool")
        contigs_table = Path("tests/data/contigs_table_1.tsv")
        imported_metacool = tmp_path / "imported.metacool"
        cool2 = tmp_path / "round2.cool"

        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(cool1),
                "--contigs",
                str(contigs_table),
                "-o",
                str(imported_metacool),
            ],
        )
        assert result.exit_code == 0, result.output

        # Flatten back to cool
        flatten(uri=imported_metacool, output=cool2)

        # Verify the output is a valid cool file
        clr2 = cooler.Cooler(str(cool2))
        assert len(clr2.bins()) > 0
        assert len(clr2.pixels()) > 0

        # Verify total pixel count is preserved
        clr1 = cooler.Cooler(str(cool1))
        assert len(clr1.pixels()) == len(clr2.pixels())

    def test_import_with_genomes_table(
        self,
        tmp_path: Path,
        runner,
        tiny_cool: Path,
        tiny_contigs_table: Path,
        tiny_genomes_table: Path,
    ):
        """Test import with optional genomes metadata table adds extra columns."""
        out = tmp_path / "with_genomes.metacool"

        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "--genomes",
                str(tiny_genomes_table),
                "-o",
                str(out),
            ],
        )

        assert result.exit_code == 0, result.output
        assert out.exists()

        mc = MetaCooler(out)
        genomes = mc.genomes
        # Extra columns from the genomes table should be present
        assert "quality" in genomes.columns
        assert "phylum" in genomes.columns
        # Values for known genomes should match the table
        mag1 = genomes.loc[genomes["name"] == "MAG_1"].iloc[0]
        assert mag1["quality"] == "high_quality"
        assert mag1["phylum"] == "Bacteroidetes"

    def test_import_invalid_cool_file(self, tmp_path: Path, runner):
        """Test that import fails gracefully for invalid cool files."""
        invalid_cool = tmp_path / "invalid.cool"
        contigs_table = Path("tests/data/contigs_table_1.tsv")
        output_metacool = tmp_path / "output.metacool"

        invalid_cool.write_text("not a cool file")

        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(invalid_cool),
                "--contigs",
                str(contigs_table),
                "-o",
                str(output_metacool),
            ],
        )

        assert result.exit_code != 0
        assert (
            "Not a valid cooler file" in str(result.exception)
            or "Error" in result.output
        )


class TestCoolRoundTripFixtures:
    """Round-trip tests using fixtures: cool -> metacool (create cool) -> cool (flatten)."""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_create_cool_from_fixture_produces_valid_metacool(
        self, tmp_path: Path, runner, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """create cool produces a valid metacool with expected genome structure."""
        out = tmp_path / "imported.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        mc = MetaCooler(out)
        assert set(mc.genomes["name"]) >= {"MAG_1", "MAG_2", "__unbinned__"}
        assert len(mc.contigs) == 4
        assert len(mc.bins()) > 0
        assert len(mc.pixels()) > 0

    def test_flatten_after_create_cool_preserves_bin_count(
        self, tmp_path: Path, runner, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """Bin count is preserved through the cool -> metacool -> cool round-trip."""
        out_metacool = tmp_path / "imported.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "-o",
                str(out_metacool),
            ],
        )
        assert result.exit_code == 0, result.output

        out_cool = tmp_path / "flattened.cool"
        result = runner.invoke(
            cli,
            ["flatten", str(out_metacool), "-o", str(out_cool)],
        )
        assert result.exit_code == 0, result.output
        assert out_cool.exists()

        original = cooler.Cooler(str(tiny_cool))
        recovered = cooler.Cooler(str(out_cool))
        assert original.shape == recovered.shape

    def test_flatten_after_create_cool_preserves_pixel_count(
        self, tmp_path: Path, runner, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """Pixel count (nnz) is preserved through the cool -> metacool -> cool round-trip."""
        out_metacool = tmp_path / "imported.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "-o",
                str(out_metacool),
            ],
        )
        assert result.exit_code == 0, result.output

        out_cool = tmp_path / "flattened.cool"
        result = runner.invoke(
            cli,
            ["flatten", str(out_metacool), "-o", str(out_cool)],
        )
        assert result.exit_code == 0, result.output

        original = cooler.Cooler(str(tiny_cool))
        recovered = cooler.Cooler(str(out_cool))
        assert original.info["nnz"] == recovered.info["nnz"]

    def test_flatten_after_create_cool_chromnames_are_contigs(
        self, tmp_path: Path, runner, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """After flatten, chromnames in the output cool match contig names."""
        out_metacool = tmp_path / "imported.metacool"
        runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(tiny_cool),
                "--contigs",
                str(tiny_contigs_table),
                "-o",
                str(out_metacool),
            ],
        )

        out_cool = tmp_path / "flattened.cool"
        runner.invoke(cli, ["flatten", str(out_metacool), "-o", str(out_cool)])

        original = cooler.Cooler(str(tiny_cool))
        recovered = cooler.Cooler(str(out_cool))
        assert set(original.chromnames) == set(recovered.chromnames)


class TestRoundTripChromLengths:
    """Tests that chrom name ↔ length mapping is preserved through cool→metacool→cool.

    These tests specifically guard against the chrom-ID-remapping bug where
    ``io.create.from_cool`` remaps bin chrom IDs but leaves bins in the original
    cool order. ``cooler.create_cooler`` then renumbers codes by first-appearance,
    causing a misalignment between the contigs table (name → expected length) and
    the actual bin data stored under each code.
    """

    DATA = Path("tests/data")

    @pytest.fixture
    def runner(self):
        from click.testing import CliRunner

        return CliRunner()

    def test_chrom_lengths_preserved_after_roundtrip_integration(self, tmp_path: Path):
        """Import assembly_1_rep1_1kb.cool → metacool → cool: every chrom keeps its length."""
        from metacooler.io.create import from_cool
        from metacooler.io.flatten import flatten

        mc_path = tmp_path / "imported.metacool"
        cool2 = tmp_path / "flat.cool"

        from_cool(
            cool_uri=self.DATA / "assembly_1_rep1_1kb.cool",
            output=mc_path,
            contig_genome_map=self.DATA / "contigs_table_1.tsv",
        )
        flatten(mc_path, cool2)

        original = cooler.Cooler(str(self.DATA / "assembly_1_rep1_1kb.cool"))
        recovered = cooler.Cooler(str(cool2))

        orig_sizes = dict(zip(original.chromnames, original.chromsizes, strict=False))
        recv_sizes = dict(zip(recovered.chromnames, recovered.chromsizes, strict=False))

        assert set(orig_sizes) == set(recv_sizes), "chrom names differ after round-trip"
        for name in orig_sizes:
            assert orig_sizes[name] == recv_sizes[name], (
                f"Length mismatch for {name}: "
                f"original={orig_sizes[name]}, recovered={recv_sizes[name]}"
            )

    def test_metacool_bins_consistent_with_contigs_table(self, tmp_path: Path):
        """In the imported metacool, each contig's max bin end matches its declared length."""
        from metacooler.io.create import from_cool

        mc_path = tmp_path / "imported.metacool"
        from_cool(
            cool_uri=self.DATA / "assembly_1_rep1_1kb.cool",
            output=mc_path,
            contig_genome_map=self.DATA / "contigs_table_1.tsv",
        )

        from metacooler.api import MetaCooler

        mc = MetaCooler(mc_path)
        contigs = mc.contigs.set_index("name")["length"]
        bins = mc.bins()

        for contig_name, grp in bins.groupby("chrom", observed=True):
            declared = int(contigs.loc[str(contig_name)])
            actual_max = int(grp["end"].max())
            assert actual_max == declared, (
                f"Bin max end {actual_max} ≠ declared length {declared} for {contig_name}"
            )

    def test_chrom_lengths_preserved_when_cool_order_differs_from_metacool_order(
        self, tmp_path: Path, tiny_cool_reversed_order: tuple
    ):
        """Round-trip preserves lengths when cool contig order ≠ genome-then-length order.

        Uses a cool where contigs are stored in ascending-length order (contig_D first),
        while genome-then-length reordering would produce descending-length order (contig_A
        first).  This directly exercises the bin-chrom-ID-remapping code path.
        """
        from metacooler.io.create import from_cool
        from metacooler.io.flatten import flatten

        cool_path, contigs_table = tiny_cool_reversed_order
        mc_path = tmp_path / "imported.metacool"
        cool2 = tmp_path / "flat.cool"

        from_cool(
            cool_uri=cool_path,
            output=mc_path,
            contig_genome_map=contigs_table,
        )
        flatten(mc_path, cool2)

        original = cooler.Cooler(str(cool_path))
        recovered = cooler.Cooler(str(cool2))

        orig_sizes = dict(zip(original.chromnames, original.chromsizes, strict=False))
        recv_sizes = dict(zip(recovered.chromnames, recovered.chromsizes, strict=False))

        assert set(orig_sizes) == set(recv_sizes), "chrom names differ"
        for name in orig_sizes:
            assert orig_sizes[name] == recv_sizes[name], (
                f"Length mismatch for {name}: "
                f"original={orig_sizes[name]}, recovered={recv_sizes[name]}"
            )

    def test_pixel_sum_per_contig_pair_preserved_after_roundtrip(
        self, tmp_path: Path, tiny_cool_reversed_order: tuple
    ):
        """Pixel counts per (chrom1, chrom2) pair are unchanged after round-trip."""
        from metacooler.io.create import from_cool
        from metacooler.io.flatten import flatten

        cool_path, contigs_table = tiny_cool_reversed_order
        mc_path = tmp_path / "imported.metacool"
        cool2 = tmp_path / "flat.cool"

        from_cool(
            cool_uri=cool_path,
            output=mc_path,
            contig_genome_map=contigs_table,
        )
        flatten(mc_path, cool2)

        original = cooler.Cooler(str(cool_path))
        recovered = cooler.Cooler(str(cool2))

        def _pair_counts(clr: cooler.Cooler) -> dict:
            bins = clr.bins()[:][["chrom"]].copy()
            pixels = clr.pixels()[:]
            pixels["chrom1"] = bins["chrom"].iloc[pixels["bin1_id"].values].values
            pixels["chrom2"] = bins["chrom"].iloc[pixels["bin2_id"].values].values
            # Normalize pair order so (A,B) and (B,A) count identically
            c1 = pixels["chrom1"].astype(str)
            c2 = pixels["chrom2"].astype(str)
            pixels["pair"] = [
                (min(a, b), max(a, b)) for a, b in zip(c1, c2, strict=False)
            ]
            return pixels.groupby("pair")["count"].sum().to_dict()

        orig_counts = _pair_counts(original)
        recv_counts = _pair_counts(recovered)

        assert orig_counts == recv_counts, (
            "Per-contig-pair pixel counts differ after round-trip"
        )

    def test_chrom_lengths_preserved_cli_roundtrip(
        self, tmp_path: Path, runner, tiny_cool_reversed_order: tuple
    ):
        """CLI round-trip (create cool + flatten) preserves chrom lengths."""
        from metacooler.cli.main import cli

        cool_path, contigs_table = tiny_cool_reversed_order
        mc_path = tmp_path / "imported.metacool"
        cool2 = tmp_path / "flat.cool"

        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(cool_path),
                "--contigs",
                str(contigs_table),
                "-o",
                str(mc_path),
            ],
        )
        assert result.exit_code == 0, result.output

        result = runner.invoke(cli, ["flatten", str(mc_path), "-o", str(cool2)])
        assert result.exit_code == 0, result.output

        original = cooler.Cooler(str(cool_path))
        recovered = cooler.Cooler(str(cool2))

        orig_sizes = dict(zip(original.chromnames, original.chromsizes, strict=False))
        recv_sizes = dict(zip(recovered.chromnames, recovered.chromsizes, strict=False))

        assert set(orig_sizes) == set(recv_sizes)
        for name in orig_sizes:
            assert orig_sizes[name] == recv_sizes[name], (
                f"CLI round-trip length mismatch for {name}: "
                f"original={orig_sizes[name]}, recovered={recv_sizes[name]}"
            )


# ===========================================================================
# create cool — strict table equality against the source cooler
# ===========================================================================


def _bins_frame(clr: cooler.Cooler):
    """Bins as ``chrom``/``start``/``end`` with ``chrom`` as plain strings."""
    import pandas as pd  # noqa: F401  (imported for the dtype cast below)

    return (
        clr.bins()[:][["chrom", "start", "end"]]
        .reset_index(drop=True)
        .astype({"chrom": str})
    )


def _coordinate_pixels(clr: cooler.Cooler):
    """Pixels re-keyed by genomic coordinate instead of bin id.

    Bin ids are only meaningful relative to a bin table, so a file whose contigs
    were reordered has legitimately different ``bin1_id``/``bin2_id`` while
    describing exactly the same contacts.  Re-keying to
    ``(contig, start)`` pairs — and normalising each pair's orientation, since
    reordering can flip which side of the upper triangle a contact lands on —
    gives a representation that is invariant to contig ordering, so the two can
    still be compared exactly.
    """
    import pandas as pd

    bins = clr.bins()[:][["chrom", "start"]].astype({"chrom": str})
    chrom = bins["chrom"].to_numpy()
    start = bins["start"].to_numpy()
    pix = clr.pixels()[:]
    side1 = list(zip(chrom[pix["bin1_id"]], start[pix["bin1_id"]], strict=True))
    side2 = list(zip(chrom[pix["bin2_id"]], start[pix["bin2_id"]], strict=True))
    ordered = [(a, b) if a <= b else (b, a) for a, b in zip(side1, side2, strict=True)]
    frame = pd.DataFrame(
        {
            "k1": [k[0] for k in ordered],
            "k2": [k[1] for k in ordered],
            "count": pix["count"].to_numpy(),
        }
    )
    return (
        frame.groupby(["k1", "k2"], as_index=False)["count"]
        .sum()
        .sort_values(["k1", "k2"])
        .reset_index(drop=True)
    )


class TestCreateCoolMatchesSource:
    """``create cool`` must reproduce its source cooler's tables exactly.

    Two distinct contracts, both checked here:

    - ``--contig-order cool`` preserves the source's contig order, so the bins
      and pixel tables must be **byte-identical** to the source — and must stay
      identical after ``flatten`` takes the file back to a plain ``.cool``.
    - The default ``genome-then-length`` deliberately regroups contigs by
      genome, so bin ids necessarily change.  The contacts themselves must not:
      re-keyed to genomic coordinates, the pixel table must still be identical.

    Existing round-trip tests only compared ``shape`` and ``nnz``, which a
    scrambled-but-same-size table would satisfy.
    """

    DATA = Path("tests/data")

    SOURCES = [
        ("fixed-1kb", "assembly_1_rep1_1kb.cool"),
        ("variable-dpnii", "assembly_1_rep1_DpnII.cool"),
    ]

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def _import(self, runner, tmp_path, source, order):
        out = tmp_path / f"imported_{order}.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "cool",
                "--cool",
                str(self.DATA / source),
                "--contigs",
                str(self.DATA / "contigs_table_1.tsv"),
                "--genomes",
                str(self.DATA / "mag_table_1.tsv"),
                "--contig-order",
                order,
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        return out

    @pytest.mark.parametrize(("label", "source"), SOURCES, ids=[s[0] for s in SOURCES])
    def test_contig_order_cool_is_table_identical(
        self, runner, tmp_path, label, source
    ):
        """``--contig-order cool`` reproduces bins and pixels exactly."""
        import pandas as pd

        out = self._import(runner, tmp_path, source, "cool")
        src = cooler.Cooler(str(self.DATA / source))
        got = cooler.Cooler(str(out))

        pd.testing.assert_frame_equal(_bins_frame(src), _bins_frame(got))
        pd.testing.assert_frame_equal(
            cast(pd.DataFrame, src.pixels()[:]).reset_index(drop=True),
            cast(pd.DataFrame, got.pixels()[:]).reset_index(drop=True),
        )
        assert got.info["sum"] == src.info["sum"]
        assert len(got.pixels()[:]) > 0

    @pytest.mark.parametrize(("label", "source"), SOURCES, ids=[s[0] for s in SOURCES])
    def test_flatten_round_trip_is_table_identical(
        self, runner, tmp_path, label, source
    ):
        """cool -> metacool -> cool under ``--contig-order cool`` is lossless."""
        import pandas as pd

        out = self._import(runner, tmp_path, source, "cool")
        flat = tmp_path / "flattened.cool"
        result = runner.invoke(cli, ["flatten", str(out), "-o", str(flat)])
        assert result.exit_code == 0, result.output

        src = cooler.Cooler(str(self.DATA / source))
        back = cooler.Cooler(str(flat))
        pd.testing.assert_frame_equal(_bins_frame(src), _bins_frame(back))
        pd.testing.assert_frame_equal(
            cast(pd.DataFrame, src.pixels()[:]).reset_index(drop=True),
            cast(pd.DataFrame, back.pixels()[:]).reset_index(drop=True),
        )

    @pytest.mark.parametrize(("label", "source"), SOURCES, ids=[s[0] for s in SOURCES])
    def test_default_ordering_preserves_contacts_not_bin_ids(
        self, runner, tmp_path, label, source
    ):
        """Regrouping by genome may renumber bins but must not move contacts."""
        import pandas as pd

        out = self._import(runner, tmp_path, source, "genome-then-length")
        src = cooler.Cooler(str(self.DATA / source))
        got = cooler.Cooler(str(out))

        # Bin ids legitimately differ...
        assert not _bins_frame(src).equals(_bins_frame(got)), (
            "genome-then-length should reorder contigs on this fixture; if it no "
            "longer does, this test has stopped exercising the reordered path"
        )
        # ...but the contacts, keyed by coordinate, must not.
        pd.testing.assert_frame_equal(_coordinate_pixels(src), _coordinate_pixels(got))
        assert got.info["sum"] == src.info["sum"]
