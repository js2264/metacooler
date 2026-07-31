"""Tests for metacooler.filters and the ``--filter`` option it backs.

The ``tiny_metacooler`` fixture is the reference throughout:

    genome        length  n_contigs  contigs
    __unbinned__     200          1  contig_D (200)
    MAG_1           1600          2  contig_A (1000), contig_B (600)
    MAG_2            400          1  contig_C (400)
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.filters import (
    FilterError,
    apply_filters,
    filter_mask,
    parse_filter,
    parse_filters,
)


@pytest.fixture
def frame() -> pd.DataFrame:
    """A stand-in genomes table with the awkward column names real files have."""
    return pd.DataFrame(
        {
            "name": ["MAG_1", "MAG_2", "MAG_3"],
            "length": [1_500_000, 400_000, 9_000],
            "Present Markers": [104, 51, np.nan],
            "GC-content": ["55.8", "42.1", "-"],
        }
    )


class TestParse:
    @pytest.mark.parametrize(
        ("expr", "table", "column", "op", "value"),
        [
            ("genomes.length > 1Mb", "genomes", "length", ">", "1Mb"),
            ("contigs.length>=10kb", "contigs", "length", ">=", "10kb"),
            ("genomes.n_contigs < 5", "genomes", "n_contigs", "<", "5"),
            ("contigs.length <= 200", "contigs", "length", "<=", "200"),
            ("genomes.name == MAG_1", "genomes", "name", "==", "MAG_1"),
            ("genomes.name != __unbinned__", "genomes", "name", "!=", "__unbinned__"),
            ("contigs.name =~ ^NODE_", "contigs", "name", "=~", "^NODE_"),
            ("contigs.name !~ ^NODE_", "contigs", "name", "!~", "^NODE_"),
            # A lone "=" is an alias for "==".
            ("genomes.name = MAG_1", "genomes", "name", "==", "MAG_1"),
            # Column names with spaces, dashes and dots need no quoting: the
            # first operator terminates the column.
            ("genomes.Present Markers >= 50", "genomes", "Present Markers", ">=", "50"),
            ("genomes.GC-content > 40", "genomes", "GC-content", ">", "40"),
            ("contigs.GC.content > 40", "contigs", "GC.content", ">", "40"),
            # Quotes around a value are stripped.
            ('genomes.name == "MAG 1"', "genomes", "name", "==", "MAG 1"),
        ],
    )
    def test_shapes(self, expr, table, column, op, value):
        filt = parse_filter(expr)
        assert (filt.table, filt.column, filt.op, filt.value) == (
            table,
            column,
            op,
            value,
        )
        assert filt.raw == expr.strip()

    def test_operator_search_is_leftmost(self):
        """A value may contain operator characters; only the first one splits."""
        filt = parse_filter("contigs.name =~ ^a>b=c")
        assert filt.op == "=~"
        assert filt.value == "^a>b=c"

    def test_missing_prefix_names_both_tables(self):
        with pytest.raises(FilterError, match="genomes.length.*contigs.length"):
            parse_filter("length > 1Mb")

    def test_unknown_table(self):
        with pytest.raises(FilterError, match="unknown table 'bins'"):
            parse_filter("bins.start > 100")

    @pytest.mark.parametrize(
        ("expr", "match"),
        [
            ("genomes.length", "no comparison operator"),
            ("genomes.length >", "nothing on the right"),
            ("> 100", "nothing on the left"),
            ("genomes. > 100", "no column"),
            ("", "empty expression"),
        ],
    )
    def test_malformed(self, expr, match):
        with pytest.raises(FilterError, match=match):
            parse_filter(expr)

    def test_bad_regex_is_rejected_at_parse_time(self):
        with pytest.raises(FilterError, match="not a valid regular expression"):
            parse_filter("contigs.name =~ [unclosed")

    def test_parse_filters_preserves_order(self):
        parsed = parse_filters(["genomes.length > 1", "contigs.length > 2"])
        assert [f.table for f in parsed] == ["genomes", "contigs"]
        assert parse_filters(None) == []


class TestUnits:
    @pytest.mark.parametrize(
        ("value", "kept"),
        [
            ("999999", 1),  # only the 1.5 Mb genome clears it
            ("1e6", 1),
            ("1_000_000", 1),
            ("1Mb", 1),
            ("1mb", 1),
            ("1M", 1),
            ("0.5Mb", 1),
            ("500kb", 1),
            ("500K", 1),
            ("10000bp", 2),
            ("1Gb", 0),
        ],
    )
    def test_size_suffixes(self, frame, value, kept):
        mask = filter_mask(
            frame, parse_filters([f"genomes.length > {value}"]), table="genomes"
        )
        assert int(mask.sum()) == kept

    def test_non_numeric_value(self, frame):
        with pytest.raises(FilterError, match="'big' is not a number"):
            filter_mask(frame, parse_filters(["genomes.length > big"]), table="genomes")


class TestEvaluate:
    def test_only_the_addressed_table_is_applied(self, frame):
        """A caller can hand the whole filter list to either table."""
        parsed = parse_filters(["contigs.length > 10", "genomes.length > 1Mb"])
        mask = filter_mask(frame, parsed, table="genomes")
        assert list(frame.loc[mask, "name"]) == ["MAG_1"]
        # No genomes filter reaches the contigs table here, so it keeps all.
        assert filter_mask(frame, parse_filters([]), table="genomes").all()

    def test_repeated_filters_are_anded(self, frame):
        parsed = parse_filters(["genomes.length > 10kb", "genomes.length < 1Mb"])
        assert list(frame.loc[filter_mask(frame, parsed, table="genomes"), "name"]) == [
            "MAG_2"
        ]

    def test_regex_is_a_search_not_a_fullmatch(self, frame):
        parsed = parse_filters(["genomes.name =~ AG_[12]"])
        assert list(apply_filters(frame, parsed, table="genomes")["name"]) == [
            "MAG_1",
            "MAG_2",
        ]

    def test_column_with_spaces(self, frame):
        parsed = parse_filters(["genomes.Present Markers >= 100"])
        assert list(apply_filters(frame, parsed, table="genomes")["name"]) == ["MAG_1"]

    def test_missing_values_never_match(self, frame):
        """NA fails every operator, negated ones included."""
        for expr in (
            "genomes.Present Markers > 0",
            "genomes.Present Markers < 1e9",
            "genomes.Present Markers != 104",
            "genomes.Present Markers == 104",
        ):
            kept = apply_filters(frame, parse_filters([expr]), table="genomes")
            assert "MAG_3" not in list(kept["name"]), expr

    def test_text_column_holding_numbers_is_coerced(self, frame, caplog):
        """MetaTOR writes '-' for a missing value; the rest is still numeric."""
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            kept = apply_filters(
                frame, parse_filters(["genomes.GC-content > 50"]), table="genomes"
            )
        assert list(kept["name"]) == ["MAG_1"]
        assert "are not numbers" in caplog.text

    def test_wholly_non_numeric_column_is_rejected_for_order_ops(self, frame):
        with pytest.raises(FilterError, match="is not numeric"):
            filter_mask(frame, parse_filters(["genomes.name > 5"]), table="genomes")

    def test_unknown_column_suggests_the_close_one(self, frame):
        with pytest.raises(FilterError, match="Did you mean 'length'"):
            filter_mask(frame, parse_filters(["genomes.lenght > 5"]), table="genomes")


class TestDumpFilter:
    def _dump(self, tmp_path, metacool, table, **kw):
        from metacooler.io.dump import dump_table

        out = tmp_path / "out.tsv"
        dump_table(metacool, table, header=True, out=out, force=True, **kw)
        return out.read_text()

    def test_genome_filter(self, tiny_metacooler, tmp_path):
        text = self._dump(
            tmp_path, tiny_metacooler, "genomes", filters=["genomes.length > 500"]
        )
        assert "MAG_1" in text
        assert "MAG_2" not in text
        assert "__unbinned__" not in text

    def test_contig_filter(self, tiny_metacooler, tmp_path):
        text = self._dump(
            tmp_path, tiny_metacooler, "contigs", filters=["contigs.length >= 400"]
        )
        for kept in ("contig_A", "contig_B", "contig_C"):
            assert kept in text
        assert "contig_D" not in text

    def test_genome_filter_narrows_the_contigs_table(self, tiny_metacooler, tmp_path):
        """A genomes.* filter selects that genome's contigs, not just its row."""
        text = self._dump(
            tmp_path, tiny_metacooler, "contigs", filters=["genomes.name == MAG_1"]
        )
        assert "contig_A" in text and "contig_B" in text
        assert "contig_C" not in text and "contig_D" not in text

    def test_virtual_genome_column(self, tiny_metacooler, tmp_path):
        text = self._dump(
            tmp_path, tiny_metacooler, "contigs", filters=["contigs.genome == MAG_2"]
        )
        assert "contig_C" in text
        assert "contig_A" not in text

    def test_names_and_filters_intersect(self, tiny_metacooler, tmp_path):
        """--genomes chooses the candidates; --filter narrows them."""
        text = self._dump(
            tmp_path,
            tiny_metacooler,
            "contigs",
            genomes=["MAG_1"],
            filters=["contigs.length > 700"],
        )
        assert "contig_A" in text
        # contig_B is in MAG_1 but too short; contig_C clears the length but is
        # not in MAG_1 — a union would have kept both.
        assert "contig_B" not in text
        assert "contig_C" not in text

    def test_pixels_are_filtered(self, tiny_metacooler, tmp_path):
        both = self._dump(tmp_path, tiny_metacooler, "pixels")
        few = self._dump(
            tmp_path, tiny_metacooler, "pixels", filters=["contigs.length >= 1kb"]
        )
        assert len(few.splitlines()) < len(both.splitlines())

    def test_empty_selection_names_the_culprit(self, tiny_metacooler, tmp_path):
        with pytest.raises(FilterError) as exc:
            self._dump(
                tmp_path,
                tiny_metacooler,
                "genomes",
                filters=["genomes.length > 1Gb", "contigs.length > 10"],
            )
        message = str(exc.value)
        assert "genomes.length > 1Gb" in message
        assert "keeps nothing" in message
        # The innocent filter is reported too, with its own count.
        assert "4 of 4 contigs" in message


class TestExtractFilter:
    def test_filter_alone_is_a_valid_selection(self, tiny_metacooler, tmp_path):
        import cooler

        from metacooler.io.extract import extract

        out = tmp_path / "e.cool"
        extract(tiny_metacooler, filters=["contigs.length >= 400"], output=out)
        assert sorted(cooler.Cooler(str(out)).chromnames) == [
            "contig_A",
            "contig_B",
            "contig_C",
        ]

    def test_split_writes_one_file_per_surviving_genome(
        self, tiny_metacooler, tmp_path
    ):
        from metacooler.io.extract import extract

        extract(
            tiny_metacooler,
            filters=["genomes.length > 300"],
            output=tmp_path / "sp",
            split=True,
        )
        written = sorted(p.name for p in tmp_path.glob("sp_*.cool"))
        assert written == ["sp_MAG_1.cool", "sp_MAG_2.cool"]

    def test_no_selector_at_all_still_raises(self, tiny_metacooler, tmp_path):
        from metacooler.io.extract import extract

        with pytest.raises(ValueError, match="--genomes, --contigs or --filter"):
            extract(tiny_metacooler, output=tmp_path / "e.cool")


class TestShowFilter:
    def test_genome_filter_stays_at_genome_level(self, tiny_metacooler, tmp_path):
        from metacooler.ops.show import show

        out = tmp_path / "s.png"
        show(tiny_metacooler, out=out, filters=["genomes.length > 300"])
        assert out.exists()

    def test_contig_filter_upgrades_the_level(self, tiny_metacooler, tmp_path, caplog):
        from metacooler.ops.show import show

        with caplog.at_level(logging.INFO, logger="metacooler"):
            show(
                tiny_metacooler,
                out=tmp_path / "s.png",
                filters=["contigs.length >= 400"],
            )
        assert "showing --level contig instead" in caplog.text

    def test_range_and_filter_are_mutually_exclusive(self, tiny_metacooler, tmp_path):
        from metacooler.ops.show import show

        with pytest.raises(ValueError, match="cannot be combined"):
            show(
                tiny_metacooler,
                out=tmp_path / "s.png",
                range1="contig_A",
                filters=["contigs.length > 1"],
            )

    def test_empty_selection_raises(self, tiny_metacooler, tmp_path):
        from metacooler.ops.show import show

        with pytest.raises(FilterError, match="Nothing passed --filter"):
            show(
                tiny_metacooler,
                out=tmp_path / "s.png",
                filters=["genomes.length > 1Gb"],
            )


class TestCreateFilter:
    def test_short_contigs_never_reach_the_bin_table(
        self, tmp_path, tiny_fasta, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        from metacooler.api import MetaCooler
        from metacooler.io.create import from_pairs

        out = tmp_path / "f.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
            filters=["contigs.length >= 400"],
            chunksize=1000,
        )
        mc = MetaCooler(out)
        assert sorted(mc.contigs["name"]) == ["contig_A", "contig_B", "contig_C"]
        # contig_D contributed no bin at all, not merely an empty genome row.
        assert "contig_D" not in set(mc.bins()["chrom"].astype(str))
        assert "__unbinned__" not in set(mc.genomes["name"])

    def test_genome_filter_drops_the_whole_genome(
        self, tmp_path, tiny_fasta, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        from metacooler.api import MetaCooler
        from metacooler.io.create import from_pairs

        out = tmp_path / "f.metacool"
        from_pairs(
            uri=out,
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map=tiny_contigs_table,
            contig_table_format="contigs",
            fasta=tiny_fasta,
            filters=["genomes.n_contigs >= 2"],
            chunksize=1000,
        )
        mc = MetaCooler(out)
        assert list(mc.genomes["name"]) == ["MAG_1"]
        assert sorted(mc.contigs["name"]) == ["contig_A", "contig_B"]

    def test_n_bins_is_reported_as_unavailable(
        self, tmp_path, tiny_fasta, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        """n_bins does not exist until the bins are tiled — say so, don't guess."""
        from metacooler.io.create import from_pairs

        with pytest.raises(FilterError, match="no bin has been tiled yet"):
            from_pairs(
                uri=tmp_path / "f.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                contig_table_format="contigs",
                fasta=tiny_fasta,
                filters=["genomes.n_bins > 3"],
                chunksize=1000,
            )

    def test_from_cool_filters_a_user_genomes_column(
        self, tmp_path, tiny_cool, tiny_contigs_table, tiny_genomes_table
    ):
        from metacooler.api import MetaCooler
        from metacooler.io.create import from_cool

        out = tmp_path / "f.metacool"
        from_cool(
            cool_uri=tiny_cool,
            output=out,
            contig_genome_map=tiny_contigs_table,
            genomes_table=tiny_genomes_table,
            filters=["genomes.quality == high_quality"],
        )
        mc = MetaCooler(out)
        assert list(mc.genomes["name"]) == ["MAG_1"]

    def test_nothing_survives(
        self, tmp_path, tiny_fasta, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        from metacooler.io.create import from_pairs

        with pytest.raises(FilterError, match="keeps nothing"):
            from_pairs(
                uri=tmp_path / "f.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                contig_table_format="contigs",
                fasta=tiny_fasta,
                filters=["contigs.length > 1Mb"],
                chunksize=1000,
            )


class TestCliOption:
    def test_bad_expression_fails_before_the_file_is_opened(self, tmp_path):
        """Parsing in the callback means a typo does not need a valid path."""
        result = CliRunner().invoke(
            cli, ["dump", "--filter", "length > 1Mb", str(tmp_path / "nope.metacool")]
        )
        assert result.exit_code != 0
        assert "must say which table" in _flat(result.output)

    def test_value_is_not_comma_split(self, tiny_metacooler):
        """Unlike --genomes, a filter value may contain a comma."""
        result = CliRunner().invoke(
            cli,
            [
                "dump",
                "-t",
                "contigs",
                "--filter",
                "contigs.name =~ ^contig_(?:A|B),?$",
                str(tiny_metacooler),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "contig_A" in result.output
        assert "contig_C" not in result.output

    def test_rejected_on_the_indexes_table(self, tiny_metacooler):
        result = CliRunner().invoke(
            cli,
            [
                "dump",
                "-t",
                "indexes",
                "--filter",
                "contigs.length > 1",
                str(tiny_metacooler),
            ],
        )
        assert result.exit_code != 0
        assert "--filter cannot be used with --table indexes" in _flat(result.output)

    def test_repeated_option_ands(self, tiny_metacooler):
        result = CliRunner().invoke(
            cli,
            [
                "dump",
                "-t",
                "contigs",
                "--filter",
                "contigs.length > 300",
                "--filter",
                "contigs.length < 800",
                str(tiny_metacooler),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "contig_B" in result.output
        assert "contig_A" not in result.output
        assert "contig_D" not in result.output
