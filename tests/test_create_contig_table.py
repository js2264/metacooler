"""``create`` refuses a contig table that names no contig of the assembly.

A contig the ``--contigs`` table does not list is collected under
``__unbinned__``, so a table matching *none* of them -- the wrong table, or the
right one read with the wrong ``--contigs-format`` -- produced a map whose only
genome was ``__unbinned__``, exited 0, and logged nothing above INFO.  The case
reproduced here is the common one: a table with a leading ID column, read with
``--contigs-format 1,2``, takes the IDs for contig names and the contig names
for genome names.

A table naming *some* contigs of the assembly is accepted, and the rest of its
names are only warned about.  Either way the table is compared with every
contig of the assembly, not with what ``--drop-unassigned`` or ``--filter``
leave of it.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from _cli_output import flat, squash
from click.testing import CliRunner, Result

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.create import from_cool, from_pairs
from metacooler.schema import validate

COMMANDS = ["pairs", "cool"]


def _table(path: Path, rows: list[tuple[str, str]]) -> Path:
    """Write a ``contigs``-format table (header + ``contig<TAB>genome`` rows)."""
    lines = ["contig_name\tgenome_name"] + [f"{c}\t{g}" for c, g in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _with_leading_id(table: Path, out: Path) -> Path:
    """*table* with an ``ID`` column in front, which ``1,2`` misreads."""
    header, *rows = table.read_text().splitlines()
    lines = [f"ID\t{header}"] + [f"{i}\t{row}" for i, row in enumerate(rows, 1)]
    out.write_text("\n".join(lines) + "\n")
    return out


@pytest.fixture
def create(tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, tiny_cool: Path):
    """Run ``create pairs`` or ``create cool`` over the tiny 4-contig assembly."""
    runner = CliRunner()

    def run(command: str, table: Path, *extra: str) -> tuple[Result, Path]:
        out = tmp_path / f"{command}.metacool"
        source = (
            ["--pairs", str(tiny_pairs), "--bins", str(tiny_bins_bed), "--no-ps"]
            if command == "pairs"
            else ["--cool", str(tiny_cool)]
        )
        args = ["create", command, *source, "--contigs", str(table), "-o", str(out)]
        return runner.invoke(cli, [*args, *extra]), out

    return run


@pytest.fixture
def no_pixel_pass(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail the test if a pair or a pixel is read: the refusal must come first."""

    def fail(*args: object, **kwargs: object) -> None:
        raise AssertionError("the pixel pass started before the table was refused")

    monkeypatch.setattr("metacooler.io.create._stream_binned_pairs", fail)
    monkeypatch.setattr("metacooler.io.create._stream_remapped_pixels", fail)


def _assert_nothing_written(out: Path) -> None:
    assert not out.exists()
    assert not out.with_name(out.name + ".partial").exists()


# ---------------------------------------------------------------------------
# A table matching no contig of the assembly
# ---------------------------------------------------------------------------


@pytest.mark.usefixtures("no_pixel_pass")
@pytest.mark.parametrize("command", COMMANDS)
class TestRefusal:
    def test_leading_id_column_read_as_contig_names(
        self, create, command: str, tmp_path: Path, tiny_contigs_table: Path
    ):
        """The reported case: exit 0 with a lone ``__unbinned__`` genome."""
        table = _with_leading_id(tiny_contigs_table, tmp_path / "with_id.tsv")

        result, out = create(command, table, "--contigs-format", "1,2")

        assert result.exit_code == 1
        message = flat(result.output)
        assert "None of the 3 contig name(s) in the contig table" in message
        assert str(table) in squash(result.output)
        # A few names from each side, quoted, next to each other.
        assert "table names: '1', '2', '3'" in message
        assert "assembly contigs: 'contig_A', 'contig_B', 'contig_C', ..." in message
        assert "no contig would be assigned to a genome" in message
        assert "--contigs-format (here '1,2')" in message
        _assert_nothing_written(out)

    def test_refused_under_drop_unassigned_too(
        self, create, command: str, tmp_path: Path
    ):
        """The refusal pre-empts the vaguer "--drop-unassigned removed all"."""
        table = _table(tmp_path / "t.tsv", [("not_a_contig", "MAG_1")])

        result, out = create(command, table, "--drop-unassigned")

        assert result.exit_code == 1
        assert "None of the 1 contig name(s)" in flat(result.output)
        assert "removed all" not in flat(result.output)
        _assert_nothing_written(out)

    def test_table_without_rows(self, create, command: str, tmp_path: Path):
        table = _table(tmp_path / "t.tsv", [])

        result, out = create(command, table)

        assert result.exit_code == 1
        assert "There is no row in the contig table" in flat(result.output)
        _assert_nothing_written(out)


def test_the_cool_file_is_named_as_the_assembly(create, tiny_cool: Path, tmp_path):
    table = _table(tmp_path / "t.tsv", [("not_a_contig", "MAG_1")])

    result, _ = create("cool", table)

    expected = f"matches a contig in the cool file {tiny_cool}"
    assert "".join(expected.split()) in squash(result.output)


# ---------------------------------------------------------------------------
# A table naming contigs beyond the assembly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("command", COMMANDS)
def test_names_beyond_the_assembly_are_ignored_with_a_warning(
    create, command: str, tmp_path: Path, tiny_contigs_table: Path
):
    """A table covering a superset of the assembly is legitimate: warn, build."""
    table = tmp_path / "superset.tsv"
    table.write_text(
        tiny_contigs_table.read_text() + "other_1\tMAG_9\nother_2\tMAG_9\n"
    )

    result, out = create(command, table)

    assert result.exit_code == 0, result.output
    output = flat(result.output)
    assert "[WARNING] 2 of the 5 name(s) in the contig table" in output
    assert "and were ignored: 'other_1', 'other_2'" in output
    validate(out)
    mc = MetaCooler(out)
    assert sorted(mc.contigs["name"]) == [f"contig_{c}" for c in "ABCD"]
    # MAG_9 has no contig of this assembly, so it is not a genome of the map.
    assert sorted(mc.genomes["name"]) == ["MAG_1", "MAG_2", "__unbinned__"]


def test_the_warning_names_at_most_ten(
    tmp_path: Path, tiny_cool: Path, tiny_contigs_table: Path, caplog
):
    table = tmp_path / "superset.tsv"
    strays = "".join(f"stray_{i:02d}\tMAG_9\n" for i in range(12))
    table.write_text(tiny_contigs_table.read_text() + strays)

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        from_cool(tiny_cool, tmp_path / "out.metacool", contig_genome_map=table)

    assert "12 of the 15 name(s)" in caplog.text
    assert "'stray_09', ..." in caplog.text
    assert "stray_10" not in caplog.text


def test_a_matching_table_draws_no_warning(
    tmp_path: Path, tiny_cool: Path, tiny_contigs_table: Path, caplog
):
    """contig_D is absent from the table: that is unbinned, not an error."""
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        from_cool(
            tiny_cool, tmp_path / "out.metacool", contig_genome_map=tiny_contigs_table
        )

    assert "are not contigs in" not in caplog.text


# ---------------------------------------------------------------------------
# --drop-unassigned and --filter narrow the contigs after the comparison
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("command", COMMANDS)
def test_a_filter_keeping_only_unlisted_contigs_is_not_a_refusal(
    create, command: str, tiny_contigs_table: Path
):
    """The table names A, B and C; the filter keeps only D (200 bp).

    Compared with the filtered set, the table would match nothing.  Compared
    with the assembly, it matches three contigs, and the map is what the user
    asked for: the unlisted contig, unbinned.
    """
    result, out = create(
        command, tiny_contigs_table, "--filter", "contigs.length < 300"
    )

    assert result.exit_code == 0, result.output
    assert "are not contigs in" not in flat(result.output)
    mc = MetaCooler(out)
    assert mc.contigs["name"].tolist() == ["contig_D"]
    assert mc.genomes["name"].tolist() == ["__unbinned__"]


@pytest.mark.parametrize("command", COMMANDS)
def test_drop_unassigned_does_not_turn_dropped_contigs_into_strays(
    create, command: str, tmp_path: Path
):
    """contig_C is listed, as unbinned, and dropped: it is still no stray."""
    table = _table(
        tmp_path / "t.tsv",
        [("contig_A", "MAG_1"), ("contig_B", "MAG_1"), ("contig_C", "NA")],
    )

    result, out = create(command, table, "--drop-unassigned")

    assert result.exit_code == 0, result.output
    assert "are not contigs in" not in flat(result.output)
    assert sorted(MetaCooler(out).contigs["name"]) == ["contig_A", "contig_B"]


# ---------------------------------------------------------------------------
# The message
# ---------------------------------------------------------------------------


def _refusal(tmp_path: Path, **kwargs: Any) -> str:
    """Run ``from_pairs`` with a table naming no contig; return its message."""
    if "contig_genome_map" not in kwargs:
        kwargs["contig_genome_map"] = _table(
            tmp_path / "unmatched.tsv", [("not_a_contig", "MAG_1")]
        )
    with pytest.raises(ValueError, match="None of the") as exc:
        from_pairs(uri=tmp_path / "out.metacool", ps=False, **kwargs)
    return str(exc.value)


@pytest.mark.parametrize("source", ["bins file", "bins table", "FASTA", "header"])
def test_the_message_names_where_the_contig_names_came_from(
    tmp_path: Path,
    tiny_pairs: Path,
    tiny_bins_bed: Path,
    tiny_fasta: Path,
    source: str,
):
    """create pairs reads its contigs from one of four places; say which."""
    headed = tmp_path / "headed.pairs"
    headed.write_text(
        "## pairs format v1.0\n"
        "#chromsize: contig_A 1000\n"
        "#chromsize: contig_B 600\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
        "r1\tcontig_A\t10\tcontig_B\t20\t+\t+\n"
    )
    kwargs, expected = {
        "bins file": (
            {"bins": tiny_bins_bed, "pairs": tiny_pairs},
            f"in the bins file {tiny_bins_bed}",
        ),
        "bins table": (
            {
                "bins": pd.read_csv(
                    tiny_bins_bed, sep="\t", names=["chrom", "start", "end"]
                ),
                "pairs": tiny_pairs,
            },
            "in the bins table",
        ),
        "FASTA": (
            {"binsize": 100, "fasta": tiny_fasta, "pairs": tiny_pairs},
            f"in the FASTA {tiny_fasta}",
        ),
        "header": ({"binsize": 100, "pairs": headed}, "in the pairs header(s)"),
    }[source]

    assert expected in _refusal(tmp_path, **kwargs)


def test_a_mapping_is_not_blamed_on_contigs_format(
    tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
):
    """A dict passed to the API has no --contigs-format to get wrong."""
    message = _refusal(
        tmp_path,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map={"not_a_contig": "MAG_1"},
    )

    assert "in the contig-to-genome mapping" in message
    assert "--contigs-format" not in message


def test_names_are_quoted_so_whitespace_shows(
    tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
):
    """Trailing spaces are invisible unquoted, and the usual cause of a miss."""
    table = _table(tmp_path / "t.tsv", [("contig_A ", "MAG_1"), ("contig_B ", "MAG_1")])

    message = _refusal(
        tmp_path, bins=tiny_bins_bed, pairs=tiny_pairs, contig_genome_map=table
    )

    assert "table names:      'contig_A ', 'contig_B '\n" in message


def test_an_empty_mapping_is_refused(
    tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
):
    with pytest.raises(ValueError, match="no row in the contig-to-genome mapping"):
        from_pairs(
            uri=tmp_path / "out.metacool",
            bins=tiny_bins_bed,
            pairs=tiny_pairs,
            contig_genome_map={},
        )
