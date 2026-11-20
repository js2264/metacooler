"""The tables as ``dump`` prints them and ``--filter`` reads them: keyed by name.

``dump -t contigs`` is ``contig, genome, length, ...`` (the genome by name, the
internal ``genome_id`` not at all), ``dump -t genomes`` is ``genome, ...``,
``dump -t clusters`` is ``round, cluster, ...`` and ``dump -t members`` is
``round, cluster, contig``: every dump joins with the others on ``contig``,
``genome`` or ``(round, cluster)``. ``contig`` and ``genome`` are therefore
reserved: a user column of either name is stored as ``user_contig`` /
``user_genome``, and one a map stored before is shown so. Filters accept the
dump's names, and a table ``dump -H`` wrote reads back without a format.
"""

from __future__ import annotations

import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.filters import FilterError
from metacooler.io._utils import (
    contigs_view,
    genomes_view,
    resolve_contig_selection,
    write_column,
)
from metacooler.io.annotate import annotate
from metacooler.io.clusters import write_clusters
from metacooler.io.create import from_pairs
from metacooler.io.dump import dump_table
from metacooler.io.reassign import reassign


def _map(directory: Path, contigs_table: str, genomes_table: str | None = None) -> Path:
    """A two-contig map: c1 in MAG_1, c2 in MAG_2, from the tables given."""
    (directory / "bins.bed").write_text("c1\t0\t1000\nc2\t0\t500\n")
    (directory / "lib.pairs").write_text(
        "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n"
        "r1\tc1\t10\tc2\t20\n"
    )
    (directory / "contigs.tsv").write_text(contigs_table)
    genomes = None
    if genomes_table is not None:
        genomes = directory / "genomes.tsv"
        genomes.write_text(genomes_table)
    path = directory / "map.metacool"
    from_pairs(
        path,
        bins=directory / "bins.bed",
        pairs=directory / "lib.pairs",
        contig_genome_map=directory / "contigs.tsv",
        genomes_table=genomes,
        ps=False,
    )
    return path


@pytest.fixture
def plain(tmp_path: Path) -> Path:
    return _map(
        tmp_path, "contig_name\tgenome_name\tgc\nc1\tMAG_1\t0.4\nc2\tMAG_2\t0.6\n"
    )


def _dump(path: Path, table: str, *args: str) -> list[list[str]]:
    result = CliRunner().invoke(cli, ["dump", str(path), "-t", table, "-H", *args])
    assert result.exit_code == 0, result.output
    return [line.split("\t") for line in result.output.strip().splitlines()]


# ---------------------------------------------------------------------------
# The dumps
# ---------------------------------------------------------------------------


def test_every_dump_is_keyed_by_name(plain: Path):
    write_clusters(
        plain,
        ["r"],
        pd.DataFrame({"round": ["r"], "name": ["K"]}),
        pd.DataFrame({"round": ["r"], "cluster": ["K"], "contig": ["c1"]}),
    )
    contigs = _dump(plain, "contigs")
    assert contigs[0] == [
        "contig",
        "genome",
        "length",
        "gc",
        "qc_intra_contig_short",
        "qc_intra_contig_long",
    ]
    assert [row[:3] for row in contigs[1:]] == [
        ["c1", "MAG_1", "1000"],
        ["c2", "MAG_2", "500"],
    ]
    genomes = _dump(plain, "genomes")
    assert genomes[0][:4] == ["genome", "length", "n_contigs", "n_bins"]
    assert [row[0] for row in genomes[1:]] == ["MAG_1", "MAG_2"]
    assert _dump(plain, "clusters")[0][:2] == ["round", "cluster"]
    assert _dump(plain, "members") == [["round", "cluster", "contig"], ["r", "K", "c1"]]
    # chroms stays cooler's own view of the contigs.
    assert _dump(plain, "chroms")[0] == ["name", "length"]


def test_the_views_keep_the_rows_and_their_ids(plain: Path):
    mc = MetaCooler(plain)
    view = contigs_view(mc.contigs, mc.genomes["name"])
    assert view.index.tolist() == mc.contigs.index.tolist()
    assert view["genome"].tolist() == ["MAG_1", "MAG_2"]
    bare = contigs_view(mc.contigs[["name", "length", "genome_id"]], mc.genomes["name"])
    assert list(bare.columns) == ["contig", "genome", "length"]
    assert list(genomes_view(mc.genomes).columns)[:2] == ["genome", "length"]


# ---------------------------------------------------------------------------
# contig and genome are reserved
# ---------------------------------------------------------------------------


def test_a_user_column_called_genome_or_contig_is_stored_as_user_(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    """A taxonomy column called genome no longer hides behind the genome a
    contig is in: create stores it as user_genome, in both tables."""
    with caplog.at_level(logging.WARNING, logger="metacooler"):
        path = _map(
            tmp_path,
            "contig_name\tgenome_name\tgenome\tcontig\nc1\tMAG_1\tE.coli\tx\n"
            "c2\tMAG_2\tB.subtilis\ty\n",
            "genome_name\tgenome\nMAG_1\tGCF_1\nMAG_2\tGCF_2\n",
        )
    assert "renamed to 'user_genome', 'user_contig'" in caplog.text
    mc = MetaCooler(path)
    assert mc.contigs["user_genome"].tolist() == ["E.coli", "B.subtilis"]
    assert mc.contigs["user_contig"].tolist() == ["x", "y"]
    assert mc.genomes["user_genome"].tolist() == ["GCF_1", "GCF_2"]
    assert "genome" not in mc.contigs.columns and "genome" not in mc.genomes.columns
    row = _dump(path, "contigs")
    assert row[0][:5] == ["contig", "genome", "length", "user_genome", "user_contig"]
    assert row[1][:5] == ["c1", "MAG_1", "1000", "E.coli", "x"]


def test_annotate_stores_them_as_create_does(plain: Path, tmp_path: Path):
    table = tmp_path / "taxonomy.tsv"
    table.write_text("contig_name\tgenome\nc1\tE.coli\nc2\tB.subtilis\n")
    annotate(plain, contigs_table=table, inplace=True)
    assert MetaCooler(plain).contigs["user_genome"].tolist() == ["E.coli", "B.subtilis"]


def _store(path: Path, table: str, column: str, values: list[str]) -> None:
    """Store a column as a map written before the name was reserved kept it."""
    with h5py.File(path, "r+") as f:
        write_column(f[table], column, np.array(values, dtype=object))


def test_a_column_a_map_stored_before_is_shown_as_user_(plain: Path):
    """dump and --filter both read a stored genome column as user_genome --
    and an existing user_genome as user_user_genome -- so the genome column
    always means the genome a contig is in."""
    _store(plain, "contigs", "genome", ["E.coli", "B.subtilis"])
    _store(plain, "contigs", "user_genome", ["e", "b"])
    _store(plain, "contigs", "contig", ["x", "y"])
    _store(plain, "genomes", "genome", ["GCF_1", "GCF_2"])
    contigs = _dump(plain, "contigs")
    header = contigs[0]
    assert header[:3] == ["contig", "genome", "length"]
    assert {"user_user_genome", "user_genome", "user_contig"} <= set(header)
    row = dict(zip(header, contigs[1], strict=True))
    assert row["genome"] == "MAG_1"
    assert row["user_user_genome"] == "E.coli"
    assert row["user_genome"] == "e"
    assert row["user_contig"] == "x"
    genomes = _dump(plain, "genomes")
    assert genomes[0][0] == "genome" and "user_genome" in genomes[0]

    mc = MetaCooler(plain)
    assert resolve_contig_selection(mc, filters=["contigs.genome == MAG_1"]) == {0}
    assert resolve_contig_selection(
        mc, filters=["contigs.user_user_genome == E.coli"]
    ) == {0}
    assert resolve_contig_selection(mc, filters=["genomes.user_genome == GCF_2"]) == {1}


# ---------------------------------------------------------------------------
# Filters read the dump's names, and their stored ones
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("expression", "kept"),
    [
        ("contigs.contig == c2", {1}),
        ("contigs.genome == MAG_1", {0}),
        ("genomes.genome == MAG_2", {1}),
        # The stored names still work.
        ("contigs.name == c2", {1}),
        ("contigs.genome_id == 0", {0}),
        ("genomes.name == MAG_1", {0}),
    ],
)
def test_filters_read_the_dumps_names(plain: Path, expression: str, kept: set[int]):
    assert resolve_contig_selection(MetaCooler(plain), filters=[expression]) == kept


def test_a_filter_on_genome_means_the_genome_a_contig_is_in(tmp_path: Path):
    """The bug this fixes: dump printed the genome's name in the genome column
    while --filter read a user column of that name."""
    path = _map(
        tmp_path,
        "contig_name\tgenome_name\tgenome\nc1\tMAG_1\tE.coli\nc2\tMAG_2\tB.subtilis\n",
    )
    mc = MetaCooler(path)
    assert resolve_contig_selection(mc, filters=["contigs.genome == MAG_1"]) == {0}
    assert resolve_contig_selection(mc, filters=["contigs.user_genome == E.coli"]) == {
        0
    }
    with pytest.raises(FilterError, match="keeps nothing"):
        resolve_contig_selection(mc, filters=["contigs.genome == E.coli"])


def test_create_filters_on_the_dumps_names(tmp_path: Path):
    """A filter means at creation what it means on the finished map."""
    (tmp_path / "bins.bed").write_text("c1\t0\t1000\nc2\t0\t500\nc3\t0\t300\n")
    (tmp_path / "lib.pairs").write_text(
        "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n"
    )
    (tmp_path / "contigs.tsv").write_text(
        "contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_2\nc3\tMAG_2\n"
    )
    out = tmp_path / "out.metacool"
    from_pairs(
        out,
        bins=tmp_path / "bins.bed",
        pairs=tmp_path / "lib.pairs",
        contig_genome_map=tmp_path / "contigs.tsv",
        ps=False,
        filters=["genomes.genome == MAG_2", "contigs.contig != c3"],
    )
    assert MetaCooler(out).contigs["name"].tolist() == ["c2"]


# ---------------------------------------------------------------------------
# A dump reads back
# ---------------------------------------------------------------------------


def test_a_dump_reads_back_as_an_assignment(
    plain: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    """dump -H names the key columns contig and genome; the default format,
    contig_name and genome_name, reads them when the table has neither."""
    dumped = tmp_path / "dumped.tsv"
    # Through the API: a CLI command would reconfigure the logger caplog reads.
    dump_table(plain, "contigs", header=True, out=dumped)
    out = tmp_path / "out.metacool"
    with caplog.at_level(logging.INFO, logger="metacooler"):
        reassign(plain, dumped, out, missing="unbinned")
    assert "as dump names them: reading those" in caplog.text
    mc = MetaCooler(out)
    assert dict(
        zip(
            mc.contigs["name"],
            mc.contigs["genome_id"].map(mc.genomes["name"]),
            strict=True,
        )
    ) == {
        "c1": "MAG_1",
        "c2": "MAG_2",
    }


def test_a_genomes_dump_and_a_metadata_table_read_back(plain: Path, tmp_path: Path):
    genomes = tmp_path / "genomes.tsv"
    genomes.write_text("genome\tnote\nMAG_1\tfirst\nMAG_2\tsecond\n")
    out = tmp_path / "out.metacool"
    reassign(plain, {"c1": "MAG_1"}, out, genomes_table=genomes)
    assert MetaCooler(out).genomes["note"].tolist() == ["first", "second"]
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("contig\tdepth\nc1\t3\nc2\t4\n")
    annotate(plain, contigs_table=metadata, inplace=True)
    assert MetaCooler(plain).contigs["depth"].tolist() == [3, 4]


def test_a_table_with_one_default_name_is_not_read_as_a_dump(
    plain: Path, tmp_path: Path
):
    """contig_name without genome_name is refused as before: the genome column
    of a dump is not looked for beside a default name, where a column called
    genome may hold something else (a taxonomy)."""
    table = tmp_path / "table.tsv"
    table.write_text("contig_name\tbin\tgenome\nc1\tMAG_9\tE.coli\n")
    with pytest.raises(ValueError, match="no column named 'genome_name'"):
        reassign(plain, table, tmp_path / "out.metacool")


def test_a_dump_reads_back_from_the_cli(plain: Path, tmp_path: Path):
    runner = CliRunner()
    dumped = tmp_path / "dumped.tsv"
    result = runner.invoke(
        cli, ["dump", str(plain), "-t", "contigs", "-H", "-o", str(dumped)]
    )
    assert result.exit_code == 0, result.output
    out = tmp_path / "out.metacool"
    result = runner.invoke(
        cli, ["reassign", str(plain), "--contigs", str(dumped), "-o", str(out)]
    )
    assert result.exit_code == 0, flat(result.output)
    # A headerless dump reads back through contig_genome (columns 1 and 2).
    bare = tmp_path / "bare.tsv"
    result = runner.invoke(cli, ["dump", str(plain), "-t", "contigs", "-o", str(bare)])
    assert result.exit_code == 0, result.output
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(plain),
            "--contigs",
            str(bare),
            "--contigs-format",
            "contig_genome",
            "-o",
            str(tmp_path / "bare.metacool"),
        ],
    )
    assert result.exit_code == 0, flat(result.output)


# ---------------------------------------------------------------------------
# A dump read back keeps none of what dump computed
# ---------------------------------------------------------------------------


def test_a_contigs_dump_reads_back_into_annotate(plain: Path, tmp_path: Path):
    """annotate given a dump stores none of the columns dump computed (length,
    the QC columns), and checks its genome column rather than storing it."""
    dumped = tmp_path / "dumped.tsv"
    dump_table(plain, "contigs", header=True, out=dumped)
    before = list(MetaCooler(plain).contigs.columns)
    annotate(plain, contigs_table=dumped, inplace=True)
    assert list(MetaCooler(plain).contigs.columns) == before
    # The genome column is checked: a dump placing c1 elsewhere is refused.
    dumped.write_text(dumped.read_text().replace("c1\tMAG_1", "c1\tMAG_2"))
    with pytest.raises(ValueError, match="does not belong to in this map"):
        annotate(plain, contigs_table=dumped, inplace=True)


def test_a_contigs_dump_reads_back_into_create(plain: Path, tmp_path: Path):
    dumped = tmp_path / "dumped.tsv"
    dump_table(plain, "contigs", header=True, out=dumped)
    out = tmp_path / "again.metacool"
    from_pairs(
        out,
        bins=tmp_path / "bins.bed",
        pairs=tmp_path / "lib.pairs",
        contig_genome_map=dumped,
        ps=False,
    )
    assert list(MetaCooler(out).contigs.columns) == list(
        MetaCooler(plain).contigs.columns
    )


def test_a_genomes_dump_reads_back(plain: Path, tmp_path: Path):
    notes = tmp_path / "notes.tsv"
    notes.write_text("genome_name\tnote\nMAG_1\tfirst\nMAG_2\tsecond\n")
    annotate(plain, genomes_table=notes, inplace=True)
    dumped = tmp_path / "genomes_dump.tsv"
    dump_table(plain, "genomes", header=True, out=dumped)
    out = tmp_path / "out.metacool"
    reassign(plain, {"c1": "MAG_1"}, out, genomes_table=dumped)
    genomes = MetaCooler(out).genomes
    assert genomes["note"].tolist() == ["first", "second"]
    assert not [c for c in genomes.columns if c.startswith("user_")]


def test_a_table_that_is_not_a_dump_keeps_a_genome_column(plain: Path, tmp_path: Path):
    """Only a dump's genome column is checked: in any other table, a column
    called genome (a taxonomy) is metadata, stored as user_genome."""
    table = tmp_path / "taxonomy.tsv"
    table.write_text("contig\tgenome\tdepth\nc1\tE.coli\t3\nc2\tB.subtilis\t4\n")
    annotate(plain, contigs_table=table, inplace=True)
    contigs = MetaCooler(plain).contigs
    assert contigs["user_genome"].tolist() == ["E.coli", "B.subtilis"]
    assert contigs["depth"].tolist() == [3, 4]


# ---------------------------------------------------------------------------
# annotate reaches a column an older map stores under a key's name
# ---------------------------------------------------------------------------


def test_annotate_drops_an_older_column_by_its_shown_name(plain: Path):
    _store(plain, "contigs", "genome", ["E.coli", "B.subtilis"])
    _store(plain, "genomes", "genome", ["GCF_1", "GCF_2"])
    with pytest.raises(
        ValueError,
        match=r"dropped by the name dump shows it under: contigs\.user_genome",
    ):
        annotate(plain, drop=["contigs.genome"], inplace=True)
    annotate(plain, drop=["contigs.user_genome", "genomes.user_genome"], inplace=True)
    mc = MetaCooler(plain)
    assert "genome" not in mc.contigs.columns and "genome" not in mc.genomes.columns


def test_annotate_replaces_an_older_column_by_its_shown_name(
    plain: Path, tmp_path: Path
):
    """A table column called genome is stored as user_genome, the name dump
    shows the older map's genome column under: it replaces that column, which
    is then stored under the new name."""
    _store(plain, "contigs", "genome", ["E.coli", "B.subtilis"])
    same = tmp_path / "same.tsv"
    same.write_text("contig_name\tgenome\nc1\tE.coli\nc2\tB.subtilis\n")
    annotate(plain, contigs_table=same, inplace=True)  # nothing changes
    assert MetaCooler(plain).contigs["genome"].tolist() == ["E.coli", "B.subtilis"]
    changed = tmp_path / "changed.tsv"
    changed.write_text("contig_name\tgenome\nc1\tE.coli\nc2\tB.sub\n")
    with pytest.raises(ValueError, match="contigs.user_genome"):
        annotate(plain, contigs_table=changed, inplace=True)
    annotate(plain, contigs_table=changed, inplace=True, replace=True)
    contigs = MetaCooler(plain).contigs
    assert "genome" not in contigs.columns
    assert contigs["user_genome"].tolist() == ["E.coli", "B.sub"]
