"""Tests for P(s) — contact probability vs genomic distance (metacooler.ops.ps)."""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metacooler.api import MetaCooler
from metacooler.io._utils import chrom_codes
from metacooler.ops.ps import (
    PS_COLUMNS,
    _distance_edges,
    contact_probability,
    plot_ps,
    ps,
)


def _brute_force_distances(mc: MetaCooler) -> np.ndarray:
    """Every intra-contig bin-pair distance, enumerated the slow way.

    Shares no code with the vectorised counter it is used to check: it walks
    ``itertools.combinations_with_replacement`` over each contig's bins, which
    is only affordable on a tiny file — the point of having one.
    """
    bins = mc.bins()
    contig_of_bin = chrom_codes(bins)
    mids = (bins["start"].to_numpy() + bins["end"].to_numpy()) // 2
    distances = []
    for contig in np.unique(contig_of_bin):
        members = np.where(contig_of_bin == contig)[0]
        for i, j in itertools.combinations_with_replacement(members, 2):
            distances.append(abs(int(mids[j]) - int(mids[i])))
    return np.array(distances, dtype=np.int64)


# ----------------------------------------------------------------------
# Distance grid
# ----------------------------------------------------------------------
def test_edges_start_at_zero_and_increase():
    edges = _distance_edges(100_000, 1000, 10)
    assert edges[0] == 0
    assert edges[1] == 1000
    assert np.all(np.diff(edges) > 0)
    assert edges[-1] > 100_000


def test_edges_survive_a_range_shorter_than_the_first_bin():
    """A file whose longest contig is one bin still yields a usable grid."""
    edges = _distance_edges(0, 1000, 10)
    assert edges.tolist() == [0, 1]


@pytest.mark.parametrize(
    ("max_distance", "min_distance", "per_decade"),
    [(-1, 1000, 10), (1000, 0, 10), (1000, 1000, 0)],
)
def test_edges_reject_nonsense(max_distance, min_distance, per_decade):
    with pytest.raises(ValueError):
        _distance_edges(max_distance, min_distance, per_decade)


# ----------------------------------------------------------------------
# The counts themselves
# ----------------------------------------------------------------------
def test_bin_pair_counts_match_brute_force(tiny_metacooler):
    """The vectorised denominator equals an explicit enumeration of every pair."""
    table = contact_probability(tiny_metacooler)
    distances = _brute_force_distances(MetaCooler(tiny_metacooler))

    for _, row in table.iterrows():
        inside = (distances >= row["dist_min"]) & (distances < row["dist_max"])
        assert row["n_bin_pairs"] == int(inside.sum()), row.to_dict()

    # And nothing is silently left out: every pair falls in a reported bin.
    assert table["n_bin_pairs"].sum() == len(distances)


def test_diagonal_is_counted(tiny_metacooler):
    """The [0, binsize) bin holds one self-pair per bin: no bin is left out."""
    mc = MetaCooler(tiny_metacooler)
    table = contact_probability(tiny_metacooler)
    first = table.iloc[0]
    assert first["dist_min"] == 0
    # Pairs at separation 0 (one per bin) plus any pair whose midpoints are
    # closer than the bin size — with a uniform grid there are none.
    assert first["n_bin_pairs"] == len(mc.bins())


def test_observed_counts_are_the_intra_contig_pixels(tiny_metacooler):
    """Every cis contact is counted once, and no trans contact leaks in."""
    mc = MetaCooler(tiny_metacooler)
    pixels = mc.pixels()
    contig_of_bin = chrom_codes(mc.bins())
    cis = pixels[
        contig_of_bin[pixels["bin1_id"].to_numpy()]
        == contig_of_bin[pixels["bin2_id"].to_numpy()]
    ]
    table = contact_probability(tiny_metacooler)
    assert table["n_contacts"].sum() == cis["count"].sum()
    assert table["n_contacts"].sum() < pixels["count"].sum()  # trans excluded


def test_avg_and_prob_are_derived_consistently(tiny_metacooler):
    table = contact_probability(tiny_metacooler)
    np.testing.assert_allclose(
        table["avg_contacts"], table["n_contacts"] / table["n_bin_pairs"]
    )
    # prob integrates to 1 over the reported range, per curve.
    width = table["dist_max"] - table["dist_min"]
    assert (table["prob"] * width).sum() == pytest.approx(1.0)


def test_columns_and_types(tiny_metacooler):
    table = contact_probability(tiny_metacooler)
    assert list(table.columns) == list(PS_COLUMNS)
    assert table["n_contacts"].dtype == np.int64
    assert table["n_bin_pairs"].dtype == np.int64
    assert (table["n_bin_pairs"] > 0).all()


def test_chunksize_does_not_change_the_answer(tiny_metacooler):
    whole = contact_probability(tiny_metacooler)
    chunked = contact_probability(tiny_metacooler, chunksize=2)
    pd.testing.assert_frame_equal(whole, chunked)


# ----------------------------------------------------------------------
# Grouping and selection
# ----------------------------------------------------------------------
def test_by_genome_partitions_the_pooled_curve(tiny_metacooler):
    """Summing the per-genome curves gives back the pooled one."""
    pooled = contact_probability(tiny_metacooler, by="all")
    per_genome = contact_probability(tiny_metacooler, by="genome")
    assert set(per_genome["group"]) == set(MetaCooler(tiny_metacooler).genomes["name"])
    assert per_genome["n_contacts"].sum() == pooled["n_contacts"].sum()
    assert per_genome["n_bin_pairs"].sum() == pooled["n_bin_pairs"].sum()


def test_genome_selection_restricts_both_counts(tiny_metacooler):
    whole = contact_probability(tiny_metacooler)
    one = contact_probability(tiny_metacooler, genomes=["MAG_1"])
    assert one["n_bin_pairs"].sum() < whole["n_bin_pairs"].sum()
    by_genome = contact_probability(tiny_metacooler, by="genome")
    mag1 = by_genome[by_genome["group"] == "MAG_1"]
    assert one["n_contacts"].sum() == mag1["n_contacts"].sum()


def test_filter_narrows_the_selection(tiny_metacooler):
    table = contact_probability(tiny_metacooler, filters=["contigs.length >= 600"])
    # contig_A (1000) and contig_B (600) survive; C (400) and D (200) do not.
    assert (
        table["n_bin_pairs"].sum()
        < contact_probability(tiny_metacooler)["n_bin_pairs"].sum()
    )


def test_unknown_grouping_is_rejected(tiny_metacooler):
    # A caller reaching the API from the CLI cannot get here (click validates
    # the choice), so the value is deliberately typed wrong on purpose.
    grouping: str = "contig"
    with pytest.raises(ValueError, match="by must be"):
        contact_probability(tiny_metacooler, by=grouping)  # ty:ignore[invalid-argument-type]


def test_empty_selection_is_rejected(tiny_metacooler):
    from metacooler.filters import FilterError

    with pytest.raises((ValueError, FilterError)):
        contact_probability(tiny_metacooler, filters=["contigs.length > 100Mb"])


def test_variable_bins_are_supported(tiny_variable_bin_metacooler):
    """Distances come from midpoints, so a fragment-level file works too."""
    table = contact_probability(tiny_variable_bin_metacooler)
    assert len(table) > 0
    assert (table["n_bin_pairs"] > 0).all()


def test_decay_is_visible_on_a_polymer_like_library(tmp_path: Path):
    """A library with exponentially distributed cis distances decays."""
    lengths = {f"c{i}": 100_000 for i in range(3)}
    bed = tmp_path / "bins.bed"
    bed.write_text(
        "".join(
            f"{c}\t{s}\t{s + 1000}\n"
            for c, ln in lengths.items()
            for s in range(0, ln, 1000)
        )
    )
    contigs = tmp_path / "contigs.tsv"
    contigs.write_text(
        "contig_name\tgenome_name\n" + "".join(f"c{i}\tMAG_1\n" for i in range(3))
    )

    rng = np.random.default_rng(7)
    lines = ["#columns: readID chrom1 pos1 chrom2 pos2\n"]
    for i in range(20_000):
        contig = f"c{rng.integers(0, 3)}"
        p1 = int(rng.integers(1, 100_000))
        p2 = int(np.clip(p1 + rng.exponential(3000) * rng.choice([-1, 1]), 1, 99_999))
        lines.append(f"r{i}\t{contig}\t{p1}\t{contig}\t{p2}\n")
    pairs = tmp_path / "decay.pairs"
    pairs.write_text("".join(lines))

    from metacooler.io.create import from_pairs

    out = tmp_path / "decay.metacool"
    from_pairs(
        uri=out,
        bins=bed,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
    )

    table = contact_probability(out)
    tail = table[table["dist_min"] >= 10_000]
    head = table[table["dist_max"] <= 10_000]
    assert head["avg_contacts"].mean() > tail["avg_contacts"].mean()
    # Contact probability at 1 kb should exceed that at 50 kb by a wide margin.
    near = table[table["dist_min"] < 2000]["avg_contacts"].max()
    far = table[table["dist_min"] > 50_000]["avg_contacts"].max()
    assert near > 10 * far


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------
def test_ps_writes_table_and_plot(tmp_path: Path, tiny_metacooler):
    table_path = tmp_path / "ps.tsv"
    plot_path = tmp_path / "ps.png"
    table = ps(tiny_metacooler, out=table_path, plot=plot_path, by="genome")
    written = pd.read_csv(table_path, sep="\t")
    assert list(written.columns) == list(PS_COLUMNS)
    assert len(written) == len(table)
    assert plot_path.stat().st_size > 0


def test_ps_refuses_to_clobber(tmp_path: Path, tiny_metacooler):
    table_path = tmp_path / "ps.tsv"
    table_path.write_text("keep me")
    with pytest.raises(FileExistsError):
        ps(tiny_metacooler, out=table_path)
    assert table_path.read_text() == "keep me"


def test_plot_rejects_an_unknown_value_column(tmp_path: Path, tiny_metacooler):
    table = contact_probability(tiny_metacooler)
    with pytest.raises(ValueError, match="not a column"):
        plot_ps(table, tmp_path / "x.png", value="nope")


def test_plot_labels_the_unbinned_genome(tmp_path: Path, tiny_metacooler, monkeypatch):
    """matplotlib hides labels starting with '_', and ``__unbinned__`` is one.

    Left to matplotlib's automatic legend, the curve of every unassigned
    contig would be drawn in its own colour with nothing to identify it.
    """
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.axes import Axes

    table = contact_probability(tiny_metacooler, by="genome")
    assert "__unbinned__" in set(table["group"])

    seen: list[list[str]] = []
    original = Axes.legend

    def spy(self, *args, **kwargs):
        if args and isinstance(args[-1], list):
            seen.append([str(label) for label in args[-1]])
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "legend", spy)
    plot_ps(table, tmp_path / "labelled.png")
    assert seen and "__unbinned__" in seen[0]


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def test_cli_ps_stdout(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(cli, ["ps", str(tiny_metacooler), "-q"])
    assert result.exit_code == 0, result.output
    header, *rows = result.output.strip().splitlines()
    assert header.split("\t") == list(PS_COLUMNS)
    assert rows


def test_cli_ps_by_genome_with_plot(tmp_path: Path, tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out = tmp_path / "ps.tsv"
    plot = tmp_path / "ps.pdf"
    result = CliRunner().invoke(
        cli,
        [
            "ps",
            str(tiny_metacooler),
            "--by",
            "genome",
            "-o",
            str(out),
            "--plot",
            str(plot),
            "--min-distance",
            "100bp",
            "-q",
        ],
    )
    assert result.exit_code == 0, result.output
    assert plot.stat().st_size > 0
    assert set(pd.read_csv(out, sep="\t")["group"]) <= set(
        MetaCooler(tiny_metacooler).genomes["name"]
    )


def test_cli_ps_rejects_a_bad_distance(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(
        cli, ["ps", str(tiny_metacooler), "--min-distance", "far", "-q"]
    )
    assert result.exit_code != 0
    assert "not a number" in result.output
