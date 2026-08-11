"""Tests for P(s) — the engine (metacooler.ps) and the op (metacooler.ops.ps)."""

from __future__ import annotations

import gzip
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat

from metacooler.api import MetaCooler
from metacooler.io._utils import chrom_codes
from metacooler.io.create import from_pairs
from metacooler.ops.ps import contact_probability, plot_ps, ps
from metacooler.ps import (
    PS_COLUMNS,
    PSCounts,
    assemble,
    bin_pair_denominators,
    distance_edges,
    pairs_distance_grid,
    position_pair_denominators,
    ps_table_path,
    tally_intra_contig,
)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _brute_force_bin_distances(mc: MetaCooler) -> np.ndarray:
    """Every intra-contig bin-pair distance, enumerated the slow way.

    Shares no code with the vectorised counter it checks: it walks
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


def _brute_force_position_pairs(lengths, lo: int, hi: int) -> int:
    """Count ordered position pairs at separations ``[lo, hi)``, by enumeration."""
    total = 0
    for length in lengths:
        for s in range(lo, hi):
            total += max(int(length) - s, 0)
    return total


@pytest.fixture
def decay_pairs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A polymer-like library: 3 contigs of 100 kb, exponential cis distances."""
    bed = tmp_path / "bins.bed"
    bed.write_text(
        "".join(
            f"c{i}\t{s}\t{s + 1000}\n"
            for i in range(3)
            for s in range(0, 100_000, 1000)
        )
    )
    contigs = tmp_path / "contigs.tsv"
    contigs.write_text("contig_name\tgenome_name\nc0\tMAG_1\nc1\tMAG_1\nc2\tMAG_2\n")
    rng = np.random.default_rng(7)
    lines = ["#columns: readID chrom1 pos1 chrom2 pos2\n"]
    for i in range(20_000):
        contig = f"c{rng.integers(0, 3)}"
        p1 = int(rng.integers(1, 100_000))
        p2 = int(np.clip(p1 + rng.exponential(3000) * rng.choice([-1, 1]), 1, 99_999))
        lines.append(f"r{i}\t{contig}\t{p1}\t{contig}\t{p2}\n")
    pairs = tmp_path / "decay.pairs"
    pairs.write_text("".join(lines))
    return pairs, bed, contigs


@pytest.fixture
def decay_metacool(tmp_path: Path, decay_pairs) -> tuple[Path, Path]:
    """The decay library as a metacool, with its create-time P(s) sidecar."""
    pairs, bed, contigs = decay_pairs
    out = tmp_path / "decay.metacool"
    from_pairs(
        uri=out,
        bins=bed,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
    )
    return out, pairs


# ----------------------------------------------------------------------
# The engine: distance grid
# ----------------------------------------------------------------------
def test_edges_start_at_zero_and_increase():
    edges = distance_edges(100_000, 1000, 10)
    assert edges[0] == 0
    assert edges[1] == 1000
    assert np.all(np.diff(edges) > 0)
    assert edges[-1] > 100_000


def test_edges_survive_a_range_shorter_than_the_first_bin():
    """A file whose longest contig is one bin still yields a usable grid."""
    assert distance_edges(0, 1000, 10).tolist() == [0, 1]


@pytest.mark.parametrize(
    ("max_distance", "min_distance", "per_decade"),
    [(-1, 1000, 10), (1000, 0, 10), (1000, 1000, 0)],
)
def test_edges_reject_nonsense(max_distance, min_distance, per_decade):
    with pytest.raises(ValueError):
        distance_edges(max_distance, min_distance, per_decade)


def test_pairs_grid_starts_at_one_bp_and_ends_at_the_longest_contig():
    """The pairs path resolves single base pairs; that is the point of it."""
    edges = pairs_distance_grid(np.array([1000, 50_000, 2000]))
    assert edges[0] == 0
    assert edges[1] == 1
    assert edges[-1] > 50_000


# ----------------------------------------------------------------------
# The engine: accumulator and denominators
# ----------------------------------------------------------------------
def test_counts_land_in_the_right_distance_bin():
    edges = np.array([0, 10, 100, 1000], dtype=np.int64)
    tally = PSCounts(n_groups=2, edges=edges)
    tally.tally(
        np.array([0, 0, 1, 1]),
        np.array([0, 50, 5, 999]),
    )
    np.testing.assert_array_equal(tally.matrix, [[1, 1, 0], [1, 0, 1]])


def test_counts_beyond_the_last_edge_are_dropped_not_folded_back():
    edges = np.array([0, 10, 100], dtype=np.int64)
    tally = PSCounts(n_groups=1, edges=edges)
    tally.tally(np.array([0, 0]), np.array([5, 100_000]))
    np.testing.assert_array_equal(tally.matrix, [[1, 0]])


def test_weights_are_pixel_counts():
    edges = np.array([0, 10, 100], dtype=np.int64)
    tally = PSCounts(n_groups=1, edges=edges)
    tally.tally(np.array([0, 0]), np.array([1, 50]), weights=np.array([7, 3]))
    np.testing.assert_array_equal(tally.matrix, [[7, 3]])


def test_tally_intra_contig_skips_trans_and_unselected():
    edges = np.array([0, 10, 100], dtype=np.int64)
    tally = PSCounts(n_groups=1, edges=edges)
    # contig 1 is outside the selection (-1); the second pair is inter-contig.
    tally_intra_contig(
        tally,
        np.array([0, -1]),
        np.array([0, 0, 1]),
        np.array([0, 1, 1]),
        np.array([5, 5, 5]),
    )
    np.testing.assert_array_equal(tally.matrix, [[1, 0]])


def test_position_pair_denominators_match_brute_force():
    """The closed form equals an explicit sum over every separation."""
    lengths = np.array([50, 120, 7])
    edges = np.array([0, 1, 10, 60, 200], dtype=np.int64)
    counts = position_pair_denominators(lengths, np.zeros(3, dtype=np.int64), 1, edges)
    expected = [
        _brute_force_position_pairs(lengths, int(edges[i]), int(edges[i + 1]))
        for i in range(len(edges) - 1)
    ]
    np.testing.assert_allclose(counts[0], expected)


def test_position_pair_denominators_split_by_group():
    lengths = np.array([100, 100])
    edges = np.array([0, 1, 50], dtype=np.int64)
    pooled = position_pair_denominators(lengths, np.zeros(2, dtype=np.int64), 1, edges)
    split = position_pair_denominators(lengths, np.array([0, 1]), 2, edges)
    np.testing.assert_allclose(split.sum(axis=0), pooled[0])
    np.testing.assert_allclose(split[0], split[1])


def test_bin_pair_denominators_match_brute_force(tiny_metacooler):
    """The vectorised bin-pair counter equals an explicit enumeration."""
    mc = MetaCooler(tiny_metacooler)
    bins = mc.bins()
    contig_of_bin = chrom_codes(bins)
    mids = (bins["start"].to_numpy() + bins["end"].to_numpy()) // 2
    binsize = mc.binsize
    assert binsize is not None  # the tiny fixture is a fixed-bin file
    edges = distance_edges(int(mids.max()), int(binsize), 10)
    fold = int(mids.max()) + int(edges[-1]) + 1
    key = contig_of_bin.astype(np.int64) * fold + mids

    counts = bin_pair_denominators(key, np.zeros(len(key), dtype=np.int64), 1, edges)
    distances = _brute_force_bin_distances(mc)
    expected = [
        int(((distances >= edges[i]) & (distances < edges[i + 1])).sum())
        for i in range(len(edges) - 1)
    ]
    np.testing.assert_allclose(counts[0], expected)
    assert counts[0].sum() == len(distances)


# ----------------------------------------------------------------------
# The engine: the output table
# ----------------------------------------------------------------------
def test_assemble_columns_and_derived_values():
    edges = np.array([0, 10, 100], dtype=np.int64)
    observed = np.array([[5.0, 20.0]])
    possible = np.array([[10.0, 40.0]])
    table = assemble(["MAG_1"], edges, observed, possible)

    assert list(table.columns) == list(PS_COLUMNS)
    assert table["genome"].tolist() == ["MAG_1", "MAG_1"]
    np.testing.assert_allclose(table["avg_contacts"], [0.5, 0.5])
    # prob integrates to 1 over the reported range: 0.5*10 + 0.5*90 = 50.
    np.testing.assert_allclose(table["prob"], [0.01, 0.01])
    width = np.array([10, 90])
    assert float((table["prob"] * width).sum()) == pytest.approx(1.0)


def test_assemble_drops_bins_nothing_could_populate():
    edges = np.array([0, 10, 100], dtype=np.int64)
    table = assemble(["a"], edges, np.array([[1.0, 0.0]]), np.array([[2.0, 0.0]]))
    assert len(table) == 1


def test_assemble_leaves_an_empty_curve_as_nan():
    edges = np.array([0, 10], dtype=np.int64)
    table = assemble(["silent"], edges, np.array([[0.0]]), np.array([[5.0]]))
    assert table["avg_contacts"].tolist() == [0.0]
    assert np.isnan(table["prob"]).all()


def test_ps_table_path_sits_beside_the_map():
    assert ps_table_path("/x/gut_T0.metacool") == Path("/x/gut_T0.ps.tsv")
    assert ps_table_path(Path("sample.v2.metacool")) == Path("sample.v2.ps.tsv")


# ----------------------------------------------------------------------
# create pairs: the sidecar table
# ----------------------------------------------------------------------
def test_create_writes_a_per_mag_sidecar(decay_metacool):
    out, _ = decay_metacool
    table = pd.read_csv(ps_table_path(out), sep="\t")
    assert list(table.columns) == list(PS_COLUMNS)
    assert set(table["genome"]) == {"MAG_1", "MAG_2"}
    assert (table["n_bin_pairs"] > 0).all()
    np.testing.assert_allclose(
        table["avg_contacts"], table["n_contacts"] / table["n_bin_pairs"]
    )


def test_sidecar_counts_every_intra_contig_pair(decay_metacool):
    """The whole library is cis, so every pair must appear in the table."""
    out, pairs = decay_metacool
    table = pd.read_csv(ps_table_path(out), sep="\t")
    n_lines = sum(
        1 for line in pairs.read_text().splitlines() if not line.startswith("#")
    )
    assert table["n_contacts"].sum() == n_lines


def test_sidecar_resolves_below_the_bin_size(decay_metacool):
    """The point of measuring on positions: distances a 1 kb grid cannot show."""
    out, _ = decay_metacool
    from_pairs_table = pd.read_csv(ps_table_path(out), sep="\t")
    from_pixels = contact_probability(out, by="genome")

    assert from_pairs_table["dist"].min() < 10  # sub-bin, from read positions
    assert from_pixels["dist"].min() >= 500  # half a 1 kb bin at best
    assert len(from_pairs_table) > len(from_pixels)


def test_no_ps_skips_the_sidecar(tmp_path: Path, decay_pairs):
    pairs, bed, contigs = decay_pairs
    out = tmp_path / "nops.metacool"
    from_pairs(uri=out, bins=bed, pairs=pairs, contig_genome_map=contigs, ps=False)
    assert out.exists()
    assert not ps_table_path(out).exists()


def test_existing_sidecar_is_refused_before_any_work(tmp_path: Path, decay_pairs):
    """The guard runs up front: an hour of binning must not end in a clash."""
    pairs, bed, contigs = decay_pairs
    out = tmp_path / "clash.metacool"
    ps_table_path(out).write_text("keep me")
    with pytest.raises(FileExistsError):
        from_pairs(uri=out, bins=bed, pairs=pairs, contig_genome_map=contigs)
    assert ps_table_path(out).read_text() == "keep me"
    assert not out.exists()  # nothing was written


def test_force_overwrites_the_sidecar(tmp_path: Path, decay_pairs):
    pairs, bed, contigs = decay_pairs
    out = tmp_path / "clash2.metacool"
    ps_table_path(out).write_text("stale")
    from_pairs(uri=out, bins=bed, pairs=pairs, contig_genome_map=contigs, force=True)
    assert pd.read_csv(ps_table_path(out), sep="\t")["n_contacts"].sum() > 0


# ----------------------------------------------------------------------
# ps --pairs
# ----------------------------------------------------------------------
def test_recomputing_from_pairs_reproduces_the_sidecar(decay_metacool):
    """Same pass, same engine, same numbers — the two must not drift."""
    out, pairs = decay_metacool
    sidecar = pd.read_csv(ps_table_path(out), sep="\t")
    recomputed = contact_probability(out, pairs=pairs, by="genome")
    pd.testing.assert_frame_equal(sidecar, recomputed)


def test_gzipped_pairs_read_identically(tmp_path: Path, decay_metacool):
    """Compression is the pairs reader's business, and it is the same reader."""
    out, pairs = decay_metacool
    gz = tmp_path / "decay.pairs.gz"
    gz.write_bytes(gzip.compress(pairs.read_bytes()))
    plain = contact_probability(out, pairs=pairs, by="genome")
    zipped = contact_probability(out, pairs=gz, by="genome")
    pd.testing.assert_frame_equal(plain, zipped)


def test_pairs_path_honours_the_genome_selection(decay_metacool):
    out, pairs = decay_metacool
    one = contact_probability(out, pairs=pairs, genomes=["MAG_2"])
    per_genome = contact_probability(out, pairs=pairs, by="genome")
    mag2 = per_genome[per_genome["genome"] == "MAG_2"]
    assert one["n_contacts"].sum() == mag2["n_contacts"].sum()
    assert one["n_contacts"].sum() < per_genome["n_contacts"].sum()


def test_pairs_path_reports_unrecognised_contigs(
    tmp_path: Path, tiny_metacooler, caplog
):
    """A name mismatch is the failure mode that silently empties the curve."""
    foreign = tmp_path / "foreign.pairs"
    foreign.write_text(
        "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(f"r{i}\tchrZ\t{i + 1}\tchrZ\t{i + 100}\n" for i in range(10))
    )
    table = contact_probability(tiny_metacooler, pairs=foreign)
    assert table["n_contacts"].sum() == 0
    assert "unrecognised contigs" in caplog.text


def test_chunksize_does_not_change_the_pairs_answer(decay_metacool):
    out, pairs = decay_metacool
    whole = contact_probability(out, pairs=pairs, by="genome")
    chunked = contact_probability(out, pairs=pairs, by="genome", chunksize=1000)
    pd.testing.assert_frame_equal(whole, chunked)


def test_pairs_curve_decays(decay_metacool):
    """An exponentially distributed library must fall off with distance."""
    out, pairs = decay_metacool
    table = contact_probability(out, pairs=pairs)
    near = table[(table["dist"] > 500) & (table["dist"] < 2000)]["avg_contacts"].mean()
    far = table[table["dist"] > 30_000]["avg_contacts"].mean()
    assert near > 10 * far


# ----------------------------------------------------------------------
# The pixel path
# ----------------------------------------------------------------------
def test_pixel_path_counts_the_intra_contig_pixels(tiny_metacooler):
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
    assert table["n_contacts"].sum() < pixels["count"].sum()


def test_pixel_diagonal_is_counted(tiny_metacooler):
    """The leading bin holds one self-pair per bin: no bin is left out."""
    mc = MetaCooler(tiny_metacooler)
    table = contact_probability(tiny_metacooler)
    assert table.iloc[0]["n_bin_pairs"] == len(mc.bins())


def test_by_genome_partitions_the_pooled_curve(tiny_metacooler):
    pooled = contact_probability(tiny_metacooler, by="all")
    per_genome = contact_probability(tiny_metacooler, by="genome")
    assert set(per_genome["genome"]) == set(MetaCooler(tiny_metacooler).genomes["name"])
    assert per_genome["n_contacts"].sum() == pooled["n_contacts"].sum()
    assert per_genome["n_bin_pairs"].sum() == pooled["n_bin_pairs"].sum()


def test_columns_and_types(tiny_metacooler):
    table = contact_probability(tiny_metacooler)
    assert list(table.columns) == list(PS_COLUMNS)
    assert table["n_contacts"].dtype == np.int64
    assert table["n_bin_pairs"].dtype == np.int64
    assert (table["n_bin_pairs"] > 0).all()


def test_pixel_chunksize_does_not_change_the_answer(tiny_metacooler):
    whole = contact_probability(tiny_metacooler)
    chunked = contact_probability(tiny_metacooler, chunksize=2)
    pd.testing.assert_frame_equal(whole, chunked)


def test_filter_narrows_the_selection(tiny_metacooler):
    narrowed = contact_probability(tiny_metacooler, filters=["contigs.length >= 600"])
    whole = contact_probability(tiny_metacooler)
    assert narrowed["n_bin_pairs"].sum() < whole["n_bin_pairs"].sum()


def test_unknown_grouping_is_rejected(tiny_metacooler):
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
    """matplotlib hides labels starting with '_', and ``__unbinned__`` is one."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.axes import Axes

    table = contact_probability(tiny_metacooler, by="genome")
    assert "__unbinned__" in set(table["genome"])

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
    # rstrip("\r"): pandas writes "\n", but a text stdout on Windows translates
    # it, so the last field of every line would otherwise carry a stray CR.
    assert header.rstrip("\r").split("\t") == list(PS_COLUMNS)
    assert rows


def test_cli_ps_from_pairs_matches_the_sidecar(tmp_path: Path, decay_metacool):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out, pairs = decay_metacool
    recomputed = tmp_path / "recomputed.tsv"
    result = CliRunner().invoke(
        cli,
        [
            "ps",
            str(out),
            "--pairs",
            str(pairs),
            "--by",
            "genome",
            "-o",
            str(recomputed),
            "-q",
        ],
    )
    assert result.exit_code == 0, result.output
    assert recomputed.read_text() == ps_table_path(out).read_text()


def test_cli_create_pairs_no_ps(tmp_path: Path, decay_pairs):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    pairs, bed, contigs = decay_pairs
    out = tmp_path / "cli_nops.metacool"
    result = CliRunner().invoke(
        cli,
        [
            "create",
            "pairs",
            "-o",
            str(out),
            "--bins",
            str(bed),
            "--pairs",
            str(pairs),
            "--contigs",
            str(contigs),
            "--no-ps",
            "-q",
        ],
    )
    assert result.exit_code == 0, result.output
    assert not ps_table_path(out).exists()


def test_cli_ps_rejects_a_bad_distance(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(
        cli, ["ps", str(tiny_metacooler), "--min-distance", "far", "-q"]
    )
    assert result.exit_code != 0
    assert "not a number" in _flat(result.output)
