"""Tests for metacooler.ops.show."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from metacooler.ops.show import (
    _format_bp,
    _join_names,
    _plural,
    show,
)


@pytest.fixture(scope="module")
def metacool(tmp_path_factory) -> Path:
    """A tiny metacool built **once** for the whole module.

    Same geometry as the ``tiny_metacooler`` fixture — MAG_1 = contig_A (1000 bp)
    + contig_B (600 bp), MAG_2 = contig_C (400 bp), ``__unbinned__`` = contig_D
    (200 bp), 100 bp bins — but module-scoped: ``show`` only ever reads its
    input, so rebuilding the file per test would pay ``from_pairs`` +
    ``cooler.create`` a few dozen times for nothing.
    """
    from metacooler.io.create import from_pairs

    tmp = tmp_path_factory.mktemp("show")
    lengths = {"contig_A": 1000, "contig_B": 600, "contig_C": 400, "contig_D": 200}

    bins_bed = tmp / "bins.bed"
    with open(bins_bed, "w") as f:
        for contig, length in lengths.items():
            for start in range(0, length, 100):
                f.write(f"{contig}\t{start}\t{min(start + 100, length)}\n")

    contigs_table = tmp / "contigs.tsv"
    with open(contigs_table, "w") as f:
        f.write("contig_name\tgenome_name\n")
        f.write("contig_A\tMAG_1\ncontig_B\tMAG_1\ncontig_C\tMAG_2\n")

    pairs = tmp / "pairs.pairs"
    rng = np.random.default_rng(42)
    names = list(lengths)
    with open(pairs, "w") as f:
        f.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
        for i in range(200):
            c1, c2 = rng.choice(names), rng.choice(names)
            p1 = rng.integers(1, lengths[c1])
            p2 = rng.integers(1, lengths[c2])
            f.write(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

    out = tmp / "show.metacool"
    from_pairs(
        uri=out,
        bins=bins_bed,
        pairs=pairs,
        contig_genome_map=contigs_table,
        contig_table_format="contigs",
        assembly_name="tiny_test",
        contig_ordering="genome-then-alpha",
        chunksize=1000,
    )
    return out


@pytest.fixture
def rendered(monkeypatch):
    """Keep the figure ``show`` renders so a test can inspect it.

    ``show`` closes its figure before returning, which is right for a
    long-running process but leaves nothing to assert on. Intercepting the close
    is the only hook that does not change what ``show`` actually draws; the
    figure is released at teardown so at most one is ever held.
    """
    figures: list = []
    # Bind the real close *before* patching: ``metacooler.ops.show.plt`` is the
    # pyplot module itself, so the patch is global, and it outlives this
    # fixture's teardown (monkeypatch unwinds later). Calling ``plt.close`` here
    # would re-enter the collector and append while we iterate — an unbounded
    # loop, not a leak that shows up as a test failure.
    real_close = plt.close
    monkeypatch.setattr(plt, "close", figures.append)
    yield figures
    for fig in list(figures):
        real_close(fig)


def _panels(figures: list) -> list:
    """Return the heatmap axes of the single figure in *figures*."""
    assert len(figures) == 1, f"expected one figure, got {len(figures)}"
    return [ax for ax in figures[0].axes if ax.images]


class TestShow:
    def test_show_genome(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="genome", out=out, dpi=50)
        assert out.exists()

    def test_show_contig(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="contig", out=out, dpi=50)
        assert out.exists()

    def test_show_bin_requires_a_selection(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="--genomes, --contigs, or --range"):
            show(metacool, level="bin", out=out, dpi=50)

    def test_show_bin_genome(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="bin", genomes=["MAG_1"], out=out, dpi=50)
        assert out.exists()

    def test_show_invalid_level(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="Unknown level"):
            show(metacool, level="invalid", out=out, dpi=50)  # type: ignore

    def test_show_invalid_transform(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        with pytest.raises(ValueError, match="Unknown transform"):
            show(metacool, level="genome", transform="invalid", out=out, dpi=50)  # type: ignore

    def test_show_range1_contig(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="bin", range1="contig_A", out=out, dpi=50)
        assert out.exists()

    def test_show_range1_genome(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="bin", range1="MAG_1", out=out, dpi=50)
        assert out.exists()

    def test_show_range1_with_coords(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="bin", range1="contig_A:0-500", out=out, dpi=50)
        assert out.exists()

    def test_show_range1_range2_inter_genome(self, metacool, tmp_path):
        out = tmp_path / "show.png"
        show(metacool, level="bin", range1="MAG_1", range2="MAG_2", out=out, dpi=50)
        assert out.exists()


class TestGenomesSelection:
    """``--genomes`` puts every name on both axes, whatever the count."""

    def test_two_genomes_are_symmetric_like_any_other_count(
        self, metacool, tmp_path, rendered
    ):
        """Two names are not a special case: both go on both axes."""
        show(
            metacool,
            level="contig",
            genomes=["MAG_1", "MAG_2"],
            out=tmp_path / "pair.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        # MAG_1 has two contigs, MAG_2 one: 3x3, holding both intra blocks and
        # the inter block — not a 2x1 rows-vs-columns slice.
        mat = ax.images[0].get_array()
        assert mat.shape == (3, 3)
        assert np.allclose(mat, mat.T)
        assert ax.get_title() == "2 genomes: MAG_1, MAG_2 — contig level"
        assert ax.get_xlabel() == ax.get_ylabel() == "2 genomes (3 contigs)"

    def test_two_genomes_at_bin_level_are_symmetric_too(
        self, metacool, tmp_path, rendered
    ):
        show(
            metacool,
            level="bin",
            genomes=["MAG_1", "MAG_2"],
            out=tmp_path / "pair_bin.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        # 16 bins (MAG_1) + 4 (MAG_2) on both axes.
        mat = ax.images[0].get_array()
        assert mat.shape == (20, 20)
        assert np.allclose(mat, mat.T)

    def test_single_name_is_the_intra_genome_view(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="contig",
            genomes=["MAG_1"],
            out=tmp_path / "one.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert ax.images[0].get_array().shape == (2, 2)
        assert ax.get_title() == "MAG_1 — contig level"

    def test_three_genomes_are_one_symmetric_matrix(self, metacool, tmp_path, rendered):
        """All selected genomes against all of them, not one subplot each."""
        names = ["MAG_1", "MAG_2", "__unbinned__"]
        show(
            metacool,
            level="contig",
            genomes=names,
            out=tmp_path / "multi.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        # MAG_1 has 2 contigs, MAG_2 and __unbinned__ one each.
        mat = ax.images[0].get_array()
        assert mat.shape == (4, 4)
        assert np.allclose(mat, mat.T)
        assert ax.get_title() == "3 genomes: MAG_1, MAG_2, __unbinned__ — contig level"
        assert ax.get_xlabel() == "3 genomes (4 contigs)"
        assert ax.get_ylabel() == "3 genomes (4 contigs)"

    def test_the_block_order_follows_the_given_order(
        self, metacool, tmp_path, rendered
    ):
        """The user controls the layout: blocks appear in the order requested."""
        show(
            metacool,
            level="contig",
            genomes=["MAG_2", "__unbinned__", "MAG_1"],
            out=tmp_path / "order.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert [t.get_text() for t in ax.get_yticklabels()] == [
            "contig_C",  # MAG_2
            "contig_D",  # __unbinned__
            "contig_A",  # MAG_1
            "contig_B",
        ]

    def test_genome_blocks_get_a_heavier_boundary_line(
        self, metacool, tmp_path, rendered
    ):
        """At bin level a genome block ends where a contig does, but not vice versa."""
        show(
            metacool,
            level="bin",
            genomes=["MAG_1", "MAG_2", "__unbinned__"],
            out=tmp_path / "seps.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        by_width: dict[float, set[float]] = {}
        for line in ax.lines:
            # Keep the axhlines: they span the full x axis in axes fractions, so
            # their y data is the constant row position. axvlines are the mirror.
            if tuple(line.get_xdata()) != (0, 1):
                continue
            by_width.setdefault(line.get_linewidth(), set()).update(line.get_ydata())
        light, heavy = (by_width[w] for w in sorted(by_width))
        # MAG_1 = contig_A (10 bins) + contig_B (6); MAG_2 = 4; __unbinned__ = 2.
        assert heavy == {15.5, 19.5}, "genome blocks"
        assert light == {9.5}, "the contig boundary inside MAG_1, drawn once"

    def test_three_genomes_at_bin_level(self, metacool, tmp_path, rendered):
        """Bins of every selected genome, concatenated on both axes."""
        show(
            metacool,
            level="bin",
            genomes=["MAG_1", "MAG_2", "__unbinned__"],
            out=tmp_path / "multi_bin.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        # 100 bp bins over 1000 + 600 (MAG_1) + 400 (MAG_2) + 200 bp.
        mat = ax.images[0].get_array()
        assert mat.shape == (22, 22)
        assert np.allclose(mat, mat.T)
        # Too many contigs to name; the ticks name the genome blocks instead.
        assert [t.get_text() for t in ax.get_yticklabels()] == [
            "MAG_1",
            "MAG_2",
            "__unbinned__",
        ]
        assert ax.get_xlabel() == "3 genomes — 4 contigs (22 bins)"

    def test_genome_level_multi_is_the_n_by_n_block(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="genome",
            genomes=["MAG_1", "MAG_2", "__unbinned__"],
            out=tmp_path / "gmulti.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        mat = ax.images[0].get_array()
        assert mat.shape == (3, 3)
        assert np.allclose(mat, mat.T)
        assert ax.get_title() == "3 genomes: MAG_1, MAG_2, __unbinned__ — genome level"
        # The genomes *are* the rows here, so the label says it once.
        assert ax.get_ylabel() == "3 genomes"

    def test_duplicate_names_are_rejected(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="more than"):
            show(
                metacool,
                level="contig",
                genomes=["MAG_1", "MAG_2", "MAG_1"],
                out=tmp_path / "dup.png",
                dpi=50,
            )

    def test_empty_selection_is_rejected(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="empty"):
            show(
                metacool,
                level="contig",
                genomes=[],
                out=tmp_path / "empty.png",
                dpi=50,
            )

    def test_unknown_name_lists_the_available_genomes(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="Unknown genome"):
            show(
                metacool,
                level="contig",
                genomes=["MAG_1", "NOPE"],
                out=tmp_path / "bad.png",
                dpi=50,
            )


class TestContigsSelection:
    """``--contigs`` is ``--genomes`` in contig space: same symmetric rule."""

    def test_contigs_at_contig_level(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="contig",
            contigs=["contig_C", "contig_A"],
            out=tmp_path / "c.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        mat = ax.images[0].get_array()
        assert mat.shape == (2, 2)
        assert np.allclose(mat, mat.T)
        # One row per selected contig, in the order given.
        assert [t.get_text() for t in ax.get_yticklabels()] == ["contig_C", "contig_A"]
        assert ax.get_title() == "2 contigs: contig_C, contig_A — contig level"
        # "2 contigs (2 contigs)" would say it twice.
        assert ax.get_ylabel() == "2 contigs"

    def test_contigs_at_bin_level(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="bin",
            contigs=["contig_A", "contig_C"],
            out=tmp_path / "cb.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        mat = ax.images[0].get_array()
        assert mat.shape == (14, 14)  # contig_A 10 bins + contig_C 4
        assert np.allclose(mat, mat.T)
        assert [t.get_text() for t in ax.get_yticklabels()] == ["contig_A", "contig_C"]
        assert ax.get_ylabel() == "2 contigs (14 bins)"

    def test_one_contig_at_bin_level_is_ticked_in_bp(
        self, metacool, tmp_path, rendered
    ):
        show(
            metacool,
            level="bin",
            contigs=["contig_C"],
            out=tmp_path / "c1.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert ax.get_title() == "contig_C — bin level"
        assert ax.get_ylabel() == "contig_C — position (bp)"
        assert ax.get_yticklabels()[0].get_text() == "0 bp"

    def test_genome_level_is_upgraded_to_contig(self, metacool, tmp_path, rendered):
        """A contig has no row of its own at genome level, so the level moves."""
        show(
            metacool,
            level="genome",
            contigs=["contig_A", "contig_C"],
            out=tmp_path / "cg.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert ax.get_title().endswith("— contig level")

    def test_contigs_and_genomes_are_mutually_exclusive(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="mutually exclusive"):
            show(
                metacool,
                level="contig",
                genomes=["MAG_1"],
                contigs=["contig_C"],
                out=tmp_path / "x.png",
                dpi=50,
            )

    def test_unknown_contig_is_reported_as_a_contig(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="Unknown contig: 'NOPE'"):
            show(
                metacool,
                level="contig",
                contigs=["NOPE"],
                out=tmp_path / "x.png",
                dpi=50,
            )

    def test_contigs_cannot_be_combined_with_a_range(self, metacool, tmp_path):
        with pytest.raises(ValueError, match="cannot be combined"):
            show(
                metacool,
                level="bin",
                contigs=["contig_A"],
                range1="contig_C",
                out=tmp_path / "x.png",
                dpi=50,
            )


class TestTitlesAndLabels:
    """Every figure names what it shows and at which resolution."""

    @pytest.mark.parametrize(
        ("level", "expected_axis"), [("genome", "3 genomes"), ("contig", "4 contigs")]
    )
    def test_unscoped_title_and_labels(
        self, metacool, tmp_path, rendered, level, expected_axis
    ):
        show(metacool, level=level, out=tmp_path / f"{level}.png", dpi=50)
        (ax,) = _panels(rendered)
        assert ax.get_title() == f"All genomes — {level} level"
        assert ax.get_xlabel() == expected_axis
        assert ax.get_ylabel() == expected_axis

    def test_scoped_label_names_the_genome_and_counts_rows(
        self, metacool, tmp_path, rendered
    ):
        show(
            metacool, level="contig", genomes=["MAG_2"], out=tmp_path / "s.png", dpi=50
        )
        (ax,) = _panels(rendered)
        # Singular, not "1 contigs".
        assert ax.get_ylabel() == "MAG_2 (1 contig)"

    def test_bin_level_title_names_the_range(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="bin",
            range1="contig_A:200-500",
            out=tmp_path / "r.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert ax.get_title() == "contig_A:200-500 — bin level"

    def test_asymmetric_range_title_names_both_axes(self, metacool, tmp_path, rendered):
        show(
            metacool,
            level="bin",
            range1="contig_A",
            range2="contig_C",
            out=tmp_path / "r2.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        assert ax.get_title() == "contig_A vs contig_C — bin level"


class TestBinAxisTicks:
    """At bin level the axes carry coordinates, not bin indices."""

    def test_range_ticks_are_bp_positions_from_the_range_start(
        self, metacool, tmp_path, rendered
    ):
        """A bin index starts at 0 even when the requested range does not."""
        show(
            metacool,
            level="bin",
            range1="contig_A:200-500",
            out=tmp_path / "r.png",
            dpi=50,
        )
        (ax,) = _panels(rendered)
        labels = [t.get_text() for t in ax.get_xticklabels()]
        # 100 bp bins, so the selection is bins [200,300) … [400,500).
        assert labels == ["200 bp", "300 bp", "400 bp"]
        assert labels == [t.get_text() for t in ax.get_yticklabels()]
        assert ax.get_xlabel() == "contig_A:200-500 — position (bp)"

    def test_multi_contig_axis_is_labelled_by_contig(
        self, metacool, tmp_path, rendered
    ):
        """bp positions restart per contig, so a multi-contig axis names them."""
        show(metacool, level="bin", genomes=["MAG_1"], out=tmp_path / "g.png", dpi=50)
        (ax,) = _panels(rendered)
        assert [t.get_text() for t in ax.get_xticklabels()] == ["contig_A", "contig_B"]
        assert ax.get_xlabel() == "MAG_1 — 2 contigs (16 bins)"


class TestLabelHelpers:
    """The label vocabulary, unit-tested without rendering anything."""

    @pytest.mark.parametrize(
        ("n", "noun", "expected"),
        [(0, "bins", "0 bins"), (1, "bins", "1 bin"), (2, "contigs", "2 contigs")],
    )
    def test_plural(self, n, noun, expected):
        assert _plural(n, noun) == expected

    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            (0, "0 bp"),
            (999, "999 bp"),
            (1_000, "1 kb"),
            (120_000, "120 kb"),
            (1_500_000, "1.5 Mb"),
            (1_234_567, "1.23 Mb"),
        ],
    )
    def test_format_bp(self, pos, expected):
        assert _format_bp(pos) == expected

    def test_join_names_truncates(self):
        assert _join_names(["a", "b"]) == "a, b"
        assert _join_names(list("abcdefg"), limit=3) == "a, b, c (+4 more)"
