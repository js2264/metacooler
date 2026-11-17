"""Tests for contact-quality statistics (metacooler.qc / ops.qc / cmd_qc)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat

import metacooler.qc as qc_engine
from metacooler.api import MetaCooler
from metacooler.io.create import from_cool, from_pairs
from metacooler.ops.qc import qc
from metacooler.qc import (
    QC_CATEGORIES,
    QC_CONTIG_COLUMNS,
    QC_GENOME_COLUMNS,
    QC_METADATA_KEY,
    BinMaps,
    ContactQC,
    strip_qc_columns,
    strip_qc_metadata,
)


def _assert_contig_columns_add_up(path: Path | str) -> None:
    """The /contigs QC columns are the block's intra-contig totals, per contig."""
    mc = MetaCooler(path)
    block = mc.metadata[QC_METADATA_KEY]
    contigs = mc.contigs
    genomes = mc.genomes
    for column, key in zip(
        QC_CONTIG_COLUMNS, ("intra_contig_short", "intra_contig_long"), strict=True
    ):
        assert contigs[column].dtype == np.int64
        assert int(contigs[column].sum()) == block[key], column
        per_genome = np.bincount(
            contigs["genome_id"], weights=contigs[column], minlength=len(genomes)
        )
        np.testing.assert_array_equal(genomes[column].to_numpy(), per_genome)


# ----------------------------------------------------------------------
# Fixtures: a file whose category counts are known exactly by construction
# ----------------------------------------------------------------------
@pytest.fixture
def known_pairs(tmp_path: Path) -> tuple[Path, Path, Path, dict[str, int]]:
    """Pairs whose contact categories are counted by hand.

    Four contigs of 10 kb in 1 kb bins, two MAGs:
    ``MAG_1 = {c1, c2}``, ``MAG_2 = {c3}``, ``c4`` unassigned.
    """
    lengths = {"c1": 10_000, "c2": 10_000, "c3": 10_000, "c4": 10_000}

    bins = tmp_path / "bins.bed"
    bins.write_text(
        "".join(
            f"{c}\t{s}\t{s + 1000}\n"
            for c, ln in lengths.items()
            for s in range(0, ln, 1000)
        )
    )
    contigs = tmp_path / "contigs.tsv"
    contigs.write_text(
        "contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_1\nc3\tMAG_2\n",
    )

    rows = [
        # (chrom1, pos1, chrom2, pos2), one line each
        # 3 intra-contig short (< 10 kb default threshold is everything here,
        # so use an explicit small threshold in the tests that care)
        ("c1", 1_000, "c1", 1_500),  # same contig, 500 bp apart
        ("c1", 2_000, "c1", 2_400),  # same contig, 400 bp
        ("c2", 1_000, "c2", 9_000),  # same contig, 8000 bp
        # 2 inter-contig within MAG_1
        ("c1", 1_000, "c2", 5_000),
        ("c2", 3_000, "c1", 7_000),
        # 4 inter-genome
        ("c1", 1_000, "c3", 1_000),
        ("c2", 2_000, "c3", 2_000),
        ("c3", 3_000, "c4", 3_000),
        ("c1", 4_000, "c4", 4_000),
    ]
    pairs = tmp_path / "known.pairs"
    pairs.write_text(
        "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(
            f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\n" for i, (c1, p1, c2, p2) in enumerate(rows)
        )
    )
    expected = {
        "intra_contig_short": 2,  # 500 bp and 400 bp, at a 1 kb threshold
        "intra_contig_long": 1,  # 8000 bp
        "inter_contig_intra_genome": 2,
        "inter_genome": 4,
    }
    return pairs, bins, contigs, expected


@pytest.fixture
def known_metacool(tmp_path: Path, known_pairs) -> tuple[Path, dict[str, int]]:
    """A metacool built from :func:`known_pairs`, QC'd at a 1 kb threshold."""
    pairs, bins, contigs, expected = known_pairs
    out = tmp_path / "known.metacool"
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        qc_distance_threshold=1000,
    )
    return out, expected


# ----------------------------------------------------------------------
# The engine
# ----------------------------------------------------------------------
def test_categories_partition_every_contact():
    """Every contact lands in exactly one category, so the parts sum to the whole."""
    tally = ContactQC(n_genomes=2, distance_threshold=1000)
    rng = np.random.default_rng(0)
    n = 500
    contig1 = rng.integers(0, 4, n)
    contig2 = rng.integers(0, 4, n)
    genome_of_contig = np.array([0, 0, 1, 1])
    tally.tally(
        contig1,
        contig2,
        genome_of_contig[contig1],
        genome_of_contig[contig2],
        rng.integers(0, 5000, n),
    )
    summary = tally.summary()
    assert summary["n_contacts"] == n
    assert sum(summary[c] for c in QC_CATEGORIES) == n
    assert summary["intra_genome"] + summary["inter_genome"] == n
    assert summary["intra_contig"] + summary["inter_contig"] == n


def test_inter_genome_contact_is_credited_to_both_genomes():
    """A contact between two MAGs counts for each of them, but once overall."""
    tally = ContactQC(n_genomes=2, distance_threshold=1000)
    tally.tally(
        np.array([0]), np.array([1]), np.array([0]), np.array([1]), np.array([0])
    )
    assert tally.summary()["inter_genome"] == 1
    table = tally.genome_table(["A", "B"])
    assert table["qc_inter_genome"].tolist() == [1, 1]
    # Neither genome kept any contact to itself, so both ratios are 0.
    assert table["qc_intra_genome_ratio"].tolist() == [0.0, 0.0]


def test_empty_tally_reports_none_not_zero():
    """A file with no contact has no measured ratio — not a ratio of zero."""
    summary = ContactQC(n_genomes=1, distance_threshold=1000).summary()
    assert summary["n_contacts"] == 0
    assert summary["trans_contig_ratio"] is None
    assert summary["intra_genome_fraction"] is None


def test_weights_are_summed_exactly():
    """Pixel counts are weights, not one contact per row."""
    tally = ContactQC(n_genomes=1, distance_threshold=1000)
    tally.tally(
        np.array([0, 0]),
        np.array([0, 0]),
        np.array([0, 0]),
        np.array([0, 0]),
        np.array([10, 5000]),
        weights=np.array([7, 3], dtype=np.int64),
    )
    summary = tally.summary()
    assert summary["intra_contig_short"] == 7
    assert summary["intra_contig_long"] == 3


def test_genome_table_rejects_mismatched_names():
    with pytest.raises(ValueError, match="expected 2 genome"):
        ContactQC(n_genomes=2).genome_table(["only_one"])


def test_negative_threshold_rejected():
    with pytest.raises(ValueError, match="distance_threshold"):
        ContactQC(n_genomes=1, distance_threshold=-1)


# ----------------------------------------------------------------------
# The fused, blocked tally: the definition's numbers, at any block size
# ----------------------------------------------------------------------
def _by_definition(contig1, contig2, genome1, genome2, distance, weights, unbinned):
    """The four categories spelled out, one masked selection each."""
    w = np.ones(len(contig1), dtype=np.int64) if weights is None else weights
    same_contig = contig1 == contig2
    same_genome = (genome1 == genome2) & (same_contig | (genome1 != unbinned))
    masks = [
        same_contig & (distance < 1000),
        same_contig & (distance >= 1000),
        ~same_contig & same_genome,
        ~same_genome,
    ]
    totals = np.array([w[m].sum() for m in masks])

    def credit(genome: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return np.bincount(genome[mask], weights=w[mask].astype(float), minlength=6)

    per_genome = np.stack([credit(genome1, m) for m in masks], axis=1)
    per_genome[:, 3] += credit(genome2, masks[3])
    return totals, per_genome


def _random_contacts(rng: np.random.Generator, n: int) -> tuple[np.ndarray, ...]:
    """*n* contacts over 30 contigs in 6 genomes, half of them intra-contig."""
    genome_of = rng.integers(0, 6, 30)
    contig1 = rng.integers(0, 30, n)
    contig2 = np.where(rng.random(n) < 0.5, contig1, rng.integers(0, 30, n))
    distance = rng.integers(0, 5000, n)
    return contig1, contig2, genome_of[contig1], genome_of[contig2], distance


class TestBlockedTally:
    @pytest.mark.parametrize("block", [1, 7, 1 << 18])
    @pytest.mark.parametrize("unbinned", [None, 0])
    @pytest.mark.parametrize("weights", [None, "int32", "int64", "float64"])
    def test_matches_the_definition(self, monkeypatch, block, unbinned, weights):
        monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", block)
        rng = np.random.default_rng(block)
        contacts = _random_contacts(rng, 2_000)
        w = None
        if weights == "float64":
            w = rng.random(2_000) * 9
        elif weights is not None:
            w = rng.integers(1, 1_000, 2_000).astype(weights)
        tally = ContactQC(6, distance_threshold=1000, unbinned_genome=unbinned)
        tally.tally(*contacts, weights=w)
        totals, per_genome = _by_definition(
            *contacts, w, -1 if unbinned is None else unbinned
        )
        # Exact, floats included: within a (genome, category) cell the fused
        # bincount adds the same values in the same order as a masked pass.
        np.testing.assert_array_equal(tally.counts, totals)
        np.testing.assert_array_equal(tally.per_genome, per_genome)

    def test_float_bits_do_not_depend_on_the_block(self, monkeypatch):
        """Once float, every chunk is one pass: float addition does not associate."""
        rng = np.random.default_rng(3)
        chunks = [
            (_random_contacts(rng, 1_500), w)
            for w in (
                rng.integers(1, 100, 1_500),
                rng.random(1_500) * 7,
                rng.integers(1, 100, 1_500),  # integers onto float accumulators
                None,
                (rng.random(1_500) * 3).astype(np.float32),
            )
        ]
        tallies = []
        for block in (1, 1 << 18):
            monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", block)
            tally = ContactQC(6, distance_threshold=1000, unbinned_genome=2)
            for contacts, w in chunks:
                tally.tally(*contacts, weights=w)
            tallies.append(tally)
        blocked, whole = tallies
        assert blocked.counts.dtype == np.float64
        assert blocked.counts.tobytes() == whole.counts.tobytes()
        assert blocked.per_genome.tobytes() == whole.per_genome.tobytes()

    def test_a_block_never_undercuts_the_genome_table(self, monkeypatch):
        monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", 2)
        tally = ContactQC(n_genomes=10)
        assert [len(range(100)[s]) for s in tally._blocks(100)] == [40, 40, 20]
        assert list(tally._blocks(0)) == []
        tally._promote_for(np.array([0.5]))
        assert list(tally._blocks(100)) == [slice(0, 100)]

    def test_an_empty_float_chunk_promotes_nothing(self):
        """As before: a chunk with no pixel leaves the counts integers."""
        tally = ContactQC(n_genomes=2)
        none = np.array([], dtype=np.int64)
        maps = BinMaps(none, none, none)
        tally.tally_pixels(none, none, np.array([], dtype=np.float64), maps)
        assert tally.counts.dtype == np.int64
        assert tally.summary()["n_contacts"] == 0

    def test_pixels_are_looked_up_block_by_block(self, monkeypatch):
        monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", 5)
        rng = np.random.default_rng(4)
        # 8 contigs of four 1 kb bins, in 3 genomes.
        contig = np.repeat(np.arange(8), 4)
        start = np.tile(np.arange(0, 4000, 1000), 8)
        genome_of_contig = np.array([0, 0, 1, 1, 2, 2, 0, 1])
        maps = BinMaps(contig, genome_of_contig[contig], start + 500)
        bin1 = rng.integers(0, 32, 500).astype(np.int32)
        bin2 = rng.integers(0, 32, 500)
        counts = rng.integers(1, 9, 500).astype(np.int32)
        by_pixels = ContactQC(3, distance_threshold=1500, unbinned_genome=2)
        by_pixels.tally_pixels(bin1, bin2, counts, maps)
        direct = ContactQC(3, distance_threshold=1500, unbinned_genome=2)
        direct.tally(
            maps.contig[bin1],
            maps.contig[bin2],
            maps.genome[bin1],
            maps.genome[bin2],
            np.abs(maps.mid[bin2] - maps.mid[bin1]),
            counts,
        )
        np.testing.assert_array_equal(by_pixels.counts, direct.counts)
        np.testing.assert_array_equal(by_pixels.per_genome, direct.per_genome)


# ----------------------------------------------------------------------
# Per contig: the two categories a reassignment cannot change
# ----------------------------------------------------------------------
def _per_contig_by_definition(contig1, contig2, distance, weights, n_contigs):
    """Short/long intra-contig counts per contig, one masked selection each."""
    w = np.ones(len(contig1)) if weights is None else weights.astype(float)
    same = contig1 == contig2
    long_range = distance >= 1000
    return np.column_stack(
        [
            np.bincount(contig1[m], weights=w[m], minlength=n_contigs)
            for m in (same & ~long_range, same & long_range)
        ]
    )


class TestPerContigTally:
    @pytest.mark.parametrize("block", [1, 7, 1 << 18])
    @pytest.mark.parametrize("weights", [None, "int32", "float64"])
    def test_matches_the_definition(self, monkeypatch, block, weights):
        monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", block)
        rng = np.random.default_rng(block)
        contacts = _random_contacts(rng, 2_000)
        w = None
        if weights == "float64":
            w = rng.random(2_000) * 9
        elif weights is not None:
            w = rng.integers(1, 1_000, 2_000).astype(weights)
        tally = ContactQC(6, distance_threshold=1000, n_contigs=30)
        tally.tally(*contacts, weights=w)

        expected = _per_contig_by_definition(
            contacts[0], contacts[1], contacts[4], w, 30
        )
        per_contig = tally.per_contig
        assert per_contig is not None
        if weights == "float64":
            np.testing.assert_allclose(per_contig, expected, rtol=1e-12)
        else:
            assert per_contig.dtype == np.int64
            np.testing.assert_array_equal(per_contig, expected)
        # The per-genome cells of both categories are the per-contig counts,
        # summed over each genome's contigs.
        genome_of = np.zeros(30, dtype=np.int64)
        genome_of[contacts[0]] = contacts[2]
        for k in range(2):
            np.testing.assert_allclose(
                tally.per_genome[:, k],
                np.bincount(genome_of, weights=per_contig[:, k], minlength=6),
            )

    def test_the_queue_is_summed_in_batches(self, monkeypatch):
        """Queued records never outgrow the table by more than one block."""
        monkeypatch.setattr(qc_engine, "_TALLY_BLOCK", 4)
        rng = np.random.default_rng(5)
        tally = ContactQC(6, distance_threshold=1000, n_contigs=30)
        tally.tally(*_random_contacts(rng, 2_000))
        assert 0 < tally._contig_queued < 60
        assert tally.per_contig is not None
        assert tally._contig_queued == 0

    def test_a_float_chunk_widens_the_table(self):
        rng = np.random.default_rng(6)
        chunks = [
            (_random_contacts(rng, 500), rng.integers(1, 9, 500)),
            (_random_contacts(rng, 500), rng.random(500) * 3),
        ]
        tally = ContactQC(6, distance_threshold=1000, n_contigs=30)
        expected = np.zeros((30, 2))
        for contacts, w in chunks:
            tally.tally(*contacts, weights=w)
            expected += _per_contig_by_definition(
                contacts[0], contacts[1], contacts[4], w, 30
            )
        assert tally.per_contig is not None
        assert tally.per_contig.dtype == np.float64
        np.testing.assert_allclose(tally.per_contig, expected, rtol=1e-12)

    def test_without_n_contigs_nothing_is_counted_per_contig(self):
        tally = ContactQC(6, distance_threshold=1000)
        tally.tally(*_random_contacts(np.random.default_rng(7), 100))
        assert tally.per_contig is None
        with pytest.raises(ValueError, match="does not count per contig"):
            tally.contig_table(["a"])

    def test_contig_table(self):
        tally = ContactQC(1, distance_threshold=1000, n_contigs=2)
        tally.tally(
            np.array([0, 0, 1]),
            np.array([0, 0, 0]),
            np.zeros(3, dtype=np.int64),
            np.zeros(3, dtype=np.int64),
            np.array([10, 5000, 0]),
        )
        table = tally.contig_table(["c1", "c2"])
        assert list(table.columns) == ["name", *QC_CONTIG_COLUMNS]
        assert table[list(QC_CONTIG_COLUMNS)].to_numpy().tolist() == [[1, 1], [0, 0]]
        with pytest.raises(ValueError, match="expected 2 contig name"):
            tally.contig_table(["c1"])

    def test_set_intra_contig_rebuilds_both_categories(self):
        tally = ContactQC(3, distance_threshold=1000, unbinned_genome=0)
        # Two inter-contig contacts, and one intra-contig one to be replaced.
        tally.tally(
            np.array([0, 1, 3]),
            np.array([2, 2, 3]),
            np.array([0, 0, 1]),
            np.array([2, 2, 1]),
            np.array([0, 0, 0]),
        )
        inter = tally.per_genome[:, 2:].copy()
        per_contig = np.array([[5, 1], [0, 2], [3, 0], [1, 1]])
        tally.set_intra_contig(per_contig, np.array([0, 0, 2, 1]))

        summary = tally.summary()
        assert summary["intra_contig_short"] == 9
        assert summary["intra_contig_long"] == 4
        assert summary["inter_genome"] == 2
        assert summary["n_contacts"] == 15
        assert tally.per_genome[:, 0].tolist() == [5, 1, 3]
        assert tally.per_genome[:, 1].tolist() == [3, 1, 0]
        np.testing.assert_array_equal(tally.per_genome[:, 2:], inter)
        assert tally.per_contig is not None
        np.testing.assert_array_equal(tally.per_contig, per_contig)
        assert tally.per_genome.dtype == np.int64

    def test_set_intra_contig_takes_float_counts(self):
        tally = ContactQC(2, distance_threshold=1000)
        tally.set_intra_contig(np.array([[0.5, 1.0], [2.0, 0.25]]), np.array([1, 1]))
        assert tally.counts.dtype == np.float64
        assert tally.per_genome[:, :2].tolist() == [[0.0, 0.0], [2.5, 1.25]]

    @pytest.mark.parametrize(
        ("per_contig", "genome", "match"),
        [
            (np.zeros((3, 2)), np.zeros(2, dtype=np.int64), "one \\(short, long\\)"),
            (np.zeros((2, 3)), np.zeros(2, dtype=np.int64), "one \\(short, long\\)"),
            (np.zeros((2, 2)), np.array([0, 2]), "genome ids must lie in"),
            (np.zeros((2, 2)), np.array([-1, 0]), "genome ids must lie in"),
        ],
    )
    def test_set_intra_contig_refuses_misaligned_input(self, per_contig, genome, match):
        tally = ContactQC(2, distance_threshold=1000)
        with pytest.raises(ValueError, match=match):
            tally.set_intra_contig(per_contig, genome)


# ----------------------------------------------------------------------
# create pairs: the exact, position-based tally
# ----------------------------------------------------------------------
def test_from_pairs_counts_match_hand_count(known_metacool):
    path, expected = known_metacool
    block = MetaCooler(path).metadata[QC_METADATA_KEY]
    assert block["source"] == "pairs"
    assert block["distance_threshold"] == 1000
    for category, count in expected.items():
        assert block[category] == count, category
    assert block["n_contacts"] == sum(expected.values())
    assert block["trans_contig_ratio"] == pytest.approx(6 / 9, abs=1e-6)
    assert block["n_pairs_read"] == 9
    assert block["n_pairs_discarded"] == 0


def test_from_pairs_writes_genome_columns(known_metacool):
    path, _ = known_metacool
    genomes = MetaCooler(path).genomes
    for column in QC_GENOME_COLUMNS:
        assert column in genomes.columns

    row = genomes.set_index("name")
    # MAG_1 holds c1 and c2: 3 intra-contig, 2 inter-contig cis, and 3 of the
    # 4 inter-genome contacts touch it (c1-c3, c2-c3, c1-c4).
    assert int(row.at["MAG_1", "qc_intra_contig_short"]) == 2
    assert int(row.at["MAG_1", "qc_intra_contig_long"]) == 1
    assert int(row.at["MAG_1", "qc_inter_contig_intra_genome"]) == 2
    assert int(row.at["MAG_1", "qc_inter_genome"]) == 3
    assert row.at["MAG_1", "qc_intra_genome_ratio"] == pytest.approx(5 / 8)


def test_from_pairs_writes_contig_columns(known_metacool):
    """Counted on read positions, contig by contig (see the fixture)."""
    path, _ = known_metacool
    contigs = MetaCooler(path).contigs.set_index("name")[list(QC_CONTIG_COLUMNS)]
    # c1: 500 and 400 bp apart; c2: 8000 bp apart; c3 and c4: none.
    assert contigs.loc[["c1", "c2", "c3", "c4"]].to_numpy().tolist() == [
        [2, 0],
        [0, 1],
        [0, 0],
        [0, 0],
    ]
    _assert_contig_columns_add_up(path)


def test_a_user_column_named_like_a_qc_column_is_renamed(tmp_path: Path, known_pairs):
    """Reserved, as the /genomes ones are, so store_qc cannot overwrite it."""
    pairs, bins, _, _ = known_pairs
    contigs = tmp_path / "contigs_with_qc.tsv"
    contigs.write_text(
        "contig_name\tgenome_name\tqc_intra_contig_short\n"
        "c1\tMAG_1\t99\nc2\tMAG_1\t98\nc3\tMAG_2\t97\n"
    )
    out = tmp_path / "renamed.metacool"
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        qc_distance_threshold=1000,
    )
    table = MetaCooler(out).contigs.set_index("name")
    assert int(table.at["c1", "user_qc_intra_contig_short"]) == 99
    assert int(table.at["c1", "qc_intra_contig_short"]) == 2
    _assert_contig_columns_add_up(out)


def test_qc_can_be_disabled(tmp_path: Path, known_pairs):
    pairs, bins, contigs, _ = known_pairs
    out = tmp_path / "noqc.metacool"
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        qc=False,
    )
    mc = MetaCooler(out)
    assert QC_METADATA_KEY not in mc.metadata
    assert not [c for c in mc.genomes.columns if c.startswith("qc_")]
    assert not [c for c in mc.contigs.columns if c.startswith("qc_")]


def test_qc_block_does_not_displace_user_metadata(tmp_path: Path, known_pairs):
    pairs, bins, contigs, _ = known_pairs
    out = tmp_path / "meta.metacool"
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        metadata={"commands": ["metacooler create pairs ..."], "operator": "mm"},
    )
    metadata = MetaCooler(out).metadata
    assert metadata["operator"] == "mm"
    assert metadata["commands"] == ["metacooler create pairs ..."]
    assert metadata[QC_METADATA_KEY]["source"] == "pairs"


# ----------------------------------------------------------------------
# create cool and the standalone command: the pixel-based tally
# ----------------------------------------------------------------------
def test_from_cool_tallies_from_pixels(tmp_path: Path, tiny_cool, tiny_contigs_table):
    out = tmp_path / "imported.metacool"
    from_cool(cool_uri=tiny_cool, output=out, contig_genome_map=tiny_contigs_table)
    block = MetaCooler(out).metadata[QC_METADATA_KEY]
    assert block["source"] == "pixels"
    assert block["n_contacts"] == MetaCooler(out).info["sum"]
    _assert_contig_columns_add_up(out)


def test_collect_writes_contig_columns(tmp_path: Path):
    from _cools import BSUB, ECOLI, make_cool

    from metacooler.io.collect import collect

    out = tmp_path / "collected.metacool"
    collect(
        [
            make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1),
            make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2),
        ],
        out,
        qc_distance_threshold=250,
    )
    assert MetaCooler(out).metadata[QC_METADATA_KEY]["intra_contig_long"] > 0
    _assert_contig_columns_add_up(out)


def test_pixel_and_pairs_tallies_agree_except_on_short_range(known_metacool):
    """Binning only blurs the intra-contig distance; every other total is exact."""
    path, _ = known_metacool
    pairs_block = MetaCooler(path).metadata[QC_METADATA_KEY]
    pixel_block = qc(path, distance_threshold=1000).summary

    for key in (
        "n_contacts",
        "intra_contig",
        "inter_contig",
        "inter_contig_intra_genome",
        "inter_genome",
        "intra_genome",
    ):
        assert pixel_block[key] == pairs_block[key], key
    assert pixel_block["trans_contig_ratio"] == pairs_block["trans_contig_ratio"]
    # The 400 bp and 500 bp contacts fall inside one 1 kb bin, so their binned
    # distance is 0 — still short. The blur only shows up at the boundary.
    assert (
        pixel_block["intra_contig_short"] + pixel_block["intra_contig_long"]
        == pairs_block["intra_contig_short"] + pairs_block["intra_contig_long"]
    )


def test_a_genome_with_no_contact_gets_nan_not_a_fabricated_ratio(tmp_path: Path):
    """A MAG none of the reads touched has no self-containment score."""
    bins = tmp_path / "bins.bed"
    bins.write_text(
        "".join(f"{c}\t{s}\t{s + 500}\n" for c in ("c1", "c2") for s in (0, 500))
    )
    contigs = tmp_path / "contigs.tsv"
    contigs.write_text("contig_name\tgenome_name\nc1\tMAG_1\nc2\tMAG_SILENT\n")
    pairs = tmp_path / "silent.pairs"
    pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2\nr0\tc1\t10\tc1\t400\n")

    out = tmp_path / "silent.metacool"
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
    )
    genomes = MetaCooler(out).genomes.set_index("name")
    assert genomes.at["MAG_1", "qc_intra_genome_ratio"] == 1.0
    assert np.isnan(genomes.at["MAG_SILENT", "qc_intra_genome_ratio"])
    assert int(genomes.at["MAG_SILENT", "qc_inter_genome"]) == 0


def test_qc_totals_match_the_pixel_table(tiny_metacooler):
    """The tally covers every pixel exactly once."""
    mc = MetaCooler(tiny_metacooler)
    assert qc(tiny_metacooler).summary["n_contacts"] == mc.info["sum"]


def test_qc_respects_chunking(tiny_metacooler):
    """Chunk size is a memory knob, not a result knob."""
    whole = qc(tiny_metacooler).summary
    chunked = qc(tiny_metacooler, chunksize=3).summary
    assert whole == chunked


def test_qc_store_round_trip(tiny_metacooler):
    result = qc(tiny_metacooler, distance_threshold=250, store=True)
    mc = MetaCooler(tiny_metacooler)
    assert mc.metadata[QC_METADATA_KEY] == result.summary
    genomes = mc.genomes
    for column in QC_GENOME_COLUMNS:
        assert column in genomes.columns
    np.testing.assert_allclose(
        genomes["qc_inter_genome"].to_numpy(dtype=float),
        result.per_genome["qc_inter_genome"].to_numpy(dtype=float),
    )


def test_qc_store_writes_contig_columns(tiny_metacooler):
    with h5py.File(tiny_metacooler, "r+") as f:
        for column in QC_CONTIG_COLUMNS:
            del f["contigs"][column]
    qc(tiny_metacooler, distance_threshold=250)
    assert not set(QC_CONTIG_COLUMNS) & set(MetaCooler(tiny_metacooler).contigs)

    qc(tiny_metacooler, distance_threshold=250, store=True)
    _assert_contig_columns_add_up(tiny_metacooler)


def test_store_qc_refuses_a_misaligned_contig_table(tiny_metacooler):
    """Both tables are checked before anything is written."""
    from metacooler.io._utils import store_qc

    before = MetaCooler(tiny_metacooler).metadata
    with pytest.raises(ValueError, match="/contigs has"):
        store_qc(
            tiny_metacooler,
            {"source": "pixels"},
            None,
            pd.DataFrame({column: [1] for column in QC_CONTIG_COLUMNS}),
        )
    assert MetaCooler(tiny_metacooler).metadata == before


def test_qc_store_twice_does_not_duplicate_columns(tiny_metacooler):
    qc(tiny_metacooler, store=True)
    qc(tiny_metacooler, distance_threshold=5000, store=True)
    columns = list(MetaCooler(tiny_metacooler).genomes.columns)
    assert len(columns) == len(set(columns))
    assert (
        MetaCooler(tiny_metacooler).metadata[QC_METADATA_KEY]["distance_threshold"]
        == 5000
    )


def test_qc_store_keeps_the_file_valid(tiny_metacooler):
    from metacooler.schema import validate

    qc(tiny_metacooler, store=True)
    validate(tiny_metacooler)  # raises on any structural problem
    import cooler

    cooler.Cooler(str(tiny_metacooler))  # the dual identity must survive


def test_qc_on_mmetacool_entry(tiny_mmetacool):
    uri = f"{tiny_mmetacool}::entries/entry1"
    result = qc(uri, store=True)
    assert result.summary["n_contacts"] > 0
    assert MetaCooler(uri).metadata[QC_METADATA_KEY] == result.summary


def test_qc_warns_when_bins_are_coarser_than_the_threshold(tiny_metacooler, caplog):
    qc(tiny_metacooler, distance_threshold=50)
    assert "at least the short/long distance threshold" in caplog.text


def test_qc_warns_before_replacing_a_pairs_block(known_metacool, caplog):
    path, _ = known_metacool
    qc(path, store=True)
    assert "Replacing a block measured on read positions" in caplog.text


# ----------------------------------------------------------------------
# Staleness in derived files
# ----------------------------------------------------------------------
def test_strip_helpers():
    frame = pd.DataFrame({"name": ["a"], "length": [1], "qc_inter_genome": [2]})
    assert list(strip_qc_columns(frame).columns) == ["name", "length"]
    assert strip_qc_metadata({"qc": {}, "operator": "mm"}) == {"operator": "mm"}


def test_which_columns_are_qc():
    """In /contigs only the two QC names are; in /genomes every qc_* is."""
    from metacooler.qc import is_qc_column

    assert is_qc_column("contigs", "qc_intra_contig_short")
    assert not is_qc_column("contigs", "qc_pass")
    assert is_qc_column("genomes", "qc_pass")
    assert not is_qc_column("genomes", "completeness")
    with pytest.raises(ValueError, match="'genomes' or 'contigs'"):
        is_qc_column("bins", "qc_pass")


@pytest.mark.parametrize("level", ["contig", "genome"])
def test_binnify_drops_stale_qc(
    tmp_path: Path, tiny_metacooler, level: Literal["contig", "genome"]
):
    """Genome level, /contigs holds the ligated genomes: no QC column of theirs."""
    from metacooler.ops.binnify import binnify

    qc(tiny_metacooler, store=True)
    out = tmp_path / "rebinned.metacool"
    binnify(uri=tiny_metacooler, output=out, binsize=200, level=level)
    mc = MetaCooler(out)
    assert QC_METADATA_KEY not in mc.metadata
    assert not [c for c in mc.genomes.columns if c.startswith("qc_")]
    assert not [c for c in mc.contigs.columns if c.startswith("qc_")]


def test_merge_ignores_qc_columns_when_comparing_genomes(tmp_path: Path, known_pairs):
    """Two libraries over one assembly differ in QC but are still mergeable."""
    from metacooler.io.merge import merge

    pairs, bins, contigs, _ = known_pairs
    first = tmp_path / "a.metacool"
    second = tmp_path / "b.metacool"
    from_pairs(
        uri=first,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        qc_distance_threshold=1000,
    )
    from_pairs(
        uri=second,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_ordering="genome-then-alpha",
        # A different threshold gives the same file a different qc_* split,
        # which is exactly the divergence the merge must tolerate.
        qc_distance_threshold=5000,
    )
    out = tmp_path / "merged.metacool"
    merge([first, second], out)

    mc = MetaCooler(out)
    assert QC_METADATA_KEY not in mc.metadata
    assert not [c for c in mc.genomes.columns if c.startswith("qc_")]
    assert not [c for c in mc.contigs.columns if c.startswith("qc_")]
    assert mc.info["sum"] == 2 * MetaCooler(first).info["sum"]


class TestAUserQcColumn:
    """A user /contigs column that merely starts with ``qc_`` is metadata.

    Only the two per-contig QC columns are metacooler's; ``qc_pass`` from a
    ``--contigs`` table still describes its contig in every derived file.
    """

    @pytest.fixture
    def source(self, tmp_path: Path, known_pairs) -> Path:
        pairs, bins, _, _ = known_pairs
        contigs = tmp_path / "contigs_with_qc_pass.tsv"
        contigs.write_text(
            "contig_name\tgenome_name\tqc_pass\n"
            "c1\tMAG_1\t1\nc2\tMAG_1\t0\nc3\tMAG_2\t1\n"
        )
        out = tmp_path / "qc_pass.metacool"
        from_pairs(
            uri=out,
            bins=bins,
            pairs=pairs,
            contig_genome_map=contigs,
            contig_ordering="genome-then-alpha",
            qc_distance_threshold=1000,
        )
        return out

    @staticmethod
    def _assert_kept(path: Path) -> None:
        contigs = MetaCooler(path).contigs.set_index("name")
        assert contigs["qc_pass"].dropna().to_dict() == {"c1": 1, "c2": 0, "c3": 1}
        assert not set(QC_CONTIG_COLUMNS) & set(contigs.columns)

    def test_merge_keeps_it(self, tmp_path: Path, source: Path):
        from metacooler.io.merge import merge

        copy = tmp_path / "copy.metacool"
        copy.write_bytes(source.read_bytes())
        out = tmp_path / "merged.metacool"
        merge([source, copy], out)
        self._assert_kept(out)

    def test_binnify_keeps_it(self, tmp_path: Path, source: Path):
        from metacooler.ops.binnify import binnify

        out = tmp_path / "rebinned.metacool"
        binnify(uri=source, output=out, binsize=2000)
        self._assert_kept(out)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def test_cli_qc_json(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(cli, ["qc", str(tiny_metacooler), "-q"])
    assert result.exit_code == 0, result.output
    block = json.loads(result.output)
    assert block["source"] == "pixels"
    assert block["n_contacts"] > 0


def test_cli_qc_threshold_accepts_size_suffix(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(
        cli, ["qc", str(tiny_metacooler), "--threshold", "2kb", "-q"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["distance_threshold"] == 2000


def test_cli_qc_rejects_a_bad_threshold(tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    result = CliRunner().invoke(
        cli, ["qc", str(tiny_metacooler), "--threshold", "ten-kb", "-q"]
    )
    assert result.exit_code != 0
    assert "not a number" in _flat(result.output)


def test_cli_qc_per_genome_tsv(tmp_path: Path, tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out = tmp_path / "per_genome.tsv"
    result = CliRunner().invoke(
        cli, ["qc", str(tiny_metacooler), "--per-genome", str(out), "-q"]
    )
    assert result.exit_code == 0, result.output
    table = pd.read_csv(out, sep="\t")
    assert list(table.columns) == ["name", *QC_GENOME_COLUMNS]
    assert len(table) == len(MetaCooler(tiny_metacooler).genomes)


def test_cli_qc_refuses_to_clobber(tmp_path: Path, tiny_metacooler):
    from click.testing import CliRunner

    from metacooler.cli.main import cli

    out = tmp_path / "taken.json"
    out.write_text("keep me")
    result = CliRunner().invoke(cli, ["qc", str(tiny_metacooler), "-o", str(out), "-q"])
    assert result.exit_code != 0
    assert out.read_text() == "keep me"
