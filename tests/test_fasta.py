"""``metacooler fasta``: the contigs of each genome or cluster, as FASTA.

The records are copied from the assembly byte for byte (as ``pyfastx``'s
``raw`` gives them), so a tool that wrote its FASTA files record by record --
metaTOR's ``final_clusters/*.fa`` -- gets the same bytes from here.
"""

from __future__ import annotations

import errno
import functools
import gzip
import logging
import os
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

import metacooler
from metacooler.cli.main import cli
from metacooler.filters import FilterError
from metacooler.io import fasta as fasta_module
from metacooler.io.clusters import write_clusters
from metacooler.io.create import from_pairs
from metacooler.io.fasta import LINE_WIDTH, write_fasta
from metacooler.io.mmetacool import create_mmetacool

LENGTHS = {"c1": 4000, "c2": 3000, "c3": 2000, "c4": 1000, "c5": 500}
ASSIGNMENT = {"c1": "MAG_A", "c2": "MAG_A", "c3": "MAG_B"}  # c4, c5: unbinned

# Each record of the assembly: its header line and the width of its sequence
# lines, all different, in another order than the map's (c1 ... c5).
LAYOUT = {
    "c3": (">c3 third contig, len=2000", 13),
    "c5": (">c5", 7),
    "c1": (">c1 desc one", 60),
    "xx": (">xx not in the map", 50),
    "c4": (">c4", 1000),
    "c2": (">c2", 70),
}
EXTRA_LENGTH = 300


def _sequences() -> dict[str, str]:
    rng = np.random.default_rng(7)
    lengths = {**LENGTHS, "xx": EXTRA_LENGTH}
    return {
        name: "".join(rng.choice(list("ACGTacgt"), size=n))
        for name, n in lengths.items()
    }


SEQUENCES = _sequences()


def _record(name: str, newline: str = "\n") -> str:
    header, width = LAYOUT[name]
    seq = SEQUENCES[name]
    lines = [seq[i : i + width] for i in range(0, len(seq), width)]
    return header + newline + "".join(line + newline for line in lines)


RAW = {name: _record(name) for name in LAYOUT}


def _write_assembly(path: Path, text: str | None = None) -> Path:
    path.write_text("".join(RAW.values()) if text is None else text)
    return path


def _map(path: Path, assignment: dict[str, str] | None = ASSIGNMENT) -> Path:
    pairs = path.with_suffix(".pairs")
    rng = np.random.default_rng(1)
    names = list(LENGTHS)
    rows = []
    for i in range(300):
        a, b = names[rng.integers(5)], names[rng.integers(5)]
        rows.append(
            f"r{i}\t{a}\t{rng.integers(1, LENGTHS[a])}\t{b}\t{rng.integers(1, LENGTHS[b])}\n"
        )
    pairs.write_text(
        "## pairs format v1.0\n"
        + "".join(f"#chromsize: {c} {n}\n" for c, n in LENGTHS.items())
        + "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(rows)
    )
    bins = pd.DataFrame(
        {"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())}
    )
    from_pairs(path, bins=bins, pairs=pairs, contig_genome_map=assignment, ps=False)
    return path


@pytest.fixture
def mapped(tmp_path: Path) -> tuple[Path, Path]:
    """A map of five contigs (MAG_A: c1, c2; MAG_B: c3; c4 and c5 in no
    genome), and its assembly, in another order, at other line widths."""
    work = tmp_path / "in"
    work.mkdir()
    return _map(work / "map.metacool"), _write_assembly(work / "assembly.fa")


@pytest.fixture
def clustered(mapped: tuple[Path, Path]) -> tuple[Path, Path]:
    """The map, with two rounds: overlapping (O1: c1 c2 c3, O2: c4) and final
    (F1: c2 c1, given in that order; F2: c3)."""
    path, assembly = mapped
    clusters = pd.DataFrame(
        {
            "round": ["overlapping", "overlapping", "final", "final"],
            "name": ["O1", "O2", "F1", "F2"],
            "completeness": [80.0, 5.0, 95.0, 50.0],
        }
    )
    members = pd.DataFrame(
        {
            "round": ["overlapping"] * 4 + ["final"] * 3,
            "cluster": ["O1", "O1", "O1", "O2", "F1", "F1", "F2"],
            "contig": ["c3", "c1", "c2", "c4", "c2", "c1", "c3"],
        }
    )
    write_clusters(path, ["overlapping", "final"], clusters, members)
    return path, assembly


def _bases(*names: str) -> list[str]:
    return [SEQUENCES[name] for name in names]


def _wrapped(header: str, sequence: str) -> str:
    lines = [sequence[i : i + LINE_WIDTH] for i in range(0, len(sequence), LINE_WIDTH)]
    return header + "\n" + "".join(line + "\n" for line in lines)


# ---------------------------------------------------------------------------
# Records as they are
# ---------------------------------------------------------------------------


def test_each_genome_gets_its_records_as_they_are(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    written = write_fasta(path, assembly, out)

    assert written == {"MAG_A": out / "MAG_A.fa", "MAG_B": out / "MAG_B.fa"}
    assert sorted(p.name for p in out.iterdir()) == ["MAG_A.fa", "MAG_B.fa"]
    # Byte for byte, in map order: header lines, line widths, descriptions.
    assert (out / "MAG_A.fa").read_text() == RAW["c1"] + RAW["c2"]
    assert (out / "MAG_B.fa").read_text() == RAW["c3"]


def test_nothing_is_written_next_to_the_assembly(mapped, tmp_path):
    path, assembly = mapped
    before = sorted(p.name for p in assembly.parent.iterdir())
    write_fasta(path, assembly, tmp_path / "out", stitch=True)
    assert sorted(p.name for p in assembly.parent.iterdir()) == before


def test_a_gzip_assembly_gives_the_same_files(mapped, tmp_path):
    path, assembly = mapped
    plain = write_fasta(path, assembly, tmp_path / "plain", unbinned=True)
    packed = assembly.with_name("assembly.fa.gz")
    packed.write_bytes(gzip.compress(assembly.read_bytes()))
    before = sorted(p.name for p in packed.parent.iterdir())
    from_gzip = write_fasta(path, packed, tmp_path / "gz", unbinned=True)

    assert list(plain) == list(from_gzip) == ["MAG_A", "MAG_B", "__unbinned__"]
    for name in plain:
        assert from_gzip[name].read_bytes() == plain[name].read_bytes()
    assert sorted(p.name for p in packed.parent.iterdir()) == before


def test_records_keep_their_crlf_and_blank_lines(mapped, tmp_path):
    """A record runs from its '>' to the next one: blank lines and carriage
    returns included, as pyfastx's raw gives it."""
    path, _ = mapped
    odd = {name: _record(name) for name in LAYOUT}
    odd["c1"] = _record("c1", newline="\r\n")
    odd["c2"] = _record("c2") + "\n\n"
    assembly = _write_assembly(tmp_path / "odd.fa", "".join(odd.values()))
    write_fasta(path, assembly, tmp_path / "out")
    assert (tmp_path / "out" / "MAG_A.fa").read_bytes() == (
        odd["c1"] + odd["c2"]
    ).encode()

    # Stitched, the bases alone: no carriage return, no blank line.
    write_fasta(path, assembly, tmp_path / "stitched", stitch=True, spacer=0)
    expected = _wrapped(">MAG_A", "".join(_bases("c1", "c2")))
    assert (tmp_path / "stitched" / "MAG_A.fa").read_text() == expected


def test_the_last_record_gets_a_line_break(mapped, tmp_path):
    """An assembly that ends without a line break: its last record gets one,
    so the record written after it starts on its own line."""
    path, _ = mapped
    records = [RAW[name] for name in ("c3", "c5", "xx", "c4", "c2")]
    text = "".join(records) + RAW["c1"].rstrip("\n")
    assembly = _write_assembly(tmp_path / "noeol.fa", text)
    write_fasta(path, assembly, tmp_path / "out")
    assert (tmp_path / "out" / "MAG_A.fa").read_text() == RAW["c1"] + RAW["c2"]


def test_groups_are_written_in_the_order_given(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    written = write_fasta(
        path, assembly, out, groups={"g2": ["c2", "c3", "c1"], "g1": ["c5"]}
    )
    assert list(written) == ["g2", "g1"]
    assert (out / "g2.fa").read_text() == RAW["c2"] + RAW["c3"] + RAW["c1"]
    assert (out / "g1.fa").read_text() == RAW["c5"]


def test_groups_with_the_contigs_in_none(mapped, tmp_path):
    path, assembly = mapped
    written = write_fasta(
        path, assembly, tmp_path / "out", groups={"g": ["c4"]}, unbinned=True
    )
    assert list(written) == ["g", "__unbinned__"]
    # In map order: the map holds __unbinned__ (c4, c5) first, then MAG_A, MAG_B.
    assert written["__unbinned__"].read_text() == "".join(
        RAW[name] for name in ("c5", "c1", "c2", "c3")
    )


# ---------------------------------------------------------------------------
# Stitched, single file, unbinned
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spacer", [0, 7, 100])
def test_stitched_records(mapped, tmp_path, spacer):
    path, assembly = mapped
    out = tmp_path / "out"
    write_fasta(path, assembly, out, stitch=True, spacer=spacer)
    joined = ("N" * spacer).join(_bases("c1", "c2"))
    assert (out / "MAG_A.fa").read_text() == _wrapped(">MAG_A", joined)
    assert (out / "MAG_B.fa").read_text() == _wrapped(">MAG_B", SEQUENCES["c3"])
    lines = (out / "MAG_A.fa").read_text().splitlines()[1:]
    assert {len(line) for line in lines[:-1]} == {LINE_WIDTH}


def test_unbinned_contigs_stay_as_they_are_when_stitched(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    written = write_fasta(path, assembly, out, stitch=True, unbinned=True)
    assert list(written) == ["MAG_A", "MAG_B", "__unbinned__"]
    assert (out / "__unbinned__.fa").read_text() == RAW["c4"] + RAW["c5"]


def test_one_file_names_each_record_after_its_genome(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "all.fa"
    written = write_fasta(path, assembly, out, single_file=True, unbinned=True)
    assert set(written.values()) == {out}
    expected = "".join(
        f">{genome}|{RAW[contig][1:]}"
        for genome, contig in [
            ("MAG_A", "c1"),
            ("MAG_A", "c2"),
            ("MAG_B", "c3"),
            ("__unbinned__", "c4"),
            ("__unbinned__", "c5"),
        ]
    )
    assert out.read_text() == expected
    assert out.read_text().startswith(">MAG_A|c1 desc one\n")


def test_one_file_stitched(mapped, tmp_path):
    """metaTOR's mges_binned.fa: one stitched record per cluster, 200 N."""
    path, assembly = mapped
    out = tmp_path / "all.fa"
    write_fasta(
        path, assembly, out, single_file=True, stitch=True, spacer=200, unbinned=True
    )
    expected = (
        _wrapped(">MAG_A", ("N" * 200).join(_bases("c1", "c2")))
        + _wrapped(">MAG_B", SEQUENCES["c3"])
        + f">__unbinned__|{RAW['c4'][1:]}>__unbinned__|{RAW['c5'][1:]}"
    )
    assert out.read_text() == expected


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def test_genomes_by_name(mapped, tmp_path):
    path, assembly = mapped
    written = write_fasta(path, assembly, tmp_path / "out", genomes=["MAG_B"])
    assert list(written) == ["MAG_B"]


def test_genome_and_contig_filters(mapped, tmp_path, caplog):
    path, assembly = mapped
    written = write_fasta(
        path, assembly, tmp_path / "a", filters=["genomes.length > 6kb"]
    )
    assert list(written) == ["MAG_A"]

    with caplog.at_level(logging.INFO, logger="metacooler"):
        written = write_fasta(
            path, assembly, tmp_path / "b", filters=["contigs.length > 3500"]
        )
    assert list(written) == ["MAG_A"]
    assert written["MAG_A"].read_text() == RAW["c1"]
    assert "1 genome(s) have no contig left after --filter" in caplog.text


def test_unbinned_contigs_follow_the_contig_filters_only(mapped, tmp_path):
    path, assembly = mapped
    written = write_fasta(
        path,
        assembly,
        tmp_path / "out",
        filters=["genomes.genome == MAG_B", "contigs.length >= 1000"],
        unbinned=True,
    )
    assert list(written) == ["MAG_B", "__unbinned__"]
    assert written["__unbinned__"].read_text() == RAW["c4"]


def test_a_filter_that_keeps_nothing_says_which(mapped, tmp_path):
    path, assembly = mapped
    with pytest.raises(FilterError, match="keeps nothing"):
        write_fasta(path, assembly, tmp_path / "out", filters=["genomes.length > 1Gb"])
    with pytest.raises(FilterError, match="named by --genomes"):
        write_fasta(
            path,
            assembly,
            tmp_path / "out",
            genomes=["MAG_B"],
            filters=["contigs.length > 3kb"],
        )
    assert not (tmp_path / "out").exists()


def test_the_unbinned_contigs_are_asked_for_by_flag(mapped, tmp_path):
    path, assembly = mapped
    with pytest.raises(ValueError, match="__unbinned__ is no genome"):
        write_fasta(path, assembly, tmp_path / "out", genomes=["__unbinned__"])
    with pytest.raises(KeyError, match="Not genomes of"):
        write_fasta(path, assembly, tmp_path / "out", genomes=["MAG_Z"])


def test_a_map_with_no_genome_yet(tmp_path):
    path = _map(tmp_path / "unbinned.metacool", assignment=None)
    assembly = _write_assembly(tmp_path / "assembly.fa")
    with pytest.raises(ValueError, match="has no genome to write"):
        write_fasta(path, assembly, tmp_path / "out")
    written = write_fasta(path, assembly, tmp_path / "out", unbinned=True)
    assert list(written) == ["__unbinned__"]
    assert written["__unbinned__"].read_text() == "".join(
        RAW[name] for name in ("c1", "c2", "c3", "c4", "c5")
    )


def test_a_map_with_every_contig_in_a_genome(tmp_path, caplog):
    every = {**ASSIGNMENT, "c4": "MAG_B", "c5": "MAG_B"}
    path = _map(tmp_path / "binned.metacool", assignment=every)
    assembly = _write_assembly(tmp_path / "assembly.fa")
    with caplog.at_level(logging.INFO, logger="metacooler"):
        written = write_fasta(path, assembly, tmp_path / "out", unbinned=True)
    assert list(written) == ["MAG_A", "MAG_B"]
    assert "No contig is in no genome or cluster" in caplog.text


def test_the_clusters_of_a_round(clustered, tmp_path):
    path, assembly = clustered
    out = tmp_path / "out"
    written = write_fasta(path, assembly, out, round="final", unbinned=True)
    assert list(written) == ["F1", "F2", "__unbinned__"]
    # In map order, not the order the members were given in.
    assert (out / "F1.fa").read_text() == RAW["c1"] + RAW["c2"]
    assert (out / "F2.fa").read_text() == RAW["c3"]
    assert (out / "__unbinned__.fa").read_text() == RAW["c4"] + RAW["c5"]


def test_clusters_by_name_and_by_filter(clustered, tmp_path, caplog):
    path, assembly = clustered
    by_name = write_fasta(path, assembly, tmp_path / "a", round="final", genomes=["F2"])
    assert list(by_name) == ["F2"]
    for expr in ["clusters.completeness >= 90", "clusters.n_contigs >= 2"]:
        kept = write_fasta(
            path, assembly, tmp_path / expr[9:12], round="final", filters=[expr]
        )
        assert list(kept) == ["F1"]
    with caplog.at_level(logging.INFO, logger="metacooler"):
        kept = write_fasta(
            path,
            assembly,
            tmp_path / "c",
            round="final",
            filters=["contigs.length >= 2500"],
        )
    assert list(kept) == ["F1"]
    assert "1 cluster(s) of round 'final' have no contig left" in caplog.text
    with pytest.raises(FilterError, match="keeps nothing"):
        write_fasta(
            path,
            assembly,
            tmp_path / "d",
            round="final",
            filters=["clusters.completeness > 99"],
        )


def test_a_round_refuses_what_it_does_not_read(clustered, tmp_path):
    path, assembly = clustered
    with pytest.raises(FilterError, match="filter them with 'clusters.<column>'"):
        write_fasta(
            path,
            assembly,
            tmp_path / "o",
            round="final",
            filters=["genomes.length > 1"],
        )
    with pytest.raises(FilterError, match="give the round"):
        write_fasta(path, assembly, tmp_path / "o", filters=["clusters.length > 1"])
    with pytest.raises(KeyError, match="is not a clustering round"):
        write_fasta(path, assembly, tmp_path / "o", round="recursive")
    with pytest.raises(KeyError, match="Not clusters of round 'final'"):
        write_fasta(path, assembly, tmp_path / "o", round="final", genomes=["O1"])


def test_a_cluster_without_a_contig(mapped, tmp_path):
    path, assembly = mapped
    clusters = pd.DataFrame({"round": ["r", "r"], "name": ["empty", "full"]})
    members = pd.DataFrame({"round": ["r"], "cluster": ["full"], "contig": ["c1"]})
    write_clusters(path, ["r"], clusters, members)
    assert list(write_fasta(path, assembly, tmp_path / "a", round="r")) == ["full"]
    with pytest.raises(ValueError, match="has no contig to write"):
        write_fasta(path, assembly, tmp_path / "b", round="r", genomes=["empty"])


def test_groups_that_hold_every_contig(mapped, tmp_path, caplog):
    path, assembly = mapped
    with caplog.at_level(logging.INFO, logger="metacooler"):
        written = write_fasta(
            path,
            assembly,
            tmp_path / "out",
            groups={"all": list(LENGTHS)},
            unbinned=True,
        )
    assert list(written) == ["all"]
    assert "No contig is in no genome or cluster" in caplog.text


def test_a_cluster_named_like_the_unbinned_contigs(mapped, tmp_path):
    path, assembly = mapped
    clusters = pd.DataFrame({"round": ["r"], "name": ["__unbinned__"]})
    members = pd.DataFrame(
        {"round": ["r"], "cluster": ["__unbinned__"], "contig": ["c1"]}
    )
    write_clusters(path, ["r"], clusters, members)
    assert list(write_fasta(path, assembly, tmp_path / "a", round="r")) == [
        "__unbinned__"
    ]
    with pytest.raises(ValueError, match="rename it"):
        write_fasta(path, assembly, tmp_path / "b", round="r", unbinned=True)


def test_an_entry_of_a_collection(mapped, tmp_path):
    path, assembly = mapped
    collection = tmp_path / "maps.mmetacool"
    create_mmetacool(collection, {"rep1": path})
    written = write_fasta(f"{collection}::entries/rep1", assembly, tmp_path / "out")
    assert written["MAG_A"].read_text() == RAW["c1"] + RAW["c2"]


# ---------------------------------------------------------------------------
# What --unbinned does not stand in for (review finding 1)
# ---------------------------------------------------------------------------


def test_unbinned_alone_cannot_hide_an_empty_selection(mapped, tmp_path):
    """Names or genome filters that keep no genome are refused, --unbinned or
    not: writing __unbinned__.fa alone would pass for a selection."""
    path, assembly = mapped
    with pytest.raises(FilterError, match="genomes.length > 1Gb  ->  0 of"):
        write_fasta(
            path,
            assembly,
            tmp_path / "a",
            filters=["genomes.length > 1Gb"],
            unbinned=True,
        )
    with pytest.raises(FilterError, match="1 genome\\(s\\) named by --genomes"):
        write_fasta(
            path,
            assembly,
            tmp_path / "b",
            genomes=["MAG_B"],
            filters=["contigs.length > 3kb"],
            unbinned=True,
        )
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()


def test_unbinned_alone_cannot_hide_an_empty_round(clustered, tmp_path):
    path, assembly = clustered
    with pytest.raises(FilterError, match="1 cluster\\(s\\) named by --genomes"):
        write_fasta(
            path,
            assembly,
            tmp_path / "a",
            round="final",
            genomes=["F2"],
            filters=["contigs.length >= 2500"],
            unbinned=True,
        )
    with pytest.raises(FilterError, match="clusters.completeness > 99"):
        write_fasta(
            path,
            assembly,
            tmp_path / "b",
            round="final",
            filters=["clusters.completeness > 99"],
            unbinned=True,
        )


def test_the_unbinned_genome_is_asked_for_with_unbinned(mapped, tmp_path, caplog):
    path, assembly = mapped
    alone = ["genomes.genome == __unbinned__"]
    with pytest.raises(FilterError, match="write its contigs with --unbinned"):
        write_fasta(path, assembly, tmp_path / "a", filters=alone)
    with caplog.at_level(logging.INFO, logger="metacooler"):
        written = write_fasta(
            path, assembly, tmp_path / "b", filters=alone, unbinned=True
        )
    assert list(written) == ["__unbinned__"]
    assert "0 genome(s) and __unbinned__" in caplog.text
    with pytest.raises(FilterError, match="keeps nothing"):
        write_fasta(
            path,
            assembly,
            tmp_path / "c",
            filters=[*alone, "contigs.length > 1Gb"],
            unbinned=True,
        )


def test_only_contig_filters_may_leave_the_unbinned_contigs_alone(mapped, tmp_path):
    path, assembly = mapped
    written = write_fasta(
        path, assembly, tmp_path / "a", filters=["contigs.length < 600"], unbinned=True
    )
    assert list(written) == ["__unbinned__"]
    assert written["__unbinned__"].read_text() == RAW["c5"]


def test_contig_filters_that_keep_nothing_of_a_round(clustered, tmp_path):
    path, assembly = clustered
    written = write_fasta(
        path,
        assembly,
        tmp_path / "a",
        round="final",
        filters=["contigs.length < 1500"],
        unbinned=True,
    )
    assert list(written) == ["__unbinned__"]
    with pytest.raises(FilterError, match="keeps nothing"):
        write_fasta(
            path,
            assembly,
            tmp_path / "b",
            round="final",
            filters=["contigs.length > 1Gb"],
            unbinned=True,
        )


def test_the_summary_counts_the_unbinned_contigs_apart(mapped, tmp_path, caplog):
    path, assembly = mapped
    with caplog.at_level(logging.INFO, logger="metacooler"):
        write_fasta(path, assembly, tmp_path / "out", unbinned=True)
    assert "3 FASTA file(s)" in caplog.text
    assert "2 genome(s) and __unbinned__, 5 record(s)" in caplog.text


def test_a_cluster_that_holds_no_contig_is_named(mapped, tmp_path, caplog):
    path, assembly = mapped
    clusters = pd.DataFrame({"round": ["r", "r"], "name": ["empty", "full"]})
    members = pd.DataFrame({"round": ["r"], "cluster": ["full"], "contig": ["c1"]})
    write_clusters(path, ["r"], clusters, members)
    with caplog.at_level(logging.INFO, logger="metacooler"):
        write_fasta(path, assembly, tmp_path / "a", round="r")
    assert "1 cluster(s) of round 'r' hold no contig: no file for them (empty)" in (
        caplog.text
    )
    assert "after --filter" not in caplog.text


# ---------------------------------------------------------------------------
# The order of a stitched record (review finding 2)
# ---------------------------------------------------------------------------


def test_a_stitched_cluster_comes_longest_first(mapped, tmp_path):
    """A cluster across genomes (metaTOR's mge): the map holds c4 (in no
    genome) first, then c1 and c3; stitched, c1, c3, c4; as they are, in map
    order."""
    path, assembly = mapped
    clusters = pd.DataFrame({"round": ["mge"], "name": ["MGE_1"]})
    members = pd.DataFrame(
        {"round": "mge", "cluster": "MGE_1", "contig": ["c3", "c4", "c1"]}
    )
    write_clusters(path, ["mge"], clusters, members)
    stitched = write_fasta(
        path, assembly, tmp_path / "s", round="mge", stitch=True, spacer=3
    )
    assert stitched["MGE_1"].read_text() == _wrapped(
        ">MGE_1", "NNN".join(_bases("c1", "c3", "c4"))
    )
    raw = write_fasta(path, assembly, tmp_path / "r", round="mge")
    assert raw["MGE_1"].read_text() == RAW["c4"] + RAW["c1"] + RAW["c3"]
    # A group keeps the order given, stitched too.
    given = write_fasta(
        path, assembly, tmp_path / "g", groups={"g": ["c4", "c3"]}, stitch=True
    )
    assert given["g"].read_text() == _wrapped(
        ">g", ("N" * 100).join(_bases("c4", "c3"))
    )


def test_stitched_contigs_of_one_length_keep_map_order(tmp_path):
    """In a map laid out genome-then-alpha (a, b, c, d), stitched: b and c
    (3 kb each, in map order), d, a; as they are: a, b, c, d."""
    lengths = {"a": 1000, "b": 3000, "c": 3000, "d": 2000}
    rng = np.random.default_rng(3)
    sequences = {
        n: "".join(rng.choice(list("ACGT"), size=k)) for n, k in lengths.items()
    }
    assembly = tmp_path / "abcd.fa"
    assembly.write_text("".join(f">{n}\n{q}\n" for n, q in sequences.items()))
    pairs = tmp_path / "abcd.pairs"
    pairs.write_text(
        "## pairs format v1.0\n"
        + "".join(f"#chromsize: {c} {n}\n" for c, n in lengths.items())
        + "#columns: readID chrom1 pos1 chrom2 pos2\nr0\ta\t1\tb\t1\n"
    )
    bins = pd.DataFrame(
        {"chrom": list(lengths), "start": 0, "end": list(lengths.values())}
    )
    path = tmp_path / "abcd.metacool"
    from_pairs(
        path,
        bins=bins,
        pairs=pairs,
        contig_genome_map=dict.fromkeys(lengths, "G"),
        contig_ordering="genome-then-alpha",
        ps=False,
    )
    stitched = write_fasta(path, assembly, tmp_path / "s", stitch=True, spacer=0)
    order = "".join(sequences[n] for n in "bcda")
    assert stitched["G"].read_text() == _wrapped(">G", order)
    raw = write_fasta(path, assembly, tmp_path / "r")
    assert raw["G"].read_text() == "".join(f">{n}\n{sequences[n]}\n" for n in "abcd")


# ---------------------------------------------------------------------------
# Names and lengths as the map's creator reads them (review finding 3)
# ---------------------------------------------------------------------------


def test_a_space_after_the_header_mark(mapped, tmp_path):
    path, _ = mapped
    text = "".join(raw for key, raw in RAW.items() if key != "c3")
    text += "> c3 spaced\n" + RAW["c3"].split("\n", 1)[1]
    assembly = _write_assembly(tmp_path / "spaced.fa", text)
    written = write_fasta(path, assembly, tmp_path / "out", genomes=["MAG_B"])
    assert written["MAG_B"].read_text().startswith("> c3 spaced\n")


@pytest.mark.parametrize("counted", ["with", "without"])
def test_trailing_spaces_count_as_create_counted_them(tmp_path, counted):
    """pyfaidx (create --fasta) counts a trailing space as a base: a map made
    so holds the longer length; one made from the bases alone, the shorter.
    Either is the map's contig; stitched, the bases alone."""
    bases = "ACGT" * 25  # 100 bp, in 4 lines of 25 and a space each
    text = ">s1\n" + "".join(bases[i : i + 25] + " \n" for i in range(0, 100, 25))
    assembly = tmp_path / "spaces.fa"
    assembly.write_text(text)
    length = 104 if counted == "with" else 100
    pairs = tmp_path / "s.pairs"
    pairs.write_text(
        f"## pairs format v1.0\n#chromsize: s1 {length}\n"
        "#columns: readID chrom1 pos1 chrom2 pos2\nr0\ts1\t1\ts1\t50\n"
    )
    path = tmp_path / "s.metacool"
    bins = pd.DataFrame({"chrom": ["s1"], "start": [0], "end": [length]})
    from_pairs(path, bins=bins, pairs=pairs, contig_genome_map={"s1": "G"}, ps=False)
    assert write_fasta(path, assembly, tmp_path / "raw")["G"].read_text() == text
    stitched = write_fasta(path, assembly, tmp_path / "st", stitch=True)
    assert stitched["G"].read_text() == _wrapped(">G", bases)


# ---------------------------------------------------------------------------
# Compressed and piped assemblies (review finding 6)
# ---------------------------------------------------------------------------


def _fifo(path: Path, data: bytes) -> threading.Thread:
    os.mkfifo(path)

    def feed():
        with open(path, "wb") as handle:
            handle.write(data)

    thread = threading.Thread(target=feed, daemon=True)
    thread.start()
    return thread


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_a_piped_assembly(mapped, tmp_path, compressed):
    path, assembly = mapped
    data = assembly.read_bytes()
    feeder = _fifo(tmp_path / "pipe", gzip.compress(data) if compressed else data)
    scratch = tmp_path / "scratch"
    written = write_fasta(path, tmp_path / "pipe", tmp_path / "out", temp_dir=scratch)
    feeder.join(timeout=10)
    assert written["MAG_A"].read_text() == RAW["c1"] + RAW["c2"]
    assert scratch.is_dir() and not list(scratch.iterdir())  # created, then emptied


def test_a_copy_that_cannot_be_written_names_where(mapped, tmp_path, monkeypatch):
    path, assembly = mapped
    packed = tmp_path / "assembly.fa.gz"
    packed.write_bytes(gzip.compress(assembly.read_bytes()))

    def full(source, target, length=0):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(fasta_module.shutil, "copyfileobj", full)
    with pytest.raises(OSError, match="its plain copy .*No space left") as error:
        write_fasta(path, packed, tmp_path / "out", temp_dir=tmp_path / "scratch")
    assert error.value.errno == errno.ENOSPC
    assert str(tmp_path / "scratch") in str(error.value)


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def _text_without(name: str) -> str:
    return "".join(raw for key, raw in RAW.items() if key != name)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (_text_without("c2"), "1 contig\\(s\\) to write are not in"),
        (
            _text_without("c1") + ">c1\nACGT\n",
            "another length .* c1 \\(4 bp in the FASTA, 4000 in the map\\)",
        ),
        ("".join(RAW.values()) + RAW["c1"], "named by several records .*: c1"),
        ("not a fasta\n" + "".join(RAW.values()), "is not a FASTA file"),
        ("  " + "".join(RAW.values()), "is not a FASTA file"),
        ("", "is not a FASTA file"),
    ],
    ids=["missing", "length", "twice", "text-first", "indented", "empty"],
)
def test_an_assembly_that_is_not_the_maps(mapped, tmp_path, text, message):
    path, _ = mapped
    assembly = _write_assembly(tmp_path / "other.fa", text)
    with pytest.raises(ValueError, match=message):
        write_fasta(path, assembly, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_blank_lines_before_the_first_record(mapped, tmp_path):
    path, _ = mapped
    assembly = _write_assembly(tmp_path / "blank.fa", "\n\n" + "".join(RAW.values()))
    written = write_fasta(path, assembly, tmp_path / "out", genomes=["MAG_B"])
    assert written["MAG_B"].read_text() == RAW["c3"]


def test_a_broken_gzip_assembly(mapped, tmp_path):
    path, _ = mapped
    broken = tmp_path / "broken.fa.gz"
    broken.write_bytes(gzip.compress(b">c1\nACGT\n" * 50)[:-20])
    with pytest.raises(ValueError, match="Cannot decompress"):
        write_fasta(path, broken, tmp_path / "out")


def test_every_missing_contig_is_named_or_counted(mapped, tmp_path):
    path, _ = mapped
    assembly = _write_assembly(tmp_path / "x.fa", RAW["xx"])
    with pytest.raises(ValueError, match="5 contig.* not in .*: c1, c2, c3, c4, c5\\."):
        write_fasta(path, assembly, tmp_path / "out", unbinned=True)
    assert fasta_module._shown(list("abcdefgh")) == "a, b, c, d, e and 3 more"


def test_existing_files_need_force(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    out.mkdir()
    (out / "MAG_A.fa").write_text("old\n")
    (out / "keep.txt").write_text("mine\n")
    with pytest.raises(FileExistsError, match="1 output file\\(s\\) already exist"):
        write_fasta(path, assembly, out)
    assert (out / "MAG_A.fa").read_text() == "old\n"
    assert not (out / "MAG_B.fa").exists()

    write_fasta(path, assembly, out, force=True)
    assert (out / "MAG_A.fa").read_text() == RAW["c1"] + RAW["c2"]
    assert (out / "keep.txt").read_text() == "mine\n"
    assert sorted(p.name for p in out.iterdir()) == ["MAG_A.fa", "MAG_B.fa", "keep.txt"]


def test_one_existing_file_needs_force(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "all.fa"
    out.write_text("old\n")
    with pytest.raises(FileExistsError):
        write_fasta(path, assembly, out, single_file=True)
    write_fasta(path, assembly, out, single_file=True, force=True)
    assert out.read_text().startswith(">MAG_A|c1")


def test_an_input_is_never_an_output(mapped, tmp_path):
    path, assembly = mapped
    with pytest.raises(ValueError, match="same file as the input"):
        write_fasta(path, assembly, assembly, single_file=True, force=True)
    # A genome's file that is the assembly itself.
    inside = assembly.parent / "MAG_A.fa"
    inside.write_text(assembly.read_text())
    with pytest.raises(ValueError, match="same file as the input"):
        write_fasta(path, inside, assembly.parent, force=True)


def test_an_output_directory_that_is_a_file(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    out.write_text("")
    with pytest.raises(NotADirectoryError, match="--single-file"):
        write_fasta(path, assembly, out)


def test_a_failed_write_leaves_no_file(mapped, tmp_path, monkeypatch):
    path, assembly = mapped
    out = tmp_path / "out"
    out.mkdir()
    (out / "MAG_A.fa").write_text("old\n")

    def broken(self, handle, group, prefix=None):
        handle.write(b">partial\n")
        raise KeyboardInterrupt

    monkeypatch.setattr(fasta_module._Writer, "write", broken)
    with pytest.raises(KeyboardInterrupt):
        write_fasta(path, assembly, out, force=True)
    assert sorted(p.name for p in out.iterdir()) == ["MAG_A.fa"]
    assert (out / "MAG_A.fa").read_text() == "old\n"


@pytest.mark.parametrize(
    ("groups", "kwargs", "error", "message"),
    [
        ({"g": ["c1", "cX"]}, {}, KeyError, "not contigs of"),
        ({"g": ["c1", "c1"]}, {}, ValueError, "more than once: c1"),
        ({"g": []}, {}, ValueError, "has no contig"),
        ({"g": "c1"}, {}, ValueError, "not one text"),
        ({"": ["c1"]}, {}, ValueError, "needs a name"),
        ({}, {}, ValueError, "groups is empty"),
        (["c1"], {}, ValueError, "maps the name of each group"),
        ({"a/b": ["c1"]}, {}, ValueError, "Cannot name a file after 'a/b'"),
        ({"a b": ["c1"]}, {"stitch": True}, ValueError, "first whitespace"),
        ({"__unbinned__": ["c1"]}, {"unbinned": True}, ValueError, "rename it"),
        ({"g": ["c1"]}, {"genomes": ["MAG_A"]}, ValueError, "cannot be combined"),
        ({"g": ["c1"]}, {"spacer": -1, "stitch": True}, ValueError, "0 or more"),
        ({"g": ["c1"]}, {"spacer": True, "stitch": True}, ValueError, "number of N"),
        ({"g": ["c1"]}, {"spacer": 7}, ValueError, "give stitch=True too"),
        ({"Bin": ["c1"], "bin": ["c2"]}, {}, ValueError, "where case is ignored"),
    ],
)
def test_groups_and_arguments_that_are_refused(
    mapped, tmp_path, groups, kwargs, error, message
):
    path, assembly = mapped
    with pytest.raises(error, match=message):
        write_fasta(path, assembly, tmp_path / "out", groups=groups, **kwargs)


def test_a_group_name_with_a_slash_in_one_file(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "all.fa"
    write_fasta(path, assembly, out, groups={"a/b": ["c3"]}, single_file=True)
    assert out.read_text() == f">a/b|{RAW['c3'][1:]}"


def test_the_lazy_export():
    assert metacooler.write_fasta is write_fasta


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_cli_writes_the_files(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli, ["fasta", str(path), "--fasta", str(assembly), "-o", str(out), "-q"]
    )
    assert result.exit_code == 0, result.output
    assert (out / "MAG_A.fa").read_text() == RAW["c1"] + RAW["c2"]


def test_cli_with_every_option(clustered, tmp_path):
    path, assembly = clustered
    out = tmp_path / "all.fa"
    args = [
        "fasta",
        str(path),
        "--fasta",
        str(assembly),
        "-o",
        str(out),
        "--round",
        "final",
        "--genomes",
        "F1,F2",
        "--filter",
        "clusters.completeness >= 50",
        "--unbinned",
        "--stitch",
        "--spacer",
        "200",
        "--single-file",
    ]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert out.read_text().startswith(">F1\n")
    again = CliRunner().invoke(cli, args)
    assert again.exit_code == 1
    assert "already exists" in flat(again.output)
    assert CliRunner().invoke(cli, [*args, "-f"]).exit_code == 0


def test_cli_refusals(clustered, tmp_path):
    path, assembly = clustered
    base = ["fasta", str(path), "--fasta", str(assembly), "-o", str(tmp_path / "o")]
    spacer = CliRunner().invoke(cli, [*base, "--spacer", "10"])
    assert spacer.exit_code == 2
    assert "add --stitch" in flat(spacer.output)
    table = CliRunner().invoke(cli, [*base, "--filter", "bins.start > 1"])
    assert table.exit_code == 2
    assert "unknown table 'bins'" in flat(table.output)
    missing = CliRunner().invoke(cli, [*base, "--round", "nope"])
    assert missing.exit_code == 1
    assert "is not a clustering round" in flat(missing.output)


def test_other_commands_refuse_a_clusters_filter(mapped):
    path, _ = mapped
    result = CliRunner().invoke(
        cli, ["dump", str(path), "-t", "genomes", "--filter", "clusters.length > 1"]
    )
    assert result.exit_code == 2
    assert "only 'fasta --round' reads" in flat(result.output)


def test_every_file_is_written_before_any_is_renamed(mapped, tmp_path, monkeypatch):
    """Each file is staged as a short ``.fasta-*.tmp`` beside it; the renames
    come once every file is written."""
    path, assembly = mapped
    out = tmp_path / "out"
    seen = []
    replace = os.replace

    def spy(src, dst):
        if not seen:  # every staged file exists at the first rename
            seen.extend(sorted(p.name for p in out.iterdir()))
        replace(src, dst)

    monkeypatch.setattr(fasta_module.os, "replace", spy)
    write_fasta(path, assembly, out, unbinned=True)
    assert len(seen) == 3
    assert all(name.startswith(".fasta-") and name.endswith(".tmp") for name in seen)
    assert sorted(p.name for p in out.iterdir()) == [
        "MAG_A.fa",
        "MAG_B.fa",
        "__unbinned__.fa",
    ]
    mode = (out / "MAG_A.fa").stat().st_mode & 0o777
    mask = os.umask(0)
    os.umask(mask)
    assert mode == 0o666 & ~mask


def test_a_failure_in_the_second_file_writes_none(mapped, tmp_path, monkeypatch):
    path, assembly = mapped
    out = tmp_path / "out"
    out.mkdir()
    (out / "MAG_B.fa").write_text("old\n")
    write = fasta_module._Writer.write
    calls = []

    def second_fails(self, handle, group, prefix=None):
        calls.append(group.name)
        if len(calls) == 2:
            raise OSError(28, "No space left on device")
        write(self, handle, group, prefix)

    monkeypatch.setattr(fasta_module._Writer, "write", second_fails)
    with pytest.raises(OSError, match="No space left"):
        write_fasta(path, assembly, out, force=True)
    # MAG_A was written, but not renamed: no new file, MAG_B as it was.
    assert sorted(p.name for p in out.iterdir()) == ["MAG_B.fa"]
    assert (out / "MAG_B.fa").read_text() == "old\n"


def test_a_long_name_is_staged_under_a_short_one(mapped, tmp_path):
    path, assembly = mapped
    name = "g" * 252  # 255 bytes with .fa: the longest file name
    written = write_fasta(path, assembly, tmp_path / "out", groups={name: ["c3"]})
    assert written[name].name == name + ".fa"
    assert written[name].read_text() == RAW["c3"]


def test_names_that_differ_by_case_alone_go_to_one_file(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "all.fa"
    write_fasta(
        path, assembly, out, groups={"Bin": ["c1"], "bin": ["c2"]}, single_file=True
    )
    assert out.read_text() == f">Bin|{RAW['c1'][1:]}>bin|{RAW['c2'][1:]}"


def test_one_file_in_a_new_directory(mapped, tmp_path):
    path, assembly = mapped
    out = tmp_path / "new" / "dir" / "all.fa"
    write_fasta(path, assembly, out, single_file=True)
    assert out.read_text().startswith(">MAG_A|c1 desc one\n")


def test_an_output_file_that_is_a_directory(mapped, tmp_path):
    path, assembly = mapped
    (tmp_path / "out" / "MAG_A.fa").mkdir(parents=True)
    with pytest.raises(IsADirectoryError):
        write_fasta(path, assembly, tmp_path / "out", force=True)


def test_the_api_takes_one_genome_or_filter_as_text(mapped, tmp_path):
    path, assembly = mapped
    assert list(write_fasta(path, assembly, tmp_path / "a", genomes="MAG_B")) == [
        "MAG_B"
    ]
    one = write_fasta(path, assembly, tmp_path / "b", filters="genomes.length > 6kb")
    assert list(one) == ["MAG_A"]


def test_the_fasta_filter_reads_the_clusters_table():
    from metacooler.cli.cmd_fasta import cmd_fasta
    from metacooler.filters import TABLES_WITH_CLUSTERS

    (option,) = [p for p in cmd_fasta.params if p.name == "filters"]
    callback = option.callback
    assert isinstance(callback, functools.partial)
    assert callback.keywords["tables"] == TABLES_WITH_CLUSTERS
