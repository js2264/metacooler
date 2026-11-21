"""The QC block survives a reassignment.

The guarantee: ``from_pairs`` with assignment A, then ``reassign`` to
assignment B, stores the same QC block -- every key -- and the same ``qc_*``
columns in ``/genomes`` and ``/contigs`` as ``from_pairs`` with B directly.
That holds because the two intra-contig categories do not depend on the
assignment and are carried per contig (their split was measured on read
positions, which no pixel can recover), while the two inter-contig categories
are re-measured from the pixels, which classify them exactly at any bin size.

The dataset below has all of it: intra-contig pairs below and above the
distance threshold, inter-contig pairs within and across genomes, unbinned
contigs under both assignments, and pairs ``create pairs`` discards.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Literal

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io import reassign as reassign_module
from metacooler.io.create import from_pairs
from metacooler.io.reassign import reassign
from metacooler.io.staging import copy_fresh
from metacooler.ops.qc import qc
from metacooler.qc import (
    DEFAULT_DISTANCE_THRESHOLD,
    QC_CONTIG_COLUMNS,
    QC_GENOME_COLUMNS,
    QC_METADATA_KEY,
    ContactQC,
)
from metacooler.schema import Layout, read_table, validate

# Distinct lengths, so every genome-then-length layout is fully determined.
LENGTHS = {
    "c01": 9_000,
    "c02": 8_000,
    "c03": 7_000,
    "c04": 6_000,
    "c05": 5_000,
    "c06": 4_000,
    "c07": 3_000,
    "c08": 2_500,
    "c09": 2_000,
    "c10": 1_500,
}
THRESHOLD = 2_500

# A: an earlier binning, which leaves c06..c10 unbinned.
ASSIGNMENT_A = {
    "c01": "MAG_1",
    "c02": "MAG_1",
    "c03": "MAG_2",
    "c04": "MAG_2",
    "c05": "MAG_3",
}
# B: the final binning.  Every genome changes, MAG_2 takes contigs that sit
# apart under A's layout, and c05, c08 and c10 are left unbinned.
ASSIGNMENT_B = {
    "c01": "MAG_2",
    "c02": "MAG_4",
    "c03": "MAG_1",
    "c04": "MAG_1",
    "c06": "MAG_2",
    "c07": "MAG_3",
    "c09": "MAG_3",
}

# One bin per contig, as a map built to be binned by MetaTOR; and 1 kb bins,
# finer than the threshold, so the pixels hold distances too.
BINNINGS = ("one-bin-per-contig", "1kb")
ONE_BIN_PER_CONTIG = pd.DataFrame(
    {"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())}
)

Ordering = Literal["genome-then-length", "genome-then-alpha", "fasta", "alpha"]


@pytest.fixture(scope="module")
def pairs_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Pairs covering every QC category, plus five the binning discards."""
    rng = np.random.default_rng(20261117)
    names = list(LENGTHS)
    rows = []
    for _ in range(4_000):
        c1 = names[rng.integers(len(names))]
        c2 = c1 if rng.random() < 0.6 else names[rng.integers(len(names))]
        rows.append(
            (
                c1,
                int(rng.integers(1, LENGTHS[c1])),
                c2,
                int(rng.integers(1, LENGTHS[c2])),
            )
        )
    rows += [("cX", 10, "c01", 20)] * 3  # an unknown contig
    rows += [("c01", LENGTHS["c01"] + 50, "c02", 10)] * 2  # past the contig end
    path = tmp_path_factory.mktemp("pairs") / "lib.pairs"
    path.write_text(
        "## pairs format v1.0\n"
        + "".join(f"#chromsize: {c} {n}\n" for c, n in LENGTHS.items())
        + "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(f"r{i}\t{a}\t{b}\t{c}\t{d}\n" for i, (a, b, c, d) in enumerate(rows))
    )
    return path


def _create(
    path: Path,
    pairs: Path,
    assignment: dict[str, str] | None,
    binning: str,
    ordering: Ordering = "genome-then-length",
) -> Path:
    one_bin = binning == "one-bin-per-contig"
    from_pairs(
        path,
        bins=ONE_BIN_PER_CONTIG if one_bin else None,
        binsize=None if one_bin else 1_000,
        pairs=pairs,
        contig_genome_map=assignment,
        contig_ordering=ordering,
        qc_distance_threshold=THRESHOLD,
        ps=False,
    )
    return path


def _qc_state(path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """The QC block and the qc_* columns of both tables, keyed by name."""
    mc = MetaCooler(path)
    genomes = mc.genomes.set_index("name")[list(QC_GENOME_COLUMNS)].sort_index()
    contigs = mc.contigs.set_index("name")[list(QC_CONTIG_COLUMNS)].sort_index()
    return mc.metadata[QC_METADATA_KEY], genomes, contigs


def _assert_same_qc(result: Path, expected: Path) -> None:
    block, genomes, contigs = _qc_state(result)
    want_block, want_genomes, want_contigs = _qc_state(expected)
    assert block == want_block
    pd.testing.assert_frame_equal(genomes, want_genomes)
    pd.testing.assert_frame_equal(contigs, want_contigs)


def _assert_sound(path: Path | str) -> None:
    """The output passes every check, and cooler reads it."""
    validate(path, full=True)
    clr = cooler.Cooler(str(path))
    assert len(clr.chroms()[:]) == len(LENGTHS)
    assert int(clr.pixels()[:]["count"].sum()) == clr.info["sum"]


# ---------------------------------------------------------------------------
# reassign(from_pairs(A), B) == from_pairs(B)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("binning", BINNINGS)
def test_unbinned_map_reassigned_matches_a_direct_creation(
    tmp_path: Path, pairs_file: Path, binning: str
):
    """MetaTOR's route: create with no binning yet, then reassign the final bins.

    The source is genome-grouped, so the contigs are re-laid out (and the bins
    and pixels permuted), as ``from_pairs`` with B lays them out directly.
    """
    source = _create(tmp_path / "a.metacool", pairs_file, None, binning)
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, binning)
    out = tmp_path / "out.metacool"

    reassign(source, ASSIGNMENT_B, out, missing="unbinned")

    block = MetaCooler(out).metadata[QC_METADATA_KEY]
    assert block["source"] == "pairs"
    assert block["n_pairs_discarded"] == 5
    assert block["intra_contig_short"] > 0
    assert block["intra_contig_long"] > 0
    assert block["inter_contig_intra_genome"] > 0
    assert block["inter_genome"] > 0
    assert MetaCooler(out).contigs["name"].tolist() == (
        MetaCooler(direct).contigs["name"].tolist()
    )
    _assert_same_qc(out, direct)
    _assert_sound(out)


@pytest.mark.parametrize("binning", BINNINGS)
@pytest.mark.parametrize(
    ("source_ordering", "layout", "direct_ordering"),
    [
        # Relabelled: the contigs keep their (alphabetical) rows.
        ("alpha", "auto", "alpha"),
        # Re-grouped by force: contigs, bins and pixels are permuted.
        ("alpha", "genome-then-length", "genome-then-length"),
        # Relabelled although genome-grouped: written as contig-ordering=custom.
        ("genome-then-length", "keep", "genome-then-length"),
    ],
)
def test_binned_map_reassigned_matches_a_direct_creation(
    tmp_path: Path,
    pairs_file: Path,
    binning: str,
    source_ordering: Ordering,
    layout: Layout,
    direct_ordering: Ordering,
):
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, binning, source_ordering
    )
    direct = _create(
        tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, binning, direct_ordering
    )
    out = tmp_path / "out.metacool"

    reassign(source, ASSIGNMENT_B, out, missing="unbinned", layout=layout)

    _assert_same_qc(out, direct)
    _assert_sound(out)


def test_inplace_reassignment_keeps_the_block(tmp_path: Path, pairs_file: Path):
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb", "alpha")
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, "1kb", "alpha")

    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")

    _assert_same_qc(source, direct)
    _assert_sound(source)


def test_an_entry_of_a_collection_keeps_the_block(tmp_path: Path, pairs_file: Path):
    from metacooler.io.mmetacool import create_mmetacool

    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb", "alpha")
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, "1kb", "alpha")
    collection = tmp_path / "c.mmetacool"
    create_mmetacool(collection, {"rep1": source})
    uri = f"{collection}::entries/rep1"

    reassign(uri, ASSIGNMENT_B, inplace=True, missing="unbinned")

    mc = MetaCooler(uri)
    assert mc.metadata[QC_METADATA_KEY] == MetaCooler(direct).metadata[QC_METADATA_KEY]
    validate(uri, full=True)


def test_the_cli_keeps_the_block(tmp_path: Path, pairs_file: Path):
    source = _create(tmp_path / "a.metacool", pairs_file, None, "1kb")
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, "1kb")
    table = tmp_path / "bins.tsv"
    table.write_text(
        "contig_name\tgenome_name\n"
        + "".join(f"{c}\t{g}\n" for c, g in ASSIGNMENT_B.items())
    )
    out = tmp_path / "out.metacool"

    result = CliRunner().invoke(
        cli,
        ["reassign", str(source), "--contigs", str(table), "-o", str(out)]
        + ["--missing", "unbinned"],
    )

    assert result.exit_code == 0, result.output
    _assert_same_qc(out, direct)


def test_the_block_keeps_every_key_it_does_not_derive(tmp_path: Path, pairs_file: Path):
    """``source``, ``distance_threshold`` and the pairs counts are facts about
    how the counts were measured, which a reassignment does not change."""
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb")
    with h5py.File(source, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        metadata[QC_METADATA_KEY]["note"] = "kept"
        f.attrs["metadata"] = json.dumps(metadata)
    out = tmp_path / "out.metacool"

    reassign(source, ASSIGNMENT_B, out, missing="unbinned")

    before = MetaCooler(source).metadata[QC_METADATA_KEY]
    after = MetaCooler(out).metadata[QC_METADATA_KEY]
    for key in ("source", "distance_threshold", "binsize", "n_pairs_read"):
        assert after[key] == before[key], key
    assert after["n_pairs_discarded"] == before["n_pairs_discarded"]
    assert after["note"] == "kept"
    assert after["n_contacts"] == before["n_contacts"]


# ---------------------------------------------------------------------------
# A block that cannot be carried is recomputed from the pixels
# ---------------------------------------------------------------------------


def _drop_contig_columns(path: Path) -> None:
    """Make *path* look like a map written by metacooler 0.2.0."""
    with h5py.File(path, "r+") as f:
        for column in QC_CONTIG_COLUMNS:
            del f["contigs"][column]


def _shift_one_count(path: Path) -> None:
    """Leave per-contig counts that no longer add up to the block."""
    with h5py.File(path, "r+") as f:
        column = f["contigs"][QC_CONTIG_COLUMNS[0]]
        values = column[:]
        values[0] += 1
        column[...] = values


def _edit_block(path: Path, **changes: object) -> None:
    """Set (or, given ``None``, delete) keys of *path*'s QC block."""
    with h5py.File(path, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        for key, value in changes.items():
            if value is None:
                del metadata[QC_METADATA_KEY][key]
            else:
                metadata[QC_METADATA_KEY][key] = value
        f.attrs["metadata"] = json.dumps(metadata)


def _set_version(path: Path) -> None:
    _edit_block(path, version=1)


def _drop_a_total(path: Path) -> None:
    _edit_block(path, intra_contig_long=None)


def _set_threshold_text(path: Path) -> None:
    _edit_block(path, distance_threshold="2.5 kb")


def _set_not_an_object(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        metadata[QC_METADATA_KEY] = "not a block"
        f.attrs["metadata"] = json.dumps(metadata)


@pytest.mark.parametrize(
    ("damage", "reason", "threshold"),
    [
        (_drop_contig_columns, "the source has no per-contig counts", THRESHOLD),
        (_shift_one_count, "do not add up", THRESHOLD),
        (_set_version, "version 1 block", THRESHOLD),
        (_drop_a_total, "do not add up", THRESHOLD),
        # No threshold to read back: the default one, as before.
        (
            _set_threshold_text,
            "no integer distance threshold",
            DEFAULT_DISTANCE_THRESHOLD,
        ),
        (_set_not_an_object, "not a JSON object", DEFAULT_DISTANCE_THRESHOLD),
    ],
)
@pytest.mark.parametrize("ordering", ["genome-then-length", "alpha"])
def test_a_block_that_cannot_be_carried_is_recomputed_from_the_pixels(
    tmp_path: Path,
    pairs_file: Path,
    caplog: pytest.LogCaptureFixture,
    damage,
    reason: str,
    threshold: int,
    ordering: Ordering,
):
    """Re-grouped or (``alpha``) relabelled in place of the source's rows."""
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb", ordering)
    damage(source)
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, "1kb")
    out = tmp_path / "out.metacool"

    with caplog.at_level(logging.INFO, logger="metacooler"):
        reassign(source, ASSIGNMENT_B, out, missing="unbinned")

    record = next(r for r in caplog.records if "Recomputing the source" in r.message)
    assert reason in record.getMessage()
    if damage is _set_not_an_object:
        assert record.levelname == "INFO"
    else:
        # A block measured on read positions: what goes is named, as a warning.
        assert record.levelname == "WARNING"
        assert "intra_contig_long drops to 0" in record.getMessage()
        assert "n_pairs_read" in record.getMessage()
    # What 'metacooler qc --store' measures on the same map, built directly.
    expected = qc(direct, distance_threshold=threshold)
    block, genomes, contigs = _qc_state(out)
    assert block["source"] == "pixels"
    assert block == expected.summary
    pd.testing.assert_frame_equal(
        genomes,
        expected.per_genome.set_index("name")[list(QC_GENOME_COLUMNS)].sort_index(),
    )
    # Recomputed, the block gains the per-contig counts it lacked.
    assert contigs.sum().tolist() == [
        block["intra_contig_short"],
        block["intra_contig_long"],
    ]
    _assert_sound(out)


def test_a_block_from_pixels_is_recomputed_quietly(
    tmp_path: Path, pairs_file: Path, caplog: pytest.LogCaptureFixture
):
    """Guard: recomputing a block measured on pixels loses nothing; INFO only."""
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb")
    qc(source, distance_threshold=THRESHOLD, store=True)
    _drop_contig_columns(source)

    with caplog.at_level(logging.INFO, logger="metacooler"):
        reassign(source, ASSIGNMENT_B, tmp_path / "out.metacool", missing="unbinned")

    record = next(r for r in caplog.records if "Recomputing the source" in r.message)
    assert record.levelname == "INFO"


# ---------------------------------------------------------------------------
# A run stopped during the QC pass leaves the source's QC as it was
# ---------------------------------------------------------------------------


def _interrupt_the_tally(monkeypatch: pytest.MonkeyPatch) -> None:
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(ContactQC, "tally_pixels", interrupted)


def test_an_interrupted_in_place_reassignment_leaves_the_source_as_it_was(
    tmp_path: Path, pairs_file: Path, monkeypatch: pytest.MonkeyPatch
):
    """The block is measured before anything is written, in place included.

    It used to be written after the relabelling: stopped in between, the map
    kept neither its block nor its per-contig counts, and the next run fell
    back to bin midpoints for good.
    """
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb", "alpha")
    before = source.read_bytes()

    with monkeypatch.context() as patch:
        _interrupt_the_tally(patch)
        with pytest.raises(KeyboardInterrupt):
            reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")

    assert source.read_bytes() == before
    # So the next run still carries the split measured on read positions.
    direct = _create(tmp_path / "b.metacool", pairs_file, ASSIGNMENT_B, "1kb", "alpha")
    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")
    _assert_same_qc(source, direct)


def test_an_interrupted_reassignment_writes_nothing(
    tmp_path: Path, pairs_file: Path, monkeypatch: pytest.MonkeyPatch
):
    """Nor does a run to a new output: no partial file, the old output intact."""
    source = _create(tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "1kb")
    out = tmp_path / "out.metacool"
    out.write_bytes(b"previous output")
    _interrupt_the_tally(monkeypatch)

    with pytest.raises(KeyboardInterrupt):
        reassign(source, ASSIGNMENT_B, out, missing="unbinned", force=True)

    assert out.read_bytes() == b"previous output"
    assert not out.with_name(out.name + ".partial").exists()


# ---------------------------------------------------------------------------
# The relabel path writes the tables as a direct creation does
# ---------------------------------------------------------------------------


def _contigs_table(path: Path, assignment: dict[str, str]) -> Path:
    """A contigs table of *assignment*, with an extra column, ``gc``."""
    path.write_text(
        "contig_name\tgenome_name\tgc\n"
        + "".join(
            f"{contig}\t{assignment.get(contig, 'NA')}\t0.{i}\n"
            for i, contig in enumerate(LENGTHS, start=1)
        )
    )
    return path


def _column_order(path: Path) -> dict[str, list[str]]:
    return {table: list(read_table(path, table).columns) for table in TABLES}


TABLES = ("genomes", "contigs")


@pytest.mark.parametrize("inplace", [False, True], ids=["new-file", "in-place"])
def test_relabelled_tables_keep_the_column_order_of_a_direct_creation(
    tmp_path: Path, pairs_file: Path, inplace: bool
):
    """``genome_id`` stays third, and every column where ``create`` puts it.

    The relabel path deleted ``genome_id``, and every QC column it changed, and
    wrote them again in a group listed in creation order: they moved after the
    other columns, where ``create`` and the regroup path put ``genome_id``
    third. The values were right; the tables differed only in their order.
    """
    from_pairs(
        tmp_path / "a.metacool",
        bins=ONE_BIN_PER_CONTIG,
        pairs=pairs_file,
        contig_genome_map=_contigs_table(tmp_path / "a.tsv", ASSIGNMENT_A),
        contig_ordering="fasta",
        qc_distance_threshold=THRESHOLD,
        ps=False,
    )
    from_pairs(
        tmp_path / "direct.metacool",
        bins=ONE_BIN_PER_CONTIG,
        pairs=pairs_file,
        contig_genome_map=_contigs_table(tmp_path / "b.tsv", ASSIGNMENT_B),
        contig_ordering="fasta",
        qc_distance_threshold=THRESHOLD,
        ps=False,
    )
    source = tmp_path / "a.metacool"
    out = source if inplace else tmp_path / "out.metacool"
    reassign(
        source,
        ASSIGNMENT_B,
        None if inplace else out,
        inplace=inplace,
        missing="unbinned",
    )

    order = _column_order(out)
    assert order["contigs"][:3] == ["name", "length", "genome_id"]
    assert order == _column_order(tmp_path / "direct.metacool")
    for table in TABLES:
        pd.testing.assert_frame_equal(
            read_table(out, table), read_table(tmp_path / "direct.metacool", table)
        )
    assert ".reassign-scratch" not in h5py.File(out, "r")
    _assert_sound(out)


def test_relabelling_keeps_the_columns_it_leaves_alone_as_stored(
    tmp_path: Path, pairs_file: Path
):
    """A column the relabelling does not change keeps its storage too.

    A guard: the relabel path left these columns where they were before it
    rebuilt the table, so this passes without the rebuild too; it holds the
    rebuild to keeping them rather than writing them again.
    """
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    with h5py.File(source, "r+") as f:
        # Stored otherwise than metacooler writes it: no chunks, no filter.
        f["contigs"].create_dataset("depth", data=np.arange(len(LENGTHS)) * 1.5)
    out = tmp_path / "out.metacool"
    reassign(source, ASSIGNMENT_B, out, missing="unbinned")

    with h5py.File(source, "r") as old, h5py.File(out, "r") as new:
        for column in ("name", "length", "depth"):
            stored, copied = old["contigs"][column], new["contigs"][column]
            assert copied.dtype == stored.dtype
            assert copied.chunks == stored.chunks
            assert copied.compression == stored.compression
            assert np.array_equal(copied[()], stored[()])


def test_a_failed_relabelling_in_place_leaves_the_tables_as_they_were(
    tmp_path: Path, pairs_file: Path, monkeypatch: pytest.MonkeyPatch
):
    """The tables are built aside and swapped in together, once both are built."""
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    before = {table: read_table(source, table) for table in TABLES}

    write_column = reassign_module.write_column

    def interrupted(group, name, *args, **kwargs):
        if name == "genome_id":
            raise KeyboardInterrupt
        return write_column(group, name, *args, **kwargs)

    # Stopped while writing the contigs' genome ids: once /genomes is written.
    monkeypatch.setattr(reassign_module, "write_column", interrupted)
    with pytest.raises(KeyboardInterrupt):
        reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")

    for table in TABLES:
        pd.testing.assert_frame_equal(read_table(source, table), before[table])
    assert ".reassign-scratch" not in h5py.File(source, "r")
    _assert_sound(source)


def test_relabelling_in_place_clears_a_scratch_group_left_behind(
    tmp_path: Path, pairs_file: Path
):
    """A run killed while building the tables leaves its scratch group: the next
    run removes it before building its own."""
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    with h5py.File(source, "r+") as f:
        f.create_group(".reassign-scratch/contigs")

    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")

    assert ".reassign-scratch" not in h5py.File(source, "r")
    _assert_sound(source)


def test_relabelling_in_place_links_the_columns_it_leaves_alone(
    tmp_path: Path, pairs_file: Path
):
    """The columns the relabelling leaves alone are the stored datasets, linked
    into the rebuilt table, not copies: a copy leaves the old dataset behind as
    dead space HDF5 does not reclaim, so every relabelling in place grew the
    file by its contigs' columns. A relabelling in place is now written to a
    fresh copy of the map, which holds none of the dead space of earlier
    edits, so the map holds at most the space its own edit left: no more than
    a fresh copy of it, but for the small tables it rebuilt."""
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    rng = np.random.default_rng(0)
    notes = rng.integers(ord("a"), ord("z") + 1, size=(len(LENGTHS), 100_000))
    with h5py.File(source, "r+") as f:
        # 500 kB that compression cannot shrink, stored contiguous.
        f["contigs"].create_dataset(
            "notes", data=notes.astype(np.uint8).view("S100000")[:, 0]
        )
    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")
    size = source.stat().st_size

    for assignment in (ASSIGNMENT_A, ASSIGNMENT_B):
        reassign(source, assignment, inplace=True, missing="unbinned")
        fresh = tmp_path / "fresh.metacool"
        copy_fresh(source, fresh)
        # Not 500 kB of notes left behind by the edit.
        assert source.stat().st_size - fresh.stat().st_size < 100_000
        fresh.unlink()
    assert source.stat().st_size - size < 100_000
    _assert_sound(source)


def test_the_group_attributes_of_the_tables_are_kept(tmp_path: Path, pairs_file: Path):
    """Rebuilt, ``/contigs`` and ``/genomes`` keep the attributes of the groups
    they replace, as ``annotate`` keeps them."""
    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    with h5py.File(source, "r+") as f:
        f["contigs"].attrs["provenance"] = "assembly v2"
        f["genomes"].attrs["provenance"] = "binning v1"
    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")
    with h5py.File(source, "r") as f:
        assert f["contigs"].attrs["provenance"] == "assembly v2"
        assert f["genomes"].attrs["provenance"] == "binning v1"


def test_ctrl_c_while_the_tables_are_swapped_in_is_held_until_they_are(
    tmp_path: Path, pairs_file: Path, monkeypatch: pytest.MonkeyPatch
):
    """A SIGINT that arrives while the new tables replace the old ones leaves
    the map as it was, whole: the tables are swapped in a staged copy of the
    map, which only replaces it once it validates. Swapped in the map itself,
    an interrupt delivered at once left a map whose ``/genomes`` was the new
    one and whose ``ngenomes`` the old, which no command could open (a SIGINT
    held until the swap was done fixed that, before the swap was staged)."""
    import signal

    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    expected = tmp_path / "b.metacool"
    shutil.copy2(source, expected)
    move = h5py.Group.move

    def interrupted(self, source_name, dest_name):
        monkeypatch.setattr(h5py.Group, "move", move)
        move(self, source_name, dest_name)
        signal.raise_signal(signal.SIGINT)

    monkeypatch.setattr(h5py.Group, "move", interrupted)
    with pytest.raises(KeyboardInterrupt):
        reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")

    assert signal.getsignal(signal.SIGINT) is signal.default_int_handler
    _assert_sound(source)
    assert source.read_bytes() == expected.read_bytes()
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith(".")) == []


def test_relabelling_in_place_from_a_thread(tmp_path: Path, pairs_file: Path):
    """Only the main thread may set a signal handler: a relabelling in place
    sets none, and runs from any thread (it once held SIGINT during the swap,
    and failed with ``signal only works in main thread`` there). A guard."""
    import threading

    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    errors: list[BaseException] = []

    def run() -> None:
        try:
            reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")
        except BaseException as exc:  # noqa: BLE001 - reported below
            errors.append(exc)

    thread = threading.Thread(target=run)
    thread.start()
    thread.join()
    assert errors == []
    _assert_sound(source)


def test_relabelling_in_place_under_a_handler_python_did_not_install(
    tmp_path: Path, pairs_file: Path, monkeypatch: pytest.MonkeyPatch
):
    """A host that embeds Python (R through reticulate, say) may set the
    SIGINT handler from C: ``getsignal`` then returns ``None``, which
    ``signal.signal`` cannot restore. metacooler sets no handler, so the
    host's is left alone. A guard."""
    import signal

    source = _create(
        tmp_path / "a.metacool", pairs_file, ASSIGNMENT_A, "one-bin-per-contig", "fasta"
    )
    monkeypatch.setattr(signal, "getsignal", lambda signum: None)

    def refused(*args):
        raise AssertionError("the SIGINT handler must be left alone")

    monkeypatch.setattr(signal, "signal", refused)
    reassign(source, ASSIGNMENT_B, inplace=True, missing="unbinned")
    _assert_sound(source)
