"""Tests for ``metacooler reassign`` — re-attributing contigs to genomes.

The invariant the whole feature rests on is that a reassignment is a
*relabelling*: contigs, their lengths, their bins and every contact are those
of the source file, so the contact matrix must come out identical (relabel
path) or an exact permutation of itself (re-grouping path).  Most of what is
asserted here is that, plus the two things the format makes non-trivial — the
derived ``/genomes`` table and the physical contig layout that
``genome_offset`` and ``genome_extent`` depend on.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from _cools import BSUB, ECOLI, make_cool
from click.testing import CliRunner
from numpy.testing import assert_array_equal

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.reassign import CUSTOM_ORDERING, reassign
from metacooler.qc import QC_CONTIG_COLUMNS
from metacooler.schema import Layout, validate

# The tiny fixture is built with contig-ordering="genome-then-alpha" from:
#   __unbinned__ = {contig_D}, MAG_1 = {contig_A, contig_B}, MAG_2 = {contig_C}
# so its contigs are laid out [contig_D, contig_A, contig_B, contig_C].
SOURCE_ORDER = ["contig_D", "contig_A", "contig_B", "contig_C"]


@pytest.fixture
def runner():
    return CliRunner()


def _table(path: Path, rows: dict[str, str]) -> Path:
    """Write a ``contigs``-format reassignment table."""
    lines = ["contig_name\tgenome_name\n"]
    lines += [f"{contig}\t{genome}\n" for contig, genome in rows.items()]
    path.write_text("".join(lines))
    return path


def _dense(mc: MetaCooler) -> np.ndarray:
    return mc.matrix(balance=False, sparse=True)[:, :].toarray()


def _bin_permutation(source: MetaCooler, target: MetaCooler) -> np.ndarray:
    """Return ``perm`` with ``target.bins()[i] == source.bins()[perm[i]]``."""
    src = source.bins()
    dst = target.bins()
    index = {
        key: i
        for i, key in enumerate(
            zip(src["chrom"].astype(str), src["start"], strict=True)
        )
    }
    return np.array(
        [index[key] for key in zip(dst["chrom"].astype(str), dst["start"], strict=True)]
    )


# ---------------------------------------------------------------------------
# The relabel path: nothing moves
# ---------------------------------------------------------------------------


def test_relabel_leaves_layout_and_contacts_untouched(
    tmp_path: Path, tiny_metacooler: Path
):
    """Merging MAG_2 into MAG_1 needs no re-ordering, so nothing but labels move."""
    source = MetaCooler(tiny_metacooler)
    before = _dense(source)
    out = tmp_path / "out.metacool"

    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), out)

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == SOURCE_ORDER
    assert_array_equal(_dense(result), before)
    assert result.info["ordering"] == "genome-then-alpha"
    assert result.genomes["name"].tolist() == ["__unbinned__", "MAG_1"]


def test_relabel_rebuilds_the_derived_genomes_table(
    tmp_path: Path, tiny_metacooler: Path
):
    """Every derived /genomes column follows the new assignment, not the old one."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), out)

    genomes = MetaCooler(out).genomes
    mag1 = genomes.loc[genomes["name"] == "MAG_1"].iloc[0]
    # contig_A (1000) + contig_B (600) + contig_C (400), at 100 bp bins.
    assert int(mag1["length"]) == 2000
    assert int(mag1["n_contigs"]) == 3
    assert int(mag1["n_bins"]) == 20


def test_relabel_renumbers_genome_ids_canonically(
    tmp_path: Path, tiny_metacooler: Path
):
    """A new name sorts into place and every genome_id is renumbered around it."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "AAA_new"}), out)

    result = MetaCooler(out)
    # __unbinned__ always takes index 0; the rest are alphabetical.
    assert result.genomes["name"].tolist() == ["__unbinned__", "AAA_new", "MAG_1"]
    contigs = result.contigs.set_index("name")
    assert int(contigs.at["contig_C", "genome_id"]) == 1
    assert int(contigs.at["contig_A", "genome_id"]) == 2


def test_relabel_preserves_the_pixel_table_byte_for_byte(
    tmp_path: Path, tiny_metacooler: Path
):
    """The relabel path must not rewrite /pixels, /bins or their indexes."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), out)

    with h5py.File(tiny_metacooler, "r") as src, h5py.File(out, "r") as dst:
        for group, column in (
            ("pixels", "bin1_id"),
            ("pixels", "bin2_id"),
            ("pixels", "count"),
            ("bins", "chrom"),
            ("bins", "start"),
            ("indexes", "chrom_offset"),
            ("indexes", "bin1_offset"),
        ):
            assert_array_equal(src[group][column][:], dst[group][column][:])


# ---------------------------------------------------------------------------
# The re-grouping path: contigs move, so pixels are remapped
# ---------------------------------------------------------------------------


def test_regroup_permutes_the_matrix_without_changing_a_contact(
    tmp_path: Path, tiny_metacooler: Path
):
    """Moving contig_A to MAG_2 breaks the layout, so the file is re-grouped."""
    source = MetaCooler(tiny_metacooler)
    out = tmp_path / "out.metacool"

    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == [
        "contig_D",
        "contig_B",
        "contig_A",
        "contig_C",
    ]
    assert result.info["sum"] == source.info["sum"]
    perm = _bin_permutation(source, result)
    assert_array_equal(_dense(result), _dense(source)[np.ix_(perm, perm)])


def test_regroup_keeps_genomes_contiguous(tmp_path: Path, tiny_metacooler: Path):
    """The point of re-grouping: genome_extent still answers."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    result = MetaCooler(out)
    assert result.info["ordering"] == "genome-then-alpha"
    lo, hi = result.genome_extent("MAG_2")
    # contig_A (10 bins) + contig_C (4 bins), laid out together.
    assert hi - lo == 14
    assert_array_equal(result.genome_bin_ids("MAG_2"), np.arange(lo, hi))


def test_regroup_recomputes_the_derived_genome_columns(
    tmp_path: Path, tiny_metacooler: Path
):
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    genomes = MetaCooler(out).genomes.set_index("name")
    # MAG_2 gains contig_A (1000 bp, 10 bins) on top of contig_C (400 bp, 4).
    assert int(genomes.at["MAG_2", "length"]) == 1400
    assert int(genomes.at["MAG_2", "n_contigs"]) == 2
    assert int(genomes.at["MAG_2", "n_bins"]) == 14
    assert int(genomes.at["MAG_1", "n_bins"]) == 6


def test_regroup_carries_extra_bin_columns(tmp_path: Path, tiny_metacooler: Path):
    """A balancing weight describes its bin, so it must ride the permutation."""
    source_path = tmp_path / "weighted.metacool"
    shutil.copy2(tiny_metacooler, source_path)
    source = MetaCooler(source_path)
    weight = np.arange(1, len(source.bins()) + 1, dtype=np.float64) / 10
    with h5py.File(source_path, "r+") as f:
        f["bins"].create_dataset("weight", data=weight)

    out = tmp_path / "out.metacool"
    reassign(source_path, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    result = MetaCooler(out)
    perm = _bin_permutation(MetaCooler(source_path), result)
    assert "weight" in result.bins().columns
    assert_array_equal(result.bins()["weight"].to_numpy(), weight[perm])


def test_regroup_preserves_the_declared_bin_size(tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    result = MetaCooler(out)
    assert result.binsize == MetaCooler(tiny_metacooler).binsize
    with h5py.File(out, "r") as f:
        assert f.attrs["bin-type"] == "fixed"


# ---------------------------------------------------------------------------
# layout
# ---------------------------------------------------------------------------


def test_layout_keep_downgrades_the_ordering_attribute(
    tmp_path: Path, tiny_metacooler: Path
):
    """Refusing to move contigs means the file can no longer claim to be grouped."""
    out = tmp_path / "out.metacool"
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
        out,
        layout="keep",
    )

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == SOURCE_ORDER
    assert result.info["ordering"] == CUSTOM_ORDERING
    validate(out)  # the downgrade is what keeps this passing
    with pytest.raises(ValueError, match="not contiguous in bin space"):
        result.genome_extent("MAG_2")
    # Membership-based selection is unaffected by the scattering.
    assert len(result.genome_bin_ids("MAG_2")) == 14


def test_layout_keep_leaves_the_label_alone_when_grouping_survives(
    tmp_path: Path, tiny_metacooler: Path
):
    """A reassignment that happens to preserve the layout keeps its ordering name."""
    out = tmp_path / "out.metacool"
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
        out,
        layout="keep",
    )
    assert MetaCooler(out).info["ordering"] == "genome-then-alpha"


def test_layout_can_force_a_different_strategy(tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
        out,
        layout="genome-then-length",
    )

    result = MetaCooler(out)
    assert result.info["ordering"] == "genome-then-length"
    # MAG_1 = {A (1000), B (600), C (400)} laid out longest first.
    assert result.contigs["name"].tolist() == [
        "contig_D",
        "contig_A",
        "contig_B",
        "contig_C",
    ]


@pytest.fixture
def alpha_ordered_metacooler(
    tmp_path: Path,
    tiny_fasta: Path,
    tiny_contigs_table: Path,
    tiny_bins_bed: Path,
    tiny_pairs: Path,
) -> Path:
    """A metacool whose contigs were never grouped by genome."""
    from metacooler.io.create import from_pairs

    out = tmp_path / "alpha.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=tiny_contigs_table,
        contig_table_format="contigs",
        fasta=tiny_fasta,
        contig_ordering="alpha",
        chunksize=1000,
    )
    return out


def test_auto_does_not_regroup_a_file_that_never_grouped(
    tmp_path: Path, alpha_ordered_metacooler: Path
):
    """There is no invariant to restore, so the rows (and pixels) stay put."""
    source = MetaCooler(alpha_ordered_metacooler)
    out = tmp_path / "out.metacool"

    reassign(
        alpha_ordered_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
        out,
    )

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == source.contigs["name"].tolist()
    # 'alpha' says nothing about genomes, so the reassignment cannot falsify it.
    assert result.info["ordering"] == "alpha"
    assert_array_equal(_dense(result), _dense(source))


def test_a_non_grouped_file_can_be_grouped_on_request(
    tmp_path: Path, alpha_ordered_metacooler: Path
):
    out = tmp_path / "out.metacool"
    reassign(
        alpha_ordered_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
        out,
        layout="genome-then-alpha",
    )

    result = MetaCooler(out)
    assert result.info["ordering"] == "genome-then-alpha"
    assert result.genome_extent("MAG_2")[1] > 0
    assert_array_equal(
        np.sort(result.contigs["genome_id"].to_numpy()),
        result.contigs["genome_id"].to_numpy(),
    )


@pytest.fixture
def collected_metacooler(tmp_path: Path) -> Path:
    """A ``collect``ed map: ``cool`` ordering, genomes contiguous by construction.

    Laid out [B_subtilis__chromosome, B_subtilis__plasmid1, E_coli__chromosome,
    E_coli__plasmid1, E_coli__plasmid2].
    """
    from metacooler.io.collect import collect

    ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
    bsub = make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)
    out = tmp_path / "collected.metacool"
    collect([ecoli, bsub], out)
    return out


def test_a_collected_map_is_relabelled_without_a_warning_when_genomes_stay_whole(
    tmp_path: Path, collected_metacooler: Path, caplog
):
    """E_coli's first contig borders B_subtilis's block, so both stay contiguous."""
    source = MetaCooler(collected_metacooler)
    out = tmp_path / "out.metacool"
    table = _table(tmp_path / "t.tsv", {"E_coli__chromosome": "B_subtilis"})
    with caplog.at_level("WARNING", logger="metacooler"):
        reassign(collected_metacooler, table, out, qc=False)

    assert "scattered" not in caplog.text
    result = MetaCooler(out)
    assert result.info["ordering"] == "cool"
    # 8 + 3 bins of its own, then E_coli__chromosome's 10.
    assert result.genome_extent("B_subtilis") == (0, 21)
    assert_array_equal(_dense(result), _dense(source))


def test_scattering_a_collected_map_is_warned_about(
    tmp_path: Path, collected_metacooler: Path, caplog
):
    """``cool`` does not promise grouping, so auto relabels -- but says so."""
    out = tmp_path / "out.metacool"
    table = _table(tmp_path / "t.tsv", {"E_coli__plasmid1": "B_subtilis"})
    with caplog.at_level("WARNING", logger="metacooler"):
        reassign(collected_metacooler, table, out, qc=False)

    assert "2 genome(s) are now scattered: B_subtilis, E_coli" in caplog.text
    result = MetaCooler(out)
    assert result.info["ordering"] == "cool"
    validate(out)
    with pytest.raises(ValueError, match="not contiguous in bin space"):
        result.genome_extent("E_coli")


def test_a_collected_map_can_be_regrouped_on_request(
    tmp_path: Path, collected_metacooler: Path
):
    source = MetaCooler(collected_metacooler)
    out = tmp_path / "out.metacool"
    reassign(
        collected_metacooler,
        _table(tmp_path / "t.tsv", {"E_coli__plasmid1": "B_subtilis"}),
        out,
        layout="genome-then-alpha",
        qc=False,
    )

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == [
        "B_subtilis__chromosome",
        "B_subtilis__plasmid1",
        "E_coli__plasmid1",
        "E_coli__chromosome",
        "E_coli__plasmid2",
    ]
    assert result.genome_extent("B_subtilis")[0] == 0
    assert result.genome_extent("E_coli")[1] == len(result.bins())
    assert result.info["sum"] == source.info["sum"]
    perm = _bin_permutation(source, result)
    assert_array_equal(_dense(result), _dense(source)[np.ix_(perm, perm)])


def test_unknown_layout_is_rejected(tmp_path: Path, tiny_metacooler: Path):
    with pytest.raises(ValueError, match="layout must be one of"):
        reassign(
            tiny_metacooler,
            _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
            tmp_path / "out.metacool",
            layout="whatever",  # type: ignore
        )


# ---------------------------------------------------------------------------
# --missing
# ---------------------------------------------------------------------------


def test_missing_keep_is_incremental(tmp_path: Path, tiny_metacooler: Path):
    """A table naming one contig changes that contig and nothing else."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_9"}), out)

    contigs = MetaCooler(out).contigs
    names = MetaCooler(out).genomes["name"].to_numpy()
    assignment = dict(zip(contigs["name"], names[contigs["genome_id"]], strict=True))
    assert assignment == {
        "contig_A": "MAG_1",
        "contig_B": "MAG_1",
        "contig_C": "MAG_9",
        "contig_D": "__unbinned__",
    }


def test_missing_unbinned_makes_the_table_authoritative(
    tmp_path: Path, tiny_metacooler: Path
):
    out = tmp_path / "out.metacool"
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_X"}),
        out,
        missing="unbinned",
    )

    result = MetaCooler(out)
    assert result.genomes["name"].tolist() == ["__unbinned__", "MAG_X"]
    contigs = result.contigs
    names = result.genomes["name"].to_numpy()
    assignment = dict(zip(contigs["name"], names[contigs["genome_id"]], strict=True))
    assert assignment == {
        "contig_A": "MAG_X",
        "contig_B": "__unbinned__",
        "contig_C": "__unbinned__",
        "contig_D": "__unbinned__",
    }


def test_sentinel_genome_names_unbin_a_contig(tmp_path: Path, tiny_metacooler: Path):
    """``NA`` and friends mean __unbinned__, as they do at creation time."""
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "NA"}), out)

    result = MetaCooler(out)
    contigs = result.contigs.set_index("name")
    assert int(contigs.at["contig_A", "genome_id"]) == 0
    assert result.genomes["name"].tolist()[0] == "__unbinned__"


# ---------------------------------------------------------------------------
# in place
# ---------------------------------------------------------------------------


def test_inplace_rewrites_the_source(tmp_path: Path, tiny_metacooler: Path):
    target = tmp_path / "inplace.metacool"
    shutil.copy2(tiny_metacooler, target)

    reassign(target, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), inplace=True)

    result = MetaCooler(target)
    assert result.genomes["name"].tolist() == ["__unbinned__", "MAG_1"]
    assert result.contigs["name"].tolist() == SOURCE_ORDER
    assert_array_equal(_dense(result), _dense(MetaCooler(tiny_metacooler)))


def test_inplace_re_orders_through_a_staged_map(tmp_path: Path, tiny_metacooler: Path):
    """A reassignment that needs re-ordering rewrites the whole map: in place,
    a new map is written beside it, and replaces it once it validates (it was
    refused in place)."""
    target = tmp_path / "inplace.metacool"
    shutil.copy2(tiny_metacooler, target)
    copy = tmp_path / "copy.metacool"
    table = _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"})

    reassign(target, table, copy)
    reassign(target, table, inplace=True)

    in_place, written = MetaCooler(target), MetaCooler(copy)
    pd.testing.assert_frame_equal(in_place.genomes, written.genomes)
    pd.testing.assert_frame_equal(in_place.contigs, written.contigs)
    pd.testing.assert_frame_equal(in_place.pixels()[:], written.pixels()[:])
    validate(target, full=True)


def test_inplace_refuses_to_re_order_an_entry(tmp_path: Path, tiny_mmetacool: Path):
    """An entry of a collection is not a file of its own: re-ordering it in
    place is refused, before anything is written."""
    before = tiny_mmetacool.read_bytes()
    with pytest.raises(ValueError, match="cannot re-group an entry of a collection"):
        reassign(
            f"{tiny_mmetacool}::entries/entry1",
            _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
            inplace=True,
        )
    assert tiny_mmetacool.read_bytes() == before


def test_inplace_and_output_are_mutually_exclusive(
    tmp_path: Path, tiny_metacooler: Path
):
    with pytest.raises(ValueError, match="no output path is expected"):
        reassign(
            tiny_metacooler,
            _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
            tmp_path / "out.metacool",
            inplace=True,
        )


def test_an_output_is_required(tmp_path: Path, tiny_metacooler: Path):
    with pytest.raises(ValueError, match="An output path is required"):
        reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}))


# ---------------------------------------------------------------------------
# Table handling
# ---------------------------------------------------------------------------


def test_names_absent_from_the_file_are_ignored_with_a_warning(
    tmp_path: Path, tiny_metacooler: Path, caplog
):
    out = tmp_path / "out.metacool"
    table = _table(tmp_path / "t.tsv", {"contig_C": "MAG_1", "contig_ZZZ": "MAG_1"})

    with caplog.at_level("WARNING", logger="metacooler"):
        reassign(tiny_metacooler, table, out)

    assert "contig_ZZZ" in caplog.text
    assert MetaCooler(out).contigs["name"].tolist() == SOURCE_ORDER


def test_a_table_matching_no_contig_is_refused(tmp_path: Path, tiny_metacooler: Path):
    table = _table(tmp_path / "t.tsv", {"scaffold_1": "MAG_1", "scaffold_2": "MAG_2"})
    with pytest.raises(ValueError, match="match a contig of this file"):
        reassign(tiny_metacooler, table, tmp_path / "out.metacool")


def test_an_empty_table_is_refused(tmp_path: Path, tiny_metacooler: Path):
    table = tmp_path / "t.tsv"
    table.write_text("contig_name\tgenome_name\n")
    with pytest.raises(ValueError, match="empty"):
        reassign(tiny_metacooler, table, tmp_path / "out.metacool")


@pytest.mark.parametrize("column", ["Final_cluster", "Final_bin"])
def test_metator_format_is_accepted(tmp_path: Path, tiny_metacooler: Path, column):
    """metaTOR 2.0's contig_data_final.txt (Final_cluster) and 1.x's
    (Final_bin)."""
    table = tmp_path / "contig_data_final.txt"
    table.write_text(f"Name\tSize\t{column}\ncontig_C\t400\tMetaTOR_1_0\n")
    out = tmp_path / "out.metacool"

    reassign(tiny_metacooler, table, out, contig_table_format="metator")

    assert "MetaTOR_1_0" in MetaCooler(out).genomes["name"].tolist()


def test_a_dict_mapping_is_accepted(tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, {"contig_C": "MAG_1"}, out)
    assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__", "MAG_1"]


def test_non_ascii_genome_names_are_refused(tmp_path: Path, tiny_metacooler: Path):
    table = _table(tmp_path / "t.tsv", {"contig_C": "MAG_é"})
    with pytest.raises(Exception, match=r"not ASCII|cooler could not open|ASCII"):
        reassign(tiny_metacooler, table, tmp_path / "out.metacool")
    assert not (tmp_path / "out.metacool").exists()


# ---------------------------------------------------------------------------
# Metadata carried across
# ---------------------------------------------------------------------------


@pytest.fixture
def annotated_metacooler(
    tmp_path: Path,
    tiny_fasta: Path,
    tiny_contigs_table: Path,
    tiny_bins_bed: Path,
    tiny_pairs: Path,
    tiny_genomes_table: Path,
) -> Path:
    """A metacool carrying user metadata on both /genomes and /contigs."""
    from metacooler.io.create import from_pairs

    contigs = tmp_path / "contigs_extra.tsv"
    contigs.write_text(
        "contig_name\tgenome_name\tcoverage\n"
        "contig_A\tMAG_1\t12\n"
        "contig_B\tMAG_1\t8\n"
        "contig_C\tMAG_2\t31\n"
    )
    out = tmp_path / "annotated.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=contigs,
        contig_table_format="contigs",
        genomes_table=tiny_genomes_table,
        fasta=tiny_fasta,
        contig_ordering="genome-then-alpha",
        chunksize=1000,
    )
    return out


def test_contig_metadata_survives_a_reassignment(
    tmp_path: Path, annotated_metacooler: Path
):
    """Contigs are untouched by a reassignment, so their columns are still true."""
    out = tmp_path / "out.metacool"
    reassign(
        annotated_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out
    )

    contigs = MetaCooler(out).contigs.set_index("name")
    assert int(contigs.at["contig_A", "coverage"]) == 12
    assert int(contigs.at["contig_C", "coverage"]) == 31


def test_genome_metadata_is_carried_over_for_surviving_genomes(
    tmp_path: Path, annotated_metacooler: Path
):
    out = tmp_path / "out.metacool"
    reassign(
        annotated_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out
    )

    genomes = MetaCooler(out).genomes.set_index("name")
    assert genomes.at["MAG_1", "phylum"] == "Bacteroidetes"
    assert genomes.at["MAG_2", "quality"] == "medium_quality"


def test_a_genomes_table_replaces_the_carried_metadata(
    tmp_path: Path, annotated_metacooler: Path
):
    fresh = tmp_path / "fresh_genomes.tsv"
    fresh.write_text("genome_name\tphylum\nMAG_1\tProteobacteria\n")
    out = tmp_path / "out.metacool"

    reassign(
        annotated_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
        out,
        genomes_table=fresh,
    )

    genomes = MetaCooler(out).genomes
    assert "quality" not in genomes.columns
    assert genomes.set_index("name").at["MAG_1", "phylum"] == "Proteobacteria"


def test_provenance_is_appended_not_replaced(tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    with h5py.File(tiny_metacooler, "r+") as f:
        f.attrs["metadata"] = json.dumps({"commands": ["metacooler create pairs …"]})

    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
        out,
        metadata={"commands": ["metacooler reassign …"]},
    )

    metadata = MetaCooler(out).metadata
    assert metadata["commands"] == [
        "metacooler create pairs …",
        "metacooler reassign …",
    ]
    assert metadata["reassign"][0]["n_contigs_reassigned"] == 1
    assert metadata["reassign"][0]["layout"] == "genome-then-alpha"


# ---------------------------------------------------------------------------
# A table that moves no contig
# ---------------------------------------------------------------------------

# contig_A is in MAG_1 already.
NO_MOVE = {"contig_A": "MAG_1"}
ONLY_PROVENANCE = (
    "No contig changes genome; the output differs from the input only in its "
    "provenance metadata."
)


def _nothing_moved(caplog: pytest.LogCaptureFixture) -> logging.LogRecord:
    (record,) = [
        r for r in caplog.records if r.getMessage().startswith("No contig changes")
    ]
    return record


def _drop_qc_block(path: Path) -> Path:
    with h5py.File(path, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        del metadata["qc"]
        f.attrs["metadata"] = json.dumps(metadata)
    return path


def test_a_table_moving_nothing_warns_that_only_provenance_differs(
    tmp_path: Path, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
):
    """Guard: the warning keeps its meaning when it is true."""
    source = _drop_qc_block(tiny_metacooler)
    reassign(source, _table(tmp_path / "t.tsv", NO_MOVE), tmp_path / "out.metacool")

    record = _nothing_moved(caplog)
    assert record.levelno == logging.WARNING
    assert record.getMessage() == ONLY_PROVENANCE


def test_a_table_moving_nothing_still_names_a_recomputed_qc_block(
    tmp_path: Path, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
):
    """Nothing was asked for but the table, so it warns -- truthfully.

    A map without per-contig QC counts (written by 0.2.0) has its block
    recomputed from the pixels, which the message names.
    """
    with h5py.File(tiny_metacooler, "r+") as f:
        for column in QC_CONTIG_COLUMNS:
            del f["contigs"][column]
    reassign(
        tiny_metacooler, _table(tmp_path / "t.tsv", NO_MOVE), tmp_path / "out.metacool"
    )

    record = _nothing_moved(caplog)
    assert record.levelno == logging.WARNING
    assert record.getMessage() == (
        "No contig changes genome; the output differs from the input only in its "
        "provenance metadata and its QC block (recomputed from the pixels)."
    )


def test_a_table_moving_nothing_does_not_name_a_carried_qc_block(
    tmp_path: Path, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
):
    """Carried with no contig moved, the block is the source's: not a change.

    A guard against this commit's parent, whose one message happened to be
    true here; it pins that a carried block is never named as recomputed.
    """
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", NO_MOVE), out)

    record = _nothing_moved(caplog)
    assert record.levelno == logging.WARNING
    assert record.getMessage() == ONLY_PROVENANCE
    assert MetaCooler(out).metadata["qc"] == MetaCooler(tiny_metacooler).metadata["qc"]


@pytest.mark.parametrize(
    ("change", "says"),
    [
        ("genomes", "its /genomes metadata columns (from the genomes table)"),
        ("layout", "its contig layout (re-grouped as 'genome-then-alpha')"),
        ("label", "its contig-ordering attribute (now 'genome-then-length')"),
        ("no-qc", "its QC block (dropped)"),
        ("metadata", "its metadata key(s) (sample)"),
    ],
)
def test_a_table_moving_nothing_reports_what_else_changes(
    tmp_path: Path,
    tiny_metacooler: Path,
    alpha_ordered_metacooler: Path,
    caplog: pytest.LogCaptureFixture,
    change: str,
    says: str,
):
    """Asked for on purpose, the other changes make it an INFO line naming them."""
    caplog.set_level(logging.INFO, logger="metacooler")
    source = alpha_ordered_metacooler if change == "layout" else tiny_metacooler
    genomes = tmp_path / "g.tsv"
    genomes.write_text("genome_name\tcompleteness\nMAG_1\t90\n")
    layouts: dict[str, Layout] = {
        "layout": "genome-then-alpha",
        "label": "genome-then-length",
    }
    reassign(
        source,
        _table(tmp_path / "t.tsv", NO_MOVE),
        tmp_path / "out.metacool",
        genomes_table=genomes if change == "genomes" else None,
        layout=layouts.get(change, "auto"),
        qc=change != "no-qc",
        metadata={"commands": ["metacooler reassign"], "sample": "S1"}
        if change == "metadata"
        else None,
    )

    record = _nothing_moved(caplog)
    assert record.levelno == logging.INFO
    message = record.getMessage()
    assert message.startswith("No contig changes genome; the output differs")
    assert says in message
    assert message != ONLY_PROVENANCE


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------


def test_qc_is_recomputed_for_the_new_assignment(tmp_path: Path, tiny_metacooler: Path):
    """The stale per-genome block is replaced, not carried and not left behind."""
    source = MetaCooler(tiny_metacooler)
    assert "qc" in source.metadata  # create writes one by default
    out = tmp_path / "out.metacool"

    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), out)

    result = MetaCooler(out)
    assert len(result.genomes) == 2
    assert "qc_inter_genome" in result.genomes.columns
    assert result.metadata["qc"]["n_contacts"] == source.metadata["qc"]["n_contacts"]
    # Merging two genomes can only move contacts out of the inter-genome bucket.
    assert result.metadata["qc"]["inter_genome"] < source.metadata["qc"]["inter_genome"]


def test_no_qc_drops_the_stale_block(tmp_path: Path, tiny_metacooler: Path):
    """Relabelled, the source's /contigs qc_* columns used to stay behind."""
    assert "qc_intra_contig_short" in MetaCooler(tiny_metacooler).contigs.columns
    out = tmp_path / "out.metacool"
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
        out,
        qc=False,
    )

    result = MetaCooler(out)
    assert "qc" not in result.metadata
    assert not [c for c in result.genomes.columns if str(c).startswith("qc_")]
    assert not [c for c in result.contigs.columns if str(c).startswith("qc_")]
    validate(out)


def test_qc_columns_are_not_carried_as_contig_metadata(
    tmp_path: Path, tiny_metacooler: Path, caplog: pytest.LogCaptureFixture
):
    """Carried as user metadata, the reserved names came back as user_qc_*."""
    out = tmp_path / "out.metacool"
    # contig_A joins MAG_2, which re-groups the contigs.
    reassign(tiny_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    result = MetaCooler(out)
    assert result.contigs["name"].tolist() != SOURCE_ORDER
    assert not [c for c in result.contigs.columns if str(c).startswith("user_")]
    assert "reserved column names" not in caplog.text


@pytest.fixture
def map_with_a_user_qc_column(
    tmp_path: Path, tiny_fasta: Path, tiny_bins_bed: Path, tiny_pairs: Path
) -> Path:
    """The tiny map, with a user /contigs column that merely starts with qc_."""
    from metacooler.io.create import from_pairs

    table = tmp_path / "with_qc_pass.tsv"
    table.write_text(
        "contig_name\tgenome_name\tqc_pass\n"
        "contig_A\tMAG_1\t1\ncontig_B\tMAG_1\t0\ncontig_C\tMAG_2\t1\n"
    )
    out = tmp_path / "with_qc_pass.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=table,
        fasta=tiny_fasta,
        contig_ordering="genome-then-alpha",
    )
    return out


@pytest.mark.parametrize("qc", [True, False])
@pytest.mark.parametrize("layout", ["keep", "genome-then-length", "in place"])
def test_a_user_qc_column_is_kept_as_contig_metadata(
    tmp_path: Path, map_with_a_user_qc_column: Path, layout: str, qc: bool
):
    """Only the two per-contig QC columns go with the block; qc_pass is metadata."""
    source_order = MetaCooler(map_with_a_user_qc_column).contigs["name"].tolist()
    out = tmp_path / "out.metacool"
    # contig_A joins MAG_2, which a genome-grouped layout has to re-group.
    table = _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"})
    if layout == "in place":
        reassign(map_with_a_user_qc_column, table, inplace=True, layout="keep", qc=qc)
        out = map_with_a_user_qc_column
    else:
        reassign(
            map_with_a_user_qc_column,
            table,
            out,
            layout="keep" if layout == "keep" else "genome-then-length",
            qc=qc,
        )

    order = MetaCooler(out).contigs["name"].tolist()
    assert (order != source_order) is (layout == "genome-then-length")
    contigs = MetaCooler(out).contigs.set_index("name")
    assert contigs["qc_pass"].dropna().to_dict() == {
        "contig_A": 1,
        "contig_B": 0,
        "contig_C": 1,
    }
    validate(out)


def test_no_qc_in_place_drops_the_stale_block(tmp_path: Path, tiny_metacooler: Path):
    reassign(
        tiny_metacooler,
        _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
        inplace=True,
        qc=False,
    )

    result = MetaCooler(tiny_metacooler)
    assert "qc" not in result.metadata
    assert not [c for c in result.contigs.columns if str(c).startswith("qc_")]
    validate(tiny_metacooler)


# ---------------------------------------------------------------------------
# mmetacool entries
# ---------------------------------------------------------------------------


def test_an_entry_reassigns_into_a_standalone_file(
    tmp_path: Path, tiny_mmetacool: Path
):
    from metacooler.schema import list_metacoolers

    out = tmp_path / "out.metacool"

    uri = f"{tiny_mmetacool}::{list_metacoolers(tiny_mmetacool)[0]}"
    reassign(uri, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), out)

    result = MetaCooler(out)
    assert result.genomes["name"].tolist() == ["__unbinned__", "MAG_1"]
    with h5py.File(out, "r") as f:
        assert isinstance(f.get("chroms", getlink=True), h5py.SoftLink)
    # The whole point of the soft link: cooler still opens the result.
    assert cooler.Cooler(str(out)).info["nchroms"] == 4


def test_an_entry_can_be_reassigned_in_place(tmp_path: Path, tiny_mmetacool: Path):
    from metacooler.schema import list_metacoolers

    uri = f"{tiny_mmetacool}::{list_metacoolers(tiny_mmetacool)[0]}"
    reassign(uri, _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}), inplace=True)

    assert MetaCooler(uri).genomes["name"].tolist() == ["__unbinned__", "MAG_1"]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_reassign(runner, tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__", "MAG_1"]


def test_cli_records_the_invocation(runner, tmp_path: Path, tiny_metacooler: Path):
    out = tmp_path / "out.metacool"
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "reassign" in MetaCooler(out).metadata["commands"][-1]


def test_cli_requires_an_output_or_inplace(
    runner, tmp_path: Path, tiny_metacooler: Path
):
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
        ],
    )
    assert result.exit_code != 0
    assert "--inplace" in flat(result.output)


def test_cli_rejects_output_with_inplace(runner, tmp_path: Path, tiny_metacooler: Path):
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
            "-o",
            str(tmp_path / "out.metacool"),
            "--inplace",
        ],
    )
    assert result.exit_code != 0
    assert "one or the other" in flat(result.output)


def test_cli_refuses_to_clobber_without_force(
    runner, tmp_path: Path, tiny_metacooler: Path
):
    out = tmp_path / "out.metacool"
    out.write_text("not a metacool")
    args = [
        "reassign",
        str(tiny_metacooler),
        "--contigs",
        str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
        "-o",
        str(out),
    ]

    refused = runner.invoke(cli, args)
    assert refused.exit_code != 0
    assert "--force" in flat(refused.output)

    forced = runner.invoke(cli, [*args, "--force"])
    assert forced.exit_code == 0, forced.output
    assert MetaCooler(out).genomes["name"].tolist() == ["__unbinned__", "MAG_1"]


def _fingerprint(path: Path) -> bytes:
    return path.read_bytes()


def test_an_output_that_is_the_input_is_refused_even_under_force(
    tmp_path: Path, tiny_metacooler: Path
):
    """``--force`` used to delete the source before a byte of it was read."""
    before = _fingerprint(tiny_metacooler)

    with pytest.raises(ValueError, match="same file as the input"):
        reassign(
            tiny_metacooler,
            _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"}),
            tiny_metacooler,
            force=True,
        )

    assert _fingerprint(tiny_metacooler) == before


def test_cli_refuses_an_output_that_is_the_input(
    runner, tmp_path: Path, tiny_metacooler: Path
):
    before = _fingerprint(tiny_metacooler)
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
            "-o",
            str(tiny_metacooler),
            "--force",
        ],
    )

    assert result.exit_code != 0
    assert "same file as the input" in flat(result.output)
    assert _fingerprint(tiny_metacooler) == before


def test_cli_refuses_an_output_that_is_the_contigs_table(
    runner, tmp_path: Path, tiny_metacooler: Path
):
    table = _table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})
    before = table.read_text()
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(table),
            "-o",
            str(table),
            "--force",
        ],
    )

    assert result.exit_code != 0
    assert "same file as the input" in flat(result.output)
    assert table.read_text() == before


def test_force_keeps_the_previous_output_when_the_run_fails(
    runner, tmp_path: Path, tiny_metacooler: Path
):
    """The old output is replaced on success, not deleted up front.

    A table that matches no contig fails *after* the output check; the
    finished map already at ``-o`` must survive that.
    """
    out = tmp_path / "out.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "ok.tsv", {"contig_C": "MAG_1"}), out)
    before = _fingerprint(out)

    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_metacooler),
            "--contigs",
            str(_table(tmp_path / "bad.tsv", {"no_such_contig": "MAG_1"})),
            "-o",
            str(out),
            "--force",
        ],
    )

    assert result.exit_code != 0
    assert "match a contig of this file" in flat(result.output)
    assert _fingerprint(out) == before
    validate(out)


@pytest.fixture
def empty_metacool(tmp_path: Path, tiny_contigs_table: Path, tiny_bins_bed: Path):
    """A valid metacool with zero pixels, laid out ``genome-then-alpha``."""
    from metacooler.io.create import from_pairs

    pairs = tmp_path / "empty.pairs"
    pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    out = tmp_path / "empty.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=pairs,
        contig_genome_map=tiny_contigs_table,
        contig_ordering="genome-then-alpha",
        ps=False,
    )
    return out


def test_a_zero_pixel_map_can_be_regrouped(tmp_path: Path, empty_metacool: Path):
    """The regroup path streams no pixel at all; the writer must still cope."""
    out = tmp_path / "out.metacool"
    reassign(empty_metacool, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out)

    validate(out)
    result = MetaCooler(out)
    assert result.contigs["name"].tolist() == [
        "contig_D",
        "contig_B",
        "contig_A",
        "contig_C",
    ]
    assert int(result.info["nnz"]) == 0


@pytest.fixture
def square_metacooler(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """The tiny map declared ``storage-mode="square"`` (which validate accepts)."""
    target = tmp_path / "square.metacool"
    shutil.copy2(tiny_metacooler, target)
    with h5py.File(target, "r+") as f:
        f.attrs["storage-mode"] = "square"
    validate(target)
    return target


def test_regrouping_a_square_map_is_refused(tmp_path: Path, square_metacooler: Path):
    """Remapping folds every pixel into the upper triangle.

    That is only right for a symmetric matrix, which square storage does not
    promise; a pixel and its mirror would be summed.
    """
    out = tmp_path / "out.metacool"
    with pytest.raises(ValueError, match="storage-mode 'square'"):
        reassign(
            square_metacooler, _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}), out
        )
    assert not out.exists()


def test_a_square_map_can_still_be_relabelled(tmp_path: Path, square_metacooler: Path):
    """Relabelling never touches a pixel, so the storage mode does not matter."""
    out = tmp_path / "out.metacool"
    reassign(
        square_metacooler,
        _table(tmp_path / "t.tsv", {"contig_A": "MAG_2"}),
        out,
        layout="keep",
    )

    validate(out)
    assert cooler.Cooler(str(out)).info["storage-mode"] == "square"
    assert_array_equal(_dense(MetaCooler(out)), _dense(MetaCooler(square_metacooler)))


def test_cli_rejects_a_bare_mmetacool_path(
    runner, tmp_path: Path, tiny_mmetacool: Path
):
    result = runner.invoke(
        cli,
        [
            "reassign",
            str(tiny_mmetacool),
            "--contigs",
            str(_table(tmp_path / "t.tsv", {"contig_C": "MAG_1"})),
            "-o",
            str(tmp_path / "out.metacool"),
        ],
    )
    assert result.exit_code != 0
    assert "entries" in flat(result.output)


# ---------------------------------------------------------------------------
# Round trip
# ---------------------------------------------------------------------------


def test_reassigning_back_restores_the_original_file(
    tmp_path: Path, tiny_metacooler: Path
):
    """Two reassignments that cancel out give back the file we started from."""
    source = MetaCooler(tiny_metacooler)
    first = tmp_path / "first.metacool"
    reassign(tiny_metacooler, _table(tmp_path / "a.tsv", {"contig_A": "MAG_2"}), first)

    second = tmp_path / "second.metacool"
    reassign(first, _table(tmp_path / "b.tsv", {"contig_A": "MAG_1"}), second)

    result = MetaCooler(second)
    assert result.contigs["name"].tolist() == source.contigs["name"].tolist()
    assert result.genomes["name"].tolist() == source.genomes["name"].tolist()
    assert_array_equal(
        result.contigs["genome_id"].to_numpy(),
        source.contigs["genome_id"].to_numpy(),
    )
    assert_array_equal(_dense(result), _dense(source))
