"""Tests for metacooler.io.annotate and the ``metacooler annotate`` command."""

from __future__ import annotations

import json
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.annotate import _SCRATCH, annotate
from metacooler.io.create import from_cool
from metacooler.schema import MetaCoolerFormatError, validate


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def genomes_of(uri: str | Path) -> pd.DataFrame:
    return MetaCooler(uri).genomes.set_index("name")


def contigs_of(uri: str | Path) -> pd.DataFrame:
    return MetaCooler(uri).contigs.set_index("name")


def history(uri: str | Path) -> list[dict]:
    return MetaCooler(uri).metadata.get("annotate", [])


@pytest.fixture
def checkm(tmp_path: Path) -> Path:
    """A CheckM2-like table: genome names under 'Name', one unknown genome."""
    return write(
        tmp_path / "checkm2.tsv",
        "Name\tCompleteness\tContamination\tlineage\n"
        "MAG_1\t91.5\t1.2\tBacteroidota\n"
        "MAG_2\t70\t0.5\tBacillota\n"
        "MAG_99\t50\t2\tunknown\n",
    )


@pytest.fixture
def depth(tmp_path: Path) -> Path:
    """A MetaBAT-like depth table: contig names under 'contigName'."""
    return write(
        tmp_path / "depth.tsv",
        "contigName\ttotalAvgDepth\ncontig_A\t12.5\ncontig_C\t3\n",
    )


@pytest.fixture
def mc(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """A private copy of the tiny map, safe to edit in place."""
    path = tmp_path / "edit.metacool"
    path.write_bytes(tiny_metacooler.read_bytes())
    return path


class TestAdd:
    def test_genome_columns(self, mc, checkm):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        g = genomes_of(mc)
        assert g.loc["MAG_1", "Completeness"] == 91.5
        assert g.loc["MAG_2", "Completeness"] == 70.0
        assert np.isnan(g.loc["__unbinned__", "Completeness"])
        assert g.loc["MAG_2", "lineage"] == "Bacillota"
        assert g.loc["__unbinned__", "lineage"] == ""
        assert "MAG_99" not in g.index
        validate(mc)

    def test_contig_columns(self, mc, depth):
        annotate(
            mc,
            inplace=True,
            contigs_table=depth,
            contigs_table_format="contigName",
        )
        c = contigs_of(mc)
        assert c.loc["contig_A", "totalAvgDepth"] == 12.5
        assert c.loc["contig_C", "totalAvgDepth"] == 3.0
        assert np.isnan(c.loc["contig_B", "totalAvgDepth"])
        # cooler reads /contigs as /chroms, new column included.
        chroms = cooler.Cooler(str(mc)).chroms()[:]
        assert "totalAvgDepth" in chroms.columns
        validate(mc)

    def test_new_columns_can_be_filtered_on(self, mc, checkm):
        from metacooler.io.dump import dump_table

        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        out = mc.parent / "g.tsv"
        dump_table(mc, "genomes", out=out, filters=["genomes.Completeness >= 90"])
        assert out.read_text().split()[0] == "MAG_1"

    def test_derived_names_are_prefixed(self, mc, tmp_path):
        table = write(
            tmp_path / "g.tsv",
            "genome_name\tlength\tn_bins\nMAG_1\t5\t6\nMAG_2\t7\t8\n",
        )
        annotate(mc, inplace=True, genomes_table=table)
        g = genomes_of(mc)
        assert g.loc["MAG_1", "user_length"] == 5
        assert g.loc["MAG_1", "user_n_bins"] == 6
        assert g.loc["MAG_1", "length"] != 5

    def test_columns_selects(self, mc, checkm):
        annotate(
            mc,
            inplace=True,
            genomes_table=checkm,
            genomes_table_format="Name",
            columns=["genomes.Completeness"],
        )
        g = genomes_of(mc)
        assert "Completeness" in g.columns
        assert "Contamination" not in g.columns

    def test_unknown_selected_column(self, mc, checkm):
        with pytest.raises(ValueError, match="no column 'Nope'"):
            annotate(
                mc,
                inplace=True,
                genomes_table=checkm,
                genomes_table_format="Name",
                columns=["genomes.Nope"],
            )

    def test_columns_for_a_table_not_given(self, mc, checkm):
        with pytest.raises(ValueError, match="no contigs table was given"):
            annotate(
                mc,
                inplace=True,
                genomes_table=checkm,
                genomes_table_format="Name",
                columns=["contigs.depth"],
            )

    def test_column_order_is_kept(self, mc, checkm, depth):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        before = list(MetaCooler(mc).genomes.columns)
        patch = write(mc.parent / "p.tsv", "Name\tCompleteness\nMAG_1\t10\n")
        annotate(
            mc,
            inplace=True,
            genomes_table=patch,
            genomes_table_format="Name",
            replace=True,
        )
        assert list(MetaCooler(mc).genomes.columns) == before


class TestSameMeaningAsCreate:
    def test_annotating_gives_what_create_would_have(
        self, tmp_path, tiny_cool, tiny_contigs_table
    ):
        """annotate(T) after create == create(T): same columns, values, types."""
        genomes_tsv = write(
            tmp_path / "g.tsv",
            "genome_name\tquality\tn_markers\tcompleteness\n"
            "MAG_1\thigh\t12\t0.9\nMAG_2\tlow\t7\t\n",
        )
        contigs_tsv = write(
            tmp_path / "c.tsv",
            "contig_name\tgenome_name\tgc\n"
            "contig_A\tMAG_1\t0.5\ncontig_B\tMAG_1\t0.4\ncontig_C\tMAG_2\t0.6\n",
        )
        created = tmp_path / "created.metacool"
        from_cool(
            tiny_cool,
            created,
            contig_genome_map=contigs_tsv,
            genomes_table=genomes_tsv,
            qc=False,
        )
        annotated = tmp_path / "annotated.metacool"
        from_cool(tiny_cool, annotated, contig_genome_map=tiny_contigs_table, qc=False)
        annotate(
            annotated,
            inplace=True,
            genomes_table=genomes_tsv,
            contigs_table=contigs_tsv,
        )
        with h5py.File(created, "r") as a, h5py.File(annotated, "r") as b:
            for table in ("genomes", "contigs"):
                assert list(a[table]) == list(b[table]), table
                for column in a[table]:
                    left, right = a[table][column][:], b[table][column][:]
                    assert left.dtype.kind == right.dtype.kind, (table, column)
                    if left.dtype.kind == "f":
                        np.testing.assert_array_equal(left, right)
                    else:
                        assert left.tolist() == right.tolist(), (table, column)


class TestExistingColumns:
    def test_identical_values_are_not_a_change(self, mc, checkm):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        # Nothing to write the second time: no second history entry.
        assert len(history(mc)) == 1

    def test_a_change_needs_replace(self, mc, checkm, tmp_path):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        patch = write(tmp_path / "p.tsv", "Name\tCompleteness\nMAG_1\t80\n")
        with pytest.raises(ValueError, match="genomes.Completeness: 1 row.*91.5 -> 80"):
            annotate(mc, inplace=True, genomes_table=patch, genomes_table_format="Name")

    def test_replace_keeps_unlisted_rows(self, mc, checkm, tmp_path):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        patch = write(tmp_path / "p.tsv", "Name\tCompleteness\nMAG_1\t80\n")
        annotate(
            mc,
            inplace=True,
            genomes_table=patch,
            genomes_table_format="Name",
            replace=True,
        )
        g = genomes_of(mc)
        assert g.loc["MAG_1", "Completeness"] == 80.0
        assert g.loc["MAG_2", "Completeness"] == 70.0
        assert history(mc)[-1]["replaced"] == ["genomes.Completeness"]

    def test_replace_with_missing_na_empties_unlisted_rows(self, mc, checkm, tmp_path):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        patch = write(tmp_path / "p.tsv", "Name\tCompleteness\tlineage\nMAG_1\t80\tX\n")
        annotate(
            mc,
            inplace=True,
            genomes_table=patch,
            genomes_table_format="Name",
            replace=True,
            missing="na",
        )
        g = genomes_of(mc)
        assert g.loc["MAG_1", "Completeness"] == 80.0
        assert np.isnan(g.loc["MAG_2", "Completeness"])
        assert g.loc["MAG_2", "lineage"] == ""

    def test_text_patch_of_a_numeric_column(self, mc, checkm, tmp_path):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        patch = write(tmp_path / "p.tsv", "Name\tCompleteness\nMAG_1\tN/A\n")
        annotate(
            mc,
            inplace=True,
            genomes_table=patch,
            genomes_table_format="Name",
            replace=True,
        )
        g = genomes_of(mc)
        assert g.loc["MAG_1", "Completeness"] == "N/A"
        assert g.loc["MAG_2", "Completeness"] == "70.0"


class TestDrop:
    def test_drop(self, mc, checkm, depth):
        annotate(
            mc,
            inplace=True,
            genomes_table=checkm,
            genomes_table_format="Name",
            contigs_table=depth,
            contigs_table_format="contigName",
        )
        annotate(mc, inplace=True, drop=["genomes.lineage", "contigs.totalAvgDepth"])
        assert "lineage" not in genomes_of(mc).columns
        assert "totalAvgDepth" not in contigs_of(mc).columns
        assert history(mc)[-1]["dropped"] == [
            "genomes.lineage",
            "contigs.totalAvgDepth",
        ]
        validate(mc)

    @pytest.mark.parametrize(
        "column", ["genomes.n_bins", "contigs.genome_id", "genomes.qc_inter_genome"]
    )
    def test_derived_columns_cannot_be_dropped(self, mc, column):
        with pytest.raises(ValueError, match="metacooler derives these columns"):
            annotate(mc, inplace=True, drop=[column])

    def test_absent_column(self, mc):
        with pytest.raises(ValueError, match="no such column"):
            annotate(mc, inplace=True, drop=["genomes.nope"])

    def test_dropped_and_supplied(self, mc, checkm):
        with pytest.raises(ValueError, match="both dropped and given"):
            annotate(
                mc,
                inplace=True,
                genomes_table=checkm,
                genomes_table_format="Name",
                drop=["genomes.lineage"],
            )

    def test_dropped_and_selected(self, mc, checkm):
        with pytest.raises(ValueError, match="both dropped and asked for"):
            annotate(
                mc,
                inplace=True,
                genomes_table=checkm,
                genomes_table_format="Name",
                columns=["genomes.lineage"],
                drop=["genomes.lineage"],
            )


class TestContigTableChecks:
    def test_genome_column_is_checked(self, mc, tmp_path):
        table = write(
            tmp_path / "c.tsv",
            "contig_name\tgenome_name\tcircular\ncontig_A\tMAG_2\tyes\n",
        )
        with pytest.raises(ValueError, match="does not belong to in this map"):
            annotate(mc, inplace=True, contigs_table=table)

    def test_unbinned_sentinels_match_unbinned_contigs(self, mc, tmp_path):
        table = write(
            tmp_path / "c.tsv",
            "contig_name\tgenome_name\tcircular\ncontig_D\tND\tyes\n"
            "contig_A\tMAG_1\tno\n",
        )
        annotate(mc, inplace=True, contigs_table=table)
        assert contigs_of(mc).loc["contig_D", "circular"] == "yes"

    def test_two_column_spec_names_the_genome_column(self, mc, tmp_path):
        table = write(
            tmp_path / "c.tsv",
            "Name\tFinal_bin\tGC\ncontig_A\tMAG_1\t0.4\ncontig_C\tMAG_1\t0.5\n",
        )
        with pytest.raises(ValueError, match="'contig_C' is listed under 'MAG_1'"):
            annotate(
                mc, inplace=True, contigs_table=table, contigs_table_format="metator"
            )

    def test_non_ascii_contig_values_are_refused(self, mc, tmp_path):
        table = write(tmp_path / "c.tsv", "contig_name\ttaxon\ncontig_A\tE_coli_é\n")
        with pytest.raises(MetaCoolerFormatError, match="not ASCII"):
            annotate(mc, inplace=True, contigs_table=table)
        assert "taxon" not in contigs_of(mc).columns

    def test_non_ascii_genome_values_are_fine(self, mc, tmp_path):
        table = write(tmp_path / "g.tsv", "genome_name\ttaxon\nMAG_1\tE_coli_é\n")
        annotate(mc, inplace=True, genomes_table=table)
        assert genomes_of(mc).loc["MAG_1", "taxon"] == "E_coli_é"

    def test_a_table_matching_nothing(self, mc, tmp_path):
        genomes = write(tmp_path / "g.tsv", "genome_name\tX\nfoo\t1\n")
        with pytest.raises(ValueError, match="matches no genome"):
            annotate(mc, inplace=True, genomes_table=genomes)
        contigs = write(tmp_path / "c.tsv", "contig_name\tX\nfoo\t1\n")
        with pytest.raises(ValueError, match="matches no contig"):
            annotate(mc, inplace=True, contigs_table=contigs)

    def test_a_table_without_metadata_columns_changes_nothing(self, mc, tmp_path):
        table = write(tmp_path / "c.tsv", "contig_name\ncontig_A\n")
        annotate(mc, inplace=True, contigs_table=table)
        assert history(mc) == []


class TestTargets:
    def test_output_copy_leaves_the_input_alone(self, mc, checkm, tmp_path):
        out = tmp_path / "out.metacool"
        annotate(mc, out, genomes_table=checkm, genomes_table_format="Name")
        assert "Completeness" in genomes_of(out).columns
        assert "Completeness" not in genomes_of(mc).columns
        assert history(out)[0]["added"][0] == "genomes.Completeness"

    def test_unchanged_output_copy_is_still_written(self, mc, tmp_path):
        out = tmp_path / "out.metacool"
        table = write(tmp_path / "c.tsv", "contig_name\ncontig_A\n")
        annotate(mc, out, contigs_table=table)
        assert out.exists()
        validate(out)

    def test_bare_collection_annotates_every_entry(
        self, tiny_mmetacool, checkm, tmp_path
    ):
        out = tmp_path / "out.mmetacool"
        annotate(tiny_mmetacool, out, genomes_table=checkm, genomes_table_format="Name")
        for entry in ("entry1", "entry2"):
            uri = f"{out}::entries/{entry}"
            assert genomes_of(uri).loc["MAG_1", "Completeness"] == 91.5
            validate(uri)
        assert "Completeness" not in genomes_of(f"{tiny_mmetacool}::entries/entry1")

    def test_bare_collection_in_place(self, tiny_mmetacool, depth):
        annotate(
            tiny_mmetacool,
            inplace=True,
            contigs_table=depth,
            contigs_table_format="contigName",
        )
        for entry in ("entry1", "entry2"):
            uri = f"{tiny_mmetacool}::entries/{entry}"
            assert contigs_of(uri).loc["contig_A", "totalAvgDepth"] == 12.5

    def test_entry_in_place_leaves_other_entries_alone(self, tiny_mmetacool, checkm):
        annotate(
            f"{tiny_mmetacool}::entries/entry1",
            inplace=True,
            genomes_table=checkm,
            genomes_table_format="Name",
        )
        assert "Completeness" in genomes_of(f"{tiny_mmetacool}::entries/entry1")
        assert "Completeness" not in genomes_of(f"{tiny_mmetacool}::entries/entry2")

    def test_entry_to_output_is_a_standalone_map(
        self, tiny_mmetacool, checkm, tmp_path
    ):
        out = tmp_path / "entry1.metacool"
        annotate(
            f"{tiny_mmetacool}::entries/entry1",
            out,
            genomes_table=checkm,
            genomes_table_format="Name",
        )
        assert genomes_of(out).loc["MAG_2", "Completeness"] == 70.0
        validate(out)

    def test_leftover_scratch_group_is_cleared(self, mc, checkm):
        with h5py.File(mc, "a") as f:
            f.create_group(_SCRATCH).create_dataset("junk", data=[1])
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        with h5py.File(mc, "r") as f:
            assert _SCRATCH not in f


class TestArguments:
    def test_output_and_inplace(self, mc, checkm, tmp_path):
        with pytest.raises(ValueError, match="pass one or the other"):
            annotate(mc, tmp_path / "o.metacool", inplace=True, genomes_table=checkm)

    def test_neither_output_nor_inplace(self, mc, checkm):
        with pytest.raises(ValueError, match="output path is required"):
            annotate(mc, genomes_table=checkm)

    def test_nothing_to_do(self, mc):
        with pytest.raises(ValueError, match="Nothing to do"):
            annotate(mc, inplace=True)

    def test_bad_missing(self, mc, checkm):
        with pytest.raises(ValueError, match="missing must be"):
            annotate(mc, inplace=True, genomes_table=checkm, missing="drop")  # ty: ignore[invalid-argument-type]

    def test_unqualified_column(self, mc):
        with pytest.raises(ValueError, match="'genomes.<column>'"):
            annotate(mc, inplace=True, drop=["lineage"])


class TestCli:
    def test_annotate_in_place(self, mc, checkm, depth):
        result = CliRunner().invoke(
            cli,
            [
                "annotate",
                str(mc),
                "--genomes",
                str(checkm),
                "--genomes-format",
                "Name",
                "--contigs",
                str(depth),
                "--contigs-format",
                "contigName",
                "--columns",
                "genomes.Completeness,contigs.totalAvgDepth",
                "--inplace",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Completeness" in genomes_of(mc).columns
        assert "Contamination" not in genomes_of(mc).columns
        commands = MetaCooler(mc).metadata["commands"]
        assert "annotate" in commands[-1]

    def test_replace_drop_and_output(self, mc, checkm, tmp_path):
        out = tmp_path / "out.metacool"
        runner = CliRunner()
        args = [
            "annotate",
            str(mc),
            "--genomes",
            str(checkm),
            "--genomes-format",
            "Name",
        ]
        assert runner.invoke(cli, [*args, "--inplace"]).exit_code == 0
        patch = write(tmp_path / "p.tsv", "Name\tCompleteness\nMAG_1\t1\n")
        result = runner.invoke(
            cli,
            [
                "annotate",
                str(mc),
                "--genomes",
                str(patch),
                "--genomes-format",
                "1",
                "--replace",
                "--missing",
                "na",
                "--drop",
                "genomes.lineage",
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        g = genomes_of(out)
        assert g.loc["MAG_1", "Completeness"] == 1.0
        assert np.isnan(g.loc["MAG_2", "Completeness"])
        assert "lineage" not in g.columns

    @pytest.mark.parametrize(
        ("extra", "message"),
        [
            (["--inplace", "-o", "x.metacool"], "Pass one or the other"),
            ([], "Missing option '-o'"),
        ],
    )
    def test_output_usage_errors(self, mc, checkm, extra, message):
        result = CliRunner().invoke(
            cli, ["annotate", str(mc), "--genomes", str(checkm), *extra]
        )
        assert result.exit_code != 0
        assert message in flat(result.output)

    def test_nothing_to_do(self, mc):
        result = CliRunner().invoke(cli, ["annotate", str(mc), "--inplace"])
        assert result.exit_code != 0
        assert "Nothing to do" in flat(result.output)

    def test_bad_spec_is_refused_before_reading(self, mc, depth):
        result = CliRunner().invoke(
            cli,
            [
                "annotate",
                str(mc),
                "--contigs",
                str(depth),
                "--contigs-format",
                "contigName,noheader",
                "--inplace",
            ],
        )
        assert result.exit_code != 0
        assert "no header line" in flat(result.output)

    def test_refuses_to_overwrite_an_input(self, mc, checkm):
        result = CliRunner().invoke(
            cli,
            ["annotate", str(mc), "--genomes", str(checkm), "-o", str(checkm), "-f"],
        )
        assert result.exit_code != 0


def test_provenance_records_tables_and_changes(mc, checkm, depth):
    annotate(
        mc,
        inplace=True,
        genomes_table=checkm,
        genomes_table_format="Name",
        contigs_table=depth,
        contigs_table_format="contigName",
        metadata={"commands": ["metacooler annotate ..."]},
    )
    with h5py.File(mc, "r") as f:
        meta = json.loads(f.attrs["metadata"])
    entry = meta["annotate"][-1]
    assert entry["genomes_table"] == str(checkm)
    assert entry["contigs_table"] == str(depth)
    assert entry["added"] == [
        "genomes.Completeness",
        "genomes.Contamination",
        "genomes.lineage",
        "contigs.totalAvgDepth",
    ]
    assert meta["commands"][-1] == "metacooler annotate ..."


def test_a_collected_map_takes_input_contig_names(tmp_path):
    """As at collect time: 'chromosome' plus its genome, or the stored name."""
    from _cools import BSUB, ECOLI, make_cool

    from metacooler.io.collect import collect

    ecoli = make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)
    bsub = make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)
    out = tmp_path / "collected.metacool"
    collect([ecoli, bsub], out)
    table = write(
        tmp_path / "t.tsv",
        "contig_name\tgenome_name\tcopies\n"
        "chromosome\tE_coli\t1\nB_subtilis__plasmid1\tB_subtilis\t9\n",
    )
    annotate(out, inplace=True, contigs_table=table)
    contigs = contigs_of(out)
    assert contigs.loc["E_coli__chromosome", "copies"] == 1
    assert contigs.loc["B_subtilis__plasmid1", "copies"] == 9
    assert np.isnan(contigs.loc["B_subtilis__chromosome", "copies"])

    shared = write(tmp_path / "s.tsv", "contig_name\tdepth\nchromosome\t1\n")
    with pytest.raises(ValueError, match="more than one genome"):
        annotate(out, inplace=True, contigs_table=shared)


class TestEdges:
    def test_an_empty_collection(self, tiny_mmetacool):
        with h5py.File(tiny_mmetacool, "a") as f:
            for entry in list(f["entries"]):
                del f["entries"][entry]
        with pytest.raises(ValueError, match="collection with no entry"):
            annotate(tiny_mmetacool, inplace=True, drop=["genomes.x"])

    def test_text_conflict_names_the_emptied_value(self, mc, checkm, tmp_path):
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        patch = write(tmp_path / "p.tsv", "Name\tlineage\nMAG_1\tBacteroidota\n")
        with pytest.raises(ValueError, match="'Bacillota' -> \\(empty\\)"):
            annotate(
                mc,
                inplace=True,
                genomes_table=patch,
                genomes_table_format="Name",
                missing="na",
            )

    def test_a_type_change_alone_is_described(self):
        from metacooler.io.annotate import _describe_change

        old = np.array(["1", "2"], dtype=object)
        new = np.array([1, 2])
        assert _describe_change(old, new, np.array(["a", "b"])) == (
            "its type would change"
        )

    def test_group_attributes_survive_a_rewrite(self, mc, checkm):
        with h5py.File(mc, "a") as f:
            f["genomes"].attrs["note"] = "kept"
        annotate(mc, inplace=True, genomes_table=checkm, genomes_table_format="Name")
        with h5py.File(mc, "r") as f:
            assert f["genomes"].attrs["note"] == "kept"

    def test_metadata_stored_as_bytes(self, mc, checkm):
        with h5py.File(mc, "a") as f:
            meta = json.loads(f.attrs["metadata"])
            meta["sample"] = "gut_T0"
            f.attrs["metadata"] = np.bytes_(json.dumps(meta))
        annotate(
            mc,
            inplace=True,
            genomes_table=checkm,
            genomes_table_format="Name",
            metadata={"operator": "me"},
        )
        meta = MetaCooler(mc).metadata
        assert meta["sample"] == "gut_T0"
        assert meta["operator"] == "me"
        assert meta["annotate"][-1]["added"][0] == "genomes.Completeness"
