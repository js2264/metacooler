"""Tests for ``metacooler.collect`` — single-genome cools into one metacool."""

from __future__ import annotations

import logging
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cools import BSUB, ECOLI, make_cool

from metacooler.api import MetaCooler
from metacooler.filters import FilterError
from metacooler.io.collect import (
    COLLECT_ORDERINGS,
    CONTIG_NAME_SEPARATOR,
    collect,
    collected_contig_name,
)
from metacooler.qc import QC_CONTIG_COLUMNS
from metacooler.schema import MetaCoolerFormatError, check, validate

DATA = Path(__file__).parent / "data"


@pytest.fixture
def ecoli(tmp_path: Path) -> Path:
    return make_cool(tmp_path / "E_coli.cool", ECOLI, seed=1)


@pytest.fixture
def bsub(tmp_path: Path) -> Path:
    return make_cool(tmp_path / "B_subtilis.cool", BSUB, seed=2)


@pytest.fixture
def collected(tmp_path: Path, ecoli: Path, bsub: Path) -> Path:
    out = tmp_path / "collected.metacool"
    collect([ecoli, bsub], out)
    return out


def pixels_of(path: Path | str) -> pd.DataFrame:
    return pd.DataFrame(cooler.Cooler(str(path)).pixels()[:]).reset_index(drop=True)


def joined_contacts(path: Path) -> set[tuple]:
    """Every contact as (contig1, start1, contig2, start2, count), layout-free.

    Resolved through ``bins()`` rather than ``pixels(join=True)``: on a metacool
    file ``/bins/chrom`` is a plain integer column, so cooler's join reports
    contig *ids*, which differ between orderings.
    """
    clr = cooler.Cooler(str(path))
    bins = clr.bins()[:]
    contig = bins["chrom"].astype(str).to_numpy()
    start = bins["start"].to_numpy()
    px = clr.pixels()[:]
    rows = set()
    for b1, b2, n in zip(px["bin1_id"], px["bin2_id"], px["count"], strict=True):
        # A reordering may flip which end of an inter-contig pair is bin1.
        a, b = (contig[b1], int(start[b1])), (contig[b2], int(start[b2]))
        rows.add((*min(a, b), *max(a, b), int(n)))
    return rows


# ----------------------------------------------------------------------------
# The basic shape of a collection
# ----------------------------------------------------------------------------


class TestCollectBasics:
    def test_each_input_becomes_one_genome(self, collected: Path):
        validate(collected)
        mc = MetaCooler(collected)
        genomes = mc.genomes
        # Genome rows are alphabetical, whatever the input order.
        assert genomes["name"].tolist() == ["B_subtilis", "E_coli"]
        assert genomes["length"].tolist() == [1100, 1650]
        assert genomes["n_contigs"].tolist() == [2, 3]
        assert genomes["n_bins"].tolist() == [11, 17]

    def test_every_contig_is_prefixed_with_its_genome(self, collected: Path):
        contigs = MetaCooler(collected).contigs
        assert contigs["name"].tolist() == [
            "B_subtilis__chromosome",
            "B_subtilis__plasmid1",
            "E_coli__chromosome",
            "E_coli__plasmid1",
            "E_coli__plasmid2",
        ]
        assert contigs["genome_id"].tolist() == [0, 0, 1, 1, 1]

    def test_opens_with_cooler(self, collected: Path, ecoli: Path, bsub: Path):
        clr = cooler.Cooler(str(collected))
        assert len(clr.chromnames) == 5
        a, b = cooler.Cooler(str(ecoli)), cooler.Cooler(str(bsub))
        assert clr.info["nnz"] == a.info["nnz"] + b.info["nnz"]
        assert clr.info["sum"] == a.info["sum"] + b.info["sum"]

    def test_inter_genome_blocks_are_empty(self, collected: Path):
        mc = MetaCooler(collected)
        genome_of_bin = mc.contigs["genome_id"].to_numpy()[
            mc.bins()["chrom"].cat.codes.to_numpy()
        ]
        px = pixels_of(collected)
        assert (
            genome_of_bin[px["bin1_id"].to_numpy()]
            == genome_of_bin[px["bin2_id"].to_numpy()]
        ).all()

    def test_passes_the_deep_check(self, collected: Path):
        results = check(collected, deep=True)
        assert all(r.ok for r in results), [r for r in results if not r.ok]

    def test_records_format_attributes(self, collected: Path):
        info = MetaCooler(collected).info
        assert info["ordering"] == "cool"
        assert info["bin-size"] == 100
        assert info["bin-type"] == "fixed"
        assert info["assembly"] == "unknown"

    def test_single_input_is_a_single_genome_file(self, tmp_path: Path, ecoli: Path):
        out = tmp_path / "one.metacool"
        collect([ecoli], out)
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist() == ["E_coli"]
        assert mc.contigs["name"].tolist() == [
            "E_coli__chromosome",
            "E_coli__plasmid1",
            "E_coli__plasmid2",
        ]
        assert pixels_of(out).equals(pixels_of(ecoli))

    def test_accepts_strings_and_paths(self, tmp_path: Path, ecoli: Path, bsub: Path):
        out = tmp_path / "mixed.metacool"
        collect([str(ecoli), bsub], str(out))
        validate(out)

    @pytest.mark.parametrize("as_str", [True, False])
    def test_a_lone_path_is_one_input(self, tmp_path: Path, ecoli: Path, as_str: bool):
        out = tmp_path / "lone.metacool"
        collect(str(ecoli) if as_str else ecoli, out)  # ty:ignore[invalid-argument-type]
        assert MetaCooler(out).genomes["name"].tolist() == ["E_coli"]

    def test_names_as_a_single_string_is_refused(
        self, tmp_path: Path, ecoli: Path, bsub: Path
    ):
        with pytest.raises(TypeError, match="not a single string"):
            collect([ecoli, bsub], tmp_path / "out.metacool", names="AB")


class TestRoundTrip:
    """Under the default ordering, collecting loses nothing."""

    @pytest.mark.parametrize(
        ("genome", "fixture"), [("E_coli", "ecoli"), ("B_subtilis", "bsub")]
    )
    def test_extract_gives_back_each_input(
        self, request, tmp_path: Path, collected: Path, genome: str, fixture: str
    ):
        from metacooler.io.extract import extract

        source = cooler.Cooler(str(request.getfixturevalue(fixture)))
        out = tmp_path / f"{genome}.cool"
        extract(collected, genomes=[genome], output=out)
        back = cooler.Cooler(str(out))

        prefix = f"{genome}{CONTIG_NAME_SEPARATOR}"
        assert back.chromnames == [prefix + c for c in source.chromnames]
        assert (
            back.chroms()[:]["length"].tolist() == source.chroms()[:]["length"].tolist()
        )
        assert back.bins()[["start", "end"]][:].equals(
            source.bins()[["start", "end"]][:]
        )
        assert pixels_of(out).equals(pixels_of(source.uri))

    def test_matches_create_cool_for_one_input(self, tmp_path: Path, ecoli: Path):
        """Collecting one cool is importing it with every contig in one genome."""
        from metacooler.io.create import from_cool

        table = tmp_path / "all_in_G.tsv"
        table.write_text(
            "contig_name\tgenome_name\n" + "".join(f"{c}\tG\n" for c in ECOLI)
        )
        imported = tmp_path / "imported.metacool"
        from_cool(ecoli, imported, table, contig_ordering="cool")
        collected = tmp_path / "collected.metacool"
        collect([ecoli], collected, names=["G"])

        a, b = MetaCooler(imported), MetaCooler(collected)
        core = ["name", "length", "n_contigs", "n_bins"]
        assert a.genomes[core].equals(b.genomes[core])
        assert a.contigs[["length", "genome_id"]].equals(
            b.contigs[["length", "genome_id"]]
        )
        assert b.contigs["name"].tolist() == ["G__" + n for n in a.contigs["name"]]
        assert a.bins()[["start", "end"]].equals(b.bins()[["start", "end"]])
        assert pixels_of(imported).equals(pixels_of(collected))

    def test_prefix_is_not_stacked_on_recollection(
        self, tmp_path: Path, collected: Path
    ):
        from metacooler.io.extract import extract

        extracted = tmp_path / "again" / "E_coli.cool"
        extracted.parent.mkdir()
        extract(collected, genomes=["E_coli"], output=extracted)
        out = tmp_path / "again.metacool"
        collect([extracted], out)
        assert MetaCooler(out).contigs["name"].tolist() == [
            "E_coli__chromosome",
            "E_coli__plasmid1",
            "E_coli__plasmid2",
        ]


# ----------------------------------------------------------------------------
# Contig ordering and the pixel sort
# ----------------------------------------------------------------------------


class TestOrdering:
    def test_orderings_constant_matches_cli(self):
        from metacooler.cli.cmd_collect import _ORDERINGS

        assert list(COLLECT_ORDERINGS) == _ORDERINGS

    def test_default_keeps_each_inputs_contig_order(self, collected: Path, caplog):
        names = MetaCooler(collected).contigs["name"].tolist()
        assert names[2:] == [
            "E_coli__chromosome",
            "E_coli__plasmid1",
            "E_coli__plasmid2",
        ]

    def test_default_needs_no_sort(
        self, tmp_path: Path, ecoli: Path, bsub: Path, caplog
    ):
        with caplog.at_level(logging.INFO, logger="metacooler"):
            collect([ecoli, bsub], tmp_path / "out.metacool")
        assert "without sorting" in caplog.text
        assert "Re-sorting" not in caplog.text

    def test_genome_then_length(self, tmp_path: Path, ecoli: Path, bsub: Path, caplog):
        out = tmp_path / "gtl.metacool"
        with caplog.at_level(logging.INFO, logger="metacooler"):
            collect([ecoli, bsub], out, contig_ordering="genome-then-length")
        assert "Re-sorting the pixels of 1 input(s)" in caplog.text
        mc = MetaCooler(out)
        assert mc.info["ordering"] == "genome-then-length"
        assert mc.contigs["name"].tolist() == [
            "B_subtilis__chromosome",
            "B_subtilis__plasmid1",
            "E_coli__chromosome",
            "E_coli__plasmid2",
            "E_coli__plasmid1",
        ]
        assert all(r.ok for r in check(out, deep=True))

    def test_genome_then_alpha(self, tmp_path: Path, ecoli: Path, bsub: Path):
        rev = make_cool(tmp_path / "rev.cool", {"z": 300, "a": 200}, seed=5)
        out = tmp_path / "gta.metacool"
        collect([ecoli, rev], out, contig_ordering="genome-then-alpha")
        assert MetaCooler(out).contigs["name"].tolist()[-2:] == ["rev__a", "rev__z"]
        assert all(r.ok for r in check(out, deep=True))

    @pytest.mark.parametrize("ordering", ["genome-then-length", "genome-then-alpha"])
    def test_orderings_hold_the_same_contacts(
        self, tmp_path: Path, collected: Path, ecoli: Path, bsub: Path, ordering: str
    ):
        out = tmp_path / f"{ordering}.metacool"
        collect([ecoli, bsub], out, contig_ordering=ordering)  # ty:ignore[invalid-argument-type]
        assert joined_contacts(out) == joined_contacts(collected)

    @pytest.mark.parametrize("sort_in_memory", [None, True, False])
    @pytest.mark.parametrize("chunksize", [7, 1000])
    def test_sort_strategies_and_chunking_agree(
        self, tmp_path: Path, ecoli: Path, bsub: Path, sort_in_memory, chunksize: int
    ):
        reference = tmp_path / "reference.metacool"
        collect([ecoli, bsub], reference, contig_ordering="genome-then-length")
        out = tmp_path / "out.metacool"
        collect(
            [ecoli, bsub],
            out,
            contig_ordering="genome-then-length",
            sort_in_memory=sort_in_memory,
            chunksize=chunksize,
        )
        assert pixels_of(out).equals(pixels_of(reference))
        assert all(r.ok for r in check(out, deep=True))

    def test_small_chunks_under_the_default_order(
        self, tmp_path: Path, collected: Path, ecoli: Path, bsub: Path
    ):
        out = tmp_path / "chunked.metacool"
        collect([ecoli, bsub], out, chunksize=3)
        assert pixels_of(out).equals(pixels_of(collected))

    def test_large_input_switches_to_external_sort(
        self, tmp_path: Path, ecoli: Path, bsub: Path, monkeypatch, caplog
    ):
        import metacooler.io.collect as collect_module

        reference = tmp_path / "reference.metacool"
        collect([ecoli, bsub], reference, contig_ordering="genome-then-length")
        monkeypatch.setattr(collect_module, "_SORT_IN_MEMORY_PIXEL_THRESHOLD", 10)
        out = tmp_path / "external.metacool"
        with caplog.at_level(logging.INFO, logger="metacooler"):
            collect([ecoli, bsub], out, contig_ordering="genome-then-length")
        assert "switching to cooler's external sort" in caplog.text
        assert pixels_of(out).equals(pixels_of(reference))

    def test_unknown_ordering_is_refused(self, tmp_path: Path, ecoli: Path):
        with pytest.raises(ValueError, match="Unknown contig ordering 'alpha'"):
            collect([ecoli], tmp_path / "out.metacool", contig_ordering="alpha")  # ty:ignore[invalid-argument-type]


# ----------------------------------------------------------------------------
# Genome and contig names
# ----------------------------------------------------------------------------


class TestNames:
    def test_collected_contig_name(self):
        assert collected_contig_name("E_coli", "chr1") == "E_coli__chr1"
        assert collected_contig_name("E_coli", "E_coli__chr1") == "E_coli__chr1"
        # Only the exact genome prefix is recognised.
        assert collected_contig_name("E", "E_coli__chr1") == "E__E_coli__chr1"

    def test_names_override_stems(self, tmp_path: Path, ecoli: Path, bsub: Path):
        out = tmp_path / "named.metacool"
        collect([ecoli, bsub], out, names=["Ec", "Bs"])
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist() == ["Bs", "Ec"]
        assert mc.contigs["name"].tolist()[0] == "Bs__chromosome"

    def test_names_count_must_match(self, tmp_path: Path, ecoli: Path, bsub: Path):
        with pytest.raises(ValueError, match="1 genome name.*for 2 input"):
            collect([ecoli, bsub], tmp_path / "out.metacool", names=["only_one"])

    def test_names_must_be_unique(self, tmp_path: Path, ecoli: Path, bsub: Path):
        with pytest.raises(ValueError, match="'same'"):
            collect([ecoli, bsub], tmp_path / "out.metacool", names=["same", "same"])

    def test_names_must_not_be_empty(self, tmp_path: Path, ecoli: Path, bsub: Path):
        with pytest.raises(ValueError, match="empty genome name"):
            collect([ecoli, bsub], tmp_path / "out.metacool", names=["a", " "])

    def test_clashing_stems_suggest_names(self, tmp_path: Path):
        for sub, seed in (("run1", 1), ("run2", 2)):
            (tmp_path / sub).mkdir()
            make_cool(tmp_path / sub / "sample.cool", ECOLI, seed=seed)
        with pytest.raises(ValueError, match="use --names") as exc:
            collect(
                [tmp_path / "run1" / "sample.cool", tmp_path / "run2" / "sample.cool"],
                tmp_path / "out.metacool",
            )
        assert "'sample'" in str(exc.value)

    def test_same_file_twice_under_two_names(self, tmp_path: Path, ecoli: Path):
        out = tmp_path / "twice.metacool"
        collect([ecoli, ecoli], out, names=["rep1", "rep2"])
        assert MetaCooler(out).genomes["name"].tolist() == ["rep1", "rep2"]

    def test_collected_name_clash_is_refused(self, tmp_path: Path):
        # 'A' + 'B__x' and 'A__B' + 'x' would both be collected as 'A__B__x'.
        a = make_cool(tmp_path / "A.cool", {"B__x": 300})
        ab = make_cool(tmp_path / "A__B.cool", {"x": 300})
        with pytest.raises(ValueError, match="'A__B__x'") as exc:
            collect([a, ab], tmp_path / "out.metacool")
        assert "--names" in str(exc.value)

    def test_unbinned_named_input_takes_genome_id_zero(
        self, tmp_path: Path, ecoli: Path, bsub: Path
    ):
        out = tmp_path / "unbinned.metacool"
        collect([ecoli, bsub], out, names=["E_coli", "__unbinned__"])
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist() == ["__unbinned__", "E_coli"]
        validate(out)

    def test_non_ascii_genome_name_is_refused_before_writing(
        self, tmp_path: Path, ecoli: Path, bsub: Path
    ):
        out = tmp_path / "accent.metacool"
        with pytest.raises(MetaCoolerFormatError, match="not ASCII"):
            collect([ecoli, bsub], out, names=["É_coli", "B_subtilis"])
        assert not out.exists()
        assert not out.with_name(out.name + ".partial").exists()


# ----------------------------------------------------------------------------
# What collect refuses
# ----------------------------------------------------------------------------


class TestInputs:
    def test_mcool_resolutions(self, tmp_path: Path):
        out = tmp_path / "mcool.metacool"
        collect(
            [
                DATA / "assembly_1_rep1_1kb.mcool::resolutions/2000",
                f"{DATA / 'assembly_2_rep1_1kb.mcool'}::/resolutions/2000",
            ],
            out,
        )
        mc = MetaCooler(out)
        assert mc.info["bin-size"] == 2000
        assert mc.genomes["name"].tolist() == [
            "assembly_1_rep1_1kb",
            "assembly_2_rep1_1kb",
        ]
        source = cooler.Cooler(
            f"{DATA / 'assembly_1_rep1_1kb.mcool'}::resolutions/2000"
        )
        assert int(mc.genomes["n_bins"].iloc[0]) == len(source.bins())

    def test_mcool_root_is_refused_with_its_resolutions(self, tmp_path: Path, ecoli):
        with pytest.raises(ValueError, match="multi-resolution") as exc:
            collect(
                [ecoli, DATA / "assembly_1_rep1_1kb.mcool"], tmp_path / "o.metacool"
            )
        assert "1000, 2000, 5000" in str(exc.value)
        assert "assembly_1_rep1_1kb.mcool::resolutions/1000" in str(exc.value)

    def test_metacool_is_refused(self, tmp_path: Path, ecoli: Path, tiny_metacooler):
        with pytest.raises(ValueError, match="already a metacool file") as exc:
            collect([ecoli, tiny_metacooler], tmp_path / "out.metacool")
        assert "metacooler flatten" in str(exc.value)

    def test_mmetacool_is_refused(self, tmp_path: Path, ecoli: Path, tiny_mmetacool):
        with pytest.raises(ValueError, match="mmetacool collection"):
            collect([ecoli, tiny_mmetacool], tmp_path / "out.metacool")
        with pytest.raises(ValueError, match="already a metacool file"):
            collect([f"{tiny_mmetacool}::entries/entry1"], tmp_path / "out.metacool")

    def test_not_a_cooler_is_refused(self, tmp_path: Path, ecoli: Path):
        junk = tmp_path / "junk.cool"
        junk.write_text("not hdf5")
        with pytest.raises(ValueError, match="not a cooler file"):
            collect([ecoli, junk], tmp_path / "out.metacool")

    @pytest.mark.filterwarnings("ignore:Creating a non-symmetric matrix")
    def test_square_storage_is_refused(self, tmp_path: Path, ecoli: Path):
        square = make_cool(tmp_path / "square.cool", BSUB, symmetric_upper=False)
        with pytest.raises(ValueError, match="square"):
            collect([ecoli, square], tmp_path / "out.metacool")

    def test_missing_input(self, tmp_path: Path, ecoli: Path):
        with pytest.raises(FileNotFoundError):
            collect([ecoli, tmp_path / "nope.cool"], tmp_path / "out.metacool")

    def test_missing_group_is_not_a_cooler(self, tmp_path: Path, ecoli: Path):
        with pytest.raises(ValueError, match="not a cooler file"):
            collect([f"{ecoli}::/resolutions/1000"], tmp_path / "out.metacool")

    def test_chunksize_must_be_positive(self, tmp_path: Path, ecoli: Path):
        with pytest.raises(ValueError, match="chunksize"):
            collect([ecoli], tmp_path / "out.metacool", chunksize=0)

    def test_no_inputs(self, tmp_path: Path):
        with pytest.raises(ValueError, match="at least one input"):
            collect([], tmp_path / "out.metacool")

    def test_output_may_not_be_an_input(self, tmp_path: Path, ecoli: Path, bsub: Path):
        before = ecoli.read_bytes()
        with pytest.raises(ValueError, match="never overwritten"):
            collect([ecoli, bsub], ecoli)
        assert ecoli.read_bytes() == before

    def test_mixed_resolutions_are_refused(self, tmp_path: Path, ecoli: Path):
        coarse = make_cool(tmp_path / "coarse.cool", BSUB, binsize=200)
        with pytest.raises(ValueError, match="not binned alike") as exc:
            collect([ecoli, coarse], tmp_path / "out.metacool")
        assert "100 bp" in str(exc.value) and "200 bp" in str(exc.value)

    def test_variable_bins_collect_as_variable(self, tmp_path: Path):
        def fragments(chroms: dict[str, int]) -> pd.DataFrame:
            rows = []
            for chrom, length in chroms.items():
                edges = [0, length // 5, length // 2, length]
                rows += [
                    {"chrom": chrom, "start": s, "end": e}
                    for s, e in zip(edges[:-1], edges[1:], strict=True)
                ]
            return pd.DataFrame(rows)

        a = make_cool(tmp_path / "a.cool", ECOLI, bins=fragments(ECOLI))
        b = make_cool(tmp_path / "b.cool", BSUB, bins=fragments(BSUB))
        out = tmp_path / "var.metacool"
        collect([a, b], out)
        info = MetaCooler(out).info
        assert info["bin-type"] == "variable"
        assert info["bin-size"] is None
        validate(out)

    def test_count_dtype_is_widened_not_wrapped(self, tmp_path: Path, ecoli: Path):
        big = make_cool(
            tmp_path / "big.cool", BSUB, count_dtype=np.int64, count_scale=3_000_000_000
        )
        out = tmp_path / "wide.metacool"
        collect([ecoli, big], out)
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"].dtype == np.int64
        assert int(pixels_of(out)["count"].max()) >= 3_000_000_000
        total = (
            cooler.Cooler(str(ecoli)).info["sum"] + cooler.Cooler(str(big)).info["sum"]
        )
        assert int(cooler.Cooler(str(out)).info["sum"]) == total

    def test_float_counts_stay_float(self, tmp_path: Path, ecoli: Path):
        floaty = make_cool(tmp_path / "float.cool", BSUB, count_dtype=np.float64)
        out = tmp_path / "float.metacool"
        collect([ecoli, floaty], out)
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"].dtype.kind == "f"

    def test_empty_input_contributes_bins_only(self, tmp_path: Path, ecoli: Path):
        empty = make_cool(tmp_path / "empty.cool", BSUB, n_contacts=0)
        out = tmp_path / "with_empty.metacool"
        collect([ecoli, empty], out)
        mc = MetaCooler(out)
        assert mc.genomes["n_bins"].tolist() == [17, 11]
        assert mc.info["nnz"] == cooler.Cooler(str(ecoli)).info["nnz"]

    def test_all_inputs_empty(self, tmp_path: Path):
        a = make_cool(tmp_path / "a.cool", ECOLI, n_contacts=0)
        b = make_cool(tmp_path / "b.cool", BSUB, n_contacts=0)
        out = tmp_path / "void.metacool"
        collect([a, b], out)
        assert MetaCooler(out).info["nnz"] == 0
        validate(out)

    def test_balancing_weights_are_dropped_with_a_note(
        self, tmp_path: Path, ecoli: Path, bsub: Path, caplog
    ):
        with h5py.File(ecoli, "r+") as f:
            f["bins"].create_dataset("weight", data=np.ones(len(f["bins"]["start"])))
        out = tmp_path / "weights.metacool"
        with caplog.at_level(logging.INFO, logger="metacooler"):
            collect([ecoli, bsub], out)
        assert "Balancing weights" in caplog.text
        with h5py.File(out, "r") as f:
            assert "weight" not in f["bins"]

    def test_contig_without_bins_is_left_out(
        self, tmp_path: Path, ecoli: Path, bsub: Path, caplog
    ):
        # cooler cannot write a chrom without bins, so graft one onto /chroms.
        with h5py.File(bsub, "r+") as f:
            names = [*f["chroms"]["name"][:].astype(str), "ghost"]
            lengths = [*f["chroms"]["length"][:], 500]
            del f["chroms"]["name"], f["chroms"]["length"]
            f["chroms"].create_dataset("name", data=np.array(names, dtype="S"))
            f["chroms"].create_dataset("length", data=np.array(lengths, dtype=np.int32))
            f.attrs["nchroms"] = len(names)
        out = tmp_path / "ghost.metacool"
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            collect([ecoli, bsub], out)
        assert "Left out 1 contig(s) that have no bins" in caplog.text
        assert "B_subtilis__ghost" not in MetaCooler(out).contigs["name"].tolist()
        validate(out)


# ----------------------------------------------------------------------------
# Metadata: genomes table, contigs table, metadata dict, QC
# ----------------------------------------------------------------------------


@pytest.fixture
def genomes_tsv(tmp_path: Path) -> Path:
    path = tmp_path / "species.tsv"
    path.write_text(
        "genome_name\tphylum\tgc\n"
        "E_coli\tProteobacteria\t0.51\n"
        "S_aureus\tFirmicutes\t0.33\n"
    )
    return path


class TestGenomesTable:
    def test_decorates_the_matching_genomes(
        self, tmp_path: Path, ecoli: Path, bsub: Path, genomes_tsv: Path
    ):
        out = tmp_path / "decorated.metacool"
        collect([ecoli, bsub], out, genomes_table=genomes_tsv)
        genomes = MetaCooler(out).genomes.set_index("name")
        assert genomes.loc["E_coli", "phylum"] == "Proteobacteria"
        assert genomes.loc["E_coli", "gc"] == pytest.approx(0.51)
        assert genomes.loc["B_subtilis", "phylum"] == ""
        assert "S_aureus" not in genomes.index

    def test_reports_divergence_in_both_directions_once(
        self, tmp_path: Path, ecoli: Path, bsub: Path, genomes_tsv: Path, caplog
    ):
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            collect(
                [ecoli, bsub],
                tmp_path / "out.metacool",
                genomes_table=genomes_tsv,
                filters=["contigs.length > 1"],
            )
        assert caplog.text.count("have no row in the genomes table") == 1
        assert "'B_subtilis'" in caplog.text
        assert caplog.text.count("match no input") == 1
        assert "'S_aureus'" in caplog.text
        # The create wording talks about a contigs table collect does not have.
        assert "contigs table" not in caplog.text

    def test_reserved_column_is_renamed(self, tmp_path: Path, ecoli: Path):
        table = tmp_path / "g.tsv"
        table.write_text("genome_name\tlength\nE_coli\t5Mb\n")
        out = tmp_path / "renamed.metacool"
        collect([ecoli], out, genomes_table=table)
        genomes = MetaCooler(out).genomes
        assert genomes["length"].tolist() == [1650]
        assert genomes["user_length"].tolist() == ["5Mb"]


class TestContigsTable:
    def collect_with(self, tmp_path: Path, ecoli: Path, bsub: Path, text: str):
        table = tmp_path / "contig_meta.tsv"
        table.write_text(text)
        out = tmp_path / "contig_meta.metacool"
        collect([ecoli, bsub], out, contigs_table=table)
        return MetaCooler(out).contigs.set_index("name")

    def test_original_names_with_genome_column(self, tmp_path, ecoli, bsub):
        contigs = self.collect_with(
            tmp_path,
            ecoli,
            bsub,
            "contig_name\tgenome_name\tcoverage\n"
            "chromosome\tE_coli\t10\n"
            "chromosome\tB_subtilis\t20\n"
            "plasmid2\tE_coli\t30\n",
        )
        assert contigs.loc["E_coli__chromosome", "coverage"] == 10
        assert contigs.loc["B_subtilis__chromosome", "coverage"] == 20
        assert contigs.loc["E_coli__plasmid2", "coverage"] == 30
        assert np.isnan(contigs.loc["E_coli__plasmid1", "coverage"])
        assert "genome_name" not in contigs.columns

    def test_unique_original_names_without_genome_column(self, tmp_path, ecoli, bsub):
        contigs = self.collect_with(
            tmp_path, ecoli, bsub, "contig_name\tgc\nplasmid2\t0.4\n"
        )
        assert contigs.loc["E_coli__plasmid2", "gc"] == pytest.approx(0.4)

    def test_collected_names(self, tmp_path, ecoli, bsub):
        contigs = self.collect_with(
            tmp_path,
            ecoli,
            bsub,
            "contig_name\ttaxon\nB_subtilis__plasmid1\tpBS\nE_coli__plasmid1\tpEC\n",
        )
        assert contigs.loc["B_subtilis__plasmid1", "taxon"] == "pBS"
        assert contigs.loc["E_coli__plasmid1", "taxon"] == "pEC"

    def test_shared_name_without_genome_column_is_ambiguous(
        self, tmp_path, ecoli, bsub
    ):
        with pytest.raises(ValueError, match="ambiguously") as exc:
            self.collect_with(
                tmp_path, ecoli, bsub, "contig_name\tcoverage\nchromosome\t10\n"
            )
        assert "genome_name column" in str(exc.value)

    def test_contig_under_the_wrong_genome_is_refused(self, tmp_path, ecoli, bsub):
        with pytest.raises(ValueError, match="does not come from") as exc:
            self.collect_with(
                tmp_path,
                ecoli,
                bsub,
                "contig_name\tgenome_name\tcoverage\nplasmid2\tB_subtilis\t1\n",
            )
        assert "'E_coli' (" in str(exc.value)

    def test_collected_name_under_the_wrong_genome_is_refused(
        self, tmp_path, ecoli, bsub
    ):
        with pytest.raises(ValueError, match="does not come from"):
            self.collect_with(
                tmp_path,
                ecoli,
                bsub,
                "contig_name\tgenome_name\tcoverage\nE_coli__chromosome\tB_subtilis\t1\n",
            )

    def test_unmatched_rows_warn(self, tmp_path, ecoli, bsub, caplog):
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            contigs = self.collect_with(
                tmp_path,
                ecoli,
                bsub,
                "contig_name\tgenome_name\tcoverage\n"
                "chromosome\tS_aureus\t5\n"  # a genome that was not collected
                "nowhere\tE_coli\t6\n"
                "chromosome\tE_coli\t7\n",
            )
        assert "2 of 3 row(s) match no collected contig" in caplog.text
        assert contigs.loc["E_coli__chromosome", "coverage"] == 7

    def test_repeated_rows_warn_and_last_wins(self, tmp_path, ecoli, bsub, caplog):
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            contigs = self.collect_with(
                tmp_path,
                ecoli,
                bsub,
                "contig_name\tcoverage\nplasmid2\t1\nE_coli__plasmid2\t2\n",
            )
        assert "more than one row" in caplog.text
        assert contigs.loc["E_coli__plasmid2", "coverage"] == 2

    def test_reserved_columns_are_renamed(self, tmp_path, ecoli, bsub):
        contigs = self.collect_with(
            tmp_path,
            ecoli,
            bsub,
            "contig_name\tgenome_name\tlength\tgenome_id\nchromosome\tE_coli\t999\tX\n",
        )
        assert contigs.loc["E_coli__chromosome", "length"] == 1000
        assert contigs.loc["E_coli__chromosome", "user_length"] == 999
        assert contigs.loc["E_coli__chromosome", "user_genome_id"] == "X"
        assert contigs.loc["E_coli__chromosome", "genome_id"] == 1

    def test_first_column_must_be_contig_name(self, tmp_path, ecoli, bsub):
        with pytest.raises(ValueError, match="'contig_name'"):
            self.collect_with(tmp_path, ecoli, bsub, "name\tcoverage\nchromosome\t1\n")

    def test_table_without_metadata_is_still_checked(self, tmp_path, ecoli, bsub):
        with pytest.raises(ValueError, match="does not come from"):
            self.collect_with(
                tmp_path,
                ecoli,
                bsub,
                "contig_name\tgenome_name\nplasmid2\tB_subtilis\n",
            )

    def test_non_ascii_metadata_is_refused(self, tmp_path, ecoli, bsub):
        with pytest.raises(MetaCoolerFormatError, match="not ASCII"):
            self.collect_with(
                tmp_path, ecoli, bsub, "contig_name\tnote\nplasmid2\tcafé\n"
            )

    def test_table_with_only_key_columns_adds_nothing(self, tmp_path, ecoli, bsub):
        contigs = self.collect_with(
            tmp_path,
            ecoli,
            bsub,
            "contig_name\tgenome_name\nplasmid2\tE_coli\n",
        )
        # Nothing from the table: only metacooler's own columns, QC included.
        assert list(contigs.columns) == ["length", "genome_id", *QC_CONTIG_COLUMNS]

    def test_missing_or_empty_table(self, tmp_path, ecoli, bsub):
        with pytest.raises(FileNotFoundError, match="Contigs table not found"):
            collect(
                [ecoli, bsub], tmp_path / "o.metacool", contigs_table=tmp_path / "no"
            )
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        with pytest.raises(ValueError, match="Contigs table is empty"):
            collect([ecoli, bsub], tmp_path / "o.metacool", contigs_table=empty)


class TestFilters:
    def test_contig_filter_drops_contigs_and_their_pixels(
        self, tmp_path: Path, collected: Path, ecoli: Path, bsub: Path
    ):
        out = tmp_path / "filtered.metacool"
        collect([ecoli, bsub], out, filters=["contigs.length >= 300"])
        mc = MetaCooler(out)
        assert "E_coli__plasmid1" not in mc.contigs["name"].tolist()
        # Exactly the contacts touching E_coli's 250 bp plasmid are gone.
        kept = {
            c
            for c in joined_contacts(collected)
            if "E_coli__plasmid1" not in (c[0], c[2])
        }
        assert joined_contacts(out) == kept
        assert all(r.ok for r in check(out, deep=True))

    @pytest.mark.parametrize("sort_in_memory", [True, False])
    def test_filter_inside_a_reordered_block(
        self, tmp_path: Path, ecoli: Path, bsub: Path, sort_in_memory: bool
    ):
        """Dropped bins and a re-sort in the same block, on both sort paths."""
        cool = tmp_path / "cool.metacool"
        collect([ecoli, bsub], cool, filters=["contigs.length != 1000"])
        out = tmp_path / "reordered.metacool"
        collect(
            [ecoli, bsub],
            out,
            filters=["contigs.length != 1000"],
            contig_ordering="genome-then-length",
            sort_in_memory=sort_in_memory,
            chunksize=5,
        )
        assert MetaCooler(out).contigs["name"].tolist()[-2:] == [
            "E_coli__plasmid2",
            "E_coli__plasmid1",
        ]
        assert joined_contacts(out) == joined_contacts(cool)
        assert all(r.ok for r in check(out, deep=True))

    def test_genome_filter_drops_a_whole_genome(self, tmp_path, ecoli, bsub):
        out = tmp_path / "big.metacool"
        collect([ecoli, bsub], out, filters=["genomes.length > 1200"])
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist() == ["E_coli"]
        assert mc.info["nnz"] == cooler.Cooler(str(ecoli)).info["nnz"]

    def test_n_bins_is_filterable(self, tmp_path, ecoli, bsub):
        out = tmp_path / "nbins.metacool"
        collect([ecoli, bsub], out, filters=["genomes.n_bins < 15"])
        assert MetaCooler(out).genomes["name"].tolist() == ["B_subtilis"]

    def test_filter_on_contig_metadata(self, tmp_path, ecoli, bsub):
        table = tmp_path / "cov.tsv"
        table.write_text(
            "contig_name\tgenome_name\tcoverage\n"
            "chromosome\tE_coli\t50\nchromosome\tB_subtilis\t5\n"
        )
        out = tmp_path / "cov.metacool"
        collect(
            [ecoli, bsub], out, contigs_table=table, filters=["contigs.coverage > 10"]
        )
        assert MetaCooler(out).contigs["name"].tolist() == ["E_coli__chromosome"]

    def test_filter_on_collected_names(self, tmp_path, ecoli, bsub):
        out = tmp_path / "chromosomes.metacool"
        collect([ecoli, bsub], out, filters=["contigs.name =~ __chromosome$"])
        assert MetaCooler(out).contigs["name"].tolist() == [
            "B_subtilis__chromosome",
            "E_coli__chromosome",
        ]

    def test_filter_that_keeps_nothing(self, tmp_path, ecoli, bsub):
        with pytest.raises(FilterError, match="Nothing passed"):
            collect(
                [ecoli, bsub],
                tmp_path / "out.metacool",
                filters=["genomes.length > 1Gb"],
            )


class TestMetadataAndQC:
    def test_metadata_and_assembly(self, tmp_path, ecoli, bsub):
        out = tmp_path / "meta.metacool"
        collect(
            [ecoli, bsub],
            out,
            metadata={"study": "demo", "commands": ["metacooler collect ..."]},
            assembly_name="two_species",
        )
        mc = MetaCooler(out)
        assert mc.metadata["study"] == "demo"
        assert mc.metadata["commands"] == ["metacooler collect ..."]
        assert mc.info["assembly"] == "two_species"

    def test_qc_block(self, collected: Path):
        mc = MetaCooler(collected)
        qc = mc.metadata["qc"]
        assert qc["source"] == "pixels"
        assert qc["binsize"] == 100
        assert qc["inter_genome"] == 0
        assert qc["n_contacts"] == int(mc.info["sum"])
        assert mc.genomes["qc_inter_genome"].tolist() == [0, 0]

    def test_no_qc(self, tmp_path, ecoli, bsub):
        out = tmp_path / "noqc.metacool"
        collect([ecoli, bsub], out, qc=False)
        mc = MetaCooler(out)
        assert "qc" not in mc.metadata
        assert not any(c.startswith("qc_") for c in mc.genomes.columns)

    def test_qc_threshold(self, tmp_path, ecoli, bsub):
        out = tmp_path / "qc300.metacool"
        collect([ecoli, bsub], out, qc_distance_threshold=300)
        assert MetaCooler(out).metadata["qc"]["distance_threshold"] == 300

    def test_merge_accepts_two_collections_of_the_same_species(
        self, tmp_path: Path, ecoli: Path, bsub: Path
    ):
        """Replicates collected alike are merge-compatible."""
        from metacooler.io.merge import merge

        rep2 = tmp_path / "rep2"
        rep2.mkdir()
        e2 = make_cool(rep2 / "E_coli.cool", ECOLI, seed=11)
        b2 = make_cool(rep2 / "B_subtilis.cool", BSUB, seed=12)
        one, two = tmp_path / "one.metacool", tmp_path / "two.metacool"
        collect([ecoli, bsub], one)
        collect([e2, b2], two)
        merged = tmp_path / "merged.metacool"
        merge([one, two], merged)
        total = sum(int(cooler.Cooler(str(p)).info["sum"]) for p in (one, two))
        assert int(cooler.Cooler(str(merged)).info["sum"]) == total


def test_output_staging_leaves_no_partial_on_success(tmp_path: Path, ecoli: Path):
    out = tmp_path / "staged.metacool"
    collect([ecoli], out)
    assert out.exists()
    assert not out.with_name(out.name + ".partial").exists()


def test_does_not_modify_inputs(tmp_path: Path, ecoli: Path, bsub: Path):
    before = {p: p.read_bytes() for p in (ecoli, bsub)}
    collect(
        [ecoli, bsub], tmp_path / "out.metacool", contig_ordering="genome-then-length"
    )
    assert {p: p.read_bytes() for p in (ecoli, bsub)} == before
