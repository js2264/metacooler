"""Cooler API and CLI compatibility tests for metacool and mmetacool files.

These tests verify that .metacool and .mmetacool files are fully usable with
cooler's Python API (cooler.Cooler, cooler.fileops) and CLI (cooler info,
cooler ls, cooler dump, etc.) without modification.

Coverage includes mmetacool files created via:
  - conftest tiny_mmetacool fixture (create_mmetacool API)
  - create_mmetacool / metacooler stack CLI
  - add_entry
  - remove_entry / metacooler rm CLI
  - copy_entry (same file) / metacooler cp CLI
  - copy_entry_cross_file / metacooler cp CLI (cross-file)
  - move_entry / metacooler mv CLI
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import cooler
import cooler.fileops
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.create import from_cool, from_pairs
from metacooler.io.mmetacool import (
    add_entry,
    copy_entry,
    copy_entry_cross_file,
    create_mmetacool,
    list_entries,
    move_entry,
    remove_entry,
)
from metacooler.ops.aggregate import aggregate
from metacooler.schema import COOLER_FORMAT, MetaCoolerFormatError, check, validate

# External cooler binary — intentionally subprocess, NOT CliRunner.
# These tests verify third-party cooler CLI compatibility.
_COOLER = str(Path(sys.executable).parent / "cooler")

# Static .metacool used as a source for all module-scoped mutation fixtures.
# Using a static file avoids re-creating a metacool per test.
_SOURCE_MC = Path(__file__).parent / "data" / "assembly_1_rep1_DpnII.metacool"


def _cooler_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([_COOLER, *args], capture_output=True, text=True)


def _assert_entry_cooler_compat(path: Path, entry_name: str) -> None:
    """Assert a single mmetacool entry passes cooler Python API compatibility."""
    uri = f"{path}::entries/{entry_name}"
    assert cooler.fileops.is_cooler(uri) is True
    c = cooler.Cooler(uri)
    assert c.info["format"] == COOLER_FORMAT
    assert int(c.info["format-version"]) == 3
    assert c.info["storage-mode"] == "symmetric-upper"
    assert int(c.info["nbins"]) == len(c.bins()[:])
    assert int(c.info["nchroms"]) == len(c.chroms()[:])
    bins = c.bins()[:]
    assert (bins["end"] > bins["start"]).all()
    pix = c.pixels()[:]
    assert (pix["bin1_id"] <= pix["bin2_id"]).all()
    first_chrom = c.chroms()[:]["name"].iloc[0]
    mat = c.matrix(balance=False).fetch(first_chrom)
    assert isinstance(mat, np.ndarray) and np.array_equal(mat, mat.T)


def _assert_entry_cooler_cli_compat(path: Path, entry_name: str) -> None:
    """Assert cooler CLI works on path::entries/entry_name (3 subprocess calls)."""
    uri = f"{path}::entries/{entry_name}"
    r = _cooler_cli("info", uri)
    assert r.returncode == 0, r.stderr
    assert json.loads(r.stdout)["format"] == COOLER_FORMAT
    r = _cooler_cli("dump", "-t", "bins", uri)
    assert r.returncode == 0 and len(r.stdout.strip().splitlines()) > 0
    r = _cooler_cli("dump", "-t", "pixels", uri)
    assert r.returncode == 0, r.stderr


# ===========================================================================
# .metacool — cooler Python API compatibility
# ===========================================================================


class TestMetacoolCoolerCompat:
    """cooler Python API works correctly on a .metacool file."""

    def test_fileops(self, tiny_metacooler: Path):
        assert cooler.fileops.is_cooler(str(tiny_metacooler)) is True
        assert cooler.fileops.is_cooler(f"{tiny_metacooler}::/") is True
        assert cooler.fileops.list_coolers(str(tiny_metacooler)) == ["/"]
        assert cooler.fileops.is_multires_file(str(tiny_metacooler)) is False

    def test_cooler_init_and_repr(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        assert c is not None
        assert str(tiny_metacooler.name) in repr(c)
        assert cooler.Cooler(f"{tiny_metacooler}::/") is not None

    def test_info_required_fields(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        assert c.info["format"] == COOLER_FORMAT
        assert int(c.info["format-version"]) == 3
        assert c.info["storage-mode"] == "symmetric-upper"
        assert "metacooler-format" in c.info

    def test_info_counts_consistent(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        assert int(c.info["nbins"]) == len(c.bins()[:])
        assert int(c.info["nchroms"]) == len(c.chroms()[:])
        assert int(c.info["nnz"]) == len(c.pixels()[:])

    def test_chroms(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        chroms = c.chroms()[:]
        assert {"name", "length"}.issubset(chroms.columns)
        assert set(chroms["name"]) == {"contig_A", "contig_B", "contig_C", "contig_D"}
        assert (chroms["length"] > 0).all()
        idx = chroms.set_index("name")
        assert int(idx.loc["contig_A", "length"]) == 1000
        assert int(idx.loc["contig_B", "length"]) == 600

    def test_bins(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        bins = c.bins()[:]
        assert {"chrom", "start", "end"}.issubset(bins.columns)
        assert (bins["start"] >= 0).all()
        assert (bins["end"] > bins["start"]).all()
        assert (bins["chrom"] == "contig_A").sum() == 10

    def test_pixels(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        pix = c.pixels()[:]
        assert {"bin1_id", "bin2_id", "count"}.issubset(pix.columns)
        assert (pix["bin1_id"] <= pix["bin2_id"]).all()
        assert (pix["count"] > 0).all()
        assert (pix["bin2_id"] < int(c.info["nbins"])).all()

    def test_matrix(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        mat = c.matrix(balance=False).fetch("contig_A")
        assert mat.ndim == 2 and mat.shape[0] == mat.shape[1]
        n_bins = (c.bins()[:]["chrom"] == "contig_A").sum()
        assert mat.shape == (n_bins, n_bins)
        assert np.array_equal(mat, mat.T)
        assert c.matrix(balance=False).fetch("contig_A:0-500").ndim == 2
        mat_cross = c.matrix(balance=False).fetch("contig_A", "contig_B")
        assert mat_cross.shape[0] > 0 and mat_cross.shape[1] > 0

    def test_extent(self, tiny_metacooler: Path):
        c = cooler.Cooler(str(tiny_metacooler))
        start, end = c.extent("contig_A")
        assert start < end
        assert (end - start) == (c.bins()[:]["chrom"] == "contig_A").sum()
        s2, e2 = c.extent("contig_A:0-300")
        assert s2 < e2

    def test_cooler_cp(self, tiny_metacooler: Path, tmp_path: Path):
        dst = str(tmp_path / "copy.cool")
        cooler.fileops.cp(f"{tiny_metacooler}::/", f"{dst}::/")
        assert cooler.fileops.is_cooler(f"{dst}::/")
        src = cooler.Cooler(str(tiny_metacooler))
        dst_c = cooler.Cooler(f"{dst}::/")
        assert len(dst_c.bins()[:]) == len(src.bins()[:])
        assert len(dst_c.pixels()[:]) == len(src.pixels()[:])


# ===========================================================================
# .metacool — names cooler cannot decode
# ===========================================================================


def _write_pairs_inputs(
    tmp_path: Path,
    contigs: list[str],
    genome_of: dict[str, str],
    contig_extra: dict[str, str] | None = None,
    genome_extra: dict[str, str] | None = None,
) -> dict[str, Path]:
    """Write a complete four-file input set for ``from_pairs``.

    The contig names have to agree across the fasta, the bins BED, the pairs
    and the contig→genome map, so the non-ASCII tests generate all four from
    one list rather than patching a fixture.
    """
    length = 400
    fasta = tmp_path / "t.fa"
    fasta.write_text("".join(f">{c}\n{'ACGT' * 100}\n" for c in contigs))

    header = "contig_name\tgenome_name" + ("\tnote\n" if contig_extra else "\n")
    rows = [
        f"{c}\t{genome_of[c]}" + (f"\t{contig_extra[c]}\n" if contig_extra else "\n")
        for c in contigs
    ]
    contig_map = tmp_path / "map.tsv"
    contig_map.write_text(header + "".join(rows))

    bins_bed = tmp_path / "bins.bed"
    bins_bed.write_text(
        "".join(
            f"{c}\t{s}\t{s + 100}\n" for c in contigs for s in range(0, length, 100)
        )
    )

    rng = np.random.default_rng(0)
    pair_rows = ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"]
    for i in range(50):
        c1, c2 = rng.choice(contigs), rng.choice(contigs)
        p1, p2 = rng.integers(1, length), rng.integers(1, length)
        pair_rows.append(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")
    pairs = tmp_path / "p.pairs"
    pairs.write_text("".join(pair_rows))

    inputs = {
        "fasta": fasta,
        "contig_genome_map": contig_map,
        "bins": bins_bed,
        "pairs": pairs,
    }
    if genome_extra is not None:
        genomes_table = tmp_path / "genomes.tsv"
        genomes_table.write_text(
            "genome_name\tphylum\n"
            + "".join(f"{g}\t{v}\n" for g, v in genome_extra.items())
        )
        inputs["genomes_table"] = genomes_table
    return inputs


def _create_from_pairs(tmp_path: Path, **input_kwargs) -> Path:
    out = tmp_path / "t.metacool"
    inputs = _write_pairs_inputs(tmp_path, **input_kwargs)
    from_pairs(
        uri=out,
        bins=inputs["bins"],
        pairs=inputs["pairs"],
        contig_genome_map=inputs["contig_genome_map"],
        genomes_table=inputs.get("genomes_table"),
        fasta=inputs["fasta"],
        contig_table_format="contigs",
        chunksize=1000,
    )
    return out


class TestNonAsciiNamesRefused:
    """Names cooler cannot decode are refused rather than written.

    cooler reads ``/chroms`` — a soft link to ``/contigs`` — with an ASCII
    ``.astype("U")``, so a single accented byte anywhere in that group yields a
    file that satisfies every structural rule metacooler has and still cannot
    be opened by ``cooler.Cooler()``, and therefore not by ``MetaCooler()``
    either.  That would break the invariant the whole format is built around,
    so it is refused at the point of writing.
    """

    def test_contig_name(self, tmp_path: Path):
        with pytest.raises(MetaCoolerFormatError, match=r"/contigs\.name"):
            _create_from_pairs(
                tmp_path,
                contigs=["contig_é", "contig_B"],
                genome_of={"contig_é": "MAG_1", "contig_B": "MAG_1"},
            )
        assert not (tmp_path / "t.metacool").exists()

    def test_genome_name(self, tmp_path: Path):
        # Invisible to cooler on read — /genomes is not part of the cooler view
        # — but `aggregate --level genome` writes these out as chrom names.
        with pytest.raises(MetaCoolerFormatError, match=r"/genomes\.name"):
            _create_from_pairs(
                tmp_path,
                contigs=["contig_A", "contig_B"],
                genome_of={"contig_A": "MAG_é", "contig_B": "MAG_1"},
            )

    def test_extra_contig_column(self, tmp_path: Path):
        # Not just the name column: cooler's chroms read decodes every column
        # of the group, so an accented annotation is equally fatal.
        with pytest.raises(MetaCoolerFormatError, match=r"/contigs\.note"):
            _create_from_pairs(
                tmp_path,
                contigs=["contig_A", "contig_B"],
                genome_of={"contig_A": "MAG_1", "contig_B": "MAG_2"},
                contig_extra={"contig_A": "café", "contig_B": "plain"},
            )

    def test_error_names_the_offending_value(self, tmp_path: Path):
        with pytest.raises(MetaCoolerFormatError) as excinfo:
            _create_from_pairs(
                tmp_path,
                contigs=["contig_é", "contig_B"],
                genome_of={"contig_é": "MAG_1", "contig_B": "MAG_1"},
            )
        message = str(excinfo.value)
        assert "contig_é" in message
        assert "cooler" in message

    def test_extra_genome_column_is_allowed(self, tmp_path: Path):
        # The deliberate limit of the rule: /genomes is metacooler-only
        # metadata that cooler never reads, so accented taxonomy stays legal.
        out = _create_from_pairs(
            tmp_path,
            contigs=["contig_A", "contig_B"],
            genome_of={"contig_A": "MAG_1", "contig_B": "MAG_2"},
            genome_extra={"MAG_1": "Bacteroidètes", "MAG_2": "Firmicutes"},
        )
        assert cooler.Cooler(str(out)).chroms()[:].shape[0] == 2
        genomes = MetaCooler(out).genomes
        assert set(genomes["phylum"]) == {"Bacteroidètes", "Firmicutes"}
        aggregate(out, tmp_path / "genome.cool", level="genome")
        assert cooler.Cooler(str(tmp_path / "genome.cool")).chromnames == [
            "MAG_1",
            "MAG_2",
        ]

    def test_create_cool_refuses(self, tiny_cool: Path, tmp_path: Path):
        contig_map = tmp_path / "map.tsv"
        contig_map.write_text(
            "contig_name\tgenome_name\n"
            "contig_A\tMAG_é\ncontig_B\tMAG_1\n"
            "contig_C\tMAG_2\ncontig_D\tMAG_2\n"
        )
        with pytest.raises(MetaCoolerFormatError, match=r"/genomes\.name"):
            from_cool(
                cool_uri=tiny_cool,
                output=tmp_path / "out.metacool",
                contig_genome_map=contig_map,
                contig_table_format="contigs",
            )

    def test_cli_reports_a_clean_error(self, tmp_path: Path):
        inputs = _write_pairs_inputs(
            tmp_path,
            contigs=["contig_é", "contig_B"],
            genome_of={"contig_é": "MAG_1", "contig_B": "MAG_1"},
        )
        result = CliRunner().invoke(
            cli,
            [
                "create",
                "pairs",
                "--pairs",
                str(inputs["pairs"]),
                "-o",
                str(tmp_path / "out.metacool"),
                "--fasta",
                str(inputs["fasta"]),
                "--bins",
                str(inputs["bins"]),
                "--contigs",
                str(inputs["contig_genome_map"]),
                "--contigs-format",
                "contigs",
            ],
        )
        assert result.exit_code != 0
        assert "not ASCII" in _flat(result.output)


class TestNonAsciiNamesOnDisk:
    """A file that already holds such names is diagnosed, not vouched for."""

    @staticmethod
    def _break(source: Path, dest: Path, raw: bytes = "contig_é".encode()) -> Path:
        """Copy *source* and overwrite one contig name with *raw*.

        The dataset is widened first: assigning into the original ``S8`` would
        truncate the multi-byte name rather than store it.
        """
        shutil.copy(source, dest)
        with h5py.File(dest, "a") as f:
            names = f["contigs"]["name"][:].astype("S16")
            names[0] = raw
            del f["contigs"]["name"]
            f["contigs"].create_dataset("name", data=names)
        return dest

    def test_validate_refuses(self, tiny_metacooler: Path, tmp_path: Path):
        broken = self._break(tiny_metacooler, tmp_path / "broken.metacool")
        with pytest.raises(MetaCoolerFormatError, match=r"/contigs\.name"):
            validate(broken)
        # And the same named error, rather than cooler's UnicodeDecodeError,
        # when the file is opened.
        with pytest.raises(MetaCoolerFormatError):
            MetaCooler(broken)

    def test_check_reports_it(self, tiny_metacooler: Path, tmp_path: Path):
        broken = self._break(tiny_metacooler, tmp_path / "broken.metacool")
        results = {r.name: r for r in check(broken)}
        assert results["structure"].ok is False
        assert "contig_é" in results["structure"].detail
        # The pre-existing empirical check agrees: cooler really cannot open it.
        assert results["cooler-compat"].ok is False

    def test_invalid_utf8_is_reported_not_raised(
        self, tiny_metacooler: Path, tmp_path: Path
    ):
        # A name that is not even valid UTF-8 must still produce the message:
        # the reporting path decodes with "replace" precisely so that naming
        # the offender cannot itself raise.
        broken = self._break(tiny_metacooler, tmp_path / "raw.metacool", b"contig_\xff")
        with pytest.raises(MetaCoolerFormatError, match=r"/contigs\.name"):
            validate(broken)

    def test_intact_file_still_validates(self, tiny_metacooler: Path):
        validate(tiny_metacooler)
        assert all(r.ok for r in check(tiny_metacooler))


# ===========================================================================
# .metacool — cooler CLI compatibility
# ===========================================================================


class TestMetacoolCoolerCLI:
    """cooler CLI works on a .metacool file."""

    def test_info(self, tiny_metacooler: Path):
        r = _cooler_cli("info", str(tiny_metacooler))
        assert r.returncode == 0
        assert json.loads(r.stdout)["format"] == COOLER_FORMAT

    def test_ls(self, tiny_metacooler: Path):
        r = _cooler_cli("ls", str(tiny_metacooler))
        assert r.returncode == 0
        assert r.stdout.strip().endswith("::/")

    def test_dump_bins(self, tiny_metacooler: Path):
        r = _cooler_cli("dump", "-t", "bins", str(tiny_metacooler))
        assert r.returncode == 0
        c = cooler.Cooler(str(tiny_metacooler))
        assert len(r.stdout.strip().splitlines()) == int(c.info["nbins"])

    def test_dump_chroms(self, tiny_metacooler: Path):
        r = _cooler_cli("dump", "-t", "chroms", str(tiny_metacooler))
        assert r.returncode == 0
        assert len(r.stdout.strip().splitlines()) == 4
        assert "contig_A" in r.stdout

    def test_dump_pixels(self, tiny_metacooler: Path):
        r = _cooler_cli("dump", "-t", "pixels", str(tiny_metacooler))
        assert r.returncode == 0
        assert len(r.stdout.strip().splitlines()) > 0

    def test_tree_and_attrs(self, tiny_metacooler: Path):
        r_tree = _cooler_cli("tree", str(tiny_metacooler))
        assert r_tree.returncode == 0 and "bins" in r_tree.stdout
        r_attrs = _cooler_cli("attrs", str(tiny_metacooler))
        assert r_attrs.returncode == 0 and "format" in r_attrs.stdout


# ===========================================================================
# .mmetacool — cooler Python API and CLI compatibility
# ===========================================================================


class TestMmetacoolCoolerCompat:
    """.mmetacool entries are fully cooler API compatible."""

    def test_fileops(self, tiny_mmetacool: Path):
        assert cooler.fileops.is_cooler(f"{tiny_mmetacool}::entries/entry1") is True
        assert set(cooler.fileops.list_coolers(str(tiny_mmetacool))) == {
            "/entries/entry1",
            "/entries/entry2",
        }
        with pytest.raises(KeyError):
            cooler.fileops.is_cooler(f"{tiny_mmetacool}::entries/nonexistent")

    def test_entries_match_metacooler(
        self, tiny_metacooler: Path, tiny_mmetacool: Path
    ):
        mc = cooler.Cooler(str(tiny_metacooler))
        entry = cooler.Cooler(f"{tiny_mmetacool}::entries/entry1")
        assert entry.info["format"] == COOLER_FORMAT
        assert int(entry.info["nchroms"]) == int(mc.info["nchroms"])
        assert len(entry.bins()[:]) == len(mc.bins()[:])
        assert len(entry.pixels()[:]) == len(mc.pixels()[:])

    def test_pixels_upper_triangular(self, tiny_mmetacool: Path):
        pix = cooler.Cooler(f"{tiny_mmetacool}::entries/entry1").pixels()[:]
        assert (pix["bin1_id"] <= pix["bin2_id"]).all()

    def test_matrix_and_extent_match_metacooler(
        self, tiny_metacooler: Path, tiny_mmetacool: Path
    ):
        mc = cooler.Cooler(str(tiny_metacooler))
        entry = cooler.Cooler(f"{tiny_mmetacool}::entries/entry1")
        assert np.array_equal(
            mc.matrix(balance=False).fetch("contig_A"),
            entry.matrix(balance=False).fetch("contig_A"),
        )
        assert mc.extent("contig_A") == entry.extent("contig_A")

    def test_data_identical_to_metacooler(
        self, tiny_metacooler: Path, tiny_mmetacool: Path
    ):
        mc = cooler.Cooler(str(tiny_metacooler))
        entry = cooler.Cooler(f"{tiny_mmetacool}::entries/entry1")
        pd.testing.assert_frame_equal(
            mc.chroms()[:][["name", "length"]].reset_index(drop=True),
            entry.chroms()[:][["name", "length"]].reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            mc.bins()[:][["chrom", "start", "end"]].reset_index(drop=True),
            entry.bins()[:][["chrom", "start", "end"]].reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            mc.pixels()[:].reset_index(drop=True),
            entry.pixels()[:].reset_index(drop=True),
        )


class TestMmetacoolCoolerCLI:
    """cooler CLI works on .mmetacool entries."""

    def test_ls_finds_both_entries(self, tiny_mmetacool: Path):
        r = _cooler_cli("ls", str(tiny_mmetacool))
        assert r.returncode == 0
        assert "entries/entry1" in r.stdout
        assert "entries/entry2" in r.stdout

    def test_info(self, tiny_mmetacool: Path):
        r = _cooler_cli("info", f"{tiny_mmetacool}::entries/entry1")
        assert r.returncode == 0
        assert json.loads(r.stdout)["format"] == COOLER_FORMAT

    def test_dump_bins_count_matches_api(self, tiny_mmetacool: Path):
        uri = f"{tiny_mmetacool}::entries/entry1"
        r = _cooler_cli("dump", "-t", "bins", uri)
        assert r.returncode == 0
        assert len(r.stdout.strip().splitlines()) == int(
            cooler.Cooler(uri).info["nbins"]
        )

    def test_dump_chroms_and_pixels(self, tiny_mmetacool: Path):
        uri = f"{tiny_mmetacool}::entries/entry1"
        r = _cooler_cli("dump", "-t", "chroms", uri)
        assert r.returncode == 0 and "contig_A" in r.stdout
        r = _cooler_cli("dump", "-t", "pixels", uri)
        assert r.returncode == 0 and len(r.stdout.strip().splitlines()) > 0

    def test_tree_and_attrs(self, tiny_mmetacool: Path):
        uri = f"{tiny_mmetacool}::entries/entry1"
        assert _cooler_cli("tree", uri).returncode == 0
        r = _cooler_cli("attrs", uri)
        assert r.returncode == 0 and "format" in r.stdout


# ===========================================================================
# Module-scoped mutation fixtures
# All use _SOURCE_MC directly — no per-test shutil.copy overhead.
# scope="module" ensures each fixture is created once for the entire module.
# ===========================================================================


@pytest.fixture(scope="module")
def stacked_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool built with create_mmetacool (stack API), 2 entries."""
    out = tmp_path_factory.mktemp("stacked") / "stacked.mmetacool"
    create_mmetacool(out, {"sampleA": _SOURCE_MC, "sampleB": _SOURCE_MC})
    return out


@pytest.fixture(scope="module")
def stacked_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool built with `metacooler stack` CLI, 2 entries."""
    out = tmp_path_factory.mktemp("stacked_cli") / "stacked_cli.mmetacool"
    runner = CliRunner()
    r = runner.invoke(
        cli,
        [
            "stack",
            "-o",
            str(out),
            "--names",
            "cli_a,cli_b",
            str(_SOURCE_MC),
            str(_SOURCE_MC),
        ],
    )
    assert r.exit_code == 0, r.output
    return out


@pytest.fixture(scope="module")
def added_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry added after creation."""
    out = tmp_path_factory.mktemp("added") / "added.mmetacool"
    create_mmetacool(out, {"original": _SOURCE_MC})
    add_entry(out, "added_later", _SOURCE_MC)
    return out


@pytest.fixture(scope="module")
def rmed_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with one entry removed via remove_entry."""
    out = tmp_path_factory.mktemp("rmed") / "rmed.mmetacool"
    create_mmetacool(out, {"keep_a": _SOURCE_MC, "remove_me": _SOURCE_MC})
    remove_entry(out, "remove_me")
    return out


@pytest.fixture(scope="module")
def rmed_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry removed via `metacooler rm` CLI."""
    out = tmp_path_factory.mktemp("rmed_cli") / "rmed_cli.mmetacool"
    create_mmetacool(out, {"keep": _SOURCE_MC, "to_remove": _SOURCE_MC})
    runner = CliRunner()
    r = runner.invoke(cli, ["rm", "-f", f"{out}::entries/to_remove"])
    assert r.exit_code == 0, r.output
    return out


@pytest.fixture(scope="module")
def cpd_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry copied within the same file."""
    out = tmp_path_factory.mktemp("cpd") / "cpd.mmetacool"
    create_mmetacool(out, {"original": _SOURCE_MC})
    copy_entry(out, "original", "copy_of_original")
    return out


@pytest.fixture(scope="module")
def cpd_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry copied via `metacooler cp` CLI (same file)."""
    out = tmp_path_factory.mktemp("cpd_cli") / "cpd_cli.mmetacool"
    create_mmetacool(out, {"src_entry": _SOURCE_MC})
    runner = CliRunner()
    r = runner.invoke(
        cli, ["cp", f"{out}::entries/src_entry", f"{out}::entries/dst_entry"]
    )
    assert r.exit_code == 0, r.output
    return out


@pytest.fixture(scope="module")
def cpd_cross_file_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool receiving an entry from another file via copy_entry_cross_file."""
    d = tmp_path_factory.mktemp("cross_cp")
    src_mm = d / "src.mmetacool"
    dst_mm = d / "dst.mmetacool"
    create_mmetacool(src_mm, {"source_entry": _SOURCE_MC})
    create_mmetacool(dst_mm, {"local_entry": _SOURCE_MC})
    copy_entry_cross_file(src_mm, "source_entry", dst_mm, "cross_copied")
    return dst_mm


@pytest.fixture(scope="module")
def cpd_cross_file_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool receiving a cross-file entry via `metacooler cp` CLI."""
    d = tmp_path_factory.mktemp("cross_cp_cli")
    src_mm = d / "src_cli.mmetacool"
    dst_mm = d / "dst_cli.mmetacool"
    create_mmetacool(src_mm, {"entry_to_copy": _SOURCE_MC})
    create_mmetacool(dst_mm, {"local": _SOURCE_MC})
    runner = CliRunner()
    r = runner.invoke(
        cli,
        [
            "cp",
            f"{src_mm}::entries/entry_to_copy",
            f"{dst_mm}::entries/cli_cross_copy",
        ],
    )
    assert r.exit_code == 0, r.output
    return dst_mm


@pytest.fixture(scope="module")
def mvd_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry renamed via move_entry."""
    out = tmp_path_factory.mktemp("mvd") / "mvd.mmetacool"
    create_mmetacool(out, {"old_name": _SOURCE_MC, "sibling": _SOURCE_MC})
    move_entry(out, "old_name", "new_name")
    return out


@pytest.fixture(scope="module")
def mvd_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool with an entry renamed via `metacooler mv` CLI."""
    out = tmp_path_factory.mktemp("mvd_cli") / "mvd_cli.mmetacool"
    create_mmetacool(out, {"before_rename": _SOURCE_MC, "untouched": _SOURCE_MC})
    runner = CliRunner()
    r = runner.invoke(
        cli,
        ["mv", "-f", f"{out}::entries/before_rename", f"{out}::entries/after_rename"],
    )
    assert r.exit_code == 0, r.output
    return out


@pytest.fixture(scope="module")
def compound_mmetacool(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool from a stack→add→rm→cp→mv sequence (API)."""
    out = tmp_path_factory.mktemp("compound") / "compound.mmetacool"
    create_mmetacool(
        out, {"alpha": _SOURCE_MC, "beta": _SOURCE_MC, "gamma": _SOURCE_MC}
    )
    add_entry(out, "delta", _SOURCE_MC)
    remove_entry(out, "gamma")
    copy_entry(out, "alpha", "alpha_copy")
    move_entry(out, "beta", "beta_renamed")
    return out


@pytest.fixture(scope="module")
def compound_mmetacool_cli(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """mmetacool from a stack→rm→cp→mv sequence (CLI)."""
    out = tmp_path_factory.mktemp("compound_cli") / "compound_cli.mmetacool"
    runner = CliRunner()
    r = runner.invoke(
        cli,
        [
            "stack",
            "-o",
            str(out),
            "--names",
            "alpha,beta,gamma,delta",
            str(_SOURCE_MC),
            str(_SOURCE_MC),
            str(_SOURCE_MC),
            str(_SOURCE_MC),
        ],
    )
    assert r.exit_code == 0, r.output
    r = runner.invoke(cli, ["rm", "-f", f"{out}::entries/gamma"])
    assert r.exit_code == 0, r.output
    r = runner.invoke(
        cli, ["cp", f"{out}::entries/alpha", f"{out}::entries/alpha_copy"]
    )
    assert r.exit_code == 0, r.output
    r = runner.invoke(
        cli, ["mv", "-f", f"{out}::entries/beta", f"{out}::entries/beta_renamed"]
    )
    assert r.exit_code == 0, r.output
    return out


# ===========================================================================
# stack
# ===========================================================================


class TestStackedMmetacool:
    """cooler compat for mmetacool created via create_mmetacool (stack API)."""

    def test_list_coolers(self, stacked_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(stacked_mmetacool))
        assert set(paths) == {"/entries/sampleA", "/entries/sampleB"}

    def test_sampleA_api_compat(self, stacked_mmetacool: Path):
        _assert_entry_cooler_compat(stacked_mmetacool, "sampleA")

    def test_sampleA_cli_compat(self, stacked_mmetacool: Path):
        _assert_entry_cooler_cli_compat(stacked_mmetacool, "sampleA")


class TestStackCLI:
    """cooler compat for mmetacool created via `metacooler stack`."""

    def test_list_coolers(self, stacked_mmetacool_cli: Path):
        paths = cooler.fileops.list_coolers(str(stacked_mmetacool_cli))
        assert set(paths) == {"/entries/cli_a", "/entries/cli_b"}

    def test_cli_a_compat(self, stacked_mmetacool_cli: Path):
        _assert_entry_cooler_compat(stacked_mmetacool_cli, "cli_a")
        _assert_entry_cooler_cli_compat(stacked_mmetacool_cli, "cli_a")


# ===========================================================================
# add_entry
# ===========================================================================


class TestAddedMmetacool:
    """cooler compat after add_entry."""

    def test_list_coolers(self, added_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(added_mmetacool))
        assert set(paths) == {"/entries/original", "/entries/added_later"}

    def test_added_entry_api_compat(self, added_mmetacool: Path):
        _assert_entry_cooler_compat(added_mmetacool, "added_later")

    def test_added_entry_cli_compat(self, added_mmetacool: Path):
        _assert_entry_cooler_cli_compat(added_mmetacool, "added_later")


# ===========================================================================
# remove_entry / metacooler rm
# ===========================================================================


class TestRmedMmetacool:
    """cooler compat after remove_entry."""

    def test_entries_and_absence(self, rmed_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(rmed_mmetacool))
        assert "/entries/keep_a" in paths
        assert "remove_me" not in list_entries(rmed_mmetacool)

    def test_keep_a_api_compat(self, rmed_mmetacool: Path):
        _assert_entry_cooler_compat(rmed_mmetacool, "keep_a")

    def test_keep_a_cli_compat(self, rmed_mmetacool: Path):
        _assert_entry_cooler_cli_compat(rmed_mmetacool, "keep_a")


class TestRmCLI:
    """cooler compat after `metacooler rm`."""

    def test_list_coolers(self, rmed_mmetacool_cli: Path):
        assert cooler.fileops.list_coolers(str(rmed_mmetacool_cli)) == ["/entries/keep"]

    def test_keep_compat(self, rmed_mmetacool_cli: Path):
        _assert_entry_cooler_compat(rmed_mmetacool_cli, "keep")
        _assert_entry_cooler_cli_compat(rmed_mmetacool_cli, "keep")


# ===========================================================================
# copy_entry (same file) / metacooler cp
# ===========================================================================


class TestCpdMmetacool:
    """cooler compat after copy_entry (same file)."""

    def test_list_coolers(self, cpd_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(cpd_mmetacool))
        assert set(paths) == {"/entries/original", "/entries/copy_of_original"}

    def test_copy_api_compat(self, cpd_mmetacool: Path):
        _assert_entry_cooler_compat(cpd_mmetacool, "copy_of_original")

    def test_copy_data_matches_original(self, cpd_mmetacool: Path):
        orig = cooler.Cooler(f"{cpd_mmetacool}::entries/original")
        copy = cooler.Cooler(f"{cpd_mmetacool}::entries/copy_of_original")
        pd.testing.assert_frame_equal(
            orig.pixels()[:].reset_index(drop=True),
            copy.pixels()[:].reset_index(drop=True),
        )


class TestCpCLI:
    """cooler compat after `metacooler cp` (same file)."""

    def test_list_coolers(self, cpd_mmetacool_cli: Path):
        paths = cooler.fileops.list_coolers(str(cpd_mmetacool_cli))
        assert set(paths) == {"/entries/src_entry", "/entries/dst_entry"}

    def test_dst_entry_compat(self, cpd_mmetacool_cli: Path):
        _assert_entry_cooler_compat(cpd_mmetacool_cli, "dst_entry")
        _assert_entry_cooler_cli_compat(cpd_mmetacool_cli, "dst_entry")


# ===========================================================================
# copy_entry_cross_file / metacooler cp (cross-file)
# ===========================================================================


class TestCrossFileCpMmetacool:
    """cooler compat after copy_entry_cross_file."""

    def test_list_coolers(self, cpd_cross_file_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(cpd_cross_file_mmetacool))
        assert set(paths) == {"/entries/local_entry", "/entries/cross_copied"}

    def test_cross_copied_api_compat(self, cpd_cross_file_mmetacool: Path):
        _assert_entry_cooler_compat(cpd_cross_file_mmetacool, "cross_copied")

    def test_cross_copied_data_matches_local(self, cpd_cross_file_mmetacool: Path):
        local = cooler.Cooler(f"{cpd_cross_file_mmetacool}::entries/local_entry")
        cross = cooler.Cooler(f"{cpd_cross_file_mmetacool}::entries/cross_copied")
        assert int(local.info["nbins"]) == int(cross.info["nbins"])
        assert int(local.info["nnz"]) == int(cross.info["nnz"])


class TestCrossFileCpCLI:
    """cooler compat after `metacooler cp` (cross-file)."""

    def test_list_coolers(self, cpd_cross_file_mmetacool_cli: Path):
        paths = cooler.fileops.list_coolers(str(cpd_cross_file_mmetacool_cli))
        assert set(paths) == {"/entries/local", "/entries/cli_cross_copy"}

    def test_cli_cross_copy_compat(self, cpd_cross_file_mmetacool_cli: Path):
        _assert_entry_cooler_compat(cpd_cross_file_mmetacool_cli, "cli_cross_copy")
        _assert_entry_cooler_cli_compat(cpd_cross_file_mmetacool_cli, "cli_cross_copy")


# ===========================================================================
# move_entry / metacooler mv
# ===========================================================================


class TestMvdMmetacool:
    """cooler compat after move_entry."""

    def test_entries_and_absence(self, mvd_mmetacool: Path):
        paths = cooler.fileops.list_coolers(str(mvd_mmetacool))
        assert "/entries/new_name" in paths
        assert "old_name" not in list_entries(mvd_mmetacool)

    def test_new_name_api_compat(self, mvd_mmetacool: Path):
        _assert_entry_cooler_compat(mvd_mmetacool, "new_name")

    def test_moved_entry_data_intact(self, mvd_mmetacool: Path):
        moved = cooler.Cooler(f"{mvd_mmetacool}::entries/new_name")
        sibling = cooler.Cooler(f"{mvd_mmetacool}::entries/sibling")
        assert int(moved.info["nbins"]) == int(sibling.info["nbins"])
        assert int(moved.info["nnz"]) == int(sibling.info["nnz"])

    def test_new_name_cli_compat(self, mvd_mmetacool: Path):
        _assert_entry_cooler_cli_compat(mvd_mmetacool, "new_name")


class TestMvCLI:
    """cooler compat after `metacooler mv`."""

    def test_list_coolers(self, mvd_mmetacool_cli: Path):
        paths = cooler.fileops.list_coolers(str(mvd_mmetacool_cli))
        assert set(paths) == {"/entries/after_rename", "/entries/untouched"}

    def test_after_rename_compat(self, mvd_mmetacool_cli: Path):
        _assert_entry_cooler_compat(mvd_mmetacool_cli, "after_rename")
        _assert_entry_cooler_cli_compat(mvd_mmetacool_cli, "after_rename")


# ===========================================================================
# Compound scenario: stack → add → rm → cp → mv
# ===========================================================================

_COMPOUND_EXPECTED = {"alpha", "alpha_copy", "beta_renamed", "delta"}
# Two representative entries are sufficient — all entries are copies of the
# same source and _assert_entry_cooler_compat is already comprehensive.
_COMPOUND_CHECK = ("alpha", "beta_renamed")


class TestCompoundMmetacool:
    """cooler compat for compound API scenario (stack→add→rm→cp→mv)."""

    def test_entries(self, compound_mmetacool: Path):
        entries = set(list_entries(compound_mmetacool))
        assert entries == _COMPOUND_EXPECTED
        assert "gamma" not in entries
        assert "beta" not in entries

    @pytest.mark.parametrize("entry", _COMPOUND_CHECK)
    def test_entry_api_compat(self, compound_mmetacool: Path, entry: str):
        _assert_entry_cooler_compat(compound_mmetacool, entry)

    def test_alpha_copy_data_matches_alpha(self, compound_mmetacool: Path):
        orig = cooler.Cooler(f"{compound_mmetacool}::entries/alpha")
        copy = cooler.Cooler(f"{compound_mmetacool}::entries/alpha_copy")
        pd.testing.assert_frame_equal(
            orig.pixels()[:].reset_index(drop=True),
            copy.pixels()[:].reset_index(drop=True),
        )

    @pytest.mark.parametrize("entry", _COMPOUND_CHECK)
    def test_entry_cli_compat(self, compound_mmetacool: Path, entry: str):
        _assert_entry_cooler_cli_compat(compound_mmetacool, entry)


class TestCompoundCLI:
    """cooler compat for compound CLI scenario (stack→rm→cp→mv)."""

    def test_entries(self, compound_mmetacool_cli: Path):
        entries = set(list_entries(compound_mmetacool_cli))
        assert entries == _COMPOUND_EXPECTED
        assert "gamma" not in entries

    @pytest.mark.parametrize("entry", _COMPOUND_CHECK)
    def test_entry_api_compat(self, compound_mmetacool_cli: Path, entry: str):
        _assert_entry_cooler_compat(compound_mmetacool_cli, entry)

    @pytest.mark.parametrize("entry", _COMPOUND_CHECK)
    def test_entry_cli_compat(self, compound_mmetacool_cli: Path, entry: str):
        _assert_entry_cooler_cli_compat(compound_mmetacool_cli, entry)
