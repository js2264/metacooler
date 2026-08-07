"""Comprehensive CLI tests for metacooler, including failure/stress tests.

Covers every CLI command (create, dump, info, tree, attrs, aggregate,
extract, show, merge, metadata) with:
  - happy-path scenarios beyond those in test_cli.py
  - failure / error-path scenarios for every command
  - verbose-flag smoke tests
  - edge-case / stress scenarios

Uses Click's CliRunner so that tests work regardless of whether ``metacooler``
is present on PATH (e.g. when running ``pytest`` directly inside a venv).
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import pytest
from _cli_output import flat as _flat
from _cli_output import squash as _squash
from click.testing import CliRunner

from metacooler.cli.main import cli
from metacooler.schema import METACOOLER_FORMAT, METACOOLER_FORMAT_VERSION

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data"

_runner = CliRunner()


def _run(args: list[str]):
    """Invoke the metacooler CLI and return a Click Result.

    Relevant attributes:
      .exit_code   int
      .output      str  (combined stdout + stderr)
    """
    return _runner.invoke(cli, args, catch_exceptions=True)


def _user_metadata(output: str) -> dict:
    """Parse ``metadata list`` output, minus the blocks metacooler writes itself.

    ``create`` stores a contact-quality block under ``qc`` (see
    ``metacooler qc``), so a freshly created file's metadata is no longer
    empty.  Tests about *user* metadata compare against this.
    """
    from metacooler.qc import QC_METADATA_KEY

    data = json.loads(output)
    data.pop(QC_METADATA_KEY, None)
    return data


def _make_corrupt_hdf5(path: Path) -> None:
    """Write a plain HDF5 file that is NOT a valid metacool file."""
    with h5py.File(path, "w") as f:
        f.attrs["format"] = "HDF5::CoolerV2"
        f.create_dataset("dummy", data=[1, 2, 3])


def _make_text_file(path: Path, content: str = "not hdf5\n") -> None:
    path.write_text(content)


# ===========================================================================
# Global CLI
# ===========================================================================


class TestGlobalCLI:
    def test_help(self):
        r = _run(["--help"])
        assert r.exit_code == 0
        assert "MetaCooler" in r.output or "metacooler" in r.output.lower()

    def test_version(self):
        r = _run(["--version"])
        assert r.exit_code == 0
        assert "metacooler" in r.output.lower() or "0." in r.output

    def test_no_args_shows_help(self):
        r = _run([])
        # rich-click shows help but exits with 0 or 2 depending on version
        assert "metacooler" in r.output.lower() or "Usage" in r.output

    def test_unknown_command_fails(self):
        r = _run(["doesnotexist"])
        assert r.exit_code != 0


# ===========================================================================
# create
# ===========================================================================


class TestCreateCLI:
    # --- help ---

    def test_create_help(self):
        r = _run(["create", "--help"])
        assert r.exit_code == 0
        # Check for subcommands
        assert "cool" in r.output
        assert "pairs" in r.output

        # Check pairs subcommand has the expected options
        r_pairs = _run(["create", "pairs", "--help"])
        assert r_pairs.exit_code == 0
        for kw in ["--output", "--bins", "--pairs", "--contigs", "--fasta"]:
            assert kw in r_pairs.output

    # --- missing required options ---

    def test_create_pairs_missing_output_fails(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        r = _run(
            [
                "create",
                "pairs",
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
            ]
        )
        assert r.exit_code != 0

    def test_create_pairs_missing_bins_fails(
        self, tiny_contigs_table, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
            ]
        )
        assert r.exit_code != 0

    def test_create_pairs_missing_pairs_fails(
        self, tiny_contigs_table, tiny_bins_bed, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--contigs",
                str(tiny_contigs_table),
            ]
        )
        assert r.exit_code != 0

    def test_create_pairs_missing_contigs_fails(
        self, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
            ]
        )
        assert r.exit_code != 0

    # --- non-existent input paths ---

    def test_create_pairs_nonexistent_pairs_fails(
        self, tiny_contigs_table, tiny_bins_bed, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tmp_path / "ghost.pairs"),
                "--contigs",
                str(tiny_contigs_table),
            ]
        )
        assert r.exit_code != 0

    def test_create_pairs_nonexistent_contigs_fails(
        self, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tmp_path / "ghost.tsv"),
            ]
        )
        assert r.exit_code != 0

    def test_create_nonexistent_bins_file_fails(
        self, tiny_contigs_table, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tmp_path / "ghost.bed"),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
            ]
        )
        assert r.exit_code != 0

    # --- invalid choice values ---

    def test_create_invalid_contigs_format_fails(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--contigs-format",
                "BADFORMAT",
            ]
        )
        assert r.exit_code != 0

    def test_create_invalid_contig_order_fails(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "out.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--contig-order",
                "BADORDER",
            ]
        )
        assert r.exit_code != 0

    # --- additional happy-path options ---

    def test_create_with_assembly_name(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "named.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--assembly-name",
                "my_assembly",
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "genome-assembly" in f.attrs
            assert f.attrs["genome-assembly"] == "my_assembly"

    def test_create_with_genomes_table(
        self,
        tiny_contigs_table,
        tiny_bins_bed,
        tiny_pairs,
        tiny_genomes_table,
        tmp_path,
    ):
        out = tmp_path / "with_genomes.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--genomes",
                str(tiny_genomes_table),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_create_with_metadata_json(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        meta_path = tmp_path / "meta.json"
        meta_path.write_text(json.dumps({"project": "test_project", "version": 1.2}))
        out = tmp_path / "with_metadata.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--metadata",
                str(meta_path),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            raw = f.attrs.get("metadata", "{}")
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            stored = json.loads(raw) if isinstance(raw, str) else dict(raw)
            assert stored["project"] == "test_project"
            assert stored["version"] == "1.2"
            assert "commands" in stored

    def test_create_with_metadata_yaml(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        meta_path = tmp_path / "meta.yaml"
        meta_path.write_text("project: yaml_project\nversion: 2\n")
        out = tmp_path / "with_metadata_yaml.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--metadata",
                str(meta_path),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        with h5py.File(out, "r") as f:
            raw = f.attrs.get("metadata", "{}")
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            stored = json.loads(raw) if isinstance(raw, str) else dict(raw)
            assert stored["project"] == "yaml_project"
            assert stored["version"] == "2"
            assert "commands" in stored

    def test_create_with_metadata_invalid_fails(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        meta_path = tmp_path / "meta.json"
        meta_path.write_text("[1, 2, 3]")
        out = tmp_path / "should_fail.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--metadata",
                str(meta_path),
            ]
        )
        assert r.exit_code != 0
        assert "must contain a JSON/YAML object" in r.output

    def test_create_drop_unassigned(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "no_unbinned.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--drop-unassigned",
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        dump_r = _run(["dump", "-t", "genomes", str(out)])
        assert "__unbinned__" not in dump_r.output

    def test_create_verbose(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "verbose.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "-v",
            ]
        )
        assert r.exit_code == 0, r.output

    def test_create_small_chunksize(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        """Stress: process pairs in tiny chunks of 10."""
        out = tmp_path / "smallchunk.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--chunksize",
                "10",
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_create_binsize_as_integer(
        self, tiny_fasta, tiny_contigs_table, tiny_pairs, tmp_path
    ):
        """Bins given as integer (bp) – fasta required for contig lengths."""
        out = tmp_path / "binsize.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                "100",
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--fasta",
                str(tiny_fasta),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_create_contig_order_alpha(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "alpha.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--contig-order",
                "alpha",
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_create_contig_order_fasta(
        self, tiny_fasta, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "fasta_order.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--contig-order",
                "fasta",
                "--fasta",
                str(tiny_fasta),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_create_contig_order_mag_then_length(
        self, tiny_contigs_table, tiny_bins_bed, tiny_pairs, tmp_path
    ):
        out = tmp_path / "mbtl.metacool"
        r = _run(
            [
                "create",
                "pairs",
                "-o",
                str(out),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--contig-order",
                "genome-then-length",
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()


# ===========================================================================
# dump
# ===========================================================================


class TestDumpCLI:
    def test_dump_help(self):
        r = _run(["dump", "--help"])
        assert r.exit_code == 0
        assert "--table" in r.output

    # --- `dump … | head` is a normal end of pipeline, not a failure ---

    def test_dump_into_a_closed_pipe_reports_nothing(self):
        """``dump … | head`` used to print ``Error: [Errno 32] Broken pipe``.

        Needs a real pipe (CliRunner writes into a buffer, which never breaks)
        and an output larger than the ~64 KiB pipe buffer, or the writer
        finishes before ``head`` closes its end and no EPIPE is ever raised.
        The shipped fixture dumps ~340 kB of pixels.
        """
        import subprocess
        import sys

        source = DATA_DIR / "assembly_1_1kb.metacool"
        dump = subprocess.Popen(
            [
                sys.executable,
                "-c",
                "from metacooler.cli.main import cli; cli()",
                "dump",
                "-t",
                "pixels",
                str(source),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert dump.stdout is not None
        head = subprocess.Popen(
            ["head", "-n", "5"], stdin=dump.stdout, stdout=subprocess.PIPE
        )
        dump.stdout.close()  # only head holds the read end now
        rows = head.communicate()[0]
        stderr = dump.stderr.read() if dump.stderr else b""
        dump.wait()

        assert len(rows.splitlines()) == 5
        assert b"Broken pipe" not in stderr
        assert b"Traceback" not in stderr
        assert b"Exception ignored" not in stderr
        # 128 + SIGPIPE, i.e. what a C tool in the same position reports; a
        # bare pipeline still exits on head's status, which is 0.
        assert dump.returncode == 141

    # --- failure: non-existent input ---

    def test_dump_nonexistent_file_fails(self, tmp_path):
        r = _run(["dump", str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0

    # --- failure: invalid table name ---

    def test_dump_invalid_table_fails(self, tiny_metacooler):
        r = _run(["dump", "-t", "NOTATABLE", str(tiny_metacooler)])
        assert r.exit_code != 0

    # --- dump indexes table ---

    def test_dump_indexes(self, tiny_metacooler):
        r = _run(["dump", "-t", "indexes", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert "bin1_offset" in r.output or "position" in r.output

    # --- dump with filters ---

    def test_dump_pixels_with_genome_filter(self, tiny_metacooler):
        r = _run(["dump", "-t", "pixels", "--genomes", "MAG_1", str(tiny_metacooler)])
        assert r.exit_code == 0

    def test_dump_bins_with_genome_filter(self, tiny_metacooler):
        r = _run(["dump", "-t", "bins", "--genomes", "MAG_1", str(tiny_metacooler)])
        assert r.exit_code == 0
        lines = [ls for ls in r.output.strip().split("\n") if ls]
        assert len(lines) >= 1  # at least one bin row (no header by default)

    def test_dump_bins_with_contig_filter(self, tiny_metacooler):
        r = _run(["dump", "-t", "bins", "--contigs", "contig_A", str(tiny_metacooler)])
        assert r.exit_code == 0
        lines = [ls for ls in r.output.strip().split("\n") if ls]
        # 1000 bp / 100 bp per bin = 10 bins for contig_A (header may or may not be present)
        assert len(lines) >= 10
        assert all("contig_A" in ls for ls in lines if not ls.startswith("chrom"))

    def test_dump_contigs_with_genome_filter(self, tiny_metacooler):
        r = _run(["dump", "-t", "contigs", "--genomes", "MAG_1", str(tiny_metacooler)])
        assert r.exit_code == 0
        lines = [ls for ls in r.output.strip().split("\n") if ls]
        # contig_A + contig_B belong to MAG_1 (header may or may not be included)
        contig_names = {ls.split("\t")[0] for ls in lines if not ls.startswith("name")}
        assert "contig_A" in contig_names
        assert "contig_B" in contig_names

    # --- output to file ---

    def test_dump_to_file(self, tiny_metacooler, tmp_path):
        out = tmp_path / "chroms.tsv"
        r = _run(["dump", "-t", "chroms", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0
        assert out.exists()
        content = out.read_text()
        assert "contig_A" in content
        assert r.output.strip() == ""

    # --- verbose ---

    def test_dump_verbose(self, tiny_metacooler):
        r = _run(["dump", "-t", "chroms", "-v", str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- --header adds exactly one extra line ---

    def test_dump_header_adds_one_line(self, tiny_metacooler):
        # Default: no header (header=False)
        no_hdr = _run(["dump", "-t", "genomes", str(tiny_metacooler)])
        # With --header flag: header row included
        with_hdr = _run(["dump", "-t", "genomes", "--header", str(tiny_metacooler)])
        assert no_hdr.exit_code == 0
        assert with_hdr.exit_code == 0
        n_without = len(no_hdr.output.strip().split("\n"))
        n_with = len(with_hdr.output.strip().split("\n"))
        assert n_with == n_without + 1

    # --- all table types succeed ---

    @pytest.mark.parametrize(
        "table", ["chroms", "bins", "pixels", "contigs", "genomes", "indexes"]
    )
    def test_dump_all_table_types(self, tiny_metacooler, table):
        r = _run(["dump", "-t", table, str(tiny_metacooler)])
        assert r.exit_code == 0, f"dump -t {table} failed: {r.output}"


# ===========================================================================
# info
# ===========================================================================


class TestInfoCLI:
    def test_info_help(self):
        r = _run(["info", "--help"])
        assert r.exit_code == 0

    def test_info_all_keys_present(self, tiny_metacooler):
        r = _run(["info", str(tiny_metacooler)])
        assert r.exit_code == 0
        data = json.loads(r.output)
        for key in [
            "format",
            "format-version",
            "ngenomes",
            "ncontigs",
            "nbins",
            "nnz",
            "ordering",
        ]:
            assert key in data, f"Missing key: {key}"

    def test_info_ordering_correct(self, tiny_metacooler):
        r = _run(["info", str(tiny_metacooler)])
        data = json.loads(r.output)
        assert data["ordering"] == "genome-then-alpha"

    def test_info_ngenomes_includes_unbinned(self, tiny_metacooler):
        r = _run(["info", str(tiny_metacooler)])
        data = json.loads(r.output)
        assert data["ngenomes"] == 3  # MAG_1, MAG_2, __unbinned__

    def test_info_format_is_metacooler(self, tiny_metacooler):
        r = _run(["info", str(tiny_metacooler)])
        data = json.loads(r.output)
        assert data["format"] == METACOOLER_FORMAT

    def test_info_format_version_is_1(self, tiny_metacooler):
        r = _run(["info", str(tiny_metacooler)])
        data = json.loads(r.output)
        assert str(data["format-version"]) == str(METACOOLER_FORMAT_VERSION)

    def test_info_table_output(self, tiny_metacooler):
        r = _run(["info", "--table", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert "Genomes" in r.output or "Format" in r.output

    def test_info_verbose(self, tiny_metacooler):
        r = _run(["info", "-v", str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- failures ---

    def test_info_nonexistent_file_fails(self, tmp_path):
        r = _run(["info", str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0

    def test_info_non_metacooler_hdf5_fails(self, tmp_path):
        bad = tmp_path / "bad.h5"
        _make_corrupt_hdf5(bad)
        r = _run(["info", str(bad)])
        assert r.exit_code != 0

    def test_info_text_file_fails(self, tmp_path):
        bad = tmp_path / "notbin.txt"
        _make_text_file(bad)
        r = _run(["info", str(bad)])
        assert r.exit_code != 0


# ===========================================================================
# tree
# ===========================================================================


class TestTreeCLI:
    def test_tree_help(self):
        r = _run(["tree", "--help"])
        assert r.exit_code == 0

    def test_tree_contains_expected_groups(self, tiny_metacooler):
        r = _run(["tree", str(tiny_metacooler)])
        assert r.exit_code == 0
        for grp in ["genomes/", "contigs/", "bins/", "pixels/", "indexes/"]:
            assert grp in r.output

    def test_tree_hides_chroms(self, tiny_metacooler):
        r = _run(["tree", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert "chroms/" not in r.output

    def test_tree_hides_attributes(self, tiny_metacooler):
        r = _run(["tree", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert "@format" not in r.output
        assert "@metacooler-format" not in r.output

    def test_tree_verbose(self, tiny_metacooler):
        r = _run(["tree", "-v", str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- failures ---

    def test_tree_nonexistent_file_fails(self, tmp_path):
        r = _run(["tree", str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0

    def test_tree_text_file_fails(self, tmp_path):
        bad = tmp_path / "bad.txt"
        _make_text_file(bad)
        r = _run(["tree", str(bad)])
        assert r.exit_code != 0


# ===========================================================================
# attrs
# ===========================================================================


class TestAttrsCLI:
    def test_attrs_help(self):
        r = _run(["attrs", "--help"])
        assert r.exit_code == 0

    def test_attrs_root(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert "format" in r.output

    def test_attrs_bins_group(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler), "/bins"])
        assert r.exit_code == 0

    def test_attrs_pixels_group(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler), "/pixels"])
        assert r.exit_code == 0

    def test_attrs_genomes_group(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler), "/genomes"])
        assert r.exit_code == 0

    def test_attrs_contigs_group(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler), "/contigs"])
        assert r.exit_code == 0

    def test_attrs_verbose(self, tiny_metacooler):
        r = _run(["attrs", "-v", str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- failures ---

    def test_attrs_nonexistent_file_fails(self, tmp_path):
        r = _run(["attrs", str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0

    def test_attrs_invalid_group_path_fails(self, tiny_metacooler):
        r = _run(["attrs", str(tiny_metacooler), "/nonexistent_group"])
        assert r.exit_code != 0

    def test_attrs_text_file_fails(self, tmp_path):
        bad = tmp_path / "bad.txt"
        _make_text_file(bad)
        r = _run(["attrs", str(bad)])
        assert r.exit_code != 0


# ===========================================================================
# aggregate
# ===========================================================================


class TestAggregateCLI:
    def test_aggregate_help(self):
        r = _run(["aggregate", "--help"])
        assert r.exit_code == 0
        assert "--level" in r.output

    def test_aggregate_contig_level(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg_contig.cool"
        r = _run(["aggregate", str(tiny_metacooler), "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_aggregate_genome_level(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg_genome.cool"
        r = _run(
            ["aggregate", "--level", "genome", str(tiny_metacooler), "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_aggregate_contig_produces_valid_cooler(self, tiny_metacooler, tmp_path):
        import cooler as clr

        out = tmp_path / "agg_valid.cool"
        r = _run(["aggregate", str(tiny_metacooler), "-o", str(out)])
        assert r.exit_code == 0, r.output
        c = clr.Cooler(str(out))
        assert len(c.chromnames) > 0

    def test_aggregate_genome_produces_three_chroms(self, tiny_metacooler, tmp_path):
        import cooler as clr

        out = tmp_path / "agg_genome_valid.cool"
        r = _run(
            ["aggregate", "--level", "genome", str(tiny_metacooler), "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        c = clr.Cooler(str(out))
        assert len(c.chromnames) == 3  # MAG_1, MAG_2, __unbinned__

    def test_aggregate_stress_small_chunksize(self, tiny_metacooler, tmp_path):
        """Stress: aggregate with tiny chunksize forces many iterations."""
        out = tmp_path / "agg_stress.cool"
        r = _run(
            ["aggregate", "--chunksize", "5", str(tiny_metacooler), "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_aggregate_verbose(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg_verbose.cool"
        r = _run(["aggregate", "-v", str(tiny_metacooler), "-o", str(out)])
        assert r.exit_code == 0, r.output

    # --- failures ---

    def test_aggregate_nonexistent_input_fails(self, tmp_path):
        out = tmp_path / "agg.cool"
        r = _run(["aggregate", str(tmp_path / "ghost.metacool"), "-o", str(out)])
        assert r.exit_code != 0

    def test_aggregate_invalid_level_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "agg.cool"
        r = _run(
            ["aggregate", "--level", "BADLEVEL", str(tiny_metacooler), "-o", str(out)]
        )
        assert r.exit_code != 0

    def test_aggregate_missing_output_arg_fails(self, tiny_metacooler):
        r = _run(["aggregate", str(tiny_metacooler)])
        assert r.exit_code != 0

    def test_aggregate_missing_both_args_fails(self):
        r = _run(["aggregate"])
        assert r.exit_code != 0


# ===========================================================================
# extract
# ===========================================================================


class TestExtractCLI:
    def test_extract_help(self):
        r = _run(["extract", "--help"])
        assert r.exit_code == 0
        assert "--genomes" in r.output

    def test_extract_mag1(self, tiny_metacooler, tmp_path):
        import cooler as clr

        out = tmp_path / "mag1.cool"
        r = _run(
            ["extract", str(tiny_metacooler), "--genomes", "MAG_1", "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()
        c = clr.Cooler(str(out))
        assert set(c.chromnames) == {"contig_A", "contig_B"}

    def test_extract_mag2(self, tiny_metacooler, tmp_path):
        import cooler as clr

        out = tmp_path / "mag2.cool"
        r = _run(
            ["extract", str(tiny_metacooler), "--genomes", "MAG_2", "-o", str(out)]
        )
        assert r.exit_code == 0, r.output
        c = clr.Cooler(str(out))
        assert "contig_C" in c.chromnames

    def test_extract_unbinned(self, tiny_metacooler, tmp_path):
        out = tmp_path / "unbinned.cool"
        r = _run(
            [
                "extract",
                str(tiny_metacooler),
                "--genomes",
                "__unbinned__",
                "-o",
                str(out),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_extract_verbose(self, tiny_metacooler, tmp_path):
        out = tmp_path / "mag1_v.cool"
        r = _run(
            [
                "extract",
                "-v",
                str(tiny_metacooler),
                "--genomes",
                "MAG_1",
                "-o",
                str(out),
            ]
        )
        assert r.exit_code == 0, r.output

    def test_extract_all_genomes_produce_valid_coolers(self, tiny_metacooler, tmp_path):
        import cooler as clr

        for g in ["MAG_1", "MAG_2", "__unbinned__"]:
            out = tmp_path / f"{g}.cool"
            r = _run(["extract", str(tiny_metacooler), "--genomes", g, "-o", str(out)])
            assert r.exit_code == 0, f"extract {g} failed: {r.output}"
            c = clr.Cooler(str(out))
            assert len(c.chromnames) > 0

    # --- failures ---

    def test_extract_nonexistent_input_fails(self, tmp_path):
        out = tmp_path / "out.cool"
        r = _run(
            [
                "extract",
                str(tmp_path / "ghost.metacool"),
                "--genomes",
                "MAG_1",
                "-o",
                str(out),
            ]
        )
        assert r.exit_code != 0

    def test_extract_missing_genome_option_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.cool"
        r = _run(["extract", str(tiny_metacooler), "-o", str(out)])
        assert r.exit_code != 0

    def test_extract_missing_output_option_fails(self, tiny_metacooler):
        r = _run(["extract", str(tiny_metacooler), "--genomes", "MAG_1"])
        assert r.exit_code != 0

    def test_extract_unknown_genome_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "out.cool"
        r = _run(
            [
                "extract",
                str(tiny_metacooler),
                "--genomes",
                "NONEXISTENT",
                "-o",
                str(out),
            ]
        )
        assert r.exit_code != 0


# ===========================================================================
# show
# ===========================================================================


class TestShowCLI:
    # NOTE: ``show`` takes the metacool URI as a POSITIONAL argument, and ``-f``
    # means ``--force``.  Passing ``-f <path>`` still "works" by accident (the
    # path falls through to the positional slot) while silently enabling
    # --force, so every test here must use the positional form.
    def test_show_help(self):
        r = _run(["show", "--help"])
        assert r.exit_code == 0
        assert "--level" in r.output

    def test_show_genome_level_to_file(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_genome.png"
        r = _run(["show", "--level", "genome", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_contig_level_to_file(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_contig.png"
        r = _run(["show", "--level", "contig", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_bin_level_to_file(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_bin.png"
        # bin level needs a selection (otherwise the matrix is impossibly large)
        r = _run(
            [
                "show",
                "--level",
                "bin",
                "--genomes",
                "MAG_1",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_single_genome_view(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_intra.png"
        r = _run(["show", "--genomes", "MAG_1", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    @pytest.mark.parametrize(
        "genomes_args",
        [
            ["--genomes", "MAG_1,MAG_2"],
            ["--genomes", "MAG_1", "--genomes", "MAG_2"],
            ["--genomes", "MAG_1,", "--genomes", "MAG_2"],
        ],
        ids=["comma-list", "repeated", "mixed"],
    )
    def test_show_two_genome_view(self, tiny_metacooler, tmp_path, genomes_args):
        """``--genomes`` is a comma-separated/repeatable list, not ``nargs=2``."""
        out = tmp_path / "show_inter.png"
        r = _run(
            ["show", *genomes_args, "--dpi", "50", "-o", str(out), str(tiny_metacooler)]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_space_separated_genomes_is_rejected(self, tiny_metacooler, tmp_path):
        """The old ``--genomes A B`` form: B would land on the positional slot."""
        out = tmp_path / "show_space.png"
        r = _run(
            [
                "show",
                "--genomes",
                "MAG_1",
                "MAG_2",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 2, r.output
        assert not out.exists()

    def test_show_more_than_two_genomes(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_multi.png"
        r = _run(
            [
                "show",
                "--genomes",
                "MAG_1,MAG_2,__unbinned__",
                "--level",
                "contig",
                "--dpi",
                "50",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    @pytest.mark.parametrize(
        "contigs_args",
        [
            ["--contigs", "contig_A,contig_C"],
            ["--contigs", "contig_A", "--contigs", "contig_C"],
        ],
        ids=["comma-list", "repeated"],
    )
    def test_show_contigs_selection(self, tiny_metacooler, tmp_path, contigs_args):
        out = tmp_path / "show_contigs.png"
        r = _run(
            [
                "show",
                *contigs_args,
                "--level",
                "bin",
                "--dpi",
                "50",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_contigs_upgrades_genome_level(self, tiny_metacooler, tmp_path):
        """A contig has no row of its own at genome level, so the level moves."""
        out = tmp_path / "show_contigs_up.png"
        r = _run(
            [
                "show",
                "--contigs",
                "contig_A,contig_C",
                "--level",
                "genome",
                "--dpi",
                "50",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_genomes_and_contigs_are_mutually_exclusive(
        self, tiny_metacooler, tmp_path
    ):
        out = tmp_path / "show_both.png"
        r = _run(
            [
                "show",
                "--genomes",
                "MAG_1",
                "--contigs",
                "contig_C",
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 1, r.output
        assert "mutually exclusive" in _flat(r.output)
        assert not out.exists()

    def test_show_empty_contigs_is_rejected(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_empty_c.png"
        r = _run(["show", "--contigs", "", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 2, r.output
        assert "is empty" in _flat(r.output)

    def test_show_empty_genomes_is_rejected(self, tiny_metacooler, tmp_path):
        """Matches the shared ``--genomes ''`` policy of every other command."""
        out = tmp_path / "show_empty.png"
        r = _run(["show", "--genomes", "", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 2, r.output
        assert "is empty" in _flat(r.output)

    def test_show_transform_none(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_raw.png"
        r = _run(["show", "--transform", "none", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_custom_cmap(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_blues.png"
        r = _run(["show", "--cmap", "Blues", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_custom_dpi(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_dpi.png"
        r = _run(["show", "--dpi", "72", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_verbose(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_v.png"
        r = _run(["show", "-v", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_show_quiet(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show_q.png"
        r = _run(["show", "-q", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    # --- the -f/--force contract ---

    def test_show_refuses_to_overwrite_without_force(self, tiny_metacooler, tmp_path):
        """``-f`` on ``show`` is ``--force``, not the input file option."""
        out = tmp_path / "show_force.png"
        first = _run(["show", "--dpi", "50", "-o", str(out), str(tiny_metacooler)])
        assert first.exit_code == 0, first.output
        stamp = out.stat().st_mtime_ns

        second = _run(["show", "--dpi", "50", "-o", str(out), str(tiny_metacooler)])
        assert second.exit_code == 1, second.output
        assert "already exists" in _flat(second.output)
        assert "--force" in _flat(second.output)
        assert out.stat().st_mtime_ns == stamp  # untouched

    @pytest.mark.parametrize("flag", ["-f", "--force"])
    def test_show_force_overwrites_existing_output(
        self, tiny_metacooler, tmp_path, flag
    ):
        out = tmp_path / f"show_force{flag.strip('-')}.png"
        assert (
            _run(
                ["show", "--dpi", "50", "-o", str(out), str(tiny_metacooler)]
            ).exit_code
            == 0
        )
        r = _run(["show", "--dpi", "50", "-o", str(out), flag, str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    # --- failures ---

    def test_show_nonexistent_input_fails(self, tmp_path):
        """Missing input is rejected by MetaCoolURI at parse time (exit 2)."""
        out = tmp_path / "show.png"
        ghost = tmp_path / "ghost.metacool"
        r = _run(["show", "-o", str(out), str(ghost)])
        assert r.exit_code == 2, r.output
        flat = _flat(r.output)
        assert "Invalid value for 'METACOOL_PATH'" in flat
        # the path is a long unbreakable token, so ignore whitespace entirely
        assert f"Path'{ghost}'doesnotexist." in _squash(r.output)
        assert not out.exists()

    def test_show_invalid_level_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        r = _run(["show", "--level", "BADLEVEL", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code == 2, r.output
        flat = _flat(r.output)
        assert "Invalid value for '--level'" in flat
        assert "'BADLEVEL' is not one of 'bin', 'contig', 'genome'." in flat
        assert not out.exists()

    def test_show_invalid_transform_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "show.png"
        r = _run(
            [
                "show",
                "--transform",
                "sqrt",  # not a valid choice
                "-o",
                str(out),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 2, r.output
        flat = _flat(r.output)
        assert "Invalid value for '--transform'" in flat
        assert "'sqrt' is not one of 'log1p', 'log10', 'none'." in flat
        assert not out.exists()

    def test_show_unknown_genome_fails(self, tiny_metacooler, tmp_path):
        """An unknown genome reaches the command body and errors there (exit 1)."""
        out = tmp_path / "show.png"
        r = _run(
            ["show", "--genomes", "NONEXISTENT", "-o", str(out), str(tiny_metacooler)]
        )
        assert r.exit_code == 1, r.output
        flat = _flat(r.output)
        assert "Unknown genome: 'NONEXISTENT'" in flat
        assert "MAG_1" in flat  # the message lists the available genomes
        assert not out.exists()

    def test_show_bare_mmetacool_rejected(self, tiny_mmetacool, tmp_path):
        """A collection root is rejected by MetaCoolURI(require_entry=True)."""
        out = tmp_path / "show.png"
        r = _run(["show", "-o", str(out), str(tiny_mmetacool)])
        assert r.exit_code == 2, r.output
        flat = _flat(r.output)
        assert "Invalid value for 'METACOOL_PATH'" in flat
        assert "is an mmetacool collection" in flat
        assert not out.exists()


# ===========================================================================
# merge
# ===========================================================================


class TestMergeCLI:
    def test_merge_help(self):
        r = _run(["merge", "--help"])
        assert r.exit_code == 0

    def test_merge_two_identical_files(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_merge_pixel_counts_doubled(self, tiny_metacooler, tmp_path):
        """Merging a file with itself should double the total pixel sum."""
        from metacooler.api import MetaCooler

        out = tmp_path / "merged_double.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output

        orig = MetaCooler(tiny_metacooler)
        merged = MetaCooler(out)
        assert merged.pixels()["count"].sum() == orig.pixels()["count"].sum() * 2

    def test_merge_verbose(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged_v.metacool"
        r = _run(
            ["merge", "-o", str(out), "-v", str(tiny_metacooler), str(tiny_metacooler)]
        )
        assert r.exit_code == 0, r.output

    def test_merge_three_files(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged3.metacool"
        r = _run(
            [
                "merge",
                "-o",
                str(out),
                str(tiny_metacooler),
                str(tiny_metacooler),
                str(tiny_metacooler),
            ]
        )
        assert r.exit_code == 0, r.output
        assert out.exists()

    def test_merge_result_is_valid_metacooler(self, tiny_metacooler, tmp_path):
        from metacooler.schema import is_metacooler

        out = tmp_path / "merged_valid.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        assert is_metacooler(out)

    def test_merge_stores_command_in_metadata(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged_cmd.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler), str(tiny_metacooler)])
        assert r.exit_code == 0, r.output
        list_r = _run(["metadata", "list", str(out)])
        assert list_r.exit_code == 0
        data = json.loads(list_r.output)
        assert "commands" in data
        assert isinstance(data["commands"], list)
        assert any("merge" in cmd for cmd in data["commands"])

    # --- failures ---

    def test_merge_single_file_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "single.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler)])
        assert r.exit_code != 0

    def test_merge_zero_files_fails(self, tmp_path):
        out = tmp_path / "zero.metacool"
        r = _run(["merge", "-o", str(out)])
        assert r.exit_code != 0

    def test_merge_nonexistent_input_fails(self, tiny_metacooler, tmp_path):
        out = tmp_path / "merged.metacool"
        r = _run(
            [
                "merge",
                "-o",
                str(out),
                str(tiny_metacooler),
                str(tmp_path / "ghost.metacool"),
            ]
        )
        assert r.exit_code != 0

    def test_merge_mismatched_axes_fails(
        self, tiny_metacooler, tmp_path, tiny_contigs_table, tiny_bins_bed, tiny_pairs
    ):
        """Files with different genomes tables cannot be merged."""
        alt = tmp_path / "alt.metacool"
        r_build = _run(
            [
                "create",
                "pairs",
                "-o",
                str(alt),
                "--bins",
                str(tiny_bins_bed),
                "--pairs",
                str(tiny_pairs),
                "--contigs",
                str(tiny_contigs_table),
                "--drop-unassigned",  # different genomes table → hash mismatch
            ]
        )
        assert r_build.exit_code == 0, r_build.output

        out = tmp_path / "mismatch_merged.metacool"
        r = _run(["merge", "-o", str(out), str(tiny_metacooler), str(alt)])
        assert r.exit_code != 0


# ===========================================================================
# metadata
# ===========================================================================


class TestMetadataCLI:
    # --- help ---

    def test_metadata_help(self):
        r = _run(["metadata", "--help"])
        assert r.exit_code == 0
        for sub in ["list", "get", "add", "remove", "load"]:
            assert sub in r.output

    def test_metadata_list_help(self):
        r = _run(["metadata", "list", "--help"])
        assert r.exit_code == 0

    def test_metadata_get_help(self):
        r = _run(["metadata", "get", "--help"])
        assert r.exit_code == 0

    def test_metadata_add_help(self):
        r = _run(["metadata", "add", "--help"])
        assert r.exit_code == 0

    def test_metadata_remove_help(self):
        r = _run(["metadata", "remove", "--help"])
        assert r.exit_code == 0

    def test_metadata_load_help(self):
        r = _run(["metadata", "load", "--help"])
        assert r.exit_code == 0

    # --- list ---

    def test_metadata_list(self, tiny_metacooler):
        r = _run(["metadata", "list", str(tiny_metacooler)])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert isinstance(data, dict)
        assert "format" not in data

    def test_metadata_list_holds_only_the_qc_block_on_a_fresh_file(
        self, tiny_metacooler
    ):
        # A fresh metacool file carries no *user* metadata; the only key is the
        # contact-quality block ``create`` records itself.
        r = _run(["metadata", "list", str(tiny_metacooler)])
        assert r.exit_code == 0
        assert _user_metadata(r.output) == {}
        assert json.loads(r.output)["qc"]["source"] == "pairs"

    def test_metadata_list_verbose(self, tiny_metacooler):
        r = _run(["metadata", "list", "-v", str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- get ---

    def test_metadata_get_existing_key(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "my_key",
                "--value",
                "my_value",
            ]
        )
        r = _run(["metadata", "get", str(tiny_metacooler), "--key", "my_key"])
        assert r.exit_code == 0
        value = json.loads(r.output)
        assert value == "my_value"

    def test_metadata_get_nested_key(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "nested_key",
                "--value",
                '{"a": 1}',
            ]
        )
        r = _run(["metadata", "get", str(tiny_metacooler), "-k", "nested_key"])
        assert r.exit_code == 0
        value = json.loads(r.output)
        assert value == {"a": 1}

    def test_metadata_get_verbose(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "v_key",
                "--value",
                "v_value",
            ]
        )
        r = _run(["metadata", "get", "-v", str(tiny_metacooler), "--key", "v_key"])
        assert r.exit_code == 0

    def test_metadata_get_nonexistent_key_fails(self, tiny_metacooler):
        r = _run(
            [
                "metadata",
                "get",
                str(tiny_metacooler),
                "--key",
                "this_key_does_not_exist",
            ]
        )
        assert r.exit_code != 0

    # --- add ---

    def test_metadata_add_string_value(self, tiny_metacooler):
        r = _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "my_key",
                "--value",
                "my_value",
            ]
        )
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "my_key"])
        assert get_r.exit_code == 0
        assert "my_value" in get_r.output

    def test_metadata_add_integer_value(self, tiny_metacooler):
        r = _run(
            ["metadata", "add", str(tiny_metacooler), "-k", "int_key", "--value", "42"]
        )
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "int_key"])
        assert get_r.exit_code == 0
        assert "42" in get_r.output

    def test_metadata_add_json_number_value(self, tiny_metacooler):
        """JSON numbers are parsed and stored in the metadata attribute."""
        r = _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "json_num_key",
                "--value",
                "99",
            ]
        )
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "json_num_key"])
        assert get_r.exit_code == 0
        assert json.loads(get_r.output) == 99

    def test_metadata_add_json_dict_value(self, tiny_metacooler):
        """JSON objects can be stored via the metadata attribute."""
        r = _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "json_key",
                "--value",
                '{"a": 1}',
            ]
        )
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "json_key"])
        assert get_r.exit_code == 0
        assert json.loads(get_r.output) == {"a": 1}

    def test_metadata_add_overwrites_existing(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "ow_key",
                "--value",
                "first",
            ]
        )
        r = _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "ow_key",
                "--value",
                "second",
            ]
        )
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "ow_key"])
        assert "second" in get_r.output

    def test_metadata_add_verbose(self, tiny_metacooler):
        r = _run(
            [
                "metadata",
                "add",
                "-v",
                str(tiny_metacooler),
                "--key",
                "v_key",
                "--value",
                "v_val",
            ]
        )
        assert r.exit_code == 0

    # --- remove ---

    def test_metadata_remove_existing_key(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "rm_key",
                "--value",
                "rm_value",
            ]
        )
        r = _run(["metadata", "remove", str(tiny_metacooler), "--key", "rm_key"])
        assert r.exit_code == 0
        get_r = _run(["metadata", "get", str(tiny_metacooler), "--key", "rm_key"])
        assert get_r.exit_code != 0

    def test_metadata_remove_nonexistent_key_fails(self, tiny_metacooler):
        r = _run(
            [
                "metadata",
                "remove",
                str(tiny_metacooler),
                "--key",
                "this_does_not_exist_at_all",
            ]
        )
        assert r.exit_code != 0

    def test_metadata_remove_verbose(self, tiny_metacooler):
        _run(
            [
                "metadata",
                "add",
                str(tiny_metacooler),
                "--key",
                "rv_key",
                "--value",
                "rv_val",
            ]
        )
        r = _run(["metadata", "remove", "-v", str(tiny_metacooler), "--key", "rv_key"])
        assert r.exit_code == 0

    # --- load ---

    def test_metadata_load_json_file(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "meta.json"
        meta_file.write_text(json.dumps({"project": "test_project", "version": "2.0"}))
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0
        list_r = _run(["metadata", "list", str(tiny_metacooler)])
        assert list_r.exit_code == 0
        assert _user_metadata(list_r.output) == {
            "project": "test_project",
            "version": "2.0",
        }

    def test_metadata_load_yaml_file(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "meta.yaml"
        meta_file.write_text("lab: RSG\npipeline: metacooler\n")
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0
        list_r = _run(["metadata", "list", str(tiny_metacooler)])
        assert list_r.exit_code == 0
        assert _user_metadata(list_r.output) == {
            "lab": "RSG",
            "pipeline": "metacooler",
        }

    def test_metadata_load_yml_extension(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "meta.yml"
        meta_file.write_text("run_id: run001\n")
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0
        list_r = _run(["metadata", "list", str(tiny_metacooler)])
        assert list_r.exit_code == 0
        assert _user_metadata(list_r.output) == {"run_id": "run001"}

    def test_metadata_load_multiple_keys(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "multi.json"
        meta_file.write_text(json.dumps({"k1": "v1", "k2": "v2", "k3": 99}))
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0
        list_r = _run(["metadata", "list", str(tiny_metacooler)])
        assert list_r.exit_code == 0
        assert _user_metadata(list_r.output) == {"k1": "v1", "k2": "v2", "k3": "99"}

    def test_metadata_load_nested_dict(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "nested.json"
        meta_file.write_text(
            json.dumps({"outer": {"inner": 42, "flag": True}, "tags": [1, 2, "three"]})
        )
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0
        list_r = _run(["metadata", "list", str(tiny_metacooler)])
        assert list_r.exit_code == 0
        assert _user_metadata(list_r.output) == {
            "outer": {"inner": "42", "flag": "True"},
            "tags": ["1", "2", "three"],
        }

    def test_metadata_load_nonexistent_file_fails(self, tiny_metacooler, tmp_path):
        r = _run(
            ["metadata", "load", str(tmp_path / "ghost.json"), str(tiny_metacooler)]
        )
        assert r.exit_code != 0

    def test_metadata_load_non_dict_json_fails(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "list.json"
        meta_file.write_text(json.dumps(["item1", "item2"]))
        r = _run(["metadata", "load", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code != 0

    def test_metadata_load_verbose(self, tiny_metacooler, tmp_path):
        meta_file = tmp_path / "meta_v.json"
        meta_file.write_text(json.dumps({"verbose_key": "verbose_val"}))
        r = _run(["metadata", "load", "-v", str(meta_file), str(tiny_metacooler)])
        assert r.exit_code == 0

    # --- failures: non-existent metacool file ---

    def test_metadata_list_nonexistent_fails(self, tmp_path):
        r = _run(["metadata", "list", str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0

    def test_metadata_get_nonexistent_file_fails(self, tmp_path):
        r = _run(
            ["metadata", "get", str(tmp_path / "ghost.metacool"), "--key", "format"]
        )
        assert r.exit_code != 0

    def test_metadata_add_nonexistent_file_fails(self, tmp_path):
        r = _run(
            [
                "metadata",
                "add",
                str(tmp_path / "ghost.metacool"),
                "--key",
                "k",
                "--value",
                "v",
            ]
        )
        assert r.exit_code != 0

    def test_metadata_remove_nonexistent_file_fails(self, tmp_path):
        r = _run(["metadata", "remove", str(tmp_path / "ghost.metacool"), "--key", "k"])
        assert r.exit_code != 0

    def test_metadata_load_nonexistent_metacool_fails(self, tmp_path):
        meta_file = tmp_path / "meta.json"
        meta_file.write_text('{"key": "val"}')
        r = _run(["metadata", "load", str(meta_file), str(tmp_path / "ghost.metacool")])
        assert r.exit_code != 0
