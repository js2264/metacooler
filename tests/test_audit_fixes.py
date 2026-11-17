"""Regression tests from the audit of the 0.2.0 release.

Every test in this module pins down a defect found while re-auditing the
package after 0.2.0 was published. As in ``test_regressions.py``, each
docstring states what went wrong before the fix, and the tests assert on
observable behaviour (file contents, exit codes, raised messages) rather than
on how the fix is implemented.
"""

from __future__ import annotations

import ast
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

import metacooler.cli
from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io._utils import read_user_table
from metacooler.schema import validate


def _run(*args: object):
    return CliRunner().invoke(cli, [str(a) for a in args])


class TestQuotedFieldsInUserTables:
    """A stray '"' in a free-text column swallowed every row up to the next one.

    ``pd.read_csv`` honours quotes by default, so ``"draft`` on one row and
    ``curated"`` three rows later became one field spanning the rows between
    them. The raw field count that guards against misaligned tables counts tabs
    and could not see it: the rows vanished silently, and in ``create`` their
    contigs went to ``__unbinned__`` with exit status 0.
    """

    TABLE = (
        "contig_name\tgenome_name\tnote\n"
        'contig_A\tMAG_1\t"draft\n'
        "contig_B\tMAG_1\tok\n"
        "contig_C\tMAG_2\tok\n"
        'contig_D\tMAG_2\tcurated"\n'
    )

    def test_every_row_is_read_and_quotes_are_kept_verbatim(self, tmp_path: Path):
        f = tmp_path / "t.tsv"
        f.write_text(self.TABLE)
        df = read_user_table(f, kind="Contig table")
        assert df["contig_name"].tolist() == [
            "contig_A",
            "contig_B",
            "contig_C",
            "contig_D",
        ]
        assert df["note"].tolist() == ['"draft', "ok", "ok", 'curated"']

    def test_headerless_tables_too(self, tmp_path: Path):
        f = tmp_path / "t.tsv"
        f.write_text("".join(self.TABLE.splitlines(keepends=True)[1:]))
        df = read_user_table(f, kind="Contig table", header=False)
        assert len(df) == 4

    def test_create_assigns_every_listed_contig(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
    ):
        table = tmp_path / "t.tsv"
        table.write_text(self.TABLE)
        out = tmp_path / "o.metacool"
        result = _run(
            "create", "pairs", "--pairs", tiny_pairs, "--bins", tiny_bins_bed,
            "--contigs", table, "-o", out, "--no-ps",
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        genomes = MetaCooler(out).genomes.set_index("name")["n_contigs"]
        assert genomes.to_dict() == {"MAG_1": 2, "MAG_2": 2}
        notes = MetaCooler(out).contigs.set_index("name")["note"]
        assert notes["contig_A"] == '"draft'


class TestNetworkIntegerOverflow:
    """``network`` summed int32 counts in int32, printing negative edge weights.

    The 0.2.0 rewrite routed the sum through ``aggregate_pixels``, which keeps
    the input dtype, so a contig or genome pair totalling more than 2**31
    contacts -- the ``__unbinned__`` self-edge of a deep library, say -- wrapped
    to a negative number while ``aggregate`` reported the right total.
    """

    @pytest.fixture
    def heavy(self, tmp_path: Path) -> Path:
        bins = pd.DataFrame(
            {
                "chrom": ["c1"] * 3 + ["c2"] * 2,
                "start": [0, 10, 20, 0, 10],
                "end": [10, 20, 30, 10, 20],
            }
        )
        pixels = pd.DataFrame(
            {
                "bin1_id": [0, 1, 3],
                "bin2_id": [1, 2, 4],
                "count": np.array([1_500_000_000, 1_500_000_000, 7], dtype=np.int32),
            }
        )
        src = tmp_path / "heavy.cool"
        cooler.create_cooler(
            str(src), bins, pixels, dtypes={"count": np.int32}, ordered=True
        )
        table = tmp_path / "c.tsv"
        table.write_text("contig_name\tgenome_name\nc1\tG1\nc2\tG2\n")
        out = tmp_path / "heavy.metacool"
        result = _run(
            "create", "cool", "--cool", src, "--contigs", table, "-o", out, "--no-qc"
        )
        assert result.exit_code == 0, result.output
        assert cooler.Cooler(str(out)).pixels()[:]["count"].dtype == np.int32
        return out

    @pytest.mark.parametrize(("level", "node"), [("genome", "G1"), ("contig", "c1")])
    def test_edge_weight_is_not_wrapped(self, heavy: Path, level: str, node: str):
        result = _run("network", heavy, "--level", level)
        assert result.exit_code == 0, result.output
        rows = [line.split("\t") for line in result.stdout.splitlines()]
        edges = {(r[0], r[1]): int(r[2]) for r in rows if len(r) == 3}
        assert edges[(node, node)] == 3_000_000_000


class TestCollectOutputGuard:
    """``collect -f`` deleted its own inputs, and a previous output it then failed to replace.

    Its ``check_output`` call passed neither ``inputs=`` nor ``defer_removal``
    -- the one command added after the guard was introduced, and the one that
    missed it. ``collect A.cool B.cool --genomes g.tsv -o g.tsv -f`` deleted
    ``g.tsv`` and then failed with "Genomes table not found".
    """

    @pytest.fixture
    def cools(self, tmp_path: Path, tiny_cool: Path) -> list[Path]:
        out = []
        for name in ("A", "B"):
            dst = tmp_path / f"{name}.cool"
            dst.write_bytes(tiny_cool.read_bytes())
            out.append(dst)
        return out

    @pytest.mark.parametrize("option", ["--genomes", "--metadata"])
    def test_an_input_table_survives_being_named_as_output(
        self, tmp_path: Path, cools: list[Path], option: str
    ):
        table = tmp_path / "in.tsv"
        table.write_text("genome_name\tphylum\nA\tx\nB\ty\n")
        if option == "--metadata":
            table = tmp_path / "in.yaml"
            table.write_text("project: x\n")
        before = table.read_bytes()
        result = _run("collect", *cools, option, table, "-o", table, "-f")
        assert result.exit_code != 0
        assert "same file as the input" in flat(result.output)
        assert table.read_bytes() == before

    def test_an_input_cool_survives_being_named_as_output(self, cools: list[Path]):
        before = cools[0].read_bytes()
        result = _run("collect", *cools, "-o", cools[0], "-f")
        assert result.exit_code != 0
        assert cools[0].read_bytes() == before

    def test_a_failed_run_keeps_the_previous_output(
        self, tmp_path: Path, cools: list[Path]
    ):
        out = tmp_path / "prev.metacool"
        assert _run("collect", *cools, "-o", out).exit_code == 0
        before = out.read_bytes()
        bad = tmp_path / "bad.tsv"
        bad.write_text("no_such_column\tx\nA\t1\n")
        result = _run("collect", *cools, "--genomes", bad, "-o", out, "-f")
        assert result.exit_code != 0
        assert out.read_bytes() == before


def test_every_command_names_its_inputs_to_the_output_guard():
    """Structural guard: a new command cannot forget ``inputs=`` again.

    ``collect`` was added after the guard learned to refuse an output that is
    one of the inputs, and its bare ``check_output(output, force)`` silently
    reintroduced the input-deleting ``-f``. Every call in the CLI must say what
    the command reads.
    """
    cli_dir = Path(metacooler.cli.__file__).parent
    missing = []
    for path in sorted(cli_dir.glob("cmd_*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "check_output"
                and not any(k.arg == "inputs" for k in node.keywords)
            ):
                missing.append(f"{path.name}:{node.lineno}")
    assert not missing, f"check_output without inputs=: {missing}"


class TestUserQcColumns:
    """A user's ``qc_*`` column in ``/genomes`` was accepted, then silently lost.

    Every ``qc_*`` column of ``/genomes`` is metacooler's QC namespace -- the
    writers that derive one map from another strip it -- yet ``create``,
    ``collect`` and ``annotate`` stored a user's ``qc_pass`` there under its own
    name. The next ``reassign`` (even one changing nothing) deleted it without
    a word, and ``annotate --drop genomes.qc_pass`` refused it as derived. It is
    now renamed ``user_qc_pass`` on the way in, like any other reserved name.
    """

    @pytest.fixture
    def tables(self, tmp_path: Path, tiny_contigs_table: Path) -> tuple[Path, Path]:
        genomes = tmp_path / "g.tsv"
        genomes.write_text(
            "genome_name\tqc_pass\tphylum\nMAG_1\tyes\tA\nMAG_2\tno\tB\n"
        )
        contigs = tmp_path / "c.tsv"
        contigs.write_text(
            "contig_name\tgenome_name\tqc_flag\n"
            "contig_A\tMAG_1\tok\ncontig_B\tMAG_1\tok\ncontig_C\tMAG_2\tlow\n"
        )
        return genomes, contigs

    @pytest.fixture
    def created(
        self, tmp_path: Path, tables, tiny_pairs: Path, tiny_bins_bed: Path
    ) -> Path:
        genomes, contigs = tables
        out = tmp_path / "a.metacool"
        result = _run(
            "create", "pairs", "--pairs", tiny_pairs, "--bins", tiny_bins_bed,
            "--contigs", contigs, "--genomes", genomes, "-o", out, "--no-ps",
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        return out

    def test_create_renames_it_out_of_the_qc_namespace(self, created: Path):
        columns = set(MetaCooler(created).genomes.columns)
        assert "user_qc_pass" in columns
        assert "qc_pass" not in columns
        assert (
            MetaCooler(created).genomes.set_index("name").loc["MAG_1", "user_qc_pass"]
            == "yes"
        )

    @pytest.mark.parametrize("layout", ["auto", "keep"])
    def test_user_columns_survive_reassign(
        self, tmp_path: Path, created: Path, tables, layout: str
    ):
        _, contigs = tables
        moved = tmp_path / "moved.tsv"
        moved.write_text(
            contigs.read_text().replace("contig_C\tMAG_2", "contig_C\tMAG_1")
        )
        out = tmp_path / f"r_{layout}.metacool"
        result = _run(
            "reassign", created, "--contigs", moved, "--layout", layout, "-o", out
        )
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert "user_qc_pass" in mc.genomes.columns
        assert mc.contigs.set_index("name").loc["contig_C", "qc_flag"] == "low"

    def test_annotate_renames_it_and_can_drop_it(self, tmp_path: Path, created):
        extra = tmp_path / "x.tsv"
        extra.write_text("genome_name\tqc_checked\nMAG_1\t1\nMAG_2\t0\n")
        out = tmp_path / "ann.metacool"
        result = _run("annotate", created, "--genomes", extra, "-o", out)
        assert result.exit_code == 0, result.output
        assert "user_qc_checked" in MetaCooler(out).genomes.columns
        dropped = tmp_path / "dropped.metacool"
        result = _run(
            "annotate", out, "--drop", "genomes.user_qc_checked", "-o", dropped
        )
        assert result.exit_code == 0, result.output
        assert "user_qc_checked" not in MetaCooler(dropped).genomes.columns


class TestHdf5InvalidColumnNames:
    """A user column named ``.`` passed every check, then failed at write time.

    ``rename_key_columns`` rewrote ``/`` and integer labels but not ``.``,
    which names the HDF5 group itself: ``create`` failed only after its whole
    pairs pass, and an in-place ``annotate`` failed halfway through writing.
    """

    @pytest.mark.parametrize("bad", [".", ".."])
    def test_create_stores_it_under_a_safe_name(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, bad: str
    ):
        contigs = tmp_path / "c.tsv"
        contigs.write_text(
            f"contig_name\tgenome_name\t{bad}\n"
            "contig_A\tMAG_1\t1\ncontig_B\tMAG_1\t2\ncontig_C\tMAG_2\t3\n"
        )
        genomes = tmp_path / "g.tsv"
        genomes.write_text(f"genome_name\t{bad}\nMAG_1\tx\nMAG_2\ty\n")
        out = tmp_path / "o.metacool"
        result = _run(
            "create", "pairs", "--pairs", tiny_pairs, "--bins", tiny_bins_bed,
            "--contigs", contigs, "--genomes", genomes, "-o", out, "--no-ps",
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert "col3" in MetaCooler(out).contigs.columns
        assert "col2" in MetaCooler(out).genomes.columns


class TestAnnotateInplaceIsAllOrNothing:
    """A failure while ``annotate --inplace`` wrote its second table half-edited the map.

    Tables were rebuilt and swapped in one at a time, so a failure building
    ``/contigs`` left ``/genomes`` already edited, no history recorded, and
    the ``.annotate-scratch`` group behind.
    """

    def test_a_failed_build_leaves_the_map_untouched(
        self, tmp_path: Path, tiny_metacooler: Path, monkeypatch: pytest.MonkeyPatch
    ):
        import h5py

        import metacooler.io.annotate as annotate_mod

        target = tmp_path / "m.metacool"
        target.write_bytes(tiny_metacooler.read_bytes())
        genomes = tmp_path / "g.tsv"
        genomes.write_text("genome_name\tphylum\nMAG_1\tA\nMAG_2\tB\n")
        contigs = tmp_path / "c.tsv"
        contigs.write_text("contig_name\tdepth\ncontig_A\t3.5\ncontig_B\t1.0\n")
        real = annotate_mod.write_column

        def failing(group, name, values):
            if name == "depth":
                raise OSError("disk full")
            return real(group, name, values)

        monkeypatch.setattr(annotate_mod, "write_column", failing)
        before = {
            "genomes": list(MetaCooler(target).genomes.columns),
            "contigs": list(MetaCooler(target).contigs.columns),
            "metadata": MetaCooler(target).metadata,
        }
        with pytest.raises(OSError, match="disk full"):
            annotate_mod.annotate(
                target, genomes_table=genomes, contigs_table=contigs, inplace=True
            )
        mc = MetaCooler(target)
        assert list(mc.genomes.columns) == before["genomes"]
        assert list(mc.contigs.columns) == before["contigs"]
        assert mc.metadata == before["metadata"]
        with h5py.File(target, "r") as f:
            assert annotate_mod._SCRATCH not in f
        validate(target)
