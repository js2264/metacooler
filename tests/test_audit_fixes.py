"""Regression tests from the audit of the 0.2.0 release.

Every test in this module pins down a defect found while re-auditing the
package after 0.2.0 was published. As in ``test_regressions.py``, each
docstring states what went wrong before the fix, and the tests assert on
observable behaviour (file contents, exit codes, raised messages) rather than
on how the fix is implemented.
"""

from __future__ import annotations

from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io._utils import read_user_table


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
