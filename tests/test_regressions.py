"""Regression tests

Every test in this module pins down a defect that was found by reading the
source and fixed in the same pass. Each docstring states what *silently* went
wrong before the fix, because in almost every case the old behaviour produced a
structurally valid file (or a zero exit code) carrying the wrong numbers — the
failure mode that unit tests exist to catch.

The tests assert on observable behaviour (file contents, exit codes, log
records, raised messages) rather than on how the fix is implemented.
"""

from __future__ import annotations

import gzip
import json
import logging
import shutil
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat
from _cli_output import squash as _squash
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.create import _read_bins_bed, from_cool, from_pairs
from metacooler.io.extract import extract
from metacooler.io.flatten import flatten
from metacooler.io.genomes import parse_genomes_table
from metacooler.io.merge import merge
from metacooler.ops.binnify import binnify
from metacooler.ops.zoomify import zoomify
from metacooler.schema import CONTIGS_DTYPE, validate

DATA_DIR = Path(__file__).parent / "data"

# Larger than int32 can hold: the value every "count dtype" test below tracks
# through a copy operation.  A silent ``astype(np.int32)`` turns it negative.
BIG_COUNT = 3_000_000_000


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# Helpers (local to this module by design — conftest.py is owned elsewhere)
# ---------------------------------------------------------------------------


def _write_pairs(path: Path, rows: list[tuple[str, int | str, str, int | str]]) -> Path:
    """Write a minimal 4-field pairs file (readID, chrom1, pos1, chrom2, pos2)."""
    lines = ["#columns: readID chrom1 pos1 chrom2 pos2\n"]
    lines += [
        f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\n" for i, (c1, p1, c2, p2) in enumerate(rows)
    ]
    path.write_text("".join(lines))
    return path


def _make_metacool(out: Path, pairs: Path, bins: Path, contigs: Path, **kw) -> Path:
    """Build a tiny ``.metacool`` from the standard fixture inputs."""
    from_pairs(
        uri=out,
        bins=bins,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_table_format="contigs",
        contig_ordering=kw.pop("contig_ordering", "genome-then-alpha"),
        chunksize=1000,
        **kw,
    )
    return out


def _file_attrs(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        return dict(f.attrs)


def _forge_wide_count(path: Path, value: int = BIG_COUNT, index: int = 0) -> int:
    """Rewrite ``/pixels/count`` as int64 with one huge value, as ``binnify`` does.

    ``ops.binnify`` deliberately writes int64 counts, so a re-binned file can
    legitimately hold values beyond the int32 range.  This reproduces such a
    file cheaply, without having to sum two billion contacts.

    Returns:
        The forged value, for convenience at the call site.
    """
    with h5py.File(path, "r+") as f:
        counts = np.asarray(f["pixels"]["count"][:]).astype(np.int64)
        counts[index] = value
        del f["pixels"]["count"]
        f["pixels"].create_dataset("count", data=counts, dtype=np.int64)
        f.attrs["sum"] = int(counts.sum())
    return value


def _max_count(path: Path) -> int:
    with h5py.File(path, "r") as f:
        return int(np.asarray(f["pixels"]["count"][:]).max())


def _append_empty_genome(path: Path, name: str = "MAG_EMPTY") -> None:
    """Append a genome carrying zero contigs, keeping the file valid.

    ``from_pairs`` cannot produce such a genome (genome rows are derived from
    the contig assignments), but a file written by another tool can, and
    ``binnify --level genome`` used to turn it into a zero-length contig that
    ``schema.validate`` then rejected.
    """
    with h5py.File(path, "r+") as f:
        g = f["genomes"]
        names = [
            v.decode() if isinstance(v, bytes) else str(v) for v in g["name"][:]
        ] + [name]
        columns = {
            "name": np.array([s.encode() for s in names]),
            "length": np.append(np.asarray(g["length"][:]).astype(np.int64), 0),
            "n_contigs": np.append(np.asarray(g["n_contigs"][:]).astype(np.int32), 0),
            "n_bins": np.append(np.asarray(g["n_bins"][:]).astype(np.int64), 0),
        }
        # Extend every *other* column too — a table with one short column is a
        # different bug from the one under test, and `validate` rejects it
        # first. `create` writes qc_* columns here, and a user genomes table
        # may add arbitrary ones.
        for key in g:
            if key in columns:
                continue
            values = np.asarray(g[key][:])
            columns[key] = np.append(values, values.dtype.type(0))
        for key, data in columns.items():
            del g[key]
            g.create_dataset(key, data=data)
        offsets = np.asarray(f["indexes"]["genome_offset"][:]).astype(np.int64)
        del f["indexes"]["genome_offset"]
        f["indexes"].create_dataset(
            "genome_offset", data=np.append(offsets, offsets[-1])
        )
        f.attrs["ngenomes"] = len(names)


# ---------------------------------------------------------------------------
# #1–#5 — reading a pairs file
# ---------------------------------------------------------------------------


class TestPairsIngestion:
    def test_bgz_pairs_are_decompressed_not_read_as_binary(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """#1 — ``pa.input_stream`` has no ``.bgz`` codec.

        Before: a bgzipped pairs file was handed to the CSV parser as raw
        binary, so no contig name was ever recognised and the run produced a
        structurally valid ``.metacool`` with ``nnz=0`` — no error, no warning.
        """
        bgz = tmp_path / "pairs.pairs.bgz"
        bgz.write_bytes(gzip.compress(tiny_pairs.read_bytes()))

        plain_out = _make_metacool(
            tmp_path / "plain.metacool", tiny_pairs, tiny_bins_bed, tiny_contigs_table
        )
        bgz_out = _make_metacool(
            tmp_path / "bgz.metacool", bgz, tiny_bins_bed, tiny_contigs_table
        )

        plain, bgzd = _file_attrs(plain_out), _file_attrs(bgz_out)
        assert int(plain["nnz"]) > 0
        assert int(bgzd["nnz"]) == int(plain["nnz"])
        assert int(bgzd["sum"]) == int(plain["sum"])

    def test_pairs_on_unknown_contigs_are_counted_and_reported(
        self, tmp_path: Path, tiny_bins_bed: Path, tiny_contigs_table: Path, caplog
    ) -> None:
        """#2 — dropped pairs had no counter, no warning and no summary.

        Before: a contig-name-format mismatch discarded 100 % of the contacts
        and the run cheerfully reported ``NNZ=0``.
        """
        pairs = _write_pairs(
            tmp_path / "unknown.pairs",
            [("chr1", 10, "chr1", 20), ("chr2", 10, "chr2", 20)] * 5,
        )
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = _make_metacool(
                tmp_path / "unknown.metacool", pairs, tiny_bins_bed, tiny_contigs_table
            )
        text = caplog.text
        assert "ALL 10 pair(s) were discarded" in text
        assert "chr1" in text and "chr2" in text
        assert int(_file_attrs(out)["nnz"]) == 0

    def test_partially_dropped_pairs_are_reported_with_a_percentage(
        self, tmp_path: Path, tiny_bins_bed: Path, tiny_contigs_table: Path, caplog
    ) -> None:
        """#2 — a partial drop was equally silent; it now reports a ratio."""
        pairs = _write_pairs(
            tmp_path / "mixed.pairs",
            [("contig_A", 10, "contig_A", 20)] * 3 + [("nope", 10, "nope", 20)],
        )
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = _make_metacool(
                tmp_path / "mixed.metacool", pairs, tiny_bins_bed, tiny_contigs_table
            )
        assert "Discarded 1 of 4 pair(s)" in caplog.text
        assert int(_file_attrs(out)["nnz"]) == 1

    def test_null_positions_are_dropped_instead_of_becoming_int64_min_bins(
        self, tmp_path: Path, tiny_bins_bed: Path, tiny_contigs_table: Path, caplog
    ) -> None:
        """#3 — an empty ``pos`` field became a NaN, then ``INT64_MIN``.

        Before: pyarrow turned the empty field into a null, pandas into a
        float64 NaN, and the int64 cast into ``INT64_MIN``; the negative-position
        warning never fired because ``NaN < 0`` is False.
        """
        pairs = _write_pairs(
            tmp_path / "nullpos.pairs",
            [
                ("contig_A", 10, "contig_A", 20),
                ("contig_A", "", "contig_A", 30),
                ("contig_B", 5, "contig_B", ""),
            ],
        )
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = _make_metacool(
                tmp_path / "nullpos.metacool", pairs, tiny_bins_bed, tiny_contigs_table
            )
        assert "with a missing position" in caplog.text
        # Exactly the one complete pair survives, on an in-range bin.
        with h5py.File(out, "r") as f:
            b1 = np.asarray(f["pixels"]["bin1_id"][:])
            n_bins = int(f.attrs["nbins"])
        assert int(_file_attrs(out)["nnz"]) == 1
        assert b1.min() >= 0 and b1.max() < n_bins
        validate(out)

    def test_bins_and_binsize_are_mutually_exclusive(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """#4 — ``binsize`` was documented as ignored when ``bins`` is given.

        Before: it was passed through anyway, so every contact was binned
        *arithmetically* against a bin table that had never been tiled at that
        size — contacts scattered onto unrelated fragments, and the file still
        passed ``validate()``.
        """
        with pytest.raises(ValueError, match="mutually exclusive"):
            from_pairs(
                uri=tmp_path / "both.metacool",
                bins=tiny_bins_bed,
                binsize=1000,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
            )

    def test_fasta_ordering_without_a_reference_order_is_refused(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """#5 — ``contig-ordering='fasta'`` silently meant *alphabetical*.

        Before: with lengths coming from a BED, ``groupby("chrom")`` produced an
        alphabetical order, yet the file still recorded ``contig-ordering=
        "fasta"`` — two files could advertise the same ordering with different
        physical layouts, which ``merge`` would then scramble.  This covers the
        Python API path (the CLI path is covered in ``test_cli_create.py``).
        """
        with pytest.raises(ValueError, match="needs a reference contig order"):
            from_pairs(
                uri=tmp_path / "fasta_order.metacool",
                bins=tiny_bins_bed,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                contig_ordering="fasta",
            )

    def test_fasta_ordering_with_a_fasta_follows_the_fasta_record_order(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
        tiny_fasta: Path,
    ) -> None:
        """#5 (positive) — given a FASTA, the recorded order really is its order."""
        out = _make_metacool(
            tmp_path / "fasta_order.metacool",
            tiny_pairs,
            tiny_bins_bed,
            tiny_contigs_table,
            fasta=tiny_fasta,
            contig_ordering="fasta",
        )
        mc = MetaCooler(out)
        assert mc.contigs["name"].tolist() == [
            "contig_A",
            "contig_B",
            "contig_C",
            "contig_D",
        ]
        assert _file_attrs(out)["contig-ordering"] == "fasta"


# ---------------------------------------------------------------------------
# #6–#7 — user-supplied contig / genome tables
# ---------------------------------------------------------------------------


class TestUserTables:
    def test_colspec_pointing_at_one_column_twice_is_refused(
        self, tmp_path: Path, tiny_contigs_table: Path
    ) -> None:
        """#6 — ``--contigs-format '1,1'`` died with ``KeyError: 'contig_name'``.

        Renaming column 1 twice left the frame with duplicate labels, so
        ``df["contig_name"]`` returned a DataFrame and the bare KeyError told
        the user nothing about what was wrong with their spec.
        """
        with pytest.raises(ValueError, match="same column"):
            parse_contigs_table(tiny_contigs_table, fmt="1,1")

    def test_colspec_with_a_colliding_column_name_still_maps_correctly(
        self, tmp_path: Path, caplog
    ) -> None:
        """#6 — a pre-existing ``contig_name``/``genome_id`` label poisoned the map.

        Before: the positional rename produced *duplicate* column labels, and
        ``dict(zip(Series, DataFrame))`` then built a garbage mapping with no
        exception at all.
        """
        path = tmp_path / "collide.tsv"
        path.write_text(
            "name\tbin\tcontig_name\n"
            "contig_A\tMAG_1\tfoo\n"
            "contig_B\tMAG_1\tbar\n"
            "contig_C\tMAG_2\tbaz\n"
        )
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            mapping, df = parse_contigs_table(path, fmt="1,2")
        assert mapping == {
            "contig_A": "MAG_1",
            "contig_B": "MAG_1",
            "contig_C": "MAG_2",
        }
        assert len(set(df.columns)) == len(df.columns)  # no duplicate labels
        assert "collides" in caplog.text

    def test_reserved_column_names_are_renamed_not_fatal(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, caplog
    ) -> None:
        """#7 — an extra column called ``length``/``n_bins`` crashed creation.

        Before: the contigs merge produced ``_x``/``_y`` suffixes
        (``KeyError: "['length'] not in index"``) and the genomes path died on
        ``ValueError: cannot insert n_bins`` — even though ``docs/format.md``
        explicitly advertises extra columns.
        """
        contigs = tmp_path / "contigs_reserved.tsv"
        contigs.write_text(
            "contig_name\tgenome_name\tlength\n"
            "contig_A\tMAG_1\t111\n"
            "contig_B\tMAG_1\t222\n"
            "contig_C\tMAG_2\t333\n"
        )
        genomes = tmp_path / "genomes_reserved.tsv"
        genomes.write_text("genome_name\tlength\tn_bins\nMAG_1\t9\t8\nMAG_2\t7\t6\n")

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = _make_metacool(
                tmp_path / "reserved.metacool",
                tiny_pairs,
                tiny_bins_bed,
                contigs,
                genomes_table=genomes,
            )

        mc = MetaCooler(out)
        assert "user_length" in mc.contigs.columns
        assert "user_length" in mc.genomes.columns
        assert "user_n_bins" in mc.genomes.columns
        # metacooler's own columns keep their real meaning.
        assert mc.contigs.loc[mc.contigs["name"] == "contig_A", "length"].item() == 1000
        assert "reserved column names" in caplog.text
        validate(out)


# ---------------------------------------------------------------------------
# #8–#9 — metadata and provenance
# ---------------------------------------------------------------------------


class TestMetadataAndProvenance:
    @pytest.mark.parametrize(
        "args",
        [
            ["metadata", "add", "-k", "x", "--value", "y"],
            ["metadata", "get", "-k", "x"],
            ["metadata", "remove", "-k", "x"],
        ],
    )
    def test_metadata_subcommands_reject_an_mmetacool_root(
        self, runner: CliRunner, tiny_mmetacool: Path, args: list[str]
    ) -> None:
        """#8 — only ``metadata list`` guarded against a collection root.

        Before: ``metacooler metadata add reps.mmetacool -k x --value y`` exited
        0 and wrote the value to the collection root, where it is invisible from
        every entry.
        """
        result = runner.invoke(cli, [args[0], args[1], str(tiny_mmetacool), *args[2:]])
        assert result.exit_code != 0
        assert "mmetacool collection" in _flat(result.output)
        assert "entries/entry1" in _squash(result.output)
        # And the root was left untouched.
        with h5py.File(tiny_mmetacool, "r") as f:
            assert "metadata" not in f.attrs

    def test_metadata_load_rejects_an_mmetacool_root(
        self, runner: CliRunner, tmp_path: Path, tiny_mmetacool: Path
    ) -> None:
        """#8 — ``metadata load`` had the same hole as ``add``/``get``/``remove``."""
        meta_file = tmp_path / "meta.json"
        meta_file.write_text(json.dumps({"lab": "pasteur"}))
        result = runner.invoke(
            cli, ["metadata", "load", str(meta_file), str(tiny_mmetacool)]
        )
        assert result.exit_code != 0
        assert "mmetacool collection" in _flat(result.output)

    def test_metadata_file_cannot_hijack_the_provenance_key(
        self,
        runner: CliRunner,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """#9 — a ``--metadata`` file could overwrite ``commands``.

        Before: ``meta.update(loaded)`` let a user file silently replace the
        provenance list recording how the file was made.  ``metadata load``
        already guarded the key; ``create`` never did.
        """
        bad = tmp_path / "bad.json"
        bad.write_text(json.dumps({"commands": ["rm -rf /"]}))
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--pairs",
                str(tiny_pairs),
                "--bins",
                str(tiny_bins_bed),
                "--contigs",
                str(tiny_contigs_table),
                "--metadata",
                str(bad),
                "-o",
                str(tmp_path / "hijack.metacool"),
            ],
        )
        assert result.exit_code != 0
        assert "reserved for provenance" in result.output
        assert not (tmp_path / "hijack.metacool").exists()

    def test_legitimate_metadata_file_is_kept_alongside_provenance(
        self,
        runner: CliRunner,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """#9 (positive) — a normal metadata file survives, and so does ``commands``."""
        good = tmp_path / "good.json"
        good.write_text(json.dumps({"lab": "pasteur", "run": 7}))
        out = tmp_path / "meta_ok.metacool"
        result = runner.invoke(
            cli,
            [
                "create",
                "pairs",
                "--pairs",
                str(tiny_pairs),
                "--bins",
                str(tiny_bins_bed),
                "--contigs",
                str(tiny_contigs_table),
                "--metadata",
                str(good),
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        meta = MetaCooler(out).metadata
        assert meta["lab"] == "pasteur"
        assert meta["run"] == "7"
        assert len(meta["commands"]) == 1
        assert "create pairs" in meta["commands"][0]


# ---------------------------------------------------------------------------
# #10 — one attribute for the assembly name
# ---------------------------------------------------------------------------


class TestAssemblyName:
    """``genome-assembly`` vs the private ``assembly-name``."""

    def test_assembly_name_lands_in_the_cooler_attribute(
        self, tiny_metacooler: Path
    ) -> None:
        """#10 — metacooler wrote ``assembly-name``, cooler reads ``genome-assembly``.

        Before: the two attributes disagreed on the same file and
        ``metacooler info`` reported ``"assembly": null`` for a file that
        ``cooler info`` described perfectly well.
        """
        attrs = _file_attrs(tiny_metacooler)
        assert attrs["genome-assembly"] == "tiny_test"
        assert "assembly-name" not in attrs
        assert MetaCooler(tiny_metacooler).info["assembly"] == "tiny_test"
        assert (
            cooler.Cooler(str(tiny_metacooler)).info["genome-assembly"] == "tiny_test"
        )

    def test_shipped_fixture_no_longer_reports_a_null_assembly(self) -> None:
        """#10 — verified on disk: ``tests/data/assembly_1_1kb.metacool``."""
        path = DATA_DIR / "assembly_1_1kb.metacool"
        info = MetaCooler(path).info
        assert info["assembly"] is not None
        assert info["assembly"] == _file_attrs(path)["genome-assembly"]


# ---------------------------------------------------------------------------
# #11 — show() selection arguments
# ---------------------------------------------------------------------------


class TestShowSelection:
    """``--genomes``/``--contigs``/``--range`` never silently ignore each other.

    The original bug was ``--genome``, since removed in favour of a one-name
    ``--genomes``; these keep the invariants it was about.
    """

    def test_genome_filter_is_honoured_at_genome_level(
        self, tmp_path: Path, tiny_metacooler: Path, caplog
    ) -> None:
        """#11 — the filter was validated and then ignored.

        Before: ``show --level genome --genomes MAG_1`` rendered the *whole*
        genome-by-genome matrix and exited 0, so the user got a picture of
        something they did not ask for.  A single-genome selection is now a 1×1
        block, which the code says so out loud.
        """
        out = tmp_path / "one_genome.png"
        from metacooler.ops.show import show

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            show(tiny_metacooler, level="genome", genomes=["MAG_1"], out=out)
        assert out.exists()
        assert "single genome" in caplog.text

    def test_genomes_and_contigs_together_are_refused(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#11 — two selectors for the same axes must not be silently combined."""
        from metacooler.ops.show import show

        with pytest.raises(ValueError, match="mutually exclusive"):
            show(
                tiny_metacooler,
                level="contig",
                genomes=["MAG_1"],
                contigs=["contig_C"],
                out=tmp_path / "x.png",
            )

    def test_range_cannot_be_combined_with_genomes(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#11 — ``--range`` and ``--genomes`` select the same axes two ways."""
        from metacooler.ops.show import show

        with pytest.raises(ValueError, match="cannot be combined"):
            show(
                tiny_metacooler,
                level="bin",
                genomes=["MAG_1"],
                range1="contig_A",
                out=tmp_path / "x.png",
            )


# ---------------------------------------------------------------------------
# #12–#13 — count dtype through extract / flatten / from_cool / merge
# ---------------------------------------------------------------------------


class TestCountDtypePreservation:
    """``astype(np.int32)`` on a copied ``count`` column."""

    def test_flatten_does_not_truncate_int64_counts(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#12 — ``flatten`` hard-cast ``count`` to int32.

        Harmless at base resolution, silently wrapping for a binnified source:
        ``ops/binnify`` deliberately writes int64 counts.
        """
        src = tmp_path / "wide.metacool"
        shutil.copy(tiny_metacooler, src)
        _forge_wide_count(src)

        out = tmp_path / "flat.cool"
        flatten(src, out)
        assert _max_count(out) == BIG_COUNT

    def test_extract_does_not_truncate_int64_counts(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#12 — same hard cast in ``extract``'s pixel stream."""
        src = tmp_path / "wide.metacool"
        shutil.copy(tiny_metacooler, src)
        _forge_wide_count(src)

        out = tmp_path / "extracted.cool"
        extract(src, genomes=MetaCooler(src).genomes["name"].tolist(), output=out)
        assert _max_count(out) == BIG_COUNT

    def test_from_cool_does_not_truncate_int64_counts(
        self, tmp_path: Path, tiny_metacooler: Path, tiny_contigs_table: Path
    ) -> None:
        """#12 — ``from_cool`` cast to int32 too (and floored float counts to 0)."""
        src = tmp_path / "wide.metacool"
        shutil.copy(tiny_metacooler, src)
        _forge_wide_count(src)
        flat = tmp_path / "wide.cool"
        flatten(src, flat)

        out = tmp_path / "reimported.metacool"
        from_cool(flat, out, contig_genome_map=tiny_contigs_table)
        assert _max_count(out) == BIG_COUNT

    def test_merge_sums_counts_in_int64(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#13 — ``merge_coolers`` was called without ``dtypes={'count': int64}``.

        Before: the summed counts inherited the inputs' int32 and wrapped
        silently, so merging deep libraries corrupted exactly the bins with the
        most signal.
        """
        a, b = tmp_path / "a.metacool", tmp_path / "b.metacool"
        shutil.copy(tiny_metacooler, a)
        shutil.copy(tiny_metacooler, b)
        _forge_wide_count(a, 2_000_000_000)
        _forge_wide_count(b, 2_000_000_000)

        out = tmp_path / "merged.metacool"
        merge([a, b], out)
        assert _max_count(out) == 4_000_000_000

    def test_merge_refuses_a_resolution_mismatch(self, tmp_path: Path) -> None:
        """#13 — only ``contigs``/``genomes`` were compared, never the binning.

        Before: merging files binned at different resolutions concatenated
        pixels addressing unrelated bin IDs into a valid-looking, meaningless
        map.  The contigs here are 200 bp, so 100 bp and 150 bp bins give the
        same *number* of bins per genome — the two files agree on every table
        the old code compared and differ only in resolution.
        """
        fasta = tmp_path / "equal.fa"
        fasta.write_text("".join(f">c{i}\n{'ACGT' * 50}\n" for i in range(4)))
        pairs = _write_pairs(
            tmp_path / "equal.pairs",
            [(f"c{i % 4}", 10 + i, f"c{(i + 1) % 4}", 120 + i) for i in range(20)],
        )
        contigs = tmp_path / "equal_contigs.tsv"
        contigs.write_text(
            "contig_name\tgenome_name\nc0\tMAG_1\nc1\tMAG_1\nc2\tMAG_2\nc3\tMAG_2\n"
        )

        a, b = tmp_path / "b100.metacool", tmp_path / "b150.metacool"
        for path, size in ((a, 100), (b, 150)):
            from_pairs(
                uri=path,
                pairs=pairs,
                binsize=size,
                fasta=fasta,
                contig_genome_map=contigs,
                contig_ordering="genome-then-alpha",
                chunksize=1000,
            )
        assert MetaCooler(a).genomes["n_bins"].tolist() == (
            MetaCooler(b).genomes["n_bins"].tolist()
        )
        with pytest.raises(ValueError, match="Resolution mismatch"):
            merge([a, b], tmp_path / "merged.metacool")

    def test_merge_refuses_differing_bin_tables(
        self, tmp_path: Path, tiny_pairs: Path, tiny_contigs_table: Path
    ) -> None:
        """#13 — two variable-bin files with the same (null) bin-size still differ."""
        lengths = {"contig_A": 1000, "contig_B": 600, "contig_C": 400, "contig_D": 200}
        paths = []
        for tag, cuts in (("x", (0.25, 0.5)), ("y", (0.4, 0.7))):
            rows = []
            for contig, length in lengths.items():
                edges = [0, int(length * cuts[0]), int(length * cuts[1]), length]
                rows += [
                    (contig, edges[i], edges[i + 1]) for i in range(len(edges) - 1)
                ]
            bed = tmp_path / f"bins_{tag}.bed"
            pd.DataFrame(rows).to_csv(bed, sep="\t", header=False, index=False)
            out = tmp_path / f"var_{tag}.metacool"
            _make_metacool(out, tiny_pairs, bed, tiny_contigs_table)
            paths.append(out)

        with pytest.raises(ValueError, match="Bins table mismatch"):
            merge(paths, tmp_path / "merged.metacool")


# ---------------------------------------------------------------------------
# #26–#28, #32–#34 — binnify
# ---------------------------------------------------------------------------


class TestBinnifyGuards:
    """P2 items #26–#28 and #32–#34: where the format could lie about itself."""

    def test_variable_bin_source_refuses_a_finer_target(
        self, tmp_path: Path, tiny_variable_bin_metacooler: Path
    ) -> None:
        """#26 — the "strictly larger" guard skipped variable-bin sources.

        Before: re-binning a fragment-level file to a *finer* size was accepted;
        every source bin landed on the single new bin holding its midpoint and
        all the others stayed empty — a valid-looking, wholly wrong file.
        """
        with pytest.raises(ValueError, match="strictly larger than the effective"):
            binnify(
                tiny_variable_bin_metacooler,
                tmp_path / "finer.metacool",
                binsize=50,
            )

    def test_non_multiple_target_is_refused_with_nearest_multiples(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#27 — a non-multiple target built in a systematic checkerboard.

        Before: 1 kb → 3.5 kb alternated 3 and 4 source bins per output bin (a
        ±14 % oscillation) while the file advertised a uniform 3500 bp.
        """
        with pytest.raises(ValueError, match=r"exact multiple") as exc:
            binnify(tiny_metacooler, tmp_path / "odd.metacool", binsize=350)
        assert "Nearest valid multiples: 300 and 400" in str(exc.value)

    def test_collapsed_output_still_records_the_requested_binsize(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#28 — ``binsize`` was declared, documented and never written.

        Before: when every contig collapsed to one bin, cooler's ``get_binsize``
        returned ``None``, so the output was written ``bin-type="variable"``,
        ``bin-size=null`` — and a later ``binnify`` to a *finer* size then slid
        straight past the resolution guard (feeding #26).
        """
        coarse = tmp_path / "coarse.metacool"
        binnify(tiny_metacooler, coarse, binsize=2000)

        attrs = _file_attrs(coarse)
        assert int(attrs["bin-size"]) == 2000
        assert attrs["bin-type"] == "fixed"
        assert MetaCooler(coarse).binsize == 2000
        # One bin per contig — the case that used to lose the bin size.
        assert int(attrs["nbins"]) == int(attrs["nchroms"])

        with pytest.raises(ValueError, match="strictly larger than current binsize"):
            binnify(coarse, tmp_path / "finer.metacool", binsize=1000)

    def test_genome_level_skips_genomes_without_contigs(
        self, tmp_path: Path, tiny_metacooler: Path, caplog
    ) -> None:
        """#32 — an empty genome became a zero-length contig.

        Before: ``--level genome`` gave every genome ``n_contigs = 1``, so a
        genome with no contigs produced ``length=0``, which ``schema.validate``
        then rejected with "contig lengths must be positive" — after all the
        work was done.
        """
        src = tmp_path / "with_empty.metacool"
        shutil.copy(tiny_metacooler, src)
        _append_empty_genome(src)
        validate(src)

        out = tmp_path / "genome_level.metacool"
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            binnify(src, out, binsize=200, level="genome")

        assert "no contigs" in caplog.text
        mc = MetaCooler(out)
        assert "MAG_EMPTY" not in mc.genomes["name"].tolist()
        validate(out)

    def test_genome_level_writes_contig_length_in_the_schema_dtype(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#32 — ``length`` was written int64 while ``CONTIGS_DTYPE`` declares i4."""
        out = tmp_path / "genome_level.metacool"
        binnify(tiny_metacooler, out, binsize=200, level="genome")
        with h5py.File(out, "r") as f:
            assert f["contigs"]["length"].dtype == np.dtype(CONTIGS_DTYPE["length"])

    def test_genome_level_never_maps_a_bin_into_the_next_genome(self, caplog) -> None:
        """#33 — the ligated remap had no upper clamp.

        Before: ``genome_bin_start[gid] + ligated_pos // binsize`` was unbounded,
        so a stale ``/contigs/length`` (a length shorter than the bins that
        reference it) silently attributed that contig's pixels to the *next*
        genome's first bin — cross-genome contamination with no error.
        """
        from metacooler.ops.binnify import _tile_per_genome_ligated

        # contig_0 declares 100 bp but carries bins out to 1000 bp (stale length).
        contigs = pd.DataFrame(
            {"name": ["c0", "c1"], "length": [100, 500], "genome_id": [0, 1]}
        )
        genomes = pd.DataFrame(
            {
                "name": ["G0", "G1"],
                "length": [100, 500],
                "n_contigs": [1, 1],
                "n_bins": [10, 5],
            }
        )
        old_bins = pd.DataFrame(
            {
                "chrom": ["c0"] * 10 + ["c1"] * 5,
                "start": list(range(0, 1000, 100)) + list(range(0, 500, 100)),
                "end": list(range(100, 1100, 100)) + list(range(100, 600, 100)),
            }
        )
        codes = np.array([0] * 10 + [1] * 5, dtype=np.int64)

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            _new_bins, old_to_new, _contigs, _genomes = _tile_per_genome_ligated(
                contigs, genomes, old_bins, codes, binsize=200
            )

        # G0 owns exactly one output bin (ceil(100/200)); every one of its source
        # bins must land there, never on G1's first bin.
        g0_targets = old_to_new[:10]
        assert set(g0_targets.tolist()) == {0}
        assert old_to_new[10:].min() >= 1
        assert "fell past the end of their genome" in caplog.text

    def test_missing_contig_ordering_attribute_does_not_break_binnify(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#34 — ``mc.info.get("ordering", default)`` never returned the default.

        ``api.info`` always inserts the key (possibly ``None``), so
        ``f.attrs["contig-ordering"] = None`` reached h5py and raised a raw
        ``TypeError`` from deep inside the write.
        """
        src = tmp_path / "no_ordering.metacool"
        shutil.copy(tiny_metacooler, src)
        with h5py.File(src, "r+") as f:
            del f.attrs["contig-ordering"]

        out = tmp_path / "rebinned.metacool"
        binnify(src, out, binsize=200)
        assert _file_attrs(out)["contig-ordering"] == "genome-then-alpha"

    def test_missing_contig_ordering_attribute_does_not_break_merge(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#34 — ``merge`` had the same unguarded read of the same attribute."""
        a, b = tmp_path / "a.metacool", tmp_path / "b.metacool"
        for path in (a, b):
            shutil.copy(tiny_metacooler, path)
            with h5py.File(path, "r+") as f:
                del f.attrs["contig-ordering"]

        out = tmp_path / "merged.metacool"
        merge([a, b], out)
        assert _file_attrs(out)["contig-ordering"] == "genome-then-alpha"


# ---------------------------------------------------------------------------
# #29–#31 — zoomify
# ---------------------------------------------------------------------------


class TestZoomify:
    """P2 items #29–#31: the resolution ladder and where it is built."""

    def test_non_divisible_ladder_is_refused(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#29 — the ladder was never checked for divisibility by the base.

        Before: a 4 kb base with ``--resolutions 1000N`` produced
        5000/10000/20000/50000 — not one of which nests in the base — while the
        docstring promised "bins stay aligned to the base".
        """
        with pytest.raises(ValueError, match="not exact multiples of the base bin"):
            zoomify(
                tiny_metacooler,
                tmp_path / "pyramid.mmetacool",
                resolutions="250,500",
            )

    def test_many_levels_stop_early_without_duplicates(
        self, tmp_path: Path, tiny_metacooler: Path
    ) -> None:
        """#30 — the ``<2 bins`` test ran *after* a full re-binning pass.

        Before: ``--levels 20`` computed and threw away a dozen full passes, and
        once the bin size exceeded the longest contig every level was identical
        (one bin per contig), so the pyramid ended in duplicate entries.
        """
        from metacooler.io.mmetacool import list_entries

        out = tmp_path / "pyramid.mmetacool"
        zoomify(tiny_metacooler, out, levels=20)

        entries = sorted(int(e) for e in list_entries(out))
        assert len(entries) == len(set(entries))
        assert entries == sorted(entries)
        # Base is 100 bp and the longest contig is 1000 bp: doubling collapses
        # every contig to a single bin at 1600, so the ladder stops there.
        assert entries == [100, 200, 400, 800, 1600]
        for res in entries:
            assert MetaCooler(f"{out}::entries/{res}").binsize == res

    def test_temp_dir_is_honoured(
        self, tmp_path: Path, tiny_metacooler: Path, monkeypatch
    ) -> None:
        """#31 — intermediates always went to ``$TMPDIR``, with no way to redirect.

        Peak temporary disk was roughly a second copy of the whole pyramid, on
        whatever filesystem ``$TMPDIR`` happened to be (cooler's own ``zoomify``
        has had a ``--temp-dir`` for years).
        """
        import metacooler.ops.binnify as binnify_mod

        temp_dir = tmp_path / "scratch"
        seen: list[Path] = []
        real_binnify = binnify_mod._binnify

        def _spy(mc, output, *args, **kwargs):
            seen.append(Path(output))
            return real_binnify(mc, output, *args, **kwargs)

        monkeypatch.setattr(binnify_mod, "_binnify", _spy)

        out = tmp_path / "pyramid.mmetacool"
        zoomify(tiny_metacooler, out, levels=3, temp_dir=temp_dir)

        assert temp_dir.is_dir()
        assert seen, "expected at least one coarser level to be built"
        for path in seen:
            assert temp_dir in path.parents, f"{path} was not built under {temp_dir}"


# ---------------------------------------------------------------------------
# Format invariant
# ---------------------------------------------------------------------------


class TestCoolerCompatibilityInvariant:
    """Every ``.metacool`` this module writes must still be a plain cooler."""

    def test_every_writer_produces_a_file_cooler_can_open(
        self,
        tmp_path: Path,
        tiny_metacooler: Path,
        tiny_cool: Path,
        tiny_contigs_table: Path,
    ) -> None:
        """The central format invariant from CLAUDE.md, exercised per writer.

        ``from_pairs``, ``from_cool``, ``binnify`` (both levels) and ``merge``
        each go through their own finalisation path, so each is checked
        separately: a regression in any one of them would leave a file that
        ``MetaCooler`` still opens but ``cooler.Cooler`` does not.
        """
        produced: dict[str, Path] = {"from_pairs": tiny_metacooler}

        produced["from_cool"] = tmp_path / "imported.metacool"
        from_cool(
            tiny_cool, produced["from_cool"], contig_genome_map=tiny_contigs_table
        )

        produced["binnify_contig"] = tmp_path / "rebinned_contig.metacool"
        binnify(tiny_metacooler, produced["binnify_contig"], binsize=200)

        produced["binnify_genome"] = tmp_path / "rebinned_genome.metacool"
        binnify(
            tiny_metacooler, produced["binnify_genome"], binsize=200, level="genome"
        )

        a, b = tmp_path / "a.metacool", tmp_path / "b.metacool"
        shutil.copy(tiny_metacooler, a)
        shutil.copy(tiny_metacooler, b)
        produced["merge"] = tmp_path / "merged.metacool"
        merge([a, b], produced["merge"])

        for label, path in produced.items():
            assert cooler.fileops.is_cooler(str(path)), f"{label} is not a cooler"
            clr = cooler.Cooler(str(path))
            assert len(clr.bins()[:]) == clr.info["nbins"], label
            assert clr.chromnames, label
            validate(path)


# ===========================================================================
# Extra defects found while implementing the review (not review items).
# Both are of the same family as reviewed items, but the review missed them.
# ===========================================================================


class TestDefectsFoundDuringRemediation:
    """Bugs uncovered while landing the review's fixes, and fixed alongside."""

    def test_create_records_the_requested_binsize_when_contigs_hold_one_bin(
        self, tmp_path: Path, tiny_contigs_table: Path
    ):
        """``from_pairs(binsize=N)`` must record ``bin-size=N``, always.

        Same defect as review item #28, but in the *create* path rather than
        ``binnify`` — the review only caught the ``binnify`` half.  cooler
        infers ``bin-size`` from the bin table and yields ``None`` whenever
        every contig holds exactly one bin, so a genuinely fixed-resolution
        file was written as ``bin-type="variable"``/``bin-size=null``.  That in
        turn let a later ``binnify`` to a *finer* resolution slip past the
        "target must be coarser" guard, silently producing a wrong map.
        """
        fasta = tmp_path / "one_bin.fa"
        fasta.write_text("".join(f">contig_{c}\n{'ACGT' * 25}\n" for c in "ABCD"))
        pairs = _write_pairs(
            tmp_path / "one_bin.pairs",
            [("contig_A", 5, "contig_B", 7)] * 5,
        )
        out = tmp_path / "one_bin.metacool"
        from_pairs(
            out,
            pairs=pairs,
            contig_genome_map=tiny_contigs_table,
            fasta=fasta,
            binsize=100,
        )

        attrs = _file_attrs(out)
        assert attrs["bin-size"] == 100
        assert attrs["bin-type"] == "fixed"
        assert MetaCooler(out).binsize == 100

        # The whole point: the coarser-target guard must now engage.
        with pytest.raises(ValueError, match="strictly larger"):
            binnify(out, tmp_path / "finer.metacool", binsize=50)

    def test_merge_reports_a_resolution_mismatch_as_such(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_contigs_table: Path,
        tiny_fasta: Path,
    ):
        """A resolution mismatch must name the resolution, not the genomes table.

        ``/genomes`` carries ``n_bins``, so its hash changes whenever the
        resolution does.  With the genomes comparison ordered first, every
        differing-resolution merge reported ``Genomes table mismatch`` and sent
        the user to inspect a table that was in fact identical apart from a
        derived column.
        """
        a = tmp_path / "a.metacool"
        b = tmp_path / "b.metacool"
        for path, binsize in ((a, 100), (b, 200)):
            from_pairs(
                path,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                binsize=binsize,
                fasta=tiny_fasta,
            )

        with pytest.raises(ValueError, match="Resolution mismatch"):
            merge([a, b], tmp_path / "merged.metacool")


class TestUserTableAlignment:
    """A header/data field-count mismatch must not silently shift columns.

    Found on a real MetaTOR genome table whose header line named 11 columns
    (plus a trailing tab) while every one of its 3203 data rows carried 8
    fields. ``pandas.read_csv`` padded the *trailing* columns, so the values
    were read three columns to the left of the names they belonged to: GC
    content was stored as ``Size``, completeness as ``n_contigs``,
    contamination as ``bins`` and taxonomy as ``GC``, and the three genuinely
    trailing columns came out empty. Nothing warned, and the resulting
    ``.metacool`` was structurally valid.
    """

    def test_genomes_header_wider_than_rows_is_refused(self, tmp_path: Path) -> None:
        """Every row narrower than the header: alignment is unknowable."""
        path = tmp_path / "genomes.tsv"
        path.write_text(
            "genome_name\tSize\tn_contigs\tGC\tCompletion\n"
            "MAG_1\t45.2\t98.6\n"
            "MAG_2\t51.0\t72.4\n"
        )
        with pytest.raises(ValueError, match="misaligned"):
            parse_genomes_table(path)

    def test_contigs_header_wider_than_rows_is_refused(self, tmp_path: Path) -> None:
        """The same guard protects the contig→genome assignment itself."""
        path = tmp_path / "contigs.tsv"
        path.write_text(
            "contig_name\tgenome_name\tgc\tcoverage\nc1\tMAG_1\nc2\tMAG_2\n"
        )
        with pytest.raises(ValueError, match="misaligned"):
            parse_contigs_table(path, fmt="contigs")

    def test_rows_wider_than_header_is_refused(self, tmp_path: Path) -> None:
        """The mirror-image shift: pandas absorbs surplus fields into the index.

        ``read_csv(header=0)`` does not raise on a row with more fields than the
        header — it promotes the extra leading field into the DataFrame index.
        Every label then names the column to its right, so ``contig_name`` held
        the genome and ``genome_id`` held the coverage, and the mapping built
        from them was garbage rather than an error.
        """
        path = tmp_path / "contigs.tsv"
        path.write_text("contig_name\tgenome_name\nc1\tMAG_1\t0.52\nc2\tMAG_2\t0.48\n")
        with pytest.raises(ValueError, match="misaligned"):
            parse_contigs_table(path, fmt="contigs")

    def test_trailing_tab_on_header_is_not_a_column(self, tmp_path: Path) -> None:
        """A trailing tab is a writer artifact, not an ``Unnamed: N`` column.

        Both real MetaTOR genome tables end their header line with a tab, which
        pandas turned into an all-empty ``Unnamed: 8`` column that metacooler
        then wrote out as a junk ``/genomes`` HDF5 dataset.
        """
        path = tmp_path / "genomes.tsv"
        path.write_text("genome_name\tGC\tcompletion\t\nMAG_1\t45.2\t98.6\n")
        df = parse_genomes_table(path)
        assert list(df.columns) == ["genome_id", "GC", "completion"]

    def test_trailing_column_with_data_is_kept(self, tmp_path: Path) -> None:
        """Only *empty* unnamed trailing columns are dropped."""
        path = tmp_path / "genomes.tsv"
        path.write_text("genome_name\tGC\t\nMAG_1\t45.2\tkeep_me\n")
        df = parse_genomes_table(path)
        assert df.shape[1] == 3
        assert df.iloc[0, 2] == "keep_me"

    def test_individually_short_rows_warn_and_fill_empty(
        self, tmp_path: Path, caplog
    ) -> None:
        """Some rows short is a per-row omission, not a systematic shift.

        Alignment is unambiguous here (at least one row fills the header), so
        the table is accepted — but the short rows are reported rather than
        passed off as genuinely empty fields.
        """
        path = tmp_path / "genomes.tsv"
        path.write_text(
            "genome_name\tGC\ttaxonomy\nMAG_1\t45.2\tBacteroidetes\nMAG_2\t51.0\n"
        )
        with caplog.at_level(logging.WARNING, logger="metacooler"):
            df = parse_genomes_table(path)
        assert "fewer than 3" in caplog.text
        assert df["taxonomy"].tolist() == ["Bacteroidetes", ""]

    def test_well_formed_metator_table_is_untouched(self, tmp_path: Path) -> None:
        """The guard must not reject the tables that were always fine."""
        path = tmp_path / "contig_data_final.txt"
        path.write_text(
            "ID\tName\tSize\tFinal_bin\n"
            "1\tContig_1\t836\tMetaTOR_MGE_00001\n"
            "2\tContig_2\t512\tND\n"
        )
        mapping, df = parse_contigs_table(path, fmt="metator")
        assert mapping == {"Contig_1": "MetaTOR_MGE_00001", "Contig_2": "ND"}
        assert df["Size"].tolist() == ["836", "512"]


class TestGenomesTableDivergence:
    """``/genomes`` is derived from the contigs table; report the divergence.

    A ``--genomes`` table only contributes extra metadata columns, left-merged
    onto rows derived from the contig→genome assignments. So a bin named in the
    contigs table but absent from the genomes table still gets a row — with
    computed ``length``/``n_contigs``/``n_bins`` and blank metadata, which reads
    as metacooler having invented a genome — and a genome listed in the genomes
    table with no contigs is dropped outright. Both used to be silent: on a real
    MetaTOR dataset 66 ``MetaTOR_MGE_*`` bins appeared in ``/genomes`` that the
    MAG-only genomes table never listed, with no indication why.
    """

    @staticmethod
    def _create(tmp_path: Path, genomes_table: Path, **kwargs) -> Path:
        out = tmp_path / "out.metacool"
        from_pairs(
            out,
            pairs=kwargs.pop("pairs"),
            contig_genome_map=kwargs.pop("contig_genome_map"),
            binsize=100,
            fasta=kwargs.pop("fasta"),
            genomes_table=genomes_table,
        )
        return out

    def test_bin_missing_from_genomes_table_is_reported(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_contigs_table: Path,
        tiny_fasta: Path,
        caplog,
    ) -> None:
        """MAG_2 is assigned contigs but absent from the genomes table."""
        genomes = tmp_path / "genomes.tsv"
        genomes.write_text("genome_name\tquality\nMAG_1\thigh\n")

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = self._create(
                tmp_path,
                genomes,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                fasta=tiny_fasta,
            )

        assert "no row in the genomes table" in caplog.text
        assert "MAG_2" in caplog.text
        # __unbinned__ is metacooler's own synthetic row, not a missing entry.
        assert "__unbinned__" not in caplog.text

        mc = MetaCooler(out)
        row = mc.genomes[mc.genomes["name"] == "MAG_2"].iloc[0]
        # The row is kept, and its derived columns are real, not invented:
        # contig_C is 400 bp at binsize 100.
        assert int(row["length"]) == 400
        assert int(row["n_contigs"]) == 1
        assert int(row["n_bins"]) == 4
        assert row["quality"] == ""

    def test_genome_without_contigs_is_reported(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_contigs_table: Path,
        tiny_fasta: Path,
        caplog,
    ) -> None:
        """A genomes-table row with no contigs never reaches ``/genomes``."""
        genomes = tmp_path / "genomes.tsv"
        genomes.write_text(
            "genome_name\tquality\nMAG_1\thigh\nMAG_2\tmedium\nMAG_GHOST\tlow\n"
        )

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            out = self._create(
                tmp_path,
                genomes,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                fasta=tiny_fasta,
            )

        assert "no contig assigned to them" in caplog.text
        assert "MAG_GHOST" in caplog.text
        assert "MAG_GHOST" not in MetaCooler(out).genomes["name"].tolist()

    def test_matching_tables_stay_quiet(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_contigs_table: Path,
        tiny_fasta: Path,
        caplog,
    ) -> None:
        """No divergence, no warning — the guard must not cry wolf."""
        genomes = tmp_path / "genomes.tsv"
        genomes.write_text("genome_name\tquality\nMAG_1\thigh\nMAG_2\tmedium\n")

        with caplog.at_level(logging.WARNING, logger="metacooler"):
            self._create(
                tmp_path,
                genomes,
                pairs=tiny_pairs,
                contig_genome_map=tiny_contigs_table,
                fasta=tiny_fasta,
            )

        assert "genomes table" not in caplog.text


class TestBinsBedExtraColumns:
    """A BED with more than three columns must not shift chrom/start/end.

    ``_read_bins_bed`` passed a three-name ``names=`` list, so pandas promoted
    the surplus leading column into the index: ``chrom`` picked up the start
    coordinate, ``start`` the end, and ``end`` the fourth column. Nothing
    downstream could catch it — when ``--bins`` is a BED, contig lengths are
    derived from that same frame, so the shifted table is self-consistent and
    creation went ahead with contigs named after coordinates. It only ever
    failed by luck, when the shifted ``start > end`` happened to trip the
    coordinate validator, reporting numbers that appear nowhere in the input.
    """

    @pytest.mark.parametrize("fourth", ["5", "9999", "bin_name"])
    def test_bed4_uses_the_first_three_columns(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_contigs_table: Path,
        fourth: str,
    ) -> None:
        """A 4th column (score, coverage, or name) is read and discarded."""
        lengths = {"contig_A": 1000, "contig_B": 600, "contig_C": 400, "contig_D": 200}
        bed = tmp_path / "bins4.bed"
        bed.write_text(
            "".join(
                f"{name}\t{start}\t{min(start + 100, length)}\t{fourth}\n"
                for name, length in lengths.items()
                for start in range(0, length, 100)
            )
        )

        out = _make_metacool(
            tmp_path / "bed4.metacool", tiny_pairs, bed, tiny_contigs_table
        )

        mc = MetaCooler(out)
        # Row order is genome-grouped, so compare the mapping, not the sequence.
        assert (
            dict(zip(mc.contigs["name"], mc.contigs["length"], strict=True)) == lengths
        )
        assert mc.bins()[:]["end"].max() == 1000
        validate(out)

    def test_two_column_bed_is_refused_by_name(self, tmp_path: Path) -> None:
        """A BED missing its ``end`` column names the real problem."""
        bed = tmp_path / "bins2.bed"
        bed.write_text("contig_A\t0\ncontig_A\t100\n")
        with pytest.raises(ValueError, match="at least 3 tab-separated columns"):
            _read_bins_bed(bed)
