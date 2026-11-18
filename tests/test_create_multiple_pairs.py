"""``from_pairs`` / ``create pairs`` with several pairs files, read in one pass.

The contract: N libraries of one assembly give exactly the map, QC block,
``/genomes`` ``qc_*`` columns, P(s) table and drop report of the one file that
concatenates them (headers united).  ``merge`` of per-library maps gives the
same map and QC block (``test_merge_qc.py``), but has no P(s) table to sum.
"""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat
from click.testing import CliRunner

import metacooler.io.create as create_mod
from metacooler.cli.main import cli
from metacooler.io.create import (
    _merge_pairs_headers,
    _pairs_paths,
    _read_contig_lengths,
    from_pairs,
)
from metacooler.ops.ps import contact_probability
from metacooler.ps import ps_table_path

LENGTHS = {"contig_A": 1000, "contig_B": 600, "contig_C": 400, "contig_D": 200}
COLUMNS = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _header(lengths: dict[str, int] | None) -> str:
    lines = ["## pairs format v1.0\n"]
    lines += [f"#chromsize: {c} {n}\n" for c, n in (lengths or {}).items()]
    lines.append(COLUMNS)
    return "".join(lines)


def _body(seed: int, n: int, tag: str) -> str:
    """*n* random contacts, with an unknown contig and an out-of-range pair."""
    rng = np.random.default_rng(seed)
    names = list(LENGTHS)
    rows = []
    for i in range(n):
        c1, c2 = (str(c) for c in rng.choice(names, 2))
        p1 = int(rng.integers(1, LENGTHS[c1] + 1))
        p2 = int(rng.integers(1, LENGTHS[c2] + 1))
        rows.append(f"{tag}{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")
    rows[3] = f"{tag}_bad\tcontig_Z_{tag}\t5\tcontig_A\t10\t+\t+\n"
    rows[5] = f"{tag}_oob\tcontig_D\t9999\tcontig_A\t10\t+\t+\n"
    return "".join(rows)


def _write(path: Path, text: str) -> Path:
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.fixture
def libraries(tmp_path: Path) -> tuple[list[Path], Path]:
    """Two libraries (one gzipped) of one assembly, and their concatenation."""
    body_a, body_b = _body(1, 150, "a"), _body(2, 120, "b")
    a = _write(tmp_path / "lib_a.pairs", _header(LENGTHS) + body_a)
    b = _write(tmp_path / "lib_b.pairs.gz", _header(LENGTHS) + body_b)
    cat = _write(tmp_path / "cat.pairs", _header(LENGTHS) + body_a + body_b)
    return [a, b], cat


def _tables(path: Path) -> dict:
    """Every dataset of the file, plus its metadata minus ``commands``."""
    out: dict = {}
    with h5py.File(path, "r") as f:

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                out[name] = obj[()].tolist()

        f.visititems(visit)
        meta = json.loads(f.attrs["metadata"])
        meta.pop("commands", None)
        out["metadata"] = meta
    return out


def _ps_bytes(out: Path) -> bytes:
    return out.with_suffix(".ps.tsv").read_bytes()


def _create(out: Path, pairs, contigs: Path, **kwargs) -> Path:
    kwargs.setdefault("binsize", 100)
    from_pairs(
        uri=out,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_table_format="contigs",
        qc_distance_threshold=100,
        **kwargs,
    )
    return out


def _discarded(caplog) -> list[str]:
    return [
        r.getMessage() for r in caplog.records if r.getMessage().startswith("Discarded")
    ]


@pytest.fixture
def small_batches(monkeypatch):
    """Re-cut every file's pairs stream into 7-row batches.

    With ``chunksize=50`` the coalesced chunks then straddle the boundary
    between two files, and several batches are in flight under ``nproc > 1``.
    """
    real = create_mod._stream_pairs

    def rechunked(*args, **kwargs):
        for chunk in real(*args, **kwargs):
            for lo in range(0, len(chunk), 7):
                yield chunk.iloc[lo : lo + 7]

    monkeypatch.setattr(create_mod, "_stream_pairs", rechunked)


# ----------------------------------------------------------------------
# Several files == their concatenation
# ----------------------------------------------------------------------
class TestSeveralFilesMatchTheirConcatenation:
    @pytest.mark.parametrize("sort_in_memory", [None, True, False])
    def test_map_qc_ps_and_drop_report(
        self, tmp_path, libraries, tiny_contigs_table, caplog, sort_in_memory
    ):
        files, cat = libraries
        caplog.set_level(logging.WARNING, logger="metacooler")
        many = _create(
            tmp_path / "many.metacool",
            files,
            tiny_contigs_table,
            sort_in_memory=sort_in_memory,
        )
        reported_many = _discarded(caplog)
        caplog.clear()
        one = _create(
            tmp_path / "one.metacool",
            cat,
            tiny_contigs_table,
            sort_in_memory=sort_in_memory,
        )
        reported_one = _discarded(caplog)

        tables_many, tables_one = _tables(many), _tables(one)
        assert tables_many == tables_one
        # Spelled out: the pixels and the QC block, which a merge of
        # per-library maps gives too, and the P(s) table, which it cannot.
        for column in ("pixels/bin1_id", "pixels/bin2_id", "pixels/count"):
            assert tables_many[column] == tables_one[column]
        qc_columns = [k for k in tables_many if k.startswith("genomes/qc_")]
        assert len(qc_columns) == 5
        qc = tables_many["metadata"]["qc"]
        assert qc == tables_one["metadata"]["qc"]
        assert qc["n_pairs_read"] == 270
        assert qc["n_pairs_discarded"] == 4
        assert qc["intra_contig_short"] > 0
        assert qc["intra_contig_long"] > 0
        assert _ps_bytes(many) == _ps_bytes(one)
        assert reported_many == reported_one
        assert "Discarded 4 of 270 pair(s)" in reported_many[0]
        assert "contig_Z_a" in reported_many[0]
        assert "contig_Z_b" in reported_many[0]

    @pytest.mark.usefixtures("small_batches")
    def test_chunks_straddling_files(self, tmp_path, libraries, tiny_contigs_table):
        files, cat = libraries
        many = _create(
            tmp_path / "many.metacool", files, tiny_contigs_table, chunksize=50
        )
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table, chunksize=50)
        assert _tables(many) == _tables(one)
        assert _ps_bytes(many) == _ps_bytes(one)

    def test_bins_bed_and_fasta_order(
        self, tmp_path, libraries, tiny_contigs_table, tiny_bins_bed
    ):
        files, cat = libraries
        kwargs = {"bins": tiny_bins_bed, "binsize": None, "contig_ordering": "fasta"}
        many = _create(tmp_path / "many.metacool", files, tiny_contigs_table, **kwargs)
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table, **kwargs)
        assert _tables(many) == _tables(one)
        assert _ps_bytes(many) == _ps_bytes(one)

    def test_ps_pairs_on_the_concatenation_reproduces_the_sidecar(
        self, tmp_path, libraries, tiny_contigs_table
    ):
        files, cat = libraries
        many = _create(tmp_path / "many.metacool", files, tiny_contigs_table)
        sidecar = pd.read_csv(ps_table_path(many), sep="\t")
        recomputed = contact_probability(many, pairs=cat, by="genome")
        pd.testing.assert_frame_equal(sidecar, recomputed)

    def test_single_path_in_a_list(self, tmp_path, libraries, tiny_contigs_table):
        _, cat = libraries
        listed = _create(tmp_path / "listed.metacool", [cat], tiny_contigs_table)
        bare = _create(tmp_path / "bare.metacool", str(cat), tiny_contigs_table)
        assert _tables(listed) == _tables(bare)

    def test_header_only_file_among_others(
        self, tmp_path, libraries, tiny_contigs_table, caplog
    ):
        files, cat = libraries
        empty = _write(tmp_path / "empty.pairs", _header(LENGTHS))
        caplog.set_level(logging.WARNING, logger="metacooler")
        many = _create(
            tmp_path / "many.metacool",
            [files[0], empty, files[1]],
            tiny_contigs_table,
        )
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table)
        assert _tables(many) == _tables(one)
        assert any("has no data rows" in r.getMessage() for r in caplog.records)


# ----------------------------------------------------------------------
# #chromsize: headers across files
# ----------------------------------------------------------------------
class TestChromsizeHeaders:
    def test_conflicting_lengths_name_both_files(self, tmp_path, tiny_contigs_table):
        a = _write(tmp_path / "a.pairs", _header(LENGTHS) + _body(1, 20, "a"))
        other = {**LENGTHS, "contig_B": 650, "contig_C": 450}
        b = _write(tmp_path / "b.pairs", _header(other) + _body(2, 20, "b"))
        with pytest.raises(ValueError, match="different '#chromsize:' lengths") as exc:
            _create(tmp_path / "x.metacool", [a, b], tiny_contigs_table)
        message = str(exc.value)
        assert "2 contig(s)" in message
        assert "'contig_B': 600 in" in message
        assert str(a) in message
        assert f"but 650 in {str(b)!r}" in message
        assert not (tmp_path / "x.metacool").exists()

    def test_conflict_found_through_fasta_ordering(
        self, tmp_path, tiny_contigs_table, tiny_bins_bed
    ):
        # Bins from a BED file, order from the pairs headers: the headers are
        # read for the order alone, and are still checked.
        a = _write(tmp_path / "a.pairs", _header(LENGTHS) + _body(1, 20, "a"))
        b = _write(
            tmp_path / "b.pairs",
            _header({**LENGTHS, "contig_D": 1}) + _body(2, 20, "b"),
        )
        with pytest.raises(ValueError, match="'contig_D': 200 in"):
            _create(
                tmp_path / "x.metacool",
                [a, b],
                tiny_contigs_table,
                bins=tiny_bins_bed,
                binsize=None,
                contig_ordering="fasta",
            )

    @pytest.mark.parametrize("headed_first", [True, False])
    def test_file_without_header(self, tmp_path, tiny_contigs_table, headed_first):
        body_a, body_b = _body(1, 80, "a"), _body(2, 60, "b")
        headed = _write(tmp_path / "headed.pairs", _header(LENGTHS) + body_a)
        bare = _write(tmp_path / "bare.pairs", _header(None) + body_b)
        cat = _write(tmp_path / "cat.pairs", _header(LENGTHS) + body_a + body_b)
        files = [headed, bare] if headed_first else [bare, headed]
        many = _create(tmp_path / "many.metacool", files, tiny_contigs_table)
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table)
        assert _tables(many) == _tables(one)
        assert _ps_bytes(many) == _ps_bytes(one)

    def test_partial_headers_are_united_in_file_order(
        self, tmp_path, tiny_contigs_table
    ):
        first = {"contig_C": 400, "contig_A": 1000}
        second = {"contig_A": 1000, "contig_D": 200, "contig_B": 600}
        united = {"contig_C": 400, "contig_A": 1000, "contig_D": 200, "contig_B": 600}
        body_a, body_b = _body(1, 80, "a"), _body(2, 60, "b")
        a = _write(tmp_path / "a.pairs", _header(first) + body_a)
        b = _write(tmp_path / "b.pairs", _header(second) + body_b)
        cat = _write(tmp_path / "cat.pairs", _header(united) + body_a + body_b)
        assert _merge_pairs_headers([a, b]) == united
        assert list(_merge_pairs_headers([a, b]) or {}) == list(united)
        kwargs = {"contig_ordering": "fasta"}
        many = _create(tmp_path / "many.metacool", [a, b], tiny_contigs_table, **kwargs)
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table, **kwargs)
        tables = _tables(many)
        assert tables == _tables(one)
        assert [n.decode() for n in tables["contigs/name"]] == list(united)

    def test_no_header_anywhere_needs_a_fasta(
        self, tmp_path, tiny_contigs_table, tiny_fasta
    ):
        body_a, body_b = _body(1, 80, "a"), _body(2, 60, "b")
        a = _write(tmp_path / "a.pairs", _header(None) + body_a)
        b = _write(tmp_path / "b.pairs", _header(None) + body_b)
        assert _merge_pairs_headers([a, b]) is None
        with pytest.raises(ValueError, match="Contig lengths required"):
            _read_contig_lengths(pairs=[a, b])
        with pytest.raises(ValueError, match="Contig lengths required"):
            _create(tmp_path / "x.metacool", [a, b], tiny_contigs_table)
        cat = _write(tmp_path / "cat.pairs", _header(None) + body_a + body_b)
        many = _create(
            tmp_path / "many.metacool", [a, b], tiny_contigs_table, fasta=tiny_fasta
        )
        one = _create(
            tmp_path / "one.metacool", cat, tiny_contigs_table, fasta=tiny_fasta
        )
        assert _tables(many) == _tables(one)

    def test_identical_headers_take_the_fast_path(self, tmp_path, libraries):
        files, _ = libraries
        assert _merge_pairs_headers(files) == LENGTHS
        assert _read_contig_lengths(pairs=files) == LENGTHS


# ----------------------------------------------------------------------
# Which inputs are refused
# ----------------------------------------------------------------------
class TestPairsPaths:
    def test_single_path_forms(self, tmp_path):
        p = tmp_path / "x.pairs"
        assert _pairs_paths(p) == [p]
        assert _pairs_paths(str(p)) == [p]
        assert _pairs_paths((str(p), tmp_path / "y.pairs")) == [
            p,
            tmp_path / "y.pairs",
        ]

    def test_none_and_empty_are_refused(self):
        with pytest.raises(ValueError, match="pairs file path required"):
            _pairs_paths(None)
        with pytest.raises(ValueError, match="empty list"):
            _pairs_paths([])

    def test_same_file_twice_is_refused(self, tmp_path):
        p = tmp_path / "x.pairs"
        p.write_text(COLUMNS)
        (tmp_path / "sub").mkdir()
        alias = tmp_path / "sub" / ".." / "x.pairs"
        with pytest.raises(ValueError, match="given twice \\(also as"):
            _pairs_paths([p, tmp_path / "y.pairs", alias])

    def test_same_path_object_and_symlink_are_refused(self, tmp_path):
        p = tmp_path / "x.pairs"
        p.write_text(COLUMNS)
        # One object listed twice: detected by resolved path, not by identity.
        with pytest.raises(ValueError, match="is given twice; every one"):
            _pairs_paths([p, p])
        link = tmp_path / "link.pairs"
        link.symlink_to(p)
        with pytest.raises(ValueError, match="'.*link.pairs' is given twice"):
            _pairs_paths([p, link])

    def test_from_pairs_refuses_the_same_file_twice(
        self, tmp_path, libraries, tiny_contigs_table
    ):
        files, _ = libraries
        with pytest.raises(ValueError, match="is given twice; every one"):
            _create(tmp_path / "x.metacool", [files[0], files[0]], tiny_contigs_table)

    def test_missing_file_is_refused_before_any_is_read(
        self, tmp_path, libraries, tiny_contigs_table, monkeypatch
    ):
        files, _ = libraries

        def never(*args, **kwargs):
            raise AssertionError("a pairs file was read")

        monkeypatch.setattr(create_mod, "_stream_pairs", never)
        monkeypatch.setattr(create_mod, "_parse_pairs_header", never)
        missing = tmp_path / "lib_c.pairs.gz"
        with pytest.raises(FileNotFoundError, match="lib_c.pairs.gz"):
            _create(tmp_path / "x.metacool", [*files, missing], tiny_contigs_table)


# ----------------------------------------------------------------------
# --nproc
# ----------------------------------------------------------------------
@pytest.mark.usefixtures("small_batches")
def test_nproc_is_invisible_across_files(
    tmp_path, libraries, tiny_contigs_table, caplog
):
    files, _ = libraries
    caplog.set_level(logging.WARNING, logger="metacooler")
    tables, ps, reported = {}, {}, {}
    for nproc in (1, 2):
        caplog.clear()
        out = _create(
            tmp_path / f"n{nproc}.metacool",
            files,
            tiny_contigs_table,
            chunksize=50,
            nproc=nproc,
        )
        tables[nproc], ps[nproc] = _tables(out), _ps_bytes(out)
        reported[nproc] = _discarded(caplog)
    assert tables[1] == tables[2]
    assert ps[1] == ps[2]
    assert reported[1] == reported[2] != []


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
class TestCLI:
    def _args(self, out: Path, files: list[Path], contigs: Path) -> list[str]:
        args = ["create", "pairs", "-o", str(out), "--bins", "100"]
        for f in files:
            args += ["--pairs", str(f)]
        return [*args, "--contigs", str(contigs), "--qc-threshold", "100"]

    def test_repeated_pairs(self, tmp_path, libraries, tiny_contigs_table):
        files, cat = libraries
        out = tmp_path / "cli.metacool"
        result = CliRunner().invoke(cli, self._args(out, files, tiny_contigs_table))
        assert result.exit_code == 0, result.output
        one = _create(tmp_path / "one.metacool", cat, tiny_contigs_table)
        assert _tables(out) == _tables(one)
        assert _ps_bytes(out) == _ps_bytes(one)
        with h5py.File(out, "r") as f:
            (command,) = json.loads(f.attrs["metadata"])["commands"]
        assert f"--pairs {files[0]} --pairs {files[1]}" in command

    def test_repeated_pairs_with_nproc(self, tmp_path, libraries, tiny_contigs_table):
        files, _ = libraries
        outs = []
        for nproc in ("1", "2"):
            out = tmp_path / f"n{nproc}.metacool"
            args = [*self._args(out, files, tiny_contigs_table), "--nproc", nproc]
            result = CliRunner().invoke(cli, args)
            assert result.exit_code == 0, result.output
            outs.append(out)
        assert _tables(outs[0]) == _tables(outs[1])
        assert _ps_bytes(outs[0]) == _ps_bytes(outs[1])

    def test_same_file_twice_is_a_clean_error(
        self, tmp_path, libraries, tiny_contigs_table
    ):
        files, _ = libraries
        out = tmp_path / "x.metacool"
        args = self._args(out, [files[0], files[0]], tiny_contigs_table)
        result = CliRunner().invoke(cli, args)
        assert result.exit_code != 0
        assert "is given twice; every one" in flat(result.output)
        assert result.exception is None or isinstance(result.exception, SystemExit)

    def test_output_cannot_be_one_of_the_pairs_files(
        self, tmp_path, libraries, tiny_contigs_table
    ):
        files, _ = libraries
        before = files[1].read_bytes()
        args = [*self._args(files[1], files, tiny_contigs_table), "--force"]
        result = CliRunner().invoke(cli, args)
        assert result.exit_code != 0
        assert "same file as the input" in flat(result.output)
        assert files[1].read_bytes() == before

    def test_conflicting_headers_are_a_clean_error(self, tmp_path, tiny_contigs_table):
        a = _write(tmp_path / "a.pairs", _header(LENGTHS) + _body(1, 20, "a"))
        b = _write(
            tmp_path / "b.pairs",
            _header({**LENGTHS, "contig_A": 999}) + _body(2, 20, "b"),
        )
        out = tmp_path / "x.metacool"
        result = CliRunner().invoke(cli, self._args(out, [a, b], tiny_contigs_table))
        assert result.exit_code != 0
        assert "aligned to the same assembly" in flat(result.output)
