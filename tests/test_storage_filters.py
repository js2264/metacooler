"""Every dataset metacooler writes itself is compressed like cooler's."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import h5py
import numpy as np
import pytest

from metacooler.io._utils import H5_STORAGE, write_column
from metacooler.schema import validate

#: The datasets metacooler (not cooler) writes into every .metacool.
OWN_DATASETS = ("bins/chrom", "indexes/chrom_offset", "indexes/genome_offset")


def _assert_compressed(path: Path, *, tables: bool = True) -> None:
    with h5py.File(path, "r") as f:
        names = list(OWN_DATASETS)
        if tables:
            names += [f"genomes/{k}" for k in f["genomes"]]
            names += [f"contigs/{k}" for k in f["contigs"]]
        for name in names:
            ds = f[name]
            assert ds.compression == "gzip", name
            assert ds.compression_opts == 6, name
            assert ds.shuffle, name


def test_write_column_round_trips(tmp_path: Path):
    with h5py.File(tmp_path / "x.h5", "w") as f:
        a = write_column(f, "text", np.array(["a", "bcd"], dtype=object))
        b = write_column(f, "ints", [1, 2, 3], dtype=np.int64)
        c = write_column(f, "empty", np.array([], dtype=np.int32))
        assert a[()].tolist() == [b"a", b"bcd"]
        assert b.dtype == np.int64 and b[()].tolist() == [1, 2, 3]
        assert len(c) == 0
        for ds in (a, b, c):
            assert ds.compression == H5_STORAGE["compression"]


def test_create_compresses_its_tables(tiny_metacooler: Path):
    _assert_compressed(tiny_metacooler)
    validate(tiny_metacooler)


def test_binnify_compresses_its_tables(tiny_metacooler: Path, tmp_path: Path):
    from metacooler.ops.binnify import binnify

    out = tmp_path / "coarse.metacool"
    binnify(tiny_metacooler, out, binsize=200)
    _assert_compressed(out)


def test_merge_compresses_bins_chrom(tiny_metacooler: Path, tmp_path: Path):
    from metacooler.io.merge import merge

    out = tmp_path / "merged.metacool"
    merge([tiny_metacooler, tiny_metacooler], out)
    _assert_compressed(out)
    validate(out)


@pytest.mark.parametrize("layout", ["keep", "auto"])
def test_reassign_compresses_what_it_rewrites(
    tiny_metacooler: Path, tmp_path: Path, layout: Literal["keep", "auto"]
):
    from metacooler.io.reassign import reassign

    table = tmp_path / "t.tsv"
    # contig_C moves from MAG_2 to MAG_1: under ``auto`` the genome-grouped
    # layout survives (relabel path); ``keep`` exercises the same writer.
    table.write_text("contig_name\tgenome_name\ncontig_C\tMAG_1\n")
    out = tmp_path / f"{layout}.metacool"
    reassign(tiny_metacooler, table, out, layout=layout)
    with h5py.File(out, "r") as f:
        for name in ("contigs/genome_id", "indexes/genome_offset"):
            assert f[name].compression == "gzip", name
        for name in f["genomes"]:
            assert f["genomes"][name].compression == "gzip", name
    validate(out)


def test_write_column_keeps_encoded_bytes(tmp_path: Path):
    with h5py.File(tmp_path / "x.h5", "w") as f:
        ds = write_column(f, "names", np.array([b"k141_1", b"ab"]))
        assert ds.dtype == np.dtype("S6")
        assert ds[()].tolist() == [b"k141_1", b"ab"]


class TestFastTextChecks:
    """The vectorised paths answer exactly what the per-value paths did."""

    @pytest.mark.parametrize(
        "values, expected",
        [
            (np.array(["a", "b", None, float("nan")], dtype=object), []),
            (np.array(["a", "é", "é", "x"], dtype=object), ["é"]),
            (np.array(["a", "é".encode()], dtype=object), ["é"]),
            (np.array(["plain", "ascii"]), []),
            (np.array(["naïve"]), ["naïve"]),
            (np.array([], dtype=object), []),
            (np.array(["", ""], dtype="U"), []),
        ],
    )
    def test_non_ascii_values(self, values, expected):
        from metacooler.schema import non_ascii_values

        assert non_ascii_values(values) == expected

    @pytest.mark.parametrize(
        "values, dup",
        [
            (np.array([b"a", b"b", b"c"]), False),
            (np.array([b"a", b"b", b"a"]), True),
            (np.array([b"a"]), False),
            (np.array([], dtype="S1"), False),
        ],
    )
    def test_has_duplicates(self, values, dup):
        from metacooler.schema import _has_duplicates

        assert _has_duplicates(values) is dup


@pytest.mark.parametrize(
    "values, expected",
    [
        (np.array(["", ""], dtype=object), np.array([b"", b""], dtype="S1")),
        (np.array([None], dtype=object), np.array([b""], dtype="S1")),
        (
            np.array([np.float64(2.0), True, np.int64(3)], dtype=object),
            np.array([b"2.0", b"True", b"3"], dtype="S4"),
        ),
        (
            np.array(["é", "ab"], dtype=object),
            np.array(["é".encode(), b"ab"], dtype="S2"),
        ),
        (
            np.array([5, "a", 1.5, None, float("nan"), "bb"], dtype=object),
            np.array([b"5", b"a", b"1.5", b"", b"", b"bb"], dtype="S3"),
        ),
        (np.array([], dtype=object), np.array([], dtype="S1")),
        (np.array(["u", "vv"]), np.array([b"u", b"vv"], dtype="S2")),
        (np.array([b"k1", b"k22"]), np.array([b"k1", b"k22"], dtype="S3")),
    ],
)
def test_str_to_hdf5(values, expected):
    from metacooler.io._utils import _str_to_hdf5

    got = _str_to_hdf5(values)
    assert got.dtype == expected.dtype
    assert got.tolist() == expected.tolist()


class TestVectorisedBuilders:
    def test_empty_contig_set(self):
        from metacooler.io.create import _build_contig_order

        assert _build_contig_order({}, {}, ordering="genome-then-length") == []

    def test_contig_without_genome_and_no_unbinned_row(self):
        import pandas as pd

        from metacooler.io.create import _build_contigs_table

        genomes = pd.DataFrame({"name": ["G1"], "length": [10], "n_contigs": [1]})
        with pytest.raises(ValueError, match="'c2' has no genome assignment"):
            _build_contigs_table(
                ["c1", "c2"], {"c1": 5, "c2": 5}, {"c1": "G1"}, genomes
            )

    def test_zero_width_text_is_ascii(self):
        from metacooler.schema import non_ascii_values

        assert non_ascii_values(np.zeros(2, dtype="U0")) == []

    @pytest.mark.parametrize(
        "attr, value, message",
        [
            ("metacooler-format-version", 99, "only understands version"),
            ("storage-mode", "diagonal", "storage-mode"),
        ],
    )
    def test_open_refuses_bad_attributes(self, tiny_metacooler, attr, value, message):
        from metacooler.schema import MetaCoolerFormatError

        with h5py.File(tiny_metacooler, "r+") as f:
            f.attrs[attr] = value
        with pytest.raises(MetaCoolerFormatError, match=message):
            validate(tiny_metacooler, full=False)
