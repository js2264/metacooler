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
