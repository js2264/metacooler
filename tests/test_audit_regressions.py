"""Regression tests for the audit findings (see AUDIT-REPORT.md).

Each test pins a previously-silent corruption or crash so it cannot come back.
Test names carry the finding id (C1, C2, … / H1, … / M1, …).
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

from metacooler.api import MetaCooler
from metacooler.io.contigs import parse_contigs_table
from metacooler.io.create import from_pairs
from metacooler.schema import MetaCoolerFormatError, validate


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _write_fasta(path: Path, contigs: list[tuple[str, int]]) -> None:
    with open(path, "w") as f:
        for name, length in contigs:
            f.write(f">{name}\n" + "A" * length + "\n")


def _write_pairs(path: Path, lengths: dict[str, int], n: int = 400) -> None:
    rng = np.random.default_rng(0)
    # Only draw positions on contigs with a usable length.
    names = [name for name, ln in lengths.items() if ln > 1]
    rows = ["#columns: readID chrom1 pos1 chrom2 pos2\n"]
    for i in range(n):
        c1 = rng.choice(names)
        c2 = rng.choice(names)
        p1 = int(rng.integers(1, lengths[c1] + 1))
        p2 = int(rng.integers(1, lengths[c2] + 1))
        rows.append(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\n")
    path.write_text("".join(rows))


def _build(
    tmp_path: Path,
    contigs: list[tuple[str, int, str | None]],
    *,
    ordering: str = "genome-then-length",
    binsize: int = 100,
    contig_table: str | None = None,
    contig_table_format: str = "contigs",
    genomes_table: Path | None = None,
) -> Path:
    """Build a small ``.metacool`` from an in-memory contig spec.

    *contigs* is a list of ``(name, length, genome_or_None)`` — ``None`` means
    the contig is unassigned (``__unbinned__``).  *genomes_table* is an
    optional path to a ``--genomes`` metadata table.
    """
    fa = tmp_path / "a.fa"
    _write_fasta(fa, [(n, ln) for n, ln, _ in contigs])
    ct = tmp_path / "ct.tsv"
    if contig_table is None:
        with open(ct, "w") as f:
            f.write("contig_name\tgenome_name\n")
            for n, _ln, g in contigs:
                if g is not None:
                    f.write(f"{n}\t{g}\n")
    else:
        ct.write_text(contig_table)
    lengths = {n: ln for n, ln, _ in contigs}
    pairs = tmp_path / "p.pairs"
    _write_pairs(pairs, lengths)
    out = tmp_path / "o.metacool"
    from_pairs(
        uri=out,
        bins=None,
        binsize=binsize,
        pairs=str(pairs),
        contig_genome_map=str(ct),
        contig_table_format=contig_table_format,
        genomes_table=genomes_table,
        fasta=str(fa),
        contig_ordering=ordering,  # ty: ignore[invalid-argument-type]
        chunksize=37,
    )
    return out


def _genome_offset(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as f:
        return f["indexes/genome_offset"][:]


# ---------------------------------------------------------------------------
# C1 — contig_genome table with a trailing column must not wipe assignments
# ---------------------------------------------------------------------------
class TestC1ContigGenomeTrailingColumn:
    def test_three_column_contig_genome_preserves_assignments(self, tmp_path):
        # metator/awk-style export: contig, genome, size — three columns.
        table = "cA\tMAG_1\t1000\ncB\tMAG_1\t800\ncC\tMAG_2\t600\n"
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_1"), ("cC", 600, "MAG_2")],
            contig_table=table,
            contig_table_format="contig_genome",
        )
        mc = MetaCooler(out)
        # The genome assignment must survive — NOT collapse to __unbinned__.
        assert set(mc.genomes["name"]) == {"MAG_1", "MAG_2"}
        assert "__unbinned__" not in set(mc.genomes["name"])

    def test_parser_reads_first_two_of_many_columns(self, tmp_path):
        p = tmp_path / "t.tsv"
        p.write_text("cA\tMAG_1\t1000\tfoo\ncB\tMAG_2\t800\tbar\n")
        mapping, _ = parse_contigs_table(p, fmt="contig_genome")
        assert mapping == {"cA": "MAG_1", "cB": "MAG_2"}


# ---------------------------------------------------------------------------
# C2 — genome_offset must index the physical contig layout
# ---------------------------------------------------------------------------
class TestC2GenomeOffset:
    def test_unbinned_first_and_offset_is_physical_csr(self, tmp_path):
        out = _build(
            tmp_path,
            [
                ("cA", 1000, "MAG_1"),
                ("cB", 800, "MAG_1"),
                ("cC", 600, "MAG_2"),
                ("cD", 400, None),  # unbinned
            ],
        )
        validate(out)  # strict validate would reject a mis-indexed offset
        mc = MetaCooler(out)
        # __unbinned__ is genome 0 and laid out first.
        assert mc.genomes["name"].iloc[0] == "__unbinned__"
        gid = np.asarray(mc.contigs["genome_id"])
        assert gid[0] == 0
        go = _genome_offset(out)
        # genome_offset must delimit each genome's contiguous contig block.
        for i in range(len(go) - 1):
            assert (gid[go[i] : go[i + 1]] == i).all()
        assert go[0] == 0 and go[-1] == len(gid)


# ---------------------------------------------------------------------------
# C3 / H4 — binnify --level genome must not scramble contig labels or offsets
# ---------------------------------------------------------------------------
class TestBinnifyGenomeLevel:
    def test_genome_level_labels_and_offsets_correct(self, tmp_path):
        from metacooler.ops.binnify import binnify

        src = _build(
            tmp_path,
            [
                ("cA", 1000, "MAG_1"),
                ("cB", 800, "MAG_2"),
                ("cC", 600, "MAG_1"),
                ("cD", 400, None),
            ],
            binsize=100,
        )
        out = tmp_path / "binned.metacool"
        binnify(src, out, binsize=300, level="genome")
        validate(out)  # catches scrambled chrom_offset / out-of-range genome_offset
        mc = MetaCooler(out)
        contigs = mc.contigs
        bins = mc.bins()
        # Every bin must sit within the contig it is labelled with (a scrambled
        # label would put a bin's end past its labelled contig's length).
        lengths = contigs["length"].to_numpy()
        chrom_codes = bins["chrom"].cat.codes.to_numpy()
        assert (bins["end"].to_numpy() <= lengths[chrom_codes]).all()
        # genome_offset stays a valid index into the contigs table.
        go = _genome_offset(out)
        assert go[-1] == len(contigs)
        assert go.max() <= len(contigs)

    def test_genome_level_conserves_counts(self, tmp_path):
        from metacooler.ops.binnify import binnify

        src = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2"), ("cC", 600, None)],
        )
        with h5py.File(src) as f:
            total_before = int(f["pixels/count"][:].sum())
        out = tmp_path / "b.metacool"
        binnify(src, out, binsize=300, level="genome")
        with h5py.File(out) as f:
            total_after = int(f["pixels/count"][:].sum())
        assert total_after == total_before


# ---------------------------------------------------------------------------
# H1 — genome ranges must be correct or refuse on non-contiguous orderings
# ---------------------------------------------------------------------------
class TestH1GenomeRange:
    def test_grouped_ordering_genome_extent_is_intra_genome(self, tmp_path):
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2"), ("cC", 600, "MAG_1")],
            ordering="genome-then-length",
        )
        mc = MetaCooler(out)
        lo, hi = mc.genome_extent("MAG_1")
        bins = mc.bins()
        chrom_codes = bins["chrom"].cat.codes.to_numpy()
        gid = mc.contigs["genome_id"].to_numpy()
        # Every bin in the returned range belongs to MAG_1.
        mag1 = int(mc.genomes.index[mc.genomes["name"] == "MAG_1"][0])
        assert (gid[chrom_codes[lo:hi]] == mag1).all()

    def test_non_contiguous_ordering_raises(self, tmp_path):
        # fasta ordering interleaves MAG_1 (cA, cC) around MAG_2 (cB).
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2"), ("cC", 600, "MAG_1")],
            ordering="fasta",
        )
        mc = MetaCooler(out)
        with pytest.raises(ValueError, match="not contiguous"):
            mc.genome_extent("MAG_1")

    def test_genome_bin_ids_correct_regardless_of_ordering(self, tmp_path):
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2"), ("cC", 600, "MAG_1")],
            ordering="fasta",
        )
        mc = MetaCooler(out)
        ids = mc.genome_bin_ids("MAG_1")
        chrom_codes = mc.bins()["chrom"].cat.codes.to_numpy()
        gid = mc.contigs["genome_id"].to_numpy()
        mag1 = int(mc.genomes.index[mc.genomes["name"] == "MAG_1"][0])
        assert (gid[chrom_codes[ids]] == mag1).all()
        # And it captures ALL such bins.
        assert ids.size == int((gid[chrom_codes] == mag1).sum())


# ---------------------------------------------------------------------------
# H2 — a zero/negative-length contig must not contaminate a neighbour
# ---------------------------------------------------------------------------
class TestH2ZeroLengthContig:
    def test_zero_length_contig_dropped(self, tmp_path):
        # cB has length 0 -> must be dropped, not folded into cA's last bin.
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 0, "MAG_1"), ("cC", 600, "MAG_2")],
        )
        validate(out)
        mc = MetaCooler(out)
        assert "cB" not in set(mc.contigs["name"])
        assert (mc.contigs["length"] > 0).all()

    def test_malformed_chromsize_header_negative_length(self, tmp_path):
        # A negative #chromsize used to produce a length -N contig that stole a
        # neighbour's bin; it must now be dropped.
        pairs = tmp_path / "p.pairs"
        pairs.write_text(
            "#chromsize: cA 1000\n"
            "#chromsize: cNEG -3\n"
            "#columns: readID chrom1 pos1 chrom2 pos2\n"
            "r0\tcA\t10\tcA\t900\n"
            "r1\tcNEG\t1\tcNEG\t1\n"
        )
        ct = tmp_path / "ct.tsv"
        ct.write_text("contig_name\tgenome_name\ncA\tMAG_1\ncNEG\tMAG_2\n")
        out = tmp_path / "o.metacool"
        from_pairs(
            uri=out,
            bins=None,
            binsize=100,
            pairs=str(pairs),
            contig_genome_map=str(ct),
            contig_table_format="contigs",
        )
        validate(out)
        mc = MetaCooler(out)
        assert "cNEG" not in set(mc.contigs["name"])


# ---------------------------------------------------------------------------
# H3 / M6 — duplicate contig / genome names must not explode the tables
# ---------------------------------------------------------------------------
class TestH3H6DuplicateNames:
    def test_duplicate_contig_name_no_phantom_row(self, tmp_path):
        # cA listed twice, with an extra column that triggers the merge path.
        table = (
            "contig_name\tgenome_name\tcov\n"
            "cA\tMAG_1\t10\ncA\tMAG_1\t20\ncB\tMAG_2\t5\n"
        )
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2")],
            contig_table=table,
        )
        validate(out)
        mc = MetaCooler(out)
        names = mc.contigs["name"].tolist()
        assert names.count("cA") == 1  # no phantom duplicate row
        with h5py.File(out) as f:
            assert int(f.attrs["nchroms"]) == len(mc.contigs)

    def test_duplicate_genome_name_deduped(self, tmp_path):
        gt = tmp_path / "g.tsv"
        gt.write_text("genome_name\tphylum\nMAG_1\tX\nMAG_1\tY\n")
        fa = tmp_path / "a.fa"
        _write_fasta(fa, [("cA", 1000), ("cB", 800)])
        ct = tmp_path / "ct.tsv"
        ct.write_text("contig_name\tgenome_name\ncA\tMAG_1\ncB\tMAG_1\n")
        pairs = tmp_path / "p.pairs"
        _write_pairs(pairs, {"cA": 1000, "cB": 800})
        out = tmp_path / "o.metacool"
        from_pairs(
            uri=out,
            bins=None,
            binsize=100,
            pairs=str(pairs),
            contig_genome_map=str(ct),
            genomes_table=str(gt),
            fasta=str(fa),
        )
        validate(out)
        mc = MetaCooler(out)
        assert mc.genomes["name"].tolist().count("MAG_1") == 1


# ---------------------------------------------------------------------------
# H5 — non-ASCII text round-trips wherever it is still allowed
# ---------------------------------------------------------------------------
class TestH5Unicode:
    """H5 originally asked that non-ASCII *names* round-trip; they no longer do.

    Byte-length sizing in ``_str_to_hdf5`` (the actual H5 fix) stored them
    faithfully, but cooler decodes ``/chroms`` — and so ``/contigs`` — as
    ASCII, and genome names reach cooler as chrom names via ``aggregate
    --level genome``.  A file with such a name is therefore unopenable by
    ``cooler.Cooler()``, which is the one thing the format guarantees, so both
    name columns are now refused at write time (see
    ``test_cooler_compat.TestNonAsciiNamesRefused``).  What H5 fixed still
    holds everywhere cooler does not look: the metadata columns of
    ``/genomes``.
    """

    def test_non_ascii_genome_metadata_roundtrips(self, tmp_path):
        genomes_table = tmp_path / "genomes.tsv"
        genomes_table.write_text(
            "genome_name\tphylum\nMAG_1\tBacteroidètes\nMAG_2\tΩ_phylum\n"
        )
        out = _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2")],
            genomes_table=genomes_table,
        )
        validate(out)
        mc = MetaCooler(out)  # reading /genomes must not raise UnicodeDecodeError
        assert set(mc.genomes["phylum"]) == {"Bacteroidètes", "Ω_phylum"}

    def test_non_ascii_genome_name_is_refused(self, tmp_path):
        with pytest.raises(MetaCoolerFormatError, match="not ASCII"):
            _build(tmp_path, [("cA", 1000, "café"), ("cB", 800, "Ω_genome")])

    def test_str_to_hdf5_sizes_by_byte_length(self):
        from metacooler.io._utils import _str_to_hdf5

        arr = _str_to_hdf5(np.array(["café", "Ω"], dtype=object))
        # café is 5 bytes in UTF-8; both must round-trip without truncation.
        assert arr[0].decode("utf-8") == "café"
        assert arr[1].decode("utf-8") == "Ω"


# ---------------------------------------------------------------------------
# M2 — network(level="bin", out=None) must not close stdout
# ---------------------------------------------------------------------------
class TestM2NetworkStdout:
    def test_bin_level_stdout_stays_open(self, tmp_path, capsys):
        from metacooler.io.network import network

        out = _build(tmp_path, [("cA", 500, "MAG_1"), ("cB", 400, "MAG_2")])
        network(out, level="bin", out=None)
        assert not sys.stdout.closed
        # stdout is still usable afterwards
        print("still works")
        assert "still works" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# M3 — an empty (header-only) pairs file yields a valid empty map, not a crash
# ---------------------------------------------------------------------------
class TestM3EmptyPairs:
    def test_header_only_pairs_makes_valid_empty_metacool(self, tmp_path):
        fa = tmp_path / "a.fa"
        _write_fasta(fa, [("cA", 500), ("cB", 400)])
        ct = tmp_path / "ct.tsv"
        ct.write_text("contig_name\tgenome_name\ncA\tMAG_1\ncB\tMAG_2\n")
        pairs = tmp_path / "p.pairs"
        pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2\n")  # no data rows
        out = tmp_path / "o.metacool"
        from_pairs(
            uri=out,
            bins=None,
            binsize=100,
            pairs=str(pairs),
            contig_genome_map=str(ct),
            fasta=str(fa),
        )
        validate(out)
        with h5py.File(out) as f:
            assert int(f.attrs["nnz"]) == 0

    def test_failed_create_leaves_no_partial_file(self, tmp_path):
        # A malformed bins file makes create fail; no partial output must remain.
        fa = tmp_path / "a.fa"
        _write_fasta(fa, [("cA", 500)])
        ct = tmp_path / "ct.tsv"
        ct.write_text("contig_name\tgenome_name\ncA\tMAG_1\n")
        bad_bins = tmp_path / "bins.bed"
        bad_bins.write_text("cA\t100\t50\n")  # end < start
        pairs = tmp_path / "p.pairs"
        _write_pairs(pairs, {"cA": 500})
        out = tmp_path / "o.metacool"
        with pytest.raises(ValueError):
            from_pairs(
                uri=out,
                bins=str(bad_bins),
                pairs=str(pairs),
                contig_genome_map=str(ct),
                fasta=str(fa),
            )
        assert not out.exists()


# ---------------------------------------------------------------------------
# M5 — malformed bin coordinates are rejected
# ---------------------------------------------------------------------------
class TestM5BinCoordValidation:
    @pytest.mark.parametrize(
        "row",
        ["cA\t-10\t100\n", "cA\t100\t50\n"],  # negative start, end < start
    )
    def test_bad_bins_rejected(self, tmp_path, row):
        fa = tmp_path / "a.fa"
        _write_fasta(fa, [("cA", 500)])
        ct = tmp_path / "ct.tsv"
        ct.write_text("contig_name\tgenome_name\ncA\tMAG_1\n")
        bins = tmp_path / "bins.bed"
        bins.write_text(row)
        pairs = tmp_path / "p.pairs"
        _write_pairs(pairs, {"cA": 500})
        with pytest.raises(ValueError):
            from_pairs(
                uri=tmp_path / "o.metacool",
                bins=str(bins),
                pairs=str(pairs),
                contig_genome_map=str(ct),
                fasta=str(fa),
            )


# ---------------------------------------------------------------------------
# M1 / H6 — validate() catches semantic corruption
# ---------------------------------------------------------------------------
class TestValidateCatchesCorruption:
    @pytest.fixture
    def good(self, tmp_path):
        return _build(
            tmp_path,
            [("cA", 1000, "MAG_1"), ("cB", 800, "MAG_2"), ("cC", 600, None)],
        )

    def test_missing_chroms_softlink(self, good):
        with h5py.File(good, "a") as f:
            del f["chroms"]
        with pytest.raises(MetaCoolerFormatError, match="soft link"):
            validate(good)

    def test_corrupted_genome_offset(self, good):
        with h5py.File(good, "a") as f:
            go = f["indexes/genome_offset"][:]
            go[1] += 1
            del f["indexes/genome_offset"]
            f["indexes"].create_dataset("genome_offset", data=go)
        with pytest.raises(MetaCoolerFormatError, match="genome_offset"):
            validate(good)

    def test_bins_chrom_out_of_range(self, good):
        with h5py.File(good, "a") as f:
            chrom = f["bins/chrom"][:]
            chrom[0] = 9999
            f["bins/chrom"][...] = chrom
        with pytest.raises(MetaCoolerFormatError, match="bins/chrom"):
            validate(good)

    def test_duplicate_contig_name_detected(self, good):
        with h5py.File(good, "a") as f:
            names = f["contigs/name"][:]
            names[1] = names[0]
            f["contigs/name"][...] = names
        with pytest.raises(MetaCoolerFormatError, match="contig names are not unique"):
            validate(good)

    def test_wrong_genome_n_contigs(self, good):
        with h5py.File(good, "a") as f:
            nc = f["genomes/n_contigs"][:]
            nc[0] += 5
            f["genomes/n_contigs"][...] = nc
        with pytest.raises(MetaCoolerFormatError, match="n_contigs"):
            validate(good)
