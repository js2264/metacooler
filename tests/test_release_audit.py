"""Regression tests from the pre-release audit.

Every test in this module pins down a defect found while auditing the package
for its first public release. As in ``test_regressions.py``, each docstring
states what went wrong before the fix, and the tests assert on observable
behaviour (file contents, exit codes, raised messages) rather than on how the
fix is implemented.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Literal

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner
from scipy import sparse

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io._utils import (
    check_output_path,
    chrom_codes,
    cooler_writer_module,
)
from metacooler.io.create import from_cool, from_pairs
from metacooler.io.merge import merge
from metacooler.io.mmetacool import create_mmetacool
from metacooler.ops.aggregate import (
    _aggregate_sparse,
    _write_aggregated_cool,
    aggregate,
)
from metacooler.ops.binnify import binnify
from metacooler.ops.ps import contact_probability
from metacooler.ops.range import parse_metacooler_range
from metacooler.ops.zoomify import zoomify
from metacooler.ps import bin_pair_denominators, distance_edges
from metacooler.schema import MetaCoolerFormatError, list_metacoolers, validate

Level = Literal["contig", "genome"]


def _run(*args: object):
    return CliRunner().invoke(cli, [str(a) for a in args])


def _replica(path: Path, tmp_path: Path) -> Path:
    """A byte-identical second replicate: merge refuses the *same* file twice."""
    copy = tmp_path / f"replica_{path.name}"
    shutil.copy(path, copy)
    return copy


class TestOutputNeverDestroysInput:
    """``-o <input> -f`` deleted the input, then failed with "File not found".

    The ``--force`` guard unlinked an existing output before the command ran,
    and nothing compared the output against the inputs, so naming the input as
    the output destroyed it.
    """

    @pytest.mark.parametrize(
        "command",
        [
            ["binnify", "{src}", "--binsize", "200"],
            ["merge", "{src}", "{src}"],
            ["zoomify", "{src}", "--levels", "2"],
            ["aggregate", "{src}", "--level", "genome"],
            ["flatten", "{src}"],
            ["extract", "{src}", "--genomes", "MAG_1"],
            ["dump", "-t", "contigs", "{src}"],
            ["network", "{src}"],
            ["qc", "{src}", "--table"],
        ],
        ids=lambda c: c[0],
    )
    def test_input_survives_being_named_as_output(
        self, tmp_path: Path, tiny_metacooler: Path, command: list[str]
    ):
        src = tmp_path / "input.metacool"
        shutil.copy(tiny_metacooler, src)
        before = src.read_bytes()
        args = [a.format(src=src) for a in command] + ["-o", str(src), "-f"]
        result = _run(*args)
        assert result.exit_code != 0
        assert "same file as the input" in _flat(result.output)
        assert src.read_bytes() == before

    def test_entry_uri_input_is_compared_by_file(
        self, tmp_path: Path, tiny_mmetacool: Path
    ):
        """An ``::entries/x`` URI names the same file as its bare path."""
        before = tiny_mmetacool.read_bytes()
        result = _run(
            "binnify",
            f"{tiny_mmetacool}::entries/entry1",
            "--binsize",
            "200",
            "-o",
            tiny_mmetacool,
            "-f",
        )
        assert result.exit_code != 0
        assert "same file as the input" in _flat(result.output)
        assert tiny_mmetacool.read_bytes() == before

    def test_api_guard_compares_inputs(self, tmp_path: Path):
        f = tmp_path / "x.metacool"
        f.write_bytes(b"x")
        with pytest.raises(ValueError, match="same file as the input"):
            check_output_path(f, force=True, inputs=[str(f)])
        assert f.exists()

    def test_api_guard_ignores_stdin_and_missing_inputs(self, tmp_path: Path):
        f = tmp_path / "x.metacool"
        f.write_bytes(b"x")
        check_output_path(f, force=True, inputs=["-", None, tmp_path / "nope"])
        assert not f.exists()


class TestForceKeepsPreviousOutputUntilSuccess:
    """``create --force`` deleted the previous map before doing any work.

    ``create`` stages its output under ``<name>.partial`` and renames it into
    place only once validation passes, precisely so a failed run cannot destroy
    a finished file — but the CLI's ``--force`` guard unlinked the old output
    (and its ``.ps.tsv``) up front, so a run that then failed left nothing.
    """

    def test_failed_recreate_keeps_old_map_and_ps(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ):
        out = tmp_path / "out.metacool"
        base = [
            "create",
            "pairs",
            "--pairs",
            tiny_pairs,
            "--bins",
            tiny_bins_bed,
            "-o",
            out,
        ]
        assert _run(*base, "--contigs", tiny_contigs_table).exit_code == 0
        ps = tmp_path / "out.ps.tsv"
        before_map, before_ps = out.read_bytes(), ps.read_bytes()

        bad = tmp_path / "bad.tsv"
        bad.write_text("contig_name\tgenome_name\ncontig_A\tMAG_é\n")
        result = _run(*base, "--contigs", bad, "-f")
        assert result.exit_code != 0
        assert out.read_bytes() == before_map
        assert ps.read_bytes() == before_ps
        validate(out)

    def test_successful_recreate_still_overwrites(
        self,
        tmp_path: Path,
        tiny_pairs: Path,
        tiny_bins_bed: Path,
        tiny_contigs_table: Path,
    ):
        out = tmp_path / "out.metacool"
        args = ["create", "pairs", "--pairs", tiny_pairs, "--bins", tiny_bins_bed]
        args += ["--contigs", tiny_contigs_table, "-o", out]
        assert _run(*args).exit_code == 0
        result = _run(*args, "-f", "--assembly-name", "second")
        assert result.exit_code == 0, result.output
        assert MetaCooler(out).info["assembly"] == "second"


def _scale_counts_to_float(path: Path, factor: float) -> None:
    """Rewrite ``pixels/count`` as float64 times *factor* (a rescaled map)."""
    with h5py.File(path, "a") as f:
        counts = f["pixels"]["count"][:].astype(np.float64) * factor
        del f["pixels"]["count"]
        f["pixels"].create_dataset("count", data=counts)
        f.attrs["sum"] = float(counts.sum())


def _pixel_sum(uri: str | Path) -> float:
    return float(cooler.Cooler(str(uri)).pixels()[:]["count"].sum())


class TestMerge:
    def test_float_counts_are_not_truncated(self, tmp_path: Path, tiny_metacooler):
        """``merge`` forced ``count`` to int64, flooring every float count."""
        a, b = tmp_path / "a.metacool", tmp_path / "b.metacool"
        for p in (a, b):
            shutil.copy(tiny_metacooler, p)
            _scale_counts_to_float(p, 0.3)
        out = tmp_path / "m.metacool"
        merge([a, b], out)
        assert _pixel_sum(out) == pytest.approx(_pixel_sum(a) + _pixel_sum(b))

    def test_integer_counts_widen_to_int64(self, tmp_path: Path, tiny_metacooler):
        out = tmp_path / "m.metacool"
        merge([tiny_metacooler, _replica(tiny_metacooler, tmp_path)], out)
        with h5py.File(out, "r") as f:
            assert f["pixels"]["count"].dtype == np.int64

    def test_cli_does_not_restore_the_stale_qc_block(
        self, tmp_path: Path, tiny_metacooler
    ):
        """The CLI re-merged every input's metadata, putting ``qc`` back.

        ``io.merge`` strips the block (it describes one input, not the sum), but
        the CLI then copied every input key into the output: a merged file of
        twice the contacts reported the first input's ``n_contacts``.
        """
        assert "qc" in MetaCooler(tiny_metacooler).metadata
        out = tmp_path / "m.metacool"
        result = _run(
            "merge", tiny_metacooler, _replica(tiny_metacooler, tmp_path), "-o", out
        )
        assert result.exit_code == 0, result.output
        meta = MetaCooler(out).metadata
        assert "qc" not in meta
        assert "merge" in meta["commands"][-1]

    def test_bin_size_survives_single_bin_contigs(
        self, tmp_path: Path, tiny_metacooler
    ):
        """With every contig one bin, cooler inferred ``bin-size=None``.

        The merged file then claimed variable bins, so ``zoomify`` refused it
        and ``binnify`` applied the variable-bin resolution rule.
        """
        coarse = tmp_path / "coarse.metacool"
        binnify(tiny_metacooler, coarse, binsize=2000)
        assert MetaCooler(coarse).info["bin-size"] == 2000
        out = tmp_path / "m.metacool"
        merge([coarse, _replica(coarse, tmp_path)], out)
        info = MetaCooler(out).info
        assert info["bin-size"] == 2000
        assert info["bin-type"] == "fixed"

    def test_entry_uris_are_accepted(self, tmp_path: Path, tiny_mmetacool: Path):
        """The CLI accepts entry URIs, but ``merge`` opened them as plain paths."""
        e1 = f"{tiny_mmetacool}::entries/entry1"
        e2 = f"{tiny_mmetacool}::entries/entry2"
        out = tmp_path / "m.metacool"
        result = _run("merge", e1, e2, "-o", out)
        assert result.exit_code == 0, result.output
        assert _pixel_sum(out) == _pixel_sum(e1) + _pixel_sum(e2)
        validate(out)

    def test_goes_through_the_writer_shims(
        self, tmp_path: Path, tiny_metacooler, monkeypatch
    ):
        """``merge_coolers`` ran ``cooler.create`` without the shims.

        Every other writer suppresses cooler's per-contig enum build and
        ``get_binsize`` loop; merge alone paid both, taking 65 s at 150 k
        contigs where every other command took under 10 s.
        """
        module = cooler_writer_module()
        assert module is not None
        real = module.write_bins
        seen: list[bool] = []

        def spy(*args, chrom_as_enum=True, **kwargs):
            seen.append(chrom_as_enum)
            return real(*args, chrom_as_enum=chrom_as_enum, **kwargs)

        def forbidden(bins):
            raise AssertionError("cooler's per-contig get_binsize was called")

        monkeypatch.setattr(module, "write_bins", spy)
        monkeypatch.setattr(module, "get_binsize", forbidden)
        merge(
            [tiny_metacooler, _replica(tiny_metacooler, tmp_path)],
            tmp_path / "m.metacool",
        )
        assert seen and not any(seen)

    def test_api_refuses_an_input_as_output(self, tmp_path: Path, tiny_metacooler):
        a = tmp_path / "a.metacool"
        shutil.copy(tiny_metacooler, a)
        before = a.read_bytes()
        with pytest.raises(ValueError, match="same file as the input"):
            merge([a, tiny_metacooler], a)
        assert a.read_bytes() == before


def _write_pairs(path: Path, text: str) -> None:
    """Write *text* to *path*, compressed according to its extension."""
    import pyarrow as pa

    codec = {".gz": "gzip", ".bgz": "gzip", ".bz2": "bz2", ".zst": "zstd"}
    comp = codec.get(path.suffix)
    with pa.output_stream(str(path), compression=comp) as out:
        out.write(text.encode())


class TestPairsCodecs:
    """``.zst`` pairs -- advertised in the CLI reference -- crashed on open.

    The body was read through pyarrow, which decodes every codec, but the two
    header scans opened the file with their own narrower openers: ``.zst`` died
    with ``UnicodeDecodeError`` before a pair was read, and ``.bz2`` failed
    whenever ``#chromsize`` lines were needed (an integer ``--bins``).
    """

    @pytest.mark.parametrize("suffix", [".gz", ".bgz", ".bz2", ".zst"])
    def test_every_codec_gives_the_plain_text_map(
        self, tmp_path: Path, tiny_pairs: Path, tiny_contigs_table: Path, suffix
    ):
        header = "".join(
            f"#chromsize: {c} {n}\n"
            for c, n in [
                ("contig_A", 1000),
                ("contig_B", 600),
                ("contig_C", 400),
                ("contig_D", 200),
            ]
        )
        text = header + tiny_pairs.read_text()
        plain, packed = tmp_path / "p.pairs", tmp_path / f"p.pairs{suffix}"
        plain.write_text(text)
        _write_pairs(packed, text)
        maps = []
        for src in (plain, packed):
            out = tmp_path / f"{src.name}.metacool"
            result = _run(
                "create", "pairs", "--pairs", src, "--bins", "100",
                "--contigs", tiny_contigs_table, "-o", out, "--no-ps",
            )  # fmt: skip
            assert result.exit_code == 0, result.output
            maps.append(cooler.Cooler(str(out)).pixels()[:])
        assert maps[0]["count"].sum() > 0
        assert maps[0].equals(maps[1])


class TestPositionsPastContigEnd:
    """A read past its contig's end was folded silently into the last bin.

    ``cooler cload`` refuses such a pair. metacooler binned it into the
    contig's final bin, and QC counted it as long-range intra-contig signal.
    """

    def test_dropped_counted_and_reported(
        self, tmp_path: Path, tiny_contigs_table: Path, caplog
    ):
        pairs = tmp_path / "p.pairs"
        pairs.write_text(
            "#chromsize: contig_A 1000\n#chromsize: contig_B 600\n"
            "#chromsize: contig_C 400\n#chromsize: contig_D 200\n"
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
            "r1\tcontig_A\t10\tcontig_A\t1000\t+\t+\n"  # pos == length: valid
            "r2\tcontig_A\t10\tcontig_A\t500000\t+\t+\n"  # far past the end
            "r3\tcontig_B\t601\tcontig_C\t5\t+\t+\n"  # one past the end
        )
        out = tmp_path / "o.metacool"
        from metacooler.io.create import from_pairs

        with caplog.at_level("WARNING", logger="metacooler"):
            from_pairs(
                uri=out, bins=None, binsize=100, pairs=pairs,
                contig_genome_map=tiny_contigs_table, ps=False,
            )  # fmt: skip
        pix = cooler.Cooler(str(out)).pixels()[:]
        assert int(pix["count"].sum()) == 1
        qc = MetaCooler(out).metadata["qc"]
        assert qc["n_contacts"] == 1
        assert qc["n_pairs_discarded"] == 2
        assert any("2 past the end" in r.getMessage() for r in caplog.records)


def _create(tmp_path: Path, contigs: Path, pairs: Path, bins: Path, *extra: object):
    out = tmp_path / "o.metacool"
    result = _run(
        "create", "pairs", "--pairs", pairs, "--bins", bins,
        "--contigs", contigs, "-o", out, "--no-ps", "-f", *extra,
    )  # fmt: skip
    return result, out


class TestFromCoolStorageMode:
    def test_square_cooler_is_refused(
        self, tmp_path: Path, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """A square cooler had every off-diagonal count doubled, and validated.

        Both triangles were folded into the upper one and summed.
        """
        clr = cooler.Cooler(str(tiny_cool))
        pix = clr.pixels()[:]
        lower = pix[pix.bin1_id != pix.bin2_id].rename(
            columns={"bin1_id": "bin2_id", "bin2_id": "bin1_id"}
        )
        both = pd.DataFrame(pd.concat([pix, lower], ignore_index=True))
        both = both.sort_values(by=["bin1_id", "bin2_id"])
        square = tmp_path / "square.cool"
        cooler.create_cooler(
            str(square), clr.bins()[:], both, symmetric_upper=False, ordered=True
        )
        with pytest.raises(ValueError, match="storage-mode 'square'"):
            from_cool(square, tmp_path / "o.metacool", tiny_contigs_table)


class TestUserTableColumnNames:
    def test_a_genome_id_column_does_not_hijack_the_mapping(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
    ):
        """A user column named ``genome_id`` sent every contig to __unbinned__.

        Renaming ``genome_name`` to ``genome_id`` by *name* gave the frame two
        ``genome_id`` columns; the mapping read the wrong one, and a genome
        literally called "genome_id" appeared, with no error.
        """
        table = tmp_path / "c.tsv"
        table.write_text(
            "contig_name\tgenome_name\tgenome_id\n"
            "contig_A\tMAG_1\tx1\ncontig_B\tMAG_1\tx2\ncontig_C\tMAG_2\tx3\n"
        )
        result, out = _create(tmp_path, table, tiny_pairs, tiny_bins_bed)
        assert result.exit_code == 0, result.output
        mc = MetaCooler(out)
        assert set(mc.genomes["name"]) == {"__unbinned__", "MAG_1", "MAG_2"}
        assert "genome_id_col3" in mc.contigs.columns

    def test_a_genome_id_column_in_the_genomes_table(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path, tiny_contigs_table
    ):
        """The same collision in ``--genomes`` crashed on a duplicate label."""
        genomes = tmp_path / "g.tsv"
        genomes.write_text("genome_name\tgenome_id\nMAG_1\tg1\nMAG_2\tg2\n")
        result, out = _create(
            tmp_path, tiny_contigs_table, tiny_pairs, tiny_bins_bed,
            "--genomes", genomes,
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        g = MetaCooler(out).genomes.set_index("name")
        assert g.loc["MAG_2", "genome_id_col2"] == "g2"

    def test_a_slash_in_a_column_name_keeps_the_file_readable(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
    ):
        """``reads/kb`` became a nested HDF5 group that cooler cannot open.

        ``validate()`` passed it; ``cooler.Cooler(...).chroms()`` and
        ``MetaCooler()`` then raised ``AttributeError``.
        """
        table = tmp_path / "c.tsv"
        table.write_text(
            "contig_name\tgenome_name\treads/kb\n"
            "contig_A\tMAG_1\t1.5\ncontig_B\tMAG_1\t2\ncontig_C\tMAG_2\t3\n"
        )
        result, out = _create(tmp_path, table, tiny_pairs, tiny_bins_bed)
        assert result.exit_code == 0, result.output
        assert "reads_kb" in cooler.Cooler(str(out)).chroms()[:].columns
        assert "reads_kb" in MetaCooler(out).contigs.columns

    def test_noheader_colspec_with_extra_columns(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
    ):
        """``I,J,noheader`` with a third column crashed after the whole pass.

        The extra columns kept pandas' integer labels, which h5py cannot use as
        dataset names. The docs' own ``1,5,noheader`` example triggered it.
        """
        table = tmp_path / "c.tsv"
        table.write_text(
            "contig_A\tx\ty\tz\tMAG_1\ncontig_B\tx\ty\tz\tMAG_1\n"
            "contig_C\tx\ty\tz\tMAG_2\n"
        )
        result, out = _create(
            tmp_path,
            table,
            tiny_pairs,
            tiny_bins_bed,
            "--contigs-format",
            "1,5,noheader",
        )
        assert result.exit_code == 0, result.output
        contigs = MetaCooler(out).contigs
        assert {"col2", "col3", "col4"} <= set(contigs.columns)
        assert set(MetaCooler(out).genomes["name"]) >= {"MAG_1", "MAG_2"}


class TestEmptyPixelStream:
    """An empty pixel stream crashed the external sort ("No objects to concatenate").

    The empty-stream guard only produced a writable frame for the ordered
    (in-memory) path.
    """

    def test_header_only_pairs_with_external_sort(
        self, tmp_path: Path, tiny_contigs_table: Path, tiny_bins_bed: Path
    ):
        pairs = tmp_path / "empty.pairs"
        pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
        result, out = _create(
            tmp_path, tiny_contigs_table, pairs, tiny_bins_bed,
            "--sort-strategy", "external",
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        validate(out)
        assert len(cooler.Cooler(str(out)).pixels()[:]) == 0

    def test_empty_cooler_with_external_sort(
        self, tmp_path: Path, tiny_cool: Path, tiny_contigs_table: Path
    ):
        clr = cooler.Cooler(str(tiny_cool))
        empty = tmp_path / "empty.cool"
        cooler.create_cooler(
            str(empty),
            clr.bins()[:],
            pd.DataFrame({"bin1_id": [], "bin2_id": [], "count": []}, dtype="int64"),
            ordered=True,
        )
        out = tmp_path / "o.metacool"
        from_cool(
            empty, out, contig_genome_map=tiny_contigs_table, sort_in_memory=False
        )
        validate(out)


class TestEarlyRefusals:
    def test_drop_unassigned_that_empties_the_map(
        self, tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path
    ):
        """Failed only after the whole pairs pass, inside cooler.

        The table names a contig of the assembly -- one naming none is refused
        earlier, as a wrong table -- but leaves it unassigned.
        """
        table = tmp_path / "c.tsv"
        table.write_text("contig_name\tgenome_name\ncontig_A\tNA\n")
        result, out = _create(
            tmp_path, table, tiny_pairs, tiny_bins_bed, "--drop-unassigned"
        )
        assert result.exit_code != 0
        assert "removed all 4 contig(s)" in _flat(result.output)
        assert "assigns none of them to a genome" in _flat(result.output)
        assert not out.exists()


class TestFromCoolAssembly:
    def test_source_assembly_is_kept(
        self, tmp_path: Path, tiny_cool: Path, tiny_contigs_table: Path
    ):
        """It was replaced by "imported" whenever no name was passed."""
        src = tmp_path / "named.cool"
        clr = cooler.Cooler(str(tiny_cool))
        cooler.create_cooler(
            str(src), clr.bins()[:], clr.pixels()[:], assembly="hg_meta_v2"
        )
        out = tmp_path / "o.metacool"
        from_cool(src, out, contig_genome_map=tiny_contigs_table)
        assert MetaCooler(out).info["assembly"] == "hg_meta_v2"


class TestValidateStructure:
    def test_nested_group_in_a_table_is_refused(self, tmp_path: Path, tiny_metacooler):
        f = tmp_path / "x.metacool"
        shutil.copy(tiny_metacooler, f)
        with h5py.File(f, "a") as h:
            h["contigs"].create_group("reads").create_dataset("kb", data=[1, 2, 3, 4])
        with pytest.raises(MetaCoolerFormatError, match="is a group, not a column"):
            validate(f)

    @pytest.mark.parametrize("target", ["/nowhere/contigs", "/contigs"])
    def test_chroms_link_must_resolve_to_this_maps_contigs(
        self, tmp_path: Path, tiny_mmetacool: Path, target: str
    ):
        """A dangling (or foreign) /chroms link passed; MetaCooler() then crashed."""
        with h5py.File(tiny_mmetacool, "a") as h:
            e = h["entries/entry1"]
            del e["chroms"]
            e["chroms"] = h5py.SoftLink(target)
        with pytest.raises(MetaCoolerFormatError, match="must resolve to this map"):
            validate(f"{tiny_mmetacool}::entries/entry1")


@pytest.fixture
def empty_metacool(tmp_path: Path, tiny_contigs_table: Path, tiny_bins_bed: Path):
    """A valid metacool with zero pixels (header-only pairs)."""
    pairs = tmp_path / "empty.pairs"
    pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
    out = tmp_path / "empty.metacool"
    from metacooler.io.create import from_pairs

    from_pairs(
        uri=out, bins=tiny_bins_bed, pairs=pairs,
        contig_genome_map=tiny_contigs_table, ps=False,
    )  # fmt: skip
    return out


class TestZeroPixelSources:
    """A valid zero-pixel metacool crashed binnify, zoomify and aggregate.

    Each handed cooler's unordered writer an empty stream (or one of empty
    frames), which fails with "No objects to concatenate" / IndexError.
    """

    @pytest.mark.parametrize("level", ["contig", "genome"])
    def test_binnify(self, tmp_path: Path, empty_metacool: Path, level: Level):
        out = tmp_path / "b.metacool"
        binnify(empty_metacool, out, binsize=200, level=level)
        validate(out)

    def test_zoomify(self, tmp_path: Path, empty_metacool: Path):
        out = tmp_path / "z.mmetacool"
        zoomify(empty_metacool, out, levels=2)
        assert len(list_metacoolers(out)) == 2

    @pytest.mark.parametrize("level", ["contig", "genome"])
    def test_aggregate(self, tmp_path: Path, empty_metacool: Path, level: Level):
        out = tmp_path / "a.cool"
        aggregate(empty_metacool, out, level=level)
        assert len(cooler.Cooler(str(out)).pixels()[:]) == 0


class TestAggregate:
    def test_genome_longer_than_int32_is_refused(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """A >2^31 bp genome got a negative chrom length; fetch() then failed."""
        f = tmp_path / "big.metacool"
        shutil.copy(tiny_metacooler, f)
        with h5py.File(f, "a") as h:
            lengths = h["genomes"]["length"][:]
            lengths[0] = 2_700_000_000
            h["genomes"]["length"][...] = lengths
        mc = MetaCooler(f, validate=False)
        with pytest.raises(ValueError, match="exceed the 2,147,483,647 bp"):
            _write_aggregated_cool(mc, tmp_path / "g.cool", level="genome")

    def test_sparse_aggregate_keeps_float_counts(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """``show``'s aggregate cast float counts to int64, flooring them."""
        f = tmp_path / "f.metacool"
        shutil.copy(tiny_metacooler, f)
        _scale_counts_to_float(f, 0.5)
        mc = MetaCooler(f)
        agg = _aggregate_sparse(mc, "genome")
        upper = float(sparse.triu(agg).sum())
        assert upper == pytest.approx(_pixel_sum(f))


class TestExtractSplitPrefix:
    def test_only_a_trailing_cool_is_stripped(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """``-o sample.rep1`` wrote ``sample_<g>.cool`` (Path.stem), colliding
        with ``-o sample.rep2``."""
        for rep in ("rep1", "rep2"):
            result = _run(
                "extract", tiny_metacooler, "--genomes", "MAG_1", "--split",
                "-o", tmp_path / f"sample.{rep}",
            )  # fmt: skip
            assert result.exit_code == 0, result.output
        assert (tmp_path / "sample.rep1_MAG_1.cool").exists()
        assert (tmp_path / "sample.rep2_MAG_1.cool").exists()
        result = _run(
            "extract", tiny_metacooler, "--genomes", "MAG_1", "--split",
            "-o", tmp_path / "x.cool",
        )  # fmt: skip
        assert (tmp_path / "x_MAG_1.cool").exists()


class TestCollectionInputs:
    """``stack``/``create_mmetacool``/``add_entry`` copied anything in, exit 0."""

    def test_a_plain_cool_is_refused(
        self, tmp_path: Path, tiny_metacooler: Path, tiny_cool: Path
    ):
        out = tmp_path / "s.mmetacool"
        result = _run("stack", tiny_metacooler, tiny_cool, "--names", "a,b", "-o", out)
        assert result.exit_code != 0
        assert not out.exists()

    def test_a_slash_in_an_entry_name_is_refused(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        out = tmp_path / "s.mmetacool"
        result = _run(
            "stack", tiny_metacooler, tiny_metacooler, "--names", "a,x/y", "-o", out
        )
        assert result.exit_code != 0
        assert "Invalid entry name 'x/y'" in _flat(result.output)

    def test_an_absolute_chroms_link_is_made_relative(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """An absolute ``/contigs`` link is legal standalone but dangles in an entry."""
        src = tmp_path / "abs.metacool"
        shutil.copy(tiny_metacooler, src)
        with h5py.File(src, "a") as h:
            del h["chroms"]
            h["chroms"] = h5py.SoftLink("/contigs")
        validate(src)
        out = tmp_path / "s.mmetacool"
        create_mmetacool(out, {"a": src, "b": tiny_metacooler})
        validate(f"{out}::entries/a")
        assert MetaCooler(f"{out}::entries/a").chroms()["name"].tolist()


class TestMetaCoolerApi:
    def test_root_and_repr_of_an_entry(self, tiny_mmetacool: Path):
        mc = MetaCooler(f"{tiny_mmetacool}::entries/entry1")
        assert mc.root == "/entries/entry1"
        assert "entries/entry1" in repr(mc)

    def test_info_keeps_a_fractional_sum(self, tmp_path: Path, tiny_metacooler):
        f = tmp_path / "f.metacool"
        shutil.copy(tiny_metacooler, f)
        _scale_counts_to_float(f, 0.3)
        assert MetaCooler(f).info["sum"] == pytest.approx(_pixel_sum(f))

    def test_contig_range_on_a_genome_level_file(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """Every contig is also a genome there, and ``G:1-200`` was refused."""
        g = tmp_path / "g.metacool"
        binnify(tiny_metacooler, g, binsize=200, level="genome")
        mc = MetaCooler(g)
        assert parse_metacooler_range(mc, "MAG_1:1-400") == (
            mc.contig_extent("MAG_1")[0],
            mc.contig_extent("MAG_1")[0] + 2,
        )


class TestPsBinnedCurve:
    """The binned-map P(s) used bin pairs as opportunities and full log widths.

    The diagonal bin pair counted as a full opportunity though it holds half
    the position pairs (halving the first point), and ``prob`` was normalised
    over ``hi - lo`` of every log bin although a binned map only realises
    separations at bin-distance steps -- so it grew with --bins-per-decade
    (pixel/pairs ratio 1.07 -> 1.28 -> 1.52 for bpd 5 -> 10 -> 20).
    """

    def test_weighted_denominators_match_brute_force(
        self, tiny_variable_bin_metacooler
    ):
        mc = MetaCooler(tiny_variable_bin_metacooler)
        bins = mc.bins()
        contig = chrom_codes(bins).astype(np.int64)
        start = bins["start"].to_numpy(dtype=np.int64)
        end = bins["end"].to_numpy(dtype=np.int64)
        mids, w = (start + end) // 2, end - start
        edges = distance_edges(int(mids.max()), 50, 10)
        key = contig * (int(mids.max()) + int(edges[-1]) + 1) + mids
        got = bin_pair_denominators(key, np.zeros(len(key), np.int64), 1, edges, w)[0]

        expected = np.zeros(len(edges) - 1)
        for i in range(len(bins)):
            for j in range(i, len(bins)):
                if contig[i] != contig[j]:
                    continue
                weight = w[i] * (w[i] + 1) // 2 if i == j else w[i] * w[j]
                d = abs(mids[j] - mids[i])
                k = np.searchsorted(edges, d, side="right") - 1
                if k < len(expected):
                    expected[k] += weight
        np.testing.assert_allclose(got, expected)
        # Every unordered position pair of every contig, exactly once.
        lengths = mc.contigs["length"].to_numpy(dtype=np.int64)
        assert got.sum() == (lengths * (lengths + 1) // 2).sum()

    @pytest.fixture
    def uniform(self, tmp_path: Path):
        """One contig with a constant contact rate per position pair."""
        rng = np.random.default_rng(0)
        length = 20_000
        pos = rng.integers(1, length + 1, size=(200_000, 2))
        pairs = tmp_path / "u.pairs"
        pairs.write_text(
            f"#chromsize: c1 {length}\n"
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
            + "".join(f"r\tc1\t{a}\tc1\t{b}\t+\t+\n" for a, b in pos)
        )
        table = tmp_path / "c.tsv"
        table.write_text("contig_name\tgenome_name\nc1\tM\n")
        out = tmp_path / "u.metacool"
        from_pairs(
            uri=out, bins=None, binsize=1000, pairs=pairs,
            contig_genome_map=table, ps=False, qc=False,
        )  # fmt: skip
        return out, pairs

    def test_a_uniform_map_gives_a_flat_curve(self, uniform):
        out, pairs = uniform
        pixel = contact_probability(out)["avg_contacts"].to_numpy()
        assert pixel[0] == pytest.approx(np.median(pixel[1:]), rel=0.05)
        # Same unit on both paths: contacts per position pair.
        reads = contact_probability(out, pairs=pairs)["avg_contacts"].to_numpy()
        assert np.median(pixel) == pytest.approx(np.median(reads), rel=0.02)

    def test_dist_sits_on_the_separations_a_bin_can_hold(self, uniform):
        out, pairs = uniform
        pixel = contact_probability(out)
        assert pixel["dist"].iloc[1:4].tolist() == [1000.0, 2000.0, 3000.0]
        reads = contact_probability(out, pairs=pairs)
        assert reads["dist"].iloc[1:3].tolist() == [1.0, 2.0]

    def test_prob_does_not_depend_on_bins_per_decade(self):
        data = Path(__file__).parent / "data"
        m = data / "assembly_1_rep1_1kb.metacool"
        pairs = data / "contacts_1_rep1.pairs.gz"

        def ratio(bpd: int) -> float:
            def at(t, s):
                t = t[t["dist"] > 0]
                x, y = np.log10(t["dist"]), t["prob"]
                return float(np.interp(np.log10(s), x, y))

            px = contact_probability(m, bins_per_decade=bpd)
            pr = contact_probability(m, pairs=pairs, bins_per_decade=bpd)
            return at(px, 10_000) / at(pr, 10_000)

        for bpd in (5, 10, 20):
            assert ratio(bpd) == pytest.approx(1.0, abs=0.08)


class TestQcUnbinned:
    """``__unbinned__`` was classified as one genome.

    A contact between two different unassigned contigs counted as
    inter_contig_intra_genome, inflating intra_genome_fraction by every contact
    among the unbinned contigs -- often most of a real assembly.
    """

    def test_two_unbinned_contigs_are_inter_genome(self):
        from metacooler.qc import ContactQC, unbinned_genome_id

        names = ["__unbinned__", "MAG_1"]
        tally = ContactQC(2, 1000, unbinned_genome=unbinned_genome_id(names))
        # contigs 0,1 unbinned; 2,3 in MAG_1
        c1 = np.array([0, 0, 2, 0])
        c2 = np.array([1, 0, 3, 2])
        g = np.array([0, 0, 1, 1])
        tally.tally(c1, c2, g[c1], g[c2], np.array([0, 50, 0, 0]))
        s = tally.summary()
        assert s["version"] == 2
        assert s["inter_genome"] == 2  # unbinned-unbinned, unbinned-MAG
        assert s["inter_contig_intra_genome"] == 1  # the two MAG_1 contigs
        assert s["intra_contig_short"] == 1
        assert s["intra_genome_fraction"] == pytest.approx(0.5)
        assert s["trans_contig_ratio"] == pytest.approx(0.75)  # unchanged rule

    def test_create_and_recompute_agree(self, tmp_path: Path, tiny_metacooler):
        """The pixel recompute applies the same rule as creation."""
        from metacooler.ops.qc import qc as run_qc

        created = MetaCooler(tiny_metacooler).metadata["qc"]
        again = run_qc(tiny_metacooler).summary
        for key in ("inter_genome", "inter_contig_intra_genome", "n_contacts"):
            assert created[key] == again[key]


class TestShowLog10:
    def test_log10_is_log10_of_one_plus_count(self):
        """log10(count) mapped a count of 1 to 0, the colour of an empty pixel."""
        from metacooler.ops.show import _transform_matrix

        out = _transform_matrix(np.array([[0.0, 1.0], [9.0, 99.0]]), "log10")
        np.testing.assert_allclose(out, [[0.0, np.log10(2)], [1.0, 2.0]])


class TestFilterMissingText:
    def test_negated_operators_skip_missing_text(self):
        """Missing text reads back as "" and ``!=``/``!~`` used to match it."""
        from metacooler.filters import filter_mask, parse_filter

        df = pd.DataFrame({"tax": ["Bacillota", "", "Pseudomonadota", "  "]})
        for expr, expected in [
            ("contigs.tax != Bacillota", [False, False, True, False]),
            ("contigs.tax !~ ^Bac", [False, False, True, False]),
            ("contigs.tax =~ .*", [True, False, True, False]),
        ]:
            mask = filter_mask(df, [parse_filter(expr)], table="contigs")
            assert mask.tolist() == expected, expr


class TestValidateDuplicatePixels:
    def test_a_repeated_pixel_is_refused(self, tmp_path: Path, tiny_metacooler):
        """``<=`` on bin2 let a duplicated (bin1, bin2) row pass validate and check."""
        f = tmp_path / "dup.metacool"
        shutil.copy(tiny_metacooler, f)
        with h5py.File(f, "a") as h:
            px = h["pixels"]
            cols = {k: px[k][:] for k in px}
            for k, v in cols.items():
                del px[k]
                px.create_dataset(k, data=np.insert(v, 5, v[5]))
            idx = h["indexes"]["bin1_offset"][:]
            idx[cols["bin1_id"][5] + 1 :] += 1
            h["indexes"]["bin1_offset"][...] = idx
            h.attrs["nnz"] = int(h.attrs["nnz"]) + 1
            h.attrs["sum"] = int(h.attrs["sum"]) + int(cols["count"][5])
        with pytest.raises(MetaCoolerFormatError, match="is repeated"):
            validate(f)
        result = _run("check", "--deep", f)
        assert result.exit_code != 0
        assert "duplicate pixel" in _flat(result.output)


class TestDistanceOptions:
    @pytest.mark.parametrize(
        "args",
        [
            ["qc", "{m}", "--threshold", "inf", "--table"],
            ["qc", "{m}", "--threshold", "nan", "--table"],
            ["ps", "{m}", "--max-distance", "1e30"],
            ["ps", "{m}", "--min-distance", "-5"],
        ],
        ids=lambda a: " ".join(a[2:4]),
    )
    def test_unusable_distances_are_a_usage_error(self, tiny_metacooler, args):
        """``inf``/``1e30`` escaped as OverflowError/TypeError tracebacks."""
        result = _run(*[a.format(m=tiny_metacooler) for a in args])
        assert result.exit_code == 2, result.output
        assert "not a usable distance" in _flat(result.output)


class TestReservedQcKey:
    def test_metadata_file_may_not_carry_qc(
        self, tmp_path: Path, tiny_pairs, tiny_bins_bed, tiny_contigs_table
    ):
        """create overwrote a user's ``qc`` key with the QC block, silently."""
        meta = tmp_path / "m.yaml"
        meta.write_text("qc: passed\n")
        result, _ = _create(
            tmp_path, tiny_contigs_table, tiny_pairs, tiny_bins_bed,
            "--metadata", meta,
        )  # fmt: skip
        assert result.exit_code == 2
        assert "'qc' key" in _flat(result.output)

    def test_store_tolerates_a_non_object_qc_value(
        self, tmp_path: Path, tiny_metacooler
    ):
        """``qc --store`` raised AttributeError on a string ``qc`` value."""
        f = tmp_path / "f.metacool"
        shutil.copy(tiny_metacooler, f)
        with h5py.File(f, "a") as h:
            h.attrs["metadata"] = '{"qc": "passed"}'
        result = _run("qc", f, "--store")
        assert result.exit_code == 0, result.output
        assert isinstance(MetaCooler(f).metadata["qc"], dict)


class TestUncoveredEdges:
    """Edge branches of this pass's fixes, so patch coverage stays complete."""

    def test_writer_module_is_none_for_an_unknown_layout(self, monkeypatch):
        import metacooler.io._utils as utils

        monkeypatch.setattr(
            utils, "_COOLER_WRITER_MODULES", ("cooler.no_such_module", "cooler.util")
        )
        assert utils.cooler_writer_module() is None

    def test_last_pairs_line_without_a_newline_is_read(self, tmp_path: Path):
        from metacooler.io.create import _iter_pairs_lines

        f = tmp_path / "p.pairs"
        f.write_bytes(b"#a\n#b\nr\tc\t1\tc\t2")
        assert list(_iter_pairs_lines(f, block_size=3)) == [
            "#a\n",
            "#b\n",
            "r\tc\t1\tc\t2",
        ]

    def test_blank_lines_between_header_and_body(self, tmp_path: Path):
        from metacooler.io.create import _scan_pairs_header

        with_body = tmp_path / "a.pairs"
        with_body.write_bytes(b"#h1\n#h2\n\n  \nr\tc\t1\tc\t2\n")
        assert _scan_pairs_header(with_body) == (2, True)
        header_only = tmp_path / "b.pairs"
        header_only.write_bytes(b"#h1\n\n\n")
        assert _scan_pairs_header(header_only) == (1, False)

    def test_a_chunk_entirely_past_the_contig_ends(self, tmp_path: Path):
        from metacooler.io.create import _decode_pairs_chunk, _PairsDropStats

        chunk = pd.DataFrame(
            {
                "chrom1": ["a", "a"],
                "pos1": [5, 1],
                "chrom2": ["a", "a"],
                "pos2": [999, 50],
            }
        )
        stats = _PairsDropStats()
        decoded = _decode_pairs_chunk(
            chunk, {"a": 0}, stats, contig_length=np.array([10])
        )
        assert decoded is None
        assert stats.out_of_bounds == 2
        assert stats.dropped == 2

    def test_level_chroms_rejects_an_unknown_level(self, tiny_metacooler: Path):
        from metacooler.ops.aggregate import _level_chroms

        with pytest.raises(ValueError, match="level must be 'contig' or 'genome'"):
            _level_chroms(MetaCooler(tiny_metacooler), "bin")  # ty:ignore[invalid-argument-type]
