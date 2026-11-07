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

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from _cli_output import flat as _flat
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io._utils import check_output_path, cooler_writer_module
from metacooler.io.create import from_cool
from metacooler.io.merge import merge
from metacooler.ops.binnify import binnify
from metacooler.schema import MetaCoolerFormatError, validate


def _run(*args: str):
    return CliRunner().invoke(cli, [str(a) for a in args])


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
        merge([tiny_metacooler, tiny_metacooler], out)
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
        result = _run("merge", tiny_metacooler, tiny_metacooler, "-o", out)
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
        merge([coarse, coarse], out)
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
        real = module.write_bins
        seen: list[bool] = []

        def spy(*args, chrom_as_enum=True, **kwargs):
            seen.append(chrom_as_enum)
            return real(*args, chrom_as_enum=chrom_as_enum, **kwargs)

        def forbidden(bins):
            raise AssertionError("cooler's per-contig get_binsize was called")

        monkeypatch.setattr(module, "write_bins", spy)
        monkeypatch.setattr(module, "get_binsize", forbidden)
        merge([tiny_metacooler, tiny_metacooler], tmp_path / "m.metacool")
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


def _create(tmp_path: Path, contigs: Path, pairs: Path, bins: Path, *extra: str):
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
        both = pd.concat([pix, lower]).sort_values(["bin1_id", "bin2_id"])
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
        """Failed only after the whole pairs pass, inside cooler."""
        table = tmp_path / "c.tsv"
        table.write_text("contig_name\tgenome_name\nnot_a_contig\tMAG_1\n")
        result, out = _create(
            tmp_path, table, tiny_pairs, tiny_bins_bed, "--drop-unassigned"
        )
        assert result.exit_code != 0
        assert "removed all 4 contig(s)" in _flat(result.output)


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
