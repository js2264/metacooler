"""Tests for network edge list export."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from _edges import assert_ids_name

from metacooler.api import MetaCooler
from metacooler.io.network import network

DATA = Path(__file__).parent / "data"


class TestNetwork:
    def test_network_bins_level(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "edges_bins.tsv"
        network(tiny_metacooler, level="bin", out=out)
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) > 0
        # Each line should have 3 tab-separated columns
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 3

    def test_network_contigs_level(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "edges_contigs.tsv"
        network(tiny_metacooler, level="contig", out=out)
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) > 0
        # Node names should be contig names
        first_parts = lines[0].split("\t")
        assert not first_parts[0].isdigit() or first_parts[0] == "0"

    def test_network_genomes_level(self, tmp_path: Path, tiny_metacooler: Path):
        out = tmp_path / "edges_genomes.tsv"
        network(tiny_metacooler, level="genome", out=out)
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) > 0

    def test_network_invalid_weight(self, tmp_path: Path, tiny_metacooler: Path):
        with pytest.raises(ValueError, match="not found"):
            network(tiny_metacooler, level="bin", weight="nonexistent")

    @pytest.mark.parametrize("level", ["bin", "contig", "genome"])
    def test_network_returns_the_edges_without_out(
        self, tmp_path: Path, tiny_metacooler: Path, capsys, level
    ):
        """Without *out*, nothing is printed: the edges are returned, the
        rows the file gets, in its order."""
        edges = network(tiny_metacooler, level=level)
        assert capsys.readouterr().out == ""
        assert isinstance(edges, pd.DataFrame)
        assert list(edges.columns) == ["node1", "node2", "value"]
        out = tmp_path / "edges.tsv"
        assert network(tiny_metacooler, level=level, out=out) is None
        written = pd.read_csv(out, sep="\t", header=None, names=list(edges.columns))
        if level != "bin":
            written = written.astype({"node1": str, "node2": str})
        pd.testing.assert_frame_equal(
            edges.reset_index(drop=True), written, check_dtype=False
        )

    def test_network_returned_contig_edges_sum_the_pixels(self, tiny_metacooler: Path):
        """The contig table holds one row per pair of contigs, a contig's
        contacts with itself included, and sums every pixel."""
        from metacooler.api import MetaCooler

        edges = network(tiny_metacooler, level="contig")
        assert edges is not None
        mc = MetaCooler(tiny_metacooler)
        assert int(edges["value"].sum()) == int(mc.pixels()[:]["count"].sum())
        assert not edges.duplicated(["node1", "node2"]).any()
        assert set(edges["node1"]) | set(edges["node2"]) <= set(mc.contigs["name"])

    def test_network_returns_no_edge_of_an_empty_map(self, tmp_path: Path):
        from metacooler.io.create import from_pairs

        bins = tmp_path / "bins.bed"
        bins.write_text("c1\t0\t1000\n")
        pairs = tmp_path / "empty.pairs"
        pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2\n")
        out = tmp_path / "empty.metacool"
        from_pairs(uri=out, bins=bins, pairs=pairs)
        for level in ("bin", "contig", "genome"):
            edges = network(out, level=level)  # type: ignore[arg-type]
            assert edges is not None and edges.empty
            assert list(edges.columns) == ["node1", "node2", "value"]

    @pytest.mark.parametrize("level", ["bin", "contig"])
    def test_network_writes_to_an_open_stream(
        self, tmp_path: Path, tiny_metacooler: Path, level
    ):
        """A stream gets what a file gets, header included, and stays open."""
        stream = io.StringIO()
        assert network(tiny_metacooler, level=level, out=stream, header=True) is None
        assert not stream.closed
        out = tmp_path / "edges.tsv"
        network(tiny_metacooler, level=level, out=out, header=True)
        assert stream.getvalue() == out.read_text()
        assert stream.getvalue().startswith("node1\tnode2\tvalue\n")

    def test_network_invalid_level(self, tmp_path: Path, tiny_metacooler: Path):
        with pytest.raises(ValueError, match="level must be"):
            network(uri=tiny_metacooler, level="invalid")  # type: ignore

    def test_network_contigs_node_names_match_contigs_table(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        contig_names = set(mc.contigs["name"].values)
        out = tmp_path / "edges.tsv"
        network(tiny_metacooler, level="contig", out=out)
        for line in out.read_text().strip().split("\n"):
            parts = line.split("\t")
            assert parts[0] in contig_names, f"Unknown node: {parts[0]}"
            assert parts[1] in contig_names, f"Unknown node: {parts[1]}"

    def test_network_genomes_node_names_match_genomes_table(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        genome_names = set(mc.genomes["name"].values)
        out = tmp_path / "edges.tsv"
        network(tiny_metacooler, level="genome", out=out)
        for line in out.read_text().strip().split("\n"):
            parts = line.split("\t")
            assert parts[0] in genome_names, f"Unknown node: {parts[0]}"
            assert parts[1] in genome_names, f"Unknown node: {parts[1]}"

    def test_network_upper_triangular_at_contig_level(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """At contig level, node1 id-order <= node2 id-order (no duplicated pairs)."""
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        contig_names = list(mc.contigs["name"].values)
        name_to_idx = {n: i for i, n in enumerate(contig_names)}

        out = tmp_path / "edges.tsv"
        network(tiny_metacooler, level="contig", out=out)
        seen = set()
        for line in out.read_text().strip().split("\n"):
            parts = line.split("\t")
            n1, n2 = parts[0], parts[1]
            i1 = name_to_idx[n1]
            i2 = name_to_idx[n2]
            assert i1 <= i2, f"Not upper-triangular: ({n1}, {n2})"
            pair = (n1, n2)
            assert pair not in seen, f"Duplicate edge: {pair}"
            seen.add(pair)

    def test_network_upper_triangular_at_genome_level(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """At genome level, node1 id-order <= node2 id-order (no duplicated pairs)."""
        from metacooler.api import MetaCooler

        mc = MetaCooler(tiny_metacooler)
        genome_names = list(mc.genomes["name"].values)
        name_to_idx = {n: i for i, n in enumerate(genome_names)}

        out = tmp_path / "edges.tsv"
        network(tiny_metacooler, level="genome", out=out)
        seen = set()
        for line in out.read_text().strip().split("\n"):
            parts = line.split("\t")
            n1, n2 = parts[0], parts[1]
            i1 = name_to_idx[n1]
            i2 = name_to_idx[n2]
            assert i1 <= i2, f"Not upper-triangular: ({n1}, {n2})"
            pair = (n1, n2)
            assert pair not in seen, f"Duplicate edge: {pair}"
            seen.add(pair)

    def test_network_total_weight_conserved_across_levels(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """Sum of weights should be equal at bin, contig, and genome levels."""
        bin_out = tmp_path / "bin.tsv"
        contig_out = tmp_path / "contig.tsv"
        genome_out = tmp_path / "genome.tsv"
        network(tiny_metacooler, level="bin", out=bin_out)
        network(tiny_metacooler, level="contig", out=contig_out)
        network(tiny_metacooler, level="genome", out=genome_out)

        def total_weight(path: Path) -> int:
            return sum(
                int(ln.split("\t")[2])
                for ln in path.read_text().strip().split("\n")
                if ln
            )

        assert total_weight(bin_out) == total_weight(contig_out)
        assert total_weight(contig_out) == total_weight(genome_out)

    def test_network_bin_nodes_are_integers(
        self, tmp_path: Path, tiny_metacooler: Path
    ):
        """At bin level, node1 and node2 are integer bin IDs."""
        out = tmp_path / "bin.tsv"
        network(tiny_metacooler, level="bin", out=out)
        for line in out.read_text().strip().split("\n"):
            parts = line.split("\t")
            assert parts[0].isdigit(), f"Expected integer bin ID, got: {parts[0]}"
            assert parts[1].isdigit(), f"Expected integer bin ID, got: {parts[1]}"


class TestNetworkIds:
    """``ids=True``: each contig or genome as its row in ``/contigs`` or
    ``/genomes``, the rows of the names in their order, summed without a
    name."""

    @pytest.mark.parametrize("level", ["contig", "genome"])
    @pytest.mark.parametrize("chunksize", [100, 50_000_000])
    def test_ids_name_back_to_the_names(self, level, chunksize):
        """Named through ``/contigs`` or ``/genomes``, the ids are the names:
        same rows, same order, same values, on a map whose contigs and
        genomes are not in name order (folds included at 100 pixels a chunk)."""
        path = DATA / "assembly_1_rep1_1kb.metacool"
        mc = MetaCooler(path)
        table = mc.contigs if level == "contig" else mc.genomes
        assert table["name"].tolist() != sorted(table["name"])
        by_name = network(path, level=level, chunksize=chunksize)
        by_id = network(path, level=level, chunksize=chunksize, ids=True)
        assert by_name is not None and by_id is not None
        assert len(by_id) > len(table)  # pairs, not only self-contacts
        assert_ids_name(by_id, by_name, table["name"])

    @pytest.mark.parametrize("level", ["contig", "genome"])
    def test_ids_on_the_tiny_map(self, tiny_metacooler: Path, level):
        mc = MetaCooler(tiny_metacooler)
        table = mc.contigs if level == "contig" else mc.genomes
        by_id = network(tiny_metacooler, level=level, ids=True)
        by_name = network(tiny_metacooler, level=level)
        assert by_id is not None and by_name is not None
        assert_ids_name(by_id, by_name, table["name"])
        assert (by_id["node1"] <= by_id["node2"]).all()
        assert int(by_id["value"].sum()) == int(mc.pixels()["count"].sum())

    def test_ids_change_nothing_at_bin_level(self, tiny_metacooler: Path, tmp_path):
        """Bin-level nodes are bin ids already: the same frame, the same file."""
        by_id = network(tiny_metacooler, level="bin", ids=True)
        by_name = network(tiny_metacooler, level="bin")
        assert by_id is not None and by_name is not None
        pd.testing.assert_frame_equal(by_id, by_name)
        ids_out, names_out = tmp_path / "ids.tsv", tmp_path / "names.tsv"
        network(tiny_metacooler, level="bin", ids=True, out=ids_out, header=True)
        network(tiny_metacooler, level="bin", out=names_out, header=True)
        assert ids_out.read_bytes() == names_out.read_bytes()

    @pytest.mark.parametrize("level", ["contig", "genome"])
    def test_ids_written_to_a_path_or_a_stream(self, tmp_path: Path, level):
        """``out=`` gets the ids the call returns, header included, the
        stream left open."""
        path = DATA / "assembly_1_rep1_1kb.metacool"
        by_id = network(path, level=level, ids=True)
        assert by_id is not None
        out = tmp_path / "edges.tsv"
        assert network(path, level=level, ids=True, out=out, header=True) is None
        written = pd.read_csv(out, sep="\t")
        pd.testing.assert_frame_equal(written, by_id)
        stream = io.StringIO()
        assert network(path, level=level, ids=True, out=stream, header=True) is None
        assert not stream.closed
        assert stream.getvalue() == out.read_text()
        bare = tmp_path / "bare.tsv"
        network(path, level=level, ids=True, out=bare)
        assert bare.read_text() == out.read_text().split("\n", 1)[1]

    def test_ids_refuse_what_names_refuse(self, tiny_metacooler: Path):
        with pytest.raises(ValueError, match="level must be"):
            network(tiny_metacooler, level="node", ids=True)  # type: ignore
        with pytest.raises(ValueError, match="not found"):
            network(tiny_metacooler, weight="nonexistent", ids=True)

    def test_ids_read_no_name(self, tiny_metacooler: Path, monkeypatch):
        """Without names to give, network reads neither ``/contigs`` nor
        ``/genomes`` as a table, nor the bin table: a bin's contig comes from
        ``chrom_offset``, a contig's genome from ``/contigs/genome_id``."""
        levels = ("contig", "genome")
        expected = [network(tiny_metacooler, level=level, ids=True) for level in levels]

        def refuse(self, *args, **kwargs):
            raise AssertionError("network(ids=True) read a table it does not need")

        monkeypatch.setattr(MetaCooler, "contigs", property(refuse))
        monkeypatch.setattr(MetaCooler, "genomes", property(refuse))
        monkeypatch.setattr(MetaCooler, "bins", refuse)
        for level, edges in zip(levels, expected, strict=True):
            got = network(tiny_metacooler, level=level, ids=True)
            assert got is not None and edges is not None
            pd.testing.assert_frame_equal(got, edges)
        with pytest.raises(AssertionError, match="does not need"):
            network(tiny_metacooler, level="contig")

    def test_ids_of_an_empty_map(self, tmp_path: Path):
        """A map without a pixel has no edge, its node columns ``int64``."""
        from metacooler.io.create import from_pairs

        bins = tmp_path / "bins.bed"
        bins.write_text("c1\t0\t1000\n")
        pairs = tmp_path / "empty.pairs"
        pairs.write_text("#columns: readID chrom1 pos1 chrom2 pos2\n")
        empty = tmp_path / "empty.metacool"
        from_pairs(uri=empty, bins=bins, pairs=pairs)
        for level in ("contig", "genome"):
            edges = network(empty, level=level, ids=True)
            assert edges is not None and edges.empty
            assert edges.dtypes.to_dict() == {
                "node1": np.int64,
                "node2": np.int64,
                "value": np.int64,
            }
            out = tmp_path / f"{level}.tsv"
            network(empty, level=level, ids=True, out=out, header=True)
            assert out.read_text() == "node1\tnode2\tvalue\n"
        by_id = network(empty, level="bin", ids=True)
        by_name = network(empty, level="bin")
        assert by_id is not None and by_name is not None
        pd.testing.assert_frame_equal(by_id, by_name)
