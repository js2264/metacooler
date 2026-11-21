"""Tests for network edge list export."""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import pytest

from metacooler.io.network import network


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
