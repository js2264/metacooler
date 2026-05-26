"""Tests for network edge list export."""

from __future__ import annotations

from pathlib import Path

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

    def test_network_stdout(self, tmp_path: Path, tiny_metacooler: Path, capsys):
        network(tiny_metacooler, level="contig")
        captured = capsys.readouterr()
        assert len(captured.out.strip()) > 0

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
