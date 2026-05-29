"""Contig table parsers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

import pandas as pd


def parse_contigs_table(
    path: str | Path,
    fmt: Literal["contig_genome", "metator", "contigs"] | str = "contigs",
) -> tuple[dict[str, str], pd.DataFrame]:
    """Parse a contig-to-genome assignment table.

    Args:
        path: Path to the contig table file.
        fmt: Format of the table.
            - ``contigs``: TSV with at least ``contig_name`` and ``genome_name``
                            columns. Additional columns are preserved and returned in the
                            DataFrame. Expects a header.
                        - ``contig_genome``: Two-column TSV ``contig_name<TAB>genome_name``.
                            Expects no header.
            - ``metator``: MetaTOR ``contig_data_final.txt`` with at least
                            ``Name`` and ``Final_bin`` columns. Expects a header.
                        - ``<i>,<j>``: Use column ``i`` for contig names and column ``j``
                            for genome names (1-based indexing). Expects a header.
                        - ``<i>,<j>,noheader``: Same as above but parses the file as having
                            no header.

    Returns:
        Tuple of (mapping ``contig_name → genome_name``, full DataFrame with
        all columns renamed to ``contig_name`` / ``genome_id`` internally).
        Contigs not present in the table are omitted from the
        mapping (the caller decides whether to assign them to ``__unbinned__``
        or drop them).

    Note:
        The ``contigs`` format expects input columns named ``contig_name`` and
        ``genome_name`` (as validated on read); ``genome_name`` is renamed
        internally to ``genome_id`` in the returned DataFrame.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Contig table not found: {path}")

    if fmt == "contigs":
        df = pd.read_csv(
            path,
            sep="\t",
            header=0,
            dtype=str,
            keep_default_na=False,
        )
        if df.shape[1] < 2:
            raise ValueError(
                "contigs format requires at least 2 columns: "
                "column 1 (contig_name) and column 2 (genome_name)."
            )
        if df.columns[0] != "contig_name":
            raise ValueError(
                f"contigs format expects first column to be 'contig_name', "
                f"got '{df.columns[0]}'."
            )
        if df.columns[1] != "genome_name":
            raise ValueError(
                f"contigs format expects second column to be 'genome_name', "
                f"got '{df.columns[1]}'."
            )
        df = df.rename(columns={"genome_name": "genome_id"})
        mapping = dict(zip(df["contig_name"], df["genome_id"], strict=False))
        return mapping, df

    if fmt == "contig_genome":
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["contig_name", "genome_id"],
            dtype=str,
            keep_default_na=False,
        )
        mapping = dict(zip(df["contig_name"], df["genome_id"], strict=False))
        return mapping, df

    if fmt == "metator":
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        if "Name" not in df.columns or "Final_bin" not in df.columns:
            raise ValueError("metator format requires 'Name' and 'Final_bin' columns")
        df = df.rename(columns={"Name": "contig_name", "Final_bin": "genome_id"})
        mapping = dict(zip(df["contig_name"], df["genome_id"], strict=False))
        return mapping, df

    colspec_match = re.fullmatch(r"\s*(\d+)\s*,\s*(\d+)(?:\s*,\s*(noheader))?\s*", fmt)
    if colspec_match:
        contig_col = int(colspec_match.group(1))
        genome_col = int(colspec_match.group(2))
        noheader = colspec_match.group(3) is not None

        if contig_col < 1 or genome_col < 1:
            raise ValueError("Column indices in contigs format must be >= 1")

        header: int | None = None if noheader else 0
        df = pd.read_csv(
            path,
            sep="\t",
            header=header,
            dtype=str,
            keep_default_na=False,
        )

        required_cols = max(contig_col, genome_col)
        if df.shape[1] < required_cols:
            raise ValueError(
                "contigs format column spec requires at least "
                f"{required_cols} columns in the table, got {df.shape[1]}."
            )

        contig_key = df.columns[contig_col - 1]
        genome_key = df.columns[genome_col - 1]
        df = df.rename(columns={contig_key: "contig_name", genome_key: "genome_id"})
        mapping = dict(zip(df["contig_name"], df["genome_id"], strict=False))
        return mapping, df

    raise ValueError(f"Unknown contig table format: {fmt!r}")
