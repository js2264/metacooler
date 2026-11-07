"""I/O utilities for reading genome-metadata tables.

This module provides helpers for parsing the TSV-formatted genome/MAG
metadata file that describes which genomes are present in a metagenomic
Hi-C dataset.
"""

import logging
from pathlib import Path

import pandas as pd

from metacooler.io._utils import read_user_table, rename_key_columns

logger = logging.getLogger("metacooler")


def parse_genomes_table(path: str | Path) -> pd.DataFrame:
    """Parse a genomes metadata table (with header).

    The file must be a tab-separated file **with a header row**.  The first
    column must be named ``genome_name`` (required); any additional columns
    are stored verbatim in the ``genomes`` table.

    Args:
        path: Path to the genomes TSV file.

    Returns:
        DataFrame with at least a ``genome_id`` column (renamed from
        ``genome_name``) and any extra columns.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If the file is empty, its header and data rows disagree on
            the number of fields, or its first column is not named
            ``genome_name``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Genomes table not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Genomes table is empty: {path}")
    df = read_user_table(path, kind="Genomes table")
    if df.columns[0] != "genome_name":
        raise ValueError(
            f"genomes table expects first column to be 'genome_name', "
            f"got '{df.columns[0]}'."
        )
    df = rename_key_columns(df, {0: "genome_id"}, kind="Genomes table")
    dup = df["genome_id"].duplicated()
    if dup.any():
        names = df.loc[dup, "genome_id"].unique().tolist()
        logger.warning(
            "Genomes table contains %d duplicate genome name(s): %s (last row wins).",
            len(names),
            names[:10],
        )
    return df
