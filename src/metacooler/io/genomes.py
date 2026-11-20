"""I/O utilities for reading genome-metadata tables.

This module provides helpers for parsing the TSV-formatted genome/MAG
metadata file that describes which genomes are present in a metagenomic
Hi-C dataset.
"""

import logging
from pathlib import Path

import pandas as pd

from metacooler.io._utils import drop_dump_columns, read_user_table, rename_key_columns
from metacooler.io.colspec import (
    GENOMES_FORMAT_DEFAULT,
    GENOMES_KEYS,
    dumped_table,
    parse_column_spec,
    resolve_column_spec,
)

logger = logging.getLogger("metacooler")


def parse_genomes_table(
    path: str | Path, fmt: str = GENOMES_FORMAT_DEFAULT
) -> pd.DataFrame:
    """Parse a genomes metadata table.

    One column holds the genome names; every other column is metadata, stored
    verbatim in the ``genomes`` table.

    Args:
        path: Path to the genomes TSV file.
        fmt: Which column holds the genome names, in the grammar of
            :mod:`metacooler.io.colspec`: a 1-based position or a header name,
            plus ``,noheader`` when the table has no header line (its metadata
            columns are then named ``col2``, ``col3``, ...).  Defaults to the
            column named ``genome_name``, wherever it is.

    Returns:
        DataFrame with a ``genome_id`` column (the genome names, relabelled in
        place) and every other column of the table.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If *fmt* is malformed, the file is empty, its header and
            data rows disagree on the number of fields, or the column *fmt*
            names is not in it.
    """
    spec = parse_column_spec(fmt, keys=GENOMES_KEYS, option="--genomes-format")
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Genomes table not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Genomes table is empty: {path}")
    df = read_user_table(path, kind="Genomes table", header=spec.header)
    if spec.header and dumped_table(df.columns) == "genomes":
        df = drop_dump_columns(df, "genomes", kind="Genomes table", path=path)
    (key,) = resolve_column_spec(
        spec, list(df.columns), kind="Genomes table", path=path, keys=GENOMES_KEYS
    )
    assert key is not None  # parse_column_spec refuses NA for a genomes table
    # Positional, never by name: see rename_key_columns.
    df = rename_key_columns(df, {key: "genome_id"}, kind="Genomes table")
    dup = df["genome_id"].duplicated()
    if dup.any():
        names = df.loc[dup, "genome_id"].unique().tolist()
        logger.warning(
            "Genomes table contains %d duplicate genome name(s): %s (last row wins).",
            len(names),
            names[:10],
        )
    return df
