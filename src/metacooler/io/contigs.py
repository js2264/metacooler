"""Contig table parsers."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from metacooler.io._utils import drop_dump_columns, read_user_table, rename_key_columns
from metacooler.io.colspec import (
    ASSIGNMENT_KEYS,
    CONTIGS_FORMAT_ALIASES,
    CONTIGS_FORMAT_DEFAULT,
    dumped_table,
    parse_column_spec,
    resolve_column_spec,
)

logger = logging.getLogger("metacooler")


def _check_not_empty(path: Path, kind: str) -> None:
    """Raise a clear error if *path* is empty (0 bytes).

    Args:
        path: File to check.
        kind: Human-readable table name used in the error message.

    Raises:
        ValueError: If the file is empty.
    """
    if path.stat().st_size == 0:
        raise ValueError(f"{kind} is empty: {path}")


def _warn_duplicate_contigs(df: pd.DataFrame) -> None:
    """Log a warning if the ``contig_name`` column contains duplicates.

    Duplicate contig names are tolerated (last assignment wins, matching the
    ``dict(zip)`` mapping), but they usually indicate a malformed table, so a
    warning is emitted. Downstream table builders de-duplicate before merging
    to avoid phantom contig rows.
    """
    dup = df["contig_name"].duplicated()
    if dup.any():
        names = df.loc[dup, "contig_name"].unique().tolist()
        logger.warning(
            "Contig table contains %d duplicate contig name(s): %s "
            "(last assignment wins).",
            len(names),
            names[:10],
        )


def _name_mapping(df: pd.DataFrame) -> dict[str, str]:
    """``{contig_name: genome_id}`` of a parsed table (last row wins).

    Zipped over NumPy object arrays: iterating pandas' Arrow-backed string
    columns costs ~1 us per value, four times as much, which on a
    two-million-contig table is seconds spent building one dict.
    """
    return dict(
        zip(
            df["contig_name"].to_numpy(dtype=object),
            df["genome_id"].to_numpy(dtype=object),
            strict=True,
        )
    )


def parse_contigs_table(
    path: str | Path,
    fmt: str = CONTIGS_FORMAT_DEFAULT,
) -> tuple[dict[str, str], pd.DataFrame]:
    """Parse a contig-to-genome assignment table.

    Args:
        path: Path to the contig table file.
        fmt: Which columns hold the contig and the genome names, in the grammar
            of :mod:`metacooler.io.colspec`: ``CONTIG,GENOME[,noheader]``, each
            column a 1-based position or a header name -- e.g.
            ``"Name,Final_cluster"``, ``"2,5"`` or ``"1,2,noheader"`` -- or one of
            :data:`~metacooler.io.colspec.CONTIGS_FORMAT_ALIASES`.  The default,
            ``contigs``, is ``contig_name,genome_name``: looked up by name, so
            a table given without its format is refused rather than misread.
            ``NA`` in the genome's place (``"1,NA,noheader"``,
            ``"contig_name,NA"``) says the table has no genome column: every
            contig it lists is assigned to ``__unbinned__``.

    Returns:
        Tuple of (mapping ``contig_name -> genome_name``, full DataFrame with
        the two key columns relabelled ``contig_name`` / ``genome_id`` in
        place, every other column kept -- a headerless table's as ``colN``).
        Contigs not present in the table are omitted from the mapping (the
        caller decides whether to assign them to ``__unbinned__`` or drop
        them).

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If *fmt* is malformed, the table is empty or misaligned,
            or the columns *fmt* names are not in it.
    """
    spec = parse_column_spec(
        fmt,
        keys=ASSIGNMENT_KEYS,
        aliases=CONTIGS_FORMAT_ALIASES,
        option="--contigs-format",
    )
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Contig table not found: {path}")
    _check_not_empty(path, "Contig table")

    df = read_user_table(path, kind="Contig table", header=spec.header)
    if spec.header and dumped_table(df.columns) == "contigs":
        df = drop_dump_columns(df, "contigs", kind="Contig table", path=path)
    contig_col, genome_col = resolve_column_spec(
        spec,
        list(df.columns),
        kind="Contig table",
        path=path,
        keys=ASSIGNMENT_KEYS,
        hint=f", or an alias: {', '.join(CONTIGS_FORMAT_ALIASES)}",
    )
    assert contig_col is not None  # parse_column_spec refuses NA in its place
    if genome_col is None:
        # NA: the table assigns no genome.  Give it a genome column that says
        # so on every row, and the rest of the pipeline reads an ordinary
        # assignment table; a column of the table's own that happens to be
        # called genome_id is moved aside by rename_key_columns, as ever.
        genome_col = df.shape[1]
        df = df.copy()
        df.insert(genome_col, "__no_genome__", "__unbinned__", allow_duplicates=True)
    # Positional, never by name: see rename_key_columns.
    df = rename_key_columns(
        df, {contig_col: "contig_name", genome_col: "genome_id"}, kind="Contig table"
    )
    _warn_duplicate_contigs(df)
    return _name_mapping(df), df
