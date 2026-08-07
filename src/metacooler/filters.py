"""Predicate filters over the ``/genomes`` and ``/contigs`` tables.

The CLI surface is a repeatable ``--filter`` option whose value is one
comparison, prefixed with the table it applies to::

    --filter "genomes.length > 1Mb"
    --filter "contigs.length >= 10kb" --filter "contigs.name =~ ^NODE_1"

Repeating the option ANDs the comparisons; there is deliberately no ``or``
and no parenthesised grammar.  The reason is the column names: a real
MetaTOR ``/genomes`` table carries ``Present Markers``, ``GC-content`` and
``Weighted completeness``, and a grammar with boolean connectives has to
quote or backtick those before it can tokenise them.  Here the *first*
comparison operator terminates the column name, so a column with spaces,
dashes or dots needs no quoting at all.

The table prefix is mandatory rather than inferred.  ``length`` exists in
both tables and means something different in each, and a rule that only
demands the prefix when the column happens to be ambiguous would make the
same expression legal on one file and rejected on the next.

Missing values never match.  A row whose value is NA fails every
comparison, including ``!=`` and ``!~`` — a filter selects rows it can
positively vouch for.
"""

from __future__ import annotations

import difflib
import logging
import operator
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

logger = logging.getLogger("metacooler")

#: Tables a ``--filter`` expression may address.
TABLES: tuple[str, ...] = ("genomes", "contigs")

# Longest-first at each position, so ">=" is never read as ">" against a
# value of "=1000", nor "=~" as "=" against "~^NODE".
_OPERATORS: tuple[str, ...] = (">=", "<=", "==", "!=", "=~", "!~", ">", "<", "=")
_ORDER_OPERATORS = frozenset({">", ">=", "<", "<="})
_REGEX_OPERATORS = frozenset({"=~", "!~"})
_COMPARE = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
}

# Size suffixes, because filtering by length is the whole point.  Matched
# case-insensitively: "1Mb", "1mb" and "1M" are the same threshold.
_UNIT_MULTIPLIERS = {
    "b": 1.0,
    "bp": 1.0,
    "k": 1e3,
    "kb": 1e3,
    "m": 1e6,
    "mb": 1e6,
    "g": 1e9,
    "gb": 1e9,
}
_NUMBER_WITH_UNIT_RE = re.compile(r"(?P<number>.*?)(?P<unit>[a-zA-Z]*)")

# How many column names to list when a filter names one that does not exist.
_MAX_REPORTED_COLUMNS = 25
# How many uninterpretable values to quote when a column resists coercion.
_MAX_REPORTED_BAD_VALUES = 3


class FilterError(ValueError):
    """A ``--filter`` expression is malformed or does not apply to the file."""


@dataclass(frozen=True)
class Filter:
    """One parsed ``<table>.<column> <op> <value>`` comparison.

    Attributes:
        table: ``"genomes"`` or ``"contigs"``.
        column: Column name, verbatim — may contain spaces, dashes or dots.
        op: One of ``> >= < <= == != =~ !~`` (a lone ``=`` is normalised to
            ``==``).
        value: Right-hand side as written, with surrounding quotes stripped.
        raw: The original expression, used in every error message.
    """

    table: str
    column: str
    op: str
    value: str
    raw: str

    def __str__(self) -> str:
        return self.raw


# ----------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------
def _find_operator(text: str) -> tuple[int, str] | None:
    """Return the ``(position, operator)`` of the leftmost comparison operator.

    Leftmost rather than rightmost so that a regex value containing ``>`` or
    ``=`` (``contigs.name =~ ^a>b``) stays entirely on the value side, and
    longest-first at each position so two-character operators win over their
    one-character prefixes.
    """
    for i in range(len(text)):
        for op in _OPERATORS:
            if text.startswith(op, i):
                return i, op
    return None


def _strip_quotes(value: str) -> str:
    """Drop one layer of matching surrounding quotes from a filter value."""
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value


def parse_filter(expr: str) -> Filter:
    """Parse one ``--filter`` expression.

    Args:
        expr: Expression of the form ``<table>.<column> <op> <value>``, e.g.
            ``"genomes.length > 1Mb"``.  Whitespace around each part is
            optional: ``"contigs.length>=10kb"`` parses identically.

    Returns:
        The parsed :class:`Filter`.

    Raises:
        FilterError: If the expression has no operator, names no table (or an
            unknown one), or leaves either side empty; also if a ``=~``/``!~``
            value is not a valid regular expression.
    """
    raw = expr.strip()
    if not raw:
        raise FilterError(
            "--filter was given an empty expression. Write one comparison, e.g. "
            "--filter 'genomes.length > 1Mb'."
        )

    found = _find_operator(raw)
    if found is None:
        raise FilterError(
            f"Filter {raw!r} contains no comparison operator. Write "
            f"'<table>.<column> <op> <value>', e.g. 'contigs.length >= 10kb'. "
            f"Operators: {', '.join(_OPERATORS[:-1])}."
        )
    pos, op = found
    lhs = raw[:pos].strip()
    value = _strip_quotes(raw[pos + len(op) :].strip())
    if op == "=":
        op = "=="

    if not lhs:
        raise FilterError(
            f"Filter {raw!r} has nothing on the left of {op!r}. Write "
            f"'<table>.<column> {op} <value>'."
        )
    if not value:
        raise FilterError(f"Filter {raw!r} has nothing on the right of {op!r}.")

    table, sep, column = lhs.partition(".")
    table, column = table.strip(), column.strip()
    if not sep:
        raise FilterError(
            f"Filter {raw!r} must say which table {lhs!r} belongs to: write "
            f"'genomes.{lhs}' or 'contigs.{lhs}'. Both tables have a 'length' "
            f"and a 'name', so the prefix is never optional."
        )
    if table not in TABLES:
        raise FilterError(
            f"Filter {raw!r} addresses an unknown table {table!r}. "
            f"Expected one of: {', '.join(TABLES)}."
        )
    if not column:
        raise FilterError(f"Filter {raw!r} names a table but no column.")

    if op in _REGEX_OPERATORS:
        try:
            re.compile(value)
        except re.error as exc:
            raise FilterError(
                f"Filter {raw!r}: {value!r} is not a valid regular expression ({exc})."
            ) from exc

    return Filter(table=table, column=column, op=op, value=value, raw=raw)


def parse_filters(exprs: Iterable[str] | None) -> list[Filter]:
    """Parse every ``--filter`` expression, preserving order.

    Args:
        exprs: Raw expressions, or ``None`` when the option was not given.

    Returns:
        List of parsed filters (empty when *exprs* is falsy).

    Raises:
        FilterError: On the first expression that does not parse.
    """
    return [parse_filter(e) for e in (exprs or [])]


def filters_for(filters: Sequence[Filter], table: str) -> list[Filter]:
    """Return only the filters addressing *table*."""
    return [f for f in filters if f.table == table]


def describe(filters: Sequence[Filter]) -> str:
    """Render a filter list for a log line or error message."""
    return "; ".join(f.raw for f in filters)


# ----------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------
def parse_size(text: str) -> float:
    """Parse a genomic size, honouring bp/kb/Mb/Gb suffixes.

    The one place the suffix table is interpreted, shared by ``--filter``
    values and by every option that takes a distance (``qc --threshold``,
    ``ps --min-distance``), so ``10kb`` means the same everywhere.

    Args:
        text: A plain number (``1000000``, ``1e6``) or a number with a size
            suffix (``1Mb``, ``10kb``), case-insensitive.

    Returns:
        The value in base pairs.

    Raises:
        ValueError: If *text* is neither.
    """
    try:
        return float(text)
    except ValueError:
        pass
    match = _NUMBER_WITH_UNIT_RE.fullmatch(text)
    if match:
        number, unit = match.group("number"), match.group("unit").lower()
        if number and unit in _UNIT_MULTIPLIERS:
            try:
                return float(number) * _UNIT_MULTIPLIERS[unit]
            except ValueError:
                pass
    raise ValueError(
        f"{text!r} is not a number. Write a plain number "
        f"(1000000, 1e6, 1_000_000) or add a size suffix "
        f"({', '.join(sorted(_UNIT_MULTIPLIERS))})."
    )


def _to_number(text: str, filt: Filter) -> float:
    """Parse a filter value as a number, naming the filter if it is not one.

    Raises:
        FilterError: If *text* is neither a plain number nor a number with a
            recognised size suffix.
    """
    try:
        return parse_size(text)
    except ValueError as exc:
        raise FilterError(f"Filter {filt.raw!r}: {exc}") from None


def _numeric_column(col: pd.Series, filt: Filter) -> np.ndarray:
    """Return *col* as float64, coercing a text column that holds numbers.

    Real contig tables store numeric fields as text with a placeholder for
    missing data (MetaTOR writes ``-`` for an absent shotgun coverage), so a
    column being ``object``-typed is not on its own a reason to refuse an
    order comparison.  Values that cannot be read as numbers become NaN,
    which fails every comparison, and their count is reported once.  Only a
    column where *nothing* is numeric is rejected outright.

    Raises:
        FilterError: If no value in *col* can be read as a number.
    """
    if pd.api.types.is_numeric_dtype(col) and not pd.api.types.is_bool_dtype(col):
        return col.astype(np.float64).to_numpy()

    coerced = pd.to_numeric(col, errors="coerce")
    filled = col.notna() & (col.astype(str).str.strip() != "")
    unreadable = filled & coerced.isna()
    if filled.any() and bool(unreadable.all()):
        examples = col[unreadable].astype(str).unique()[:_MAX_REPORTED_BAD_VALUES]
        raise FilterError(
            f"Filter {filt.raw!r}: column {filt.column!r} is not numeric "
            f"(dtype {col.dtype}) and none of its {int(filled.sum())} value(s) "
            f"could be read as a number (e.g. {list(examples)}). Compare it as "
            f"text with ==, != or =~ instead."
        )
    if unreadable.any():
        examples = col[unreadable].astype(str).unique()[:_MAX_REPORTED_BAD_VALUES]
        logger.warning(
            "Filter %r: %d of %d value(s) in column %r are not numbers "
            "(e.g. %s); those rows cannot satisfy the comparison and are "
            "excluded.",
            filt.raw,
            int(unreadable.sum()),
            len(col),
            filt.column,
            list(examples),
        )
    return coerced.astype(np.float64).to_numpy()


def _unknown_column_message(
    df: pd.DataFrame, filt: Filter, hints: Mapping[str, str] | None
) -> str:
    """Build the error text for a filter naming a column the table lacks."""
    if hints and filt.column in hints:
        return (
            f"Filter {filt.raw!r}: /{filt.table} has no column "
            f"{filt.column!r} at this point — {hints[filt.column]}."
        )
    columns = [str(c) for c in df.columns]
    shown = columns[:_MAX_REPORTED_COLUMNS]
    overflow = len(columns) - len(shown)
    listing = f"{shown}" + (f" ... ({overflow} more)" if overflow else "")
    close = difflib.get_close_matches(filt.column, columns, n=1)
    suggestion = f" Did you mean {close[0]!r}?" if close else ""
    return (
        f"Filter {filt.raw!r}: /{filt.table} has no column "
        f"{filt.column!r}.{suggestion} Available columns: {listing}"
    )


def _single_mask(
    df: pd.DataFrame, filt: Filter, hints: Mapping[str, str] | None
) -> np.ndarray:
    """Evaluate one filter against *df*, returning a boolean row mask."""
    if filt.column not in df.columns:
        raise FilterError(_unknown_column_message(df, filt, hints))

    col = df[filt.column]
    present = col.notna().to_numpy()

    if filt.op in _REGEX_OPERATORS:
        text = col.where(col.notna(), "").astype(str)
        hit = text.str.contains(filt.value, regex=True, na=False).to_numpy()
        return (hit if filt.op == "=~" else ~hit) & present

    if filt.op in _ORDER_OPERATORS:
        # NaN compares False under every order operator, so missing values
        # drop out without a separate mask.
        return _COMPARE[filt.op](
            _numeric_column(col, filt), _to_number(filt.value, filt)
        )

    numeric = pd.api.types.is_numeric_dtype(col) and not pd.api.types.is_bool_dtype(col)
    if numeric:
        equal = col.astype(np.float64).to_numpy() == _to_number(filt.value, filt)
    else:
        equal = col.where(col.notna(), "").astype(str).to_numpy() == filt.value
    return (equal if filt.op == "==" else ~equal) & present


def filter_mask(
    df: pd.DataFrame,
    filters: Sequence[Filter],
    *,
    table: str,
    hints: Mapping[str, str] | None = None,
) -> np.ndarray:
    """Return the boolean row mask of *df* satisfying every filter for *table*.

    Args:
        df: The table to filter.
        filters: Parsed filters.  Those addressing another table are ignored,
            so the caller can hand the whole list to both tables.
        table: Which table *df* is — ``"genomes"`` or ``"contigs"``.
        hints: Optional ``column → explanation`` map used to give a specific
            reason when a column is legitimately absent in this context (e.g.
            ``n_bins`` during ``create``, before any bin has been tiled).

    Returns:
        Boolean array of length ``len(df)``.  All-True when no filter
        addresses *table*.

    Raises:
        FilterError: If a filter names a column *df* does not have, or its
            value cannot be interpreted for that column.
    """
    mask = np.ones(len(df), dtype=bool)
    for filt in filters:
        if filt.table != table:
            continue
        mask &= _single_mask(df, filt, hints)
    return mask


def apply_filters(
    df: pd.DataFrame,
    filters: Sequence[Filter],
    *,
    table: str,
    hints: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Return the rows of *df* satisfying every filter for *table*."""
    return df[filter_mask(df, filters, table=table, hints=hints)]


def explain_empty(
    frames: Mapping[str, pd.DataFrame],
    filters: Sequence[Filter],
    *,
    context: str = "",
    hints: Mapping[str, str] | None = None,
) -> str:
    """Build the error text for a selection that every row failed.

    Evaluates each filter *on its own* so the message names the one that
    emptied the selection, rather than reporting only that the combination
    matched nothing — a mistyped threshold is otherwise indistinguishable
    from a genuinely empty intersection.

    Args:
        frames: ``table name → table`` for the tables that were filtered.
        filters: The parsed filters, in the order given.
        context: Optional trailing sentence naming what else narrowed the
            selection (e.g. an explicit ``--genomes`` list).
        hints: Passed through to :func:`filter_mask`.

    Returns:
        A multi-line message.
    """
    lines = ["Nothing passed --filter. On its own, each filter keeps:"]
    for filt in filters:
        df = frames.get(filt.table)
        if df is None:
            continue
        kept = int(filter_mask(df, [filt], table=filt.table, hints=hints).sum())
        marker = "   <-- keeps nothing" if kept == 0 else ""
        lines.append(f"  {filt.raw}  ->  {kept} of {len(df)} {filt.table}{marker}")
    if context:
        lines.append(context)
    return "\n".join(lines)
