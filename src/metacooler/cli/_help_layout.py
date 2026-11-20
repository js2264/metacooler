"""Start the help text of every panel of a help page at the same column.

rich-click draws each panel of a help page (``Options``, then one per
``COMMAND_GROUPS`` entry) as its own expanding rich ``Table``. rich hands an
expanding table's spare width to *every* column, in proportion to its content,
so a panel of short help texts pushes its name column further right than a
panel whose help wraps: no two panels started their help at the same offset.

:func:`align_help_columns` fixes the columns left of the help column of every
table to one shared width, the widest such prefix on the page, and makes the
help column the only one that takes up the rest of the line.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table

# rich-click titles each table column after its column type ("help" -> "Help").
HELP_HEADER = "Help"


def _tables(renderable: object) -> Iterator[Table]:
    """Yield the tables a panel draws (a panel with help text nests its table)."""
    from rich.containers import Renderables
    from rich.table import Table

    if isinstance(renderable, Table):
        yield renderable
    elif isinstance(renderable, Renderables):
        # rich types ``Renderables.__iter__`` as returning an Iterable.
        for inner in renderable:  # ty: ignore[not-iterable]
            yield from _tables(inner)


def _padding_width(table: Table, index: int) -> int:
    """Horizontal padding rich adds to column *index* of *table*.

    The rule of rich's ``Table._get_padding_width``, read from public
    attributes: a private method that moved would break every ``--help``.
    Only columns left of the help column are asked about, so never the last
    one, whose right padding the rule would also drop without ``pad_edge``.
    """
    _, right, _, left = table.padding
    if table.collapse_padding or (index == 0 and not table.pad_edge):
        left = 0
    return left + right


def align_help_columns(panels: Sequence[Panel], console: Console) -> None:
    """Give every table in *panels* the same offset for its help column.

    Edits the tables in place: each column left of the help column gets a
    fixed width (its content; the last one also gets whatever the table lacks
    to reach the page's widest prefix), and the help column gets the only
    ``ratio``, so an expanding table widens it alone. A table without a help
    column, or with nothing left of it, is left as rich-click built it.
    """
    from rich.measure import Measurement
    from rich.padding import Padding

    options = console.options
    layouts: list[tuple[Table, int, list[int], int]] = []
    for panel in panels:
        indent = Padding.unpack(panel.padding)[3]
        for table in _tables(panel.renderable):
            help_index = next(
                (i for i, c in enumerate(table.columns) if c.header == HELP_HEADER),
                0,
            )
            if help_index == 0:
                continue
            widths = [
                max(
                    (Measurement.get(console, options, c).maximum for c in col.cells),
                    default=0,
                )
                for col in table.columns[:help_index]
            ]
            offset = indent + sum(
                width + _padding_width(table, i) for i, width in enumerate(widths)
            )
            if table.box:
                # One rule after each column left of the help, plus the edge.
                offset += help_index + int(table.show_edge)
            layouts.append((table, help_index, widths, offset))

    if not layouts:
        return
    target = max(offset for *_, offset in layouts)
    for table, help_index, widths, offset in layouts:
        widths[-1] += target - offset
        for column, width in zip(table.columns, widths, strict=False):
            column.width = width
        table.columns[help_index].ratio = 1
