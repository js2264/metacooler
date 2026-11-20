"""Start the help text of every panel of a help page at the same column.

rich-click draws each panel of a help page (``Options``, then one per
``COMMAND_GROUPS`` entry) as its own expanding rich ``Table``. rich hands an
expanding table's spare width to *every* column, in proportion to its content,
so a panel of short help texts pushes its name column further right than a
panel whose help wraps: no two panels started their help at the same offset.

:func:`align_help_columns` fixes the columns left of the help column of every
table to one shared width, the widest such prefix on the page, and makes the
help column the only one that takes up the rest of the line.

The help column is found by position, as the last one: every rich-click theme
puts it there, and rich-click before 1.9.9 gives its columns no header to find
it by.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table


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


def align_help_columns(panels: Sequence[Panel], console: Console) -> None:
    """Give every table in *panels* the same offset for its help column.

    Edits the tables in place: each column left of the help (last) column gets
    a fixed width (its content; the one next to the help also gets whatever
    the table lacks to reach the page's widest prefix), and the help column
    gets the only ``ratio``, so an expanding table widens it alone. A table of
    one column is left as rich-click built it.

    How much padding rich adds to a column is its private
    ``Table._get_padding_width`` (rich 13 and 15 differ on it, under
    ``pad_edge`` and ``collapse_padding``, so no copy of the rule is right for
    both). Should a rich without it come, the page keeps rich-click's layout
    rather than ``--help`` failing.
    """
    from rich.measure import Measurement
    from rich.padding import Padding

    options = console.options
    layouts: list[tuple[Table, list[int], int]] = []
    for panel in panels:
        indent = Padding.unpack(panel.padding)[3]
        for table in _tables(panel.renderable):
            padding_width = getattr(table, "_get_padding_width", None)
            help_index = len(table.columns) - 1
            if padding_width is None or help_index < 1:
                continue
            widths = [
                max(
                    (Measurement.get(console, options, c).maximum for c in col.cells),
                    default=0,
                )
                for col in table.columns[:help_index]
            ]
            # The help text starts after the help column's own left padding,
            # which collapse_padding cuts (rich's _get_cells, 13 to 15 alike).
            _, right, _, left = table.padding
            if table.collapse_padding:
                left = max(0, left - right)
            offset = indent + left
            offset += sum(width + padding_width(i) for i, width in enumerate(widths))
            if table.box:
                # One rule after each column left of the help, plus the edge.
                offset += help_index + int(table.show_edge)
            layouts.append((table, widths, offset))

    if not layouts:
        return
    target = max(offset for *_, offset in layouts)
    for table, widths, offset in layouts:
        widths[-1] += target - offset
        for column, width in zip(table.columns, widths, strict=False):
            column.width = width
        table.columns[-1].ratio = 1
