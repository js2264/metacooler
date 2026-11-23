"""Every panel of a help page starts its help text at the same column.

rich-click lays out each panel (``Options``, then one per command group) as its
own expanding table, and rich spreads an expanding table's spare width over
every column in proportion to its content, so the help text of each panel
started wherever its contents put it. ``cli/_help_layout.py`` aligns them.
"""

from __future__ import annotations

import io
import re

import pytest
import rich_click.rich_click as rc
from _cli_output import _ANSI_RE
from click.testing import CliRunner
from rich.box import SQUARE, Box
from rich.console import Console
from rich.containers import Renderables
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from metacooler.cli._help_layout import align_help_columns
from metacooler.cli.main import cli


def _pages(group, prefix=()):
    """Yield ``(argv, command)`` for every help page of the CLI."""
    yield prefix, group
    for name, cmd in group.commands.items():
        if hasattr(cmd, "commands"):
            yield from _pages(cmd, (*prefix, name))
        else:
            yield (*prefix, name), cmd


def _help_columns(output: str, helps: list[str]) -> set[int]:
    """Column at which each of *helps* starts inside a panel of *output*."""
    lines = [line for line in _ANSI_RE.sub("", output).splitlines() if "│" in line]
    columns = set()
    for text in helps:
        key = " ".join(Text.from_markup(text).plain.split())[:14]
        hits = [line.index(key) for line in lines if key in line]
        if hits:  # a hidden option has help but no row
            columns.add(min(hits))
    return columns


@pytest.fixture
def width(monkeypatch):
    # rich-click reads its width once per root context, from these globals.
    monkeypatch.setattr(rc, "WIDTH", 100)
    monkeypatch.setattr(rc, "MAX_WIDTH", 100)
    return 100


@pytest.mark.parametrize(
    ("argv", "command"),
    [pytest.param(*page, id=" ".join(page[0]) or "metacooler") for page in _pages(cli)],
)
def test_every_help_page_is_aligned(width, argv, command):
    result = CliRunner().invoke(cli, [*argv, "--help"])
    assert result.exit_code == 0, result.output

    helps = [p.help for p in command.params if getattr(p, "help", None)]
    helps.append("Show this message and exit.")
    if hasattr(command, "commands"):
        helps += [c.get_short_help_str(200) for c in command.commands.values()]
    columns = _help_columns(result.output, helps)
    assert len(columns) == 1, f"help text starts at columns {sorted(columns)}"


def test_root_options_line_up_with_command_groups(width):
    result = CliRunner().invoke(cli, ["--help"])
    columns = _help_columns(
        result.output,
        ["Show the version and exit.", "Merge metacool files", "Extract a standard"],
    )
    assert len(columns) == 1


def test_extract_is_listed_on_one_line(width):
    lines = _ANSI_RE.sub("", CliRunner().invoke(cli, ["--help"]).output).splitlines()
    row = next(i for i, line in enumerate(lines) if line.startswith("│ extract "))
    help_text = lines[row].split("extract", 1)[1].rstrip(" │")
    assert help_text.strip() == (
        "Extract a standard cool file for specific genome(s)/contig(s)."
    )
    # The next line is the next command's row, not the help wrapped over.
    assert re.match(r"│ [a-z]", lines[row + 1]), lines[row + 1]


# --------------------------------------------------------------------------
# align_help_columns on hand-built rich-click-like panels
# --------------------------------------------------------------------------


def _table(
    rows,
    *,
    box: Box | None = None,
    show_edge: bool = True,
    pad_edge: bool = False,
    collapse_padding: bool = False,
    expand: bool = True,
) -> Table:
    """A table styled as rich-click styles its panels (by default).

    Its columns carry no header, as under rich-click before 1.9.9.
    """
    table = Table(
        show_header=False,
        padding=(0, 1),
        box=box,
        show_edge=show_edge,
        pad_edge=pad_edge,
        collapse_padding=collapse_padding,
        expand=expand,
    )
    for row in rows:
        table.add_row(*row)
    return table


def _options(**styles) -> Table:
    return _table(
        [("*", "--output", "-o", "Output path."), ("", "--verbose", "", "Talk.")],
        **styles,
    )


def _commands(**styles) -> Table:
    return _table(
        [("go", "Run it."), ("aggregate", "Sum contacts.")],
        **styles,
    )


def _render(panels) -> str:
    console = Console(width=70, file=io.StringIO(), color_system=None)
    align_help_columns(panels, console)
    for panel in panels:
        console.print(panel)
    file = console.file
    assert isinstance(file, io.StringIO)
    return file.getvalue()


_HELPS = ["Output path.", "Talk.", "Run it.", "Sum contacts."]


@pytest.mark.parametrize(
    "styles",
    [
        pytest.param({}, id="default"),
        pytest.param({"expand": False}, id="not-expanding"),
        pytest.param({"pad_edge": True}, id="pad-edge"),
        pytest.param({"collapse_padding": True}, id="collapse-padding"),
        pytest.param({"box": SQUARE}, id="boxed"),
        pytest.param({"box": SQUARE, "show_edge": False}, id="boxed-no-edge"),
    ],
)
def test_panels_share_one_help_column(styles):
    output = _render([Panel(_options(**styles)), Panel(_commands(**styles))])
    assert len(_help_columns(output, _HELPS)) == 1


@pytest.mark.parametrize(
    "styles",
    [
        pytest.param({"pad_edge": True}, id="pad-edge"),
        pytest.param({"collapse_padding": True}, id="collapse-padding"),
        pytest.param({"pad_edge": True, "collapse_padding": True}, id="both"),
    ],
)
def test_panels_with_different_paddings_share_one_help_column(styles):
    # rich 13 and 15 pad these differently: the padding must be rich's own.
    panels = [
        Panel(_options(), padding=(0, 1)),
        Panel(_commands(**styles), padding=(0, 4)),
    ]
    assert len(_help_columns(_render(panels), _HELPS)) == 1


def test_table_under_panel_help_text_is_aligned():
    nested = Renderables([Text("Panel help."), _commands()])
    output = _render([Panel(_options()), Panel(nested)])
    assert len(_help_columns(output, _HELPS)) == 1


def test_only_the_help_column_takes_the_spare_width():
    table = _commands()
    _render([Panel(table)])
    assert [c.width for c in table.columns] == [len("aggregate"), None]
    assert table.columns[1].ratio == 1


def test_one_column_table_is_left_alone():
    table = _table([("a",), ("b",)])
    _render([Panel(table), Panel(Text("not a table"))])
    assert [(c.width, c.ratio) for c in table.columns] == [(None, None)]


def test_rich_without_padding_rule_keeps_rich_click_layout(monkeypatch):
    monkeypatch.delattr(Table, "_get_padding_width")
    table = _commands()
    align_help_columns([Panel(table)], Console(file=io.StringIO()))
    assert [(c.width, c.ratio) for c in table.columns] == [(None, None)] * 2
