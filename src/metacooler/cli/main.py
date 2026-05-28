"""metacooler CLI root."""

from __future__ import annotations

import rich_click as click

from metacooler._version import __version__


class MetaCoolURI(click.ParamType):
    """Click parameter type for metacool/mmetacool URI arguments.

    Accepts the ``file.metacool`` or ``file.mmetacool::entries/entry`` URI
    syntax, validates that the file portion exists, and returns the raw
    URI string (so downstream code can call ``parse_uri()`` on it).
    """

    name = "METACOOL_URI"

    def convert(self, value: str, param, ctx) -> str:
        from metacooler.uri import parse_uri

        file_path, _ = parse_uri(value)
        if not file_path.exists():
            self.fail(f"Path '{file_path}' does not exist.", param, ctx)
        return str(value)


click.rich_click.TEXT_MARKUP = "rich"
click.rich_click.OPTIONS_TABLE_COLUMN_TYPES = [
    "required",
    "opt_long",
    "opt_short",
    "help",
]
click.rich_click.COMMAND_GROUPS = {
    "metacooler": [
        {
            "name": "Create and edit",
            "commands": ["create", "merge", "binnify", "zoomify"],
        },
        {
            "name": "Extract and view",
            "commands": ["extract", "flatten", "dump", "network", "aggregate", "show"],
        },
        {
            "name": "Inspect",
            "commands": ["info", "tree", "metadata", "attrs"],
        },
        {
            "name": "Manage collections",
            "commands": ["ls", "cp", "mv", "rm", "stack"],
        },
    ]
}


@click.group(
    context_settings={"help_option_names": ["-h", "--help"]},
    help="""
[b]MetaCooler[/b]: Hi-C contact maps for metagenomes.

An extension of the cooler format for multi-genome data.
""",
)
@click.version_option(version=__version__, prog_name="metacooler")
def cli() -> None:
    """metacooler CLI."""
    pass


# Import subcommands so they register
from metacooler.cli import (  # noqa: E402, F401
    cmd_aggregate,
    cmd_attrs,
    cmd_binnify,
    cmd_cp,
    cmd_create,
    cmd_create_cool,
    cmd_create_pairs,
    cmd_dump,
    cmd_extract,
    cmd_flatten,
    cmd_info,
    cmd_ls,
    cmd_merge,
    cmd_metadata,
    cmd_mv,
    cmd_network,
    cmd_rm,
    cmd_show,
    cmd_stack,
    cmd_tree,
    cmd_zoomify,
)
