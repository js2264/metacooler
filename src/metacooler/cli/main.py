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
        """Validate a metacool URI argument and return the string value.

        Args:
            value: Raw string value provided by the user on the command
                line.  Accepts plain file paths or
                ``file.mmetacool::entries/entry`` URIs.
            param: The Click ``Parameter`` object (passed through for error
                reporting).
            ctx: The Click ``Context`` object (passed through for error
                reporting).

        Returns:
            The original *value* string unchanged if the file exists.

        Raises:
            click.exceptions.BadParameter: If the file portion of the URI
                does not exist on disk.
        """

        from metacooler.uri import parse_uri

        file_path, group_path = parse_uri(value)
        if not file_path.exists():
            self.fail(f"Path '{file_path}' does not exist.", param, ctx)
        if group_path != "/":
            import h5py

            try:
                with h5py.File(file_path, "r") as f:
                    if group_path.lstrip("/") not in f:
                        self.fail(
                            f"Group '{group_path}' not found in '{file_path}'.",
                            param,
                            ctx,
                        )
            except OSError as exc:
                self.fail(f"Cannot open '{file_path}': {exc}", param, ctx)
        return str(value)


# Domain exceptions the io/ops layer raises for bad user input. These are turned
# into clean ``Error: <message>`` output instead of a Python traceback. Anything
# NOT listed here (e.g. RuntimeError, AttributeError) still surfaces as a
# traceback, since those indicate bugs rather than user mistakes.
_USER_ERRORS = (
    ValueError,
    KeyError,
    FileNotFoundError,
    FileExistsError,
    IsADirectoryError,
    PermissionError,
    NotImplementedError,
    TypeError,
)


class MetaCoolerGroup(click.RichGroup):
    """Root group that reports expected domain errors as clean CLI messages.

    Every subcommand (including nested ``create`` subcommands) is invoked
    through this group, so wrapping :meth:`invoke` gives uniform error output
    without a ``try/except`` in each ``cmd_*`` module.
    """

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except (click.ClickException, click.Abort):
            raise  # already a clean CLI error / user abort
        except _USER_ERRORS as exc:
            # KeyError stringifies with surrounding quotes; use its message.
            message = (
                str(exc.args[0]) if isinstance(exc, KeyError) and exc.args else str(exc)
            )
            raise click.ClickException(message) from exc


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
    cls=MetaCoolerGroup,
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
