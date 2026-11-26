"""metacooler CLI root."""

from __future__ import annotations

import contextlib
import os
import signal
import sys
from typing import NoReturn

import rich_click as click

from metacooler._version import __version__

# Exit status for "the reader closed the pipe", following the shell convention
# of 128 + the signal number a C program would have died from. Windows has no
# SIGPIPE, so its number is spelled out as a fallback.
_BROKEN_PIPE_EXIT_CODE = 128 + getattr(signal, "SIGPIPE", 13)


class MetaCoolURI(click.ParamType):
    """Click parameter type for metacool/mmetacool URI arguments.

    Accepts the ``file.metacool`` or ``file.mmetacool::entries/entry`` URI
    syntax, validates that the file portion exists, and returns the raw
    URI string (so downstream code can call ``parse_uri()`` on it).

    This is the single place where the "what does a bare ``.mmetacool`` path
    mean?" policy lives.  Commands that operate on one map (most of them) pass
    ``require_entry=True`` so a bare collection path is rejected up front with
    the list of available entries, instead of each command hand-rolling its own
    check — or, worse, silently operating on the collection root.

    Args:
        require_entry: If True, reject a URI that resolves to the root of an
            ``.mmetacool`` collection.
        exists: If False, do not require the file to exist (used for
            destinations).
    """

    name = "METACOOL_URI"

    def __init__(self, require_entry: bool = False, exists: bool = True) -> None:
        self.require_entry = require_entry
        self.exists = exists

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
            click.exceptions.BadParameter: If the file portion of the URI does
                not exist, the group is absent, or ``require_entry`` is set and
                the URI addresses a collection root.
        """

        from metacooler.uri import parse_uri

        file_path, group_path = parse_uri(value)
        if not file_path.exists():
            if self.exists:
                self.fail(f"Path '{file_path}' does not exist.", param, ctx)
            return str(value)
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
        elif self.require_entry:
            self._fail_on_collection_root(file_path, param, ctx)
        return str(value)

    def _fail_on_collection_root(self, file_path, param, ctx) -> None:
        """Fail with the available entry URIs if *file_path* is an mmetacool.

        Args:
            file_path: Path portion of the URI, known to exist.
            param: The Click ``Parameter`` object.
            ctx: The Click ``Context`` object.
        """
        from metacooler.schema import is_mmetacooler, list_metacoolers

        try:
            if not is_mmetacooler(file_path):
                return
        except OSError:
            return
        entries = list_metacoolers(file_path)
        listing = "\n".join(f"  {file_path}::{e}" for e in entries)
        self.fail(
            f"'{file_path}' is an mmetacool collection with {len(entries)} "
            "entries; this command needs a single map. Specify an entry with "
            f"URI syntax:\n{listing}",
            param,
            ctx,
        )


# Domain exceptions the io/ops layer raises for bad user input. These are turned
# into clean ``Error: <message>`` output instead of a Python traceback. Anything
# NOT listed here (e.g. RuntimeError, AttributeError, TypeError) still surfaces
# as a traceback, since those indicate bugs rather than user mistakes.
#
# ``OSError`` covers the whole filesystem family (FileNotFoundError,
# FileExistsError, IsADirectoryError, PermissionError) *and* h5py's "not an
# HDF5 file" error, which previously reached the user as a raw traceback in
# every command that opens a file without pre-checking it.
_USER_ERRORS = (
    ValueError,
    KeyError,
    OSError,
    NotImplementedError,
)


def _user_message(exc: BaseException) -> str:
    """The message of a user error, as the command line prints it.

    A ``KeyError`` stringifies with surrounding quotes: its message is used.
    An ``OSError`` made with an errno prints ``[Errno N]`` before its message
    (``FileNotFoundError(errno.ENOENT, "...")``, h5py's own): the message is
    used, with the file it names.
    """
    if isinstance(exc, KeyError) and exc.args:
        return str(exc.args[0])
    if isinstance(exc, OSError) and exc.errno is not None and exc.strerror:
        files = [name for name in (exc.filename, exc.filename2) if name is not None]
        return exc.strerror + (": " + " -> ".join(map(str, files)) if files else "")
    return str(exc)


def _exit_on_broken_pipe() -> NoReturn:
    """Exit quietly after the reader of our stdout went away.

    ``metacooler dump … | head`` closes the pipe as soon as ``head`` has its
    lines; that is the normal end of a pipeline, not a failure, but the
    resulting ``BrokenPipeError`` is an ``OSError`` and so used to be reported
    in the same red ``Error: [Errno 32] Broken pipe`` box as a genuinely bad
    input — which reads as "the dump failed".

    Python's interpreter shutdown would flush ``sys.stdout`` one more time and
    raise again, printing an unsuppressable ``Exception ignored in: …``, so
    stdout is pointed at the null device before exiting.  The status follows the
    shell convention for a SIGPIPE death, which is what a C tool in the same
    position would report.
    """
    with contextlib.suppress(OSError, ValueError):
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
    raise SystemExit(_BROKEN_PIPE_EXIT_CODE)


class MetaCoolerCommand(click.RichCommand):
    """Command whose help starts every panel's help text at the same column.

    rich-click lays each panel out on its own, so the help column of one
    panel lands wherever its contents put it (see ``cli/_help_layout.py``).
    This is rich-click's ``format_options`` with one extra step: align the
    rendered panels before writing them.
    """

    def format_options(self, ctx: click.Context, formatter: click.HelpFormatter):
        from rich.table import Table
        from rich_click.rich_panel import construct_panels

        from metacooler.cli._help_layout import align_help_columns

        assert isinstance(ctx, click.RichContext)
        assert isinstance(formatter, click.RichHelpFormatter)
        panels = [
            panel.render(self, ctx, formatter)
            for panel in construct_panels(self, ctx, formatter)
        ]
        # rich-click does not draw a panel whose table has no row.
        panels = [
            p
            for p in panels
            if not isinstance(p.renderable, Table) or p.renderable.rows
        ]
        align_help_columns(panels, formatter.console)
        for panel in panels:
            formatter.write(panel)


class MetaCoolerGroup(MetaCoolerCommand, click.RichGroup):
    """Root group that reports expected domain errors as clean CLI messages.

    Every subcommand (including nested ``create`` subcommands) is invoked
    through this group, so wrapping :meth:`invoke` gives uniform error output
    without a ``try/except`` in each ``cmd_*`` module.

    Its commands are :class:`MetaCoolerCommand`, and its subgroups (``create``,
    ``metadata``) are of its own class, so every help page is aligned.
    """

    command_class = MetaCoolerCommand
    group_class = type

    def make_context(self, info_name, args, parent=None, **extra):
        """Keep the command line's arguments as given, for the header of a
        command's log (``cli/logging_utils.py``): parsing consumes them."""
        given = list(args)
        ctx = super().make_context(info_name, args, parent=parent, **extra)
        if parent is None:
            from metacooler.cli.logging_utils import ARGV_KEY

            ctx.meta[ARGV_KEY] = given
        return ctx

    def invoke(self, ctx):
        try:
            result = super().invoke(ctx)
            # Flush inside the ``try`` so a pipe that closed while the last
            # buffered block was still pending is diagnosed here rather than
            # during interpreter shutdown, where it cannot be handled.
            sys.stdout.flush()
            return result
        except BrokenPipeError:
            # Must precede the _USER_ERRORS clause: BrokenPipeError is an
            # OSError, and this is not a user mistake to report.
            _exit_on_broken_pipe()
        except (click.ClickException, click.Abort):
            raise  # already a clean CLI error / user abort
        except _USER_ERRORS as exc:
            raise click.ClickException(_user_message(exc)) from exc


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
            "name": "Create",
            "commands": ["create", "collect"],
        },
        {
            "name": "Edit an existing map",
            "commands": ["reassign", "annotate", "normalize", "metadata"],
        },
        {
            "name": "Combine and re-bin",
            "commands": ["merge", "binnify", "zoomify"],
        },
        {
            "name": "Extract and view",
            "commands": [
                "extract",
                "fasta",
                "flatten",
                "dump",
                "network",
                "aggregate",
                "show",
                "ps",
            ],
        },
        {
            "name": "Inspect",
            "commands": ["info", "qc", "report", "tree", "attrs", "check"],
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
    cmd_annotate,
    cmd_attrs,
    cmd_binnify,
    cmd_check,
    cmd_collect,
    cmd_cp,
    cmd_create,
    cmd_create_cool,
    cmd_create_pairs,
    cmd_dump,
    cmd_extract,
    cmd_fasta,
    cmd_flatten,
    cmd_info,
    cmd_ls,
    cmd_merge,
    cmd_metadata,
    cmd_mv,
    cmd_network,
    cmd_normalize,
    cmd_ps,
    cmd_qc,
    cmd_reassign,
    cmd_report,
    cmd_rm,
    cmd_show,
    cmd_stack,
    cmd_tree,
    cmd_zoomify,
)
