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


class MetaCoolerGroup(click.RichGroup):
    """Root group that reports expected domain errors as clean CLI messages.

    Every subcommand (including nested ``create`` subcommands) is invoked
    through this group, so wrapping :meth:`invoke` gives uniform error output
    without a ``try/except`` in each ``cmd_*`` module.
    """

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
            "commands": [
                "extract",
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
            "commands": ["info", "qc", "tree", "metadata", "attrs", "check"],
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
    cmd_check,
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
    cmd_ps,
    cmd_qc,
    cmd_rm,
    cmd_show,
    cmd_stack,
    cmd_tree,
    cmd_zoomify,
)
