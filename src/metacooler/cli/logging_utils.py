"""Logging of the metacooler command line: to stderr, and, for a command that
writes a map or a file, to its log beside what it writes (CONVENTIONS.md § 11).

The log is ``<sample>.metacooler_<YYYYMMDDHHMMSS>.log``, in the directory of
the command's output: ``<sample>`` is the output's name, as
:func:`metacooler.uri.sample_name` gives it, and the time is when the command
started, in local time. It starts with hicstuff's header (the version, the
date, the command line, one line per parameter, then ``---``), holds a line
per record, and ends with what the command wrote and how long it took. A
command that only reads logs to stderr only.

Only the command line configures handlers: :func:`setup_logging` at the start
of every command, :func:`run_log` around the work of one that writes.
"""

from __future__ import annotations

import contextlib
import logging
import shlex
import sys
import time
from collections.abc import Generator
from datetime import datetime
from pathlib import Path

#: The logger of every metacooler module.
LOGGER_NAME = "metacooler"
#: Where the root command keeps its arguments as given, for the log's header
#: (``click.Context.meta``, which every context of a command line shares).
ARGV_KEY = "metacooler.argv"
#: A line of the log file, as hicstuff and metaTOR write theirs.
FILE_FORMATTER = logging.Formatter(
    "%(asctime)s :: %(levelname)s :: %(message)s", datefmt="%Y-%m-%d,%H:%M:%S"
)

#: Whether the command asked for debug logging (-v), which the log file gets too.
_state = {"verbose": False}


class MetacoolerFormatter(logging.Formatter):
    """Custom formatter for metacooler CLI logging.

    Renders ``[timestamp] [command] [LEVEL] message``, then appends the
    exception traceback and stack info when the record carries them — without
    that, ``logger.exception(...)`` silently degraded to a bare one-liner.
    """

    def __init__(self, command_name: str = "metacooler") -> None:
        self.command_name = command_name
        super().__init__()

    def format(self, record: logging.LogRecord) -> str:
        timestamp = datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
        line = (
            f"[{timestamp}] [{self.command_name}] "
            f"[{record.levelname}] {record.getMessage()}"
        )
        if record.exc_info and not record.exc_text:
            record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            line = f"{line}\n{record.exc_text}"
        if record.stack_info:
            line = f"{line}\n{self.formatStack(record.stack_info)}"
        return line


def setup_logging(
    verbose: bool = False,
    command_name: str = "metacooler",
    quiet: bool = False,
) -> None:
    """Set up consistent logging for metacooler CLI commands.

    Args:
        verbose: If True, log debug messages, to stderr and to the log file.
        command_name: Name of the command for the log prefix.
        quiet: If True, suppress everything below WARNING on stderr. Takes
            precedence over *verbose* there; the log file still gets
            everything from INFO up (see :func:`run_log`).
    """
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(MetacoolerFormatter(command_name))

    root = logging.getLogger(LOGGER_NAME)
    # INFO at least, whatever -q says: -q quietens stderr, never the log file.
    root.setLevel(logging.DEBUG if verbose else logging.INFO)
    # Remove existing handlers to avoid duplicates
    root.handlers = []
    root.addHandler(handler)
    # Own our own output: without this, a host application that configured the
    # root logger gets every metacooler record a second time in its own format.
    root.propagate = False
    _state["verbose"] = verbose


def duration(seconds: float) -> str:
    """A duration, in whole seconds, as hicstuff's and metaTOR's logs write it:
    ``0h 17m 31s``."""
    hours, rest = divmod(round(seconds), 3600)
    minutes, seconds = divmod(rest, 60)
    return f"{hours}h {minutes}m {seconds}s"


def log_file_name(
    output: str | Path,
    started: float,
    *,
    directory: str | Path | None = None,
    sample: str | None = None,
) -> Path:
    """The log of a command that writes *output*: in the output's directory,
    ``<sample>.metacooler_<YYYYMMDDHHMMSS>.log``, ``<sample>`` being the
    output's name (:func:`metacooler.uri.sample_name`) and the time when the
    command started (*started*, seconds since the epoch), in local time.

    Args:
        output: The map, file, directory or prefix the command writes, or the
            map it edits in place (a URI of an entry names the entry).
        started: When the command started, as :func:`time.time` gives it.
        directory: The log's directory, when it is not the output's.
        sample: The sample's name, when it is not the output's.

    Returns:
        The log file's path.
    """
    from metacooler.uri import parse_uri, sample_name

    file_path, _ = parse_uri(str(output))
    stamp = time.strftime("%Y%m%d%H%M%S", time.localtime(started))
    sample = sample_name(output) if sample is None else sample
    name = f"{sample}.metacooler_{stamp}.log" if sample else f"metacooler_{stamp}.log"
    return (file_path.parent if directory is None else Path(directory)) / name


def _header(started: float) -> list[str]:
    """The header of a log, as hicstuff writes it: the version, the date, the
    command line, and the value of each parameter of the running command (as
    given, or by default), by its long option's name, then ``---``."""
    import click

    from metacooler._version import __version__

    ctx = click.get_current_context(silent=True)
    argv = ctx.meta.get(ARGV_KEY, sys.argv[1:]) if ctx is not None else sys.argv[1:]
    lines = [
        f"## metacooler: v{__version__} log file",
        f"## date: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(started))}",
        f"## cmd: {shlex.join(['metacooler', *map(str, argv)])}",
    ]
    if ctx is not None:
        for param in ctx.command.params:
            if param.name is None or param.name not in ctx.params:
                continue
            long = [opt[2:] for opt in param.opts if opt.startswith("--")]
            value = ctx.params[param.name]
            if isinstance(value, (list, tuple)):
                value = ",".join(map(str, value))
            lines.append(f"## {long[0] if long else param.name}: {value}")
    lines.append("---")
    return lines


def _emit(handler: logging.Handler, level: int, message: str) -> None:
    """Write one record to *handler* alone: to the log file, not to stderr."""
    record = logging.getLogger(LOGGER_NAME).makeRecord(
        LOGGER_NAME, level, __file__, 0, message, (), None
    )
    handler.handle(record)


@contextlib.contextmanager
def run_log(
    output: str | Path | None,
    what: str,
    *,
    directory: str | Path | None = None,
    sample: str | None = None,
    mkdir: bool = False,
) -> Generator[Path | None]:
    """The log of a command that writes *output*, around its work.

    Starts the log file (see :func:`log_file_name`) with hicstuff's header,
    sends it every record of the ``metacooler`` logger from INFO up (from
    DEBUG with -v; -q quietens stderr only), and ends it with *what*, how long
    the work took and where it went (``Map written after 0h 0m 2s:
    out.metacool``), or with the error the work raised. The header and the
    last line go to the file only: stderr shows what it always showed.

    Nothing is logged to a file when *output* is None (the command writes to
    stdout), nor when the output's directory does not exist (the command then
    fails as it would have), nor when the log cannot be opened there (said on
    stderr).

    Args:
        output: The map, file, directory or prefix the command writes, or the
            map it edits in place; None when it writes nothing but stdout.
        what: What the command writes, for the last line: ``Map written``.
        directory: The log's directory, when it is not the output's: a
            report's, which is its output.
        sample: The sample's name, when it is not the output's: the map's,
            for its report.
        mkdir: Create the log's directory first, as the command does.

    Yields:
        The log file's path, or None when there is none.
    """
    if output is None:
        yield None
        return
    started, clock = time.time(), time.monotonic()
    path = log_file_name(output, started, directory=directory, sample=sample)
    logger = logging.getLogger(LOGGER_NAME)
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    if not path.parent.is_dir():
        yield None
        return
    try:
        handler = logging.FileHandler(path, mode="a", encoding="utf-8")
    except OSError as error:
        logger.warning("No log file for this run: %s", error)
        yield None
        return
    handler.setLevel(logging.DEBUG if _state["verbose"] else logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    for line in _header(started):
        _emit(handler, logging.INFO, line)
    handler.setFormatter(FILE_FORMATTER)
    logger.addHandler(handler)
    try:
        yield path
    except BaseException as error:
        took = duration(time.monotonic() - clock)
        # An exit with status 0 (click's ctx.exit(), sys.exit(0)) is no failure.
        if getattr(error, "exit_code", getattr(error, "code", 1)) in (0, None):
            _emit(handler, logging.INFO, f"{what} after {took}: {output}")
        else:
            reason = str(error) or type(error).__name__
            _emit(handler, logging.ERROR, f"Failed after {took}: {reason}")
        raise
    else:
        took = duration(time.monotonic() - clock)
        _emit(handler, logging.INFO, f"{what} after {took}: {output}")
    finally:
        logger.removeHandler(handler)
        handler.close()
