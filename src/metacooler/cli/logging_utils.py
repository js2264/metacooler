"""Logging utilities for metacooler CLI."""

from __future__ import annotations

import logging
import sys
from datetime import datetime


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
        verbose: If True, set log level to DEBUG.
        command_name: Name of the command for the log prefix.
        quiet: If True, suppress everything below WARNING.  Takes precedence
            over *verbose*.
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

    root = logging.getLogger("metacooler")
    root.setLevel(level)
    # Remove existing handlers to avoid duplicates
    root.handlers = []
    root.addHandler(handler)
    # Own our own output: without this, a host application that configured the
    # root logger gets every metacooler record a second time in its own format.
    root.propagate = False
