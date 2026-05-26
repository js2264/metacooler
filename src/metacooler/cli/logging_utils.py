"""Logging utilities for metacooler CLI."""

from __future__ import annotations

import logging
import sys
from datetime import datetime


class MetacoolerFormatter(logging.Formatter):
    """Custom formatter for metacooler CLI logging."""

    def __init__(self, command_name: str = "metacooler") -> None:
        self.command_name = command_name
        super().__init__()

    def format(self, record: logging.LogRecord) -> str:
        timestamp = datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
        return f"[{timestamp}] [{self.command_name}] [{record.levelname}] {record.getMessage()}"


def setup_logging(verbose: bool = False, command_name: str = "metacooler") -> None:
    """Set up consistent logging for metacooler CLI commands.

    Args:
        verbose: If True, set log level to DEBUG. Otherwise INFO.
        command_name: Name of the command for the log prefix.
    """
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(MetacoolerFormatter(command_name))

    root = logging.getLogger("metacooler")
    root.setLevel(level)
    # Remove existing handlers to avoid duplicates
    root.handlers = []
    root.addHandler(handler)
