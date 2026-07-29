"""Shared Click option callbacks and output guards for the metacooler CLI.

Single home for the idioms that were previously copy-pasted across
``cmd_*.py`` modules: comma-separated multi-value parsing and the
"refuse to clobber an existing output" check.
"""

from __future__ import annotations

from pathlib import Path

import rich_click as click


def parse_comma_list(ctx, param, value):
    """Parse a comma-separated or repeated multi-value option into a flat list.

    Accepts both ``--genomes a,b`` and ``--genomes a --genomes b``.  An
    explicitly supplied but wholly empty value (``--genomes ''``) is rejected
    rather than silently disabling the filter, which used to turn an unset
    shell variable into "dump everything".

    Args:
        ctx: Click context (unused, part of the callback signature).
        param: The Click ``Parameter`` being converted.
        value: Tuple of raw strings from a ``multiple=True`` option.

    Returns:
        A flat list of stripped, non-empty tokens.  ``[]`` when the option was
        not given at all.

    Raises:
        click.BadParameter: If the option was given but yields no tokens.
    """
    if not value:
        return []
    result: list[str] = []
    for v in value:
        result.extend(v.split(","))
    parsed = [x.strip() for x in result if x.strip()]
    if not parsed:
        name = param.opts[0] if param.opts else param.name
        raise click.BadParameter(
            f"{name} was given but is empty. Omit the option to disable the filter."
        )
    return parsed


def check_output(path: str | Path | None, force: bool = False) -> None:
    """Refuse to overwrite an existing output file unless ``--force`` was given.

    Thin CLI-facing wrapper over :func:`metacooler.io._utils.check_output_path`
    that reports the failure as a clean ``Error:`` line instead of a traceback.

    Args:
        path: Output path about to be written.  ``None`` is a no-op (the
            command writes to stdout).
        force: Whether the user passed ``-f/--force``.

    Raises:
        click.ClickException: If *path* exists and *force* is False.
    """
    if path is None:
        return
    from metacooler.io._utils import check_output_path

    try:
        check_output_path(path, force)
    except OSError as exc:
        raise click.ClickException(str(exc)) from exc


def force_option(f):
    """Attach the standard ``-f/--force`` flag to a command."""
    return click.option(
        "-f",
        "--force",
        is_flag=True,
        help="Overwrite output file if it already exists.",
    )(f)


def verbose_option(f):
    """Attach the standard ``-v/--verbose`` flag to a command."""
    return click.option(
        "-v",
        "--verbose",
        is_flag=True,
        help="Increase verbosity (debug logging).",
    )(f)


def quiet_option(f):
    """Attach the standard ``-q/--quiet`` flag to a command."""
    return click.option(
        "-q",
        "--quiet",
        is_flag=True,
        help="Suppress informational logging (warnings and errors still shown).",
    )(f)
