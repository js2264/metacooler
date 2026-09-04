"""Shared Click option callbacks and output guards for the metacooler CLI.

Single home for the idioms that were previously copy-pasted across
``cmd_*.py`` modules: comma-separated multi-value parsing and the
"refuse to clobber an existing output" check.
"""

from __future__ import annotations

from pathlib import Path

import rich_click as click

# Help text for the ``--contigs-format`` option, shared by every command that
# reads a contig-to-genome assignment table (``create pairs``, ``create cool``,
# ``reassign``).  One copy, so the four formats cannot drift apart between the
# command that writes a file and the command that re-labels it.
CONTIGS_FORMAT_HELP = """
    Format of the contig assignment table. Four formats are supported: [dim]\\[default: contigs][/dim]\n
    - [yellow]contigs[/yellow]: TSV [white][bold]with a header row[/white][/bold]. Column 1 must be named `contig_name` and column 2 `genome_name`. Additional columns are stored verbatim.\n
    - [yellow]contig_genome[/yellow]: Two-column TSV (`contig<TAB>genome`) [white][bold]with no header[/white][/bold]. First value is contig name, second is genome assignment.\n
    - [yellow]metator[/yellow]: MetaTOR output table (ie `contig_data_final.txt`). Expects a header with columns `Name` (contig) and `Final_bin` (genome). Other columns ignored.\n
    - [yellow]I,J(,noheader)[/yellow]: Column-index spec (1-based). Use column `I` for contigs and `J` for genomes. Header expected by default. Append `,noheader` (e.g., `1,5,noheader`) to treat the first row as data.
    """


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


def parse_filter_list(ctx, param, value):
    """Validate every ``--filter`` expression, passing the raw strings through.

    Deliberately *not* comma-split like :func:`parse_comma_list`: a filter
    value is free text (a regex, a genome name) and may legitimately contain a
    comma.  Repeat the option to AND several comparisons.

    Parsing here rather than in the io/ops layer means a typo is reported
    before the file is opened, and reported against the option that caused it.

    Args:
        ctx: Click context (unused, part of the callback signature).
        param: The Click ``Parameter`` being converted.
        value: Tuple of raw strings from a ``multiple=True`` option.

    Returns:
        The expressions as a list, unchanged.  ``[]`` when the option was not
        given.

    Raises:
        click.BadParameter: If any expression does not parse.
    """
    from metacooler.filters import FilterError, parse_filter

    exprs = [v for v in (value or ())]
    for expr in exprs:
        try:
            parse_filter(expr)
        except FilterError as exc:
            raise click.BadParameter(str(exc)) from exc
    return exprs


def filter_option(f):
    """Attach the standard repeatable ``--filter`` option to a command."""
    return click.option(
        "--filter",
        "filters",
        type=str,
        multiple=True,
        callback=parse_filter_list,
        metavar="EXPR",
        help=(
            "Keep only rows matching [cyan]<table>.<column> <op> <value>[/cyan], "
            "where [cyan]<table>[/cyan] is [yellow]genomes[/yellow] or "
            "[yellow]contigs[/yellow] — e.g. "
            "[cyan]--filter 'genomes.length > 1Mb'[/cyan]. Operators: "
            "[cyan]> >= < <= == != =~ !~[/cyan] ([cyan]=~[/cyan] is a regex "
            "match). Numbers accept [cyan]bp/kb/Mb/Gb[/cyan] suffixes. Repeat "
            "the option to require several conditions at once."
        ),
    )(f)


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
