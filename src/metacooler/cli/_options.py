"""Shared Click option callbacks and output guards for the metacooler CLI.

Single home for the idioms that were previously copy-pasted across
``cmd_*.py`` modules: comma-separated multi-value parsing and the
"refuse to clobber an existing output" check.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import rich_click as click

from metacooler.io.colspec import (
    ALIAS_FALLBACKS,
    ASSIGNMENT_KEYS,
    CONTIG_METADATA_FORMAT_DEFAULT,
    CONTIGS_FORMAT_ALIASES,
    CONTIGS_FORMAT_DEFAULT,
    GENOMES_FORMAT_DEFAULT,
    GENOMES_KEYS,
    parse_column_spec,
)

# Help texts for the column-spec options (grammar in metacooler.io.colspec).
# One copy each, so the commands that read the same kind of table cannot drift
# apart.  An *assignment* table (``create pairs``, ``create cool``,
# ``reassign``) names two columns; a per-contig *metadata* table (``collect``,
# ``annotate``) and a genomes table name one.
CONTIGS_FORMAT_HELP = (
    "Which columns of the [cyan]--contigs[/cyan] table hold the contig and the "
    "genome names: [cyan]CONTIG,GENOME[/cyan], each a 1-based position or a "
    "header name — e.g. [cyan]2,5[/cyan] or [cyan]Name,Final_cluster[/cyan] — plus "
    "[cyan],noheader[/cyan] when the table has no header line. "
    "[cyan]NA[/cyan] as the genome column ([cyan]1,NA,noheader[/cyan], "
    "[cyan]contig_name,NA[/cyan]) says the table assigns no genome: "
    "[cyan]create[/cyan] bins nothing yet (every contig in __unbinned__), "
    "[cyan]reassign[/cyan] unbins the contigs it lists. "
    f"[dim]\\[default: {CONTIGS_FORMAT_DEFAULT}][/dim]\n\n"
    "Aliases for the tables of common binners: "
    + ", ".join(
        f"[yellow]{name}[/yellow] ({spec}"
        + "".join(
            f"; {fallback} if the table has no {column}"
            for column, fallback in ALIAS_FALLBACKS.get(name, {}).items()
        )
        + ")"
        for name, spec in CONTIGS_FORMAT_ALIASES.items()
    )
    + "."
)
CONTIG_METADATA_FORMAT_HELP = (
    "Which column of the [cyan]--contigs[/cyan] table holds the contig names: a "
    "1-based position or a header name — e.g. [cyan]contigName[/cyan] or "
    "[cyan]2[/cyan] — plus [cyan],noheader[/cyan] when the table has no header "
    f"line. [dim]\\[default: {CONTIG_METADATA_FORMAT_DEFAULT}][/dim]"
)
GENOMES_FORMAT_HELP = (
    "Which column of the [cyan]--genomes[/cyan] table holds the genome names: a "
    "1-based position or a header name — e.g. [cyan]Name[/cyan] (CheckM2) or "
    "[cyan]user_genome[/cyan] (GTDB-Tk) — plus [cyan],noheader[/cyan] when the "
    f"table has no header line. [dim]\\[default: {GENOMES_FORMAT_DEFAULT}][/dim]"
)


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


def check_output(
    path: str | Path | None,
    force: bool = False,
    *,
    inputs: Iterable[str | Path | None] = (),
    defer_removal: bool = False,
) -> None:
    """Refuse to overwrite an existing output file unless ``--force`` was given.

    Thin CLI-facing wrapper over :func:`metacooler.io._utils.check_output_path`
    that reports the failure as a clean ``Error:`` line instead of a traceback.

    Args:
        path: Output path about to be written.  ``None`` is a no-op (the
            command writes to stdout).
        force: Whether the user passed ``-f/--force``.
        inputs: Every path or URI the command reads.  An output that is one of
            them is refused even under ``--force``.
        defer_removal: The command stages its output and renames it into
            place on success, so leave an existing file alone until then.

    Raises:
        click.ClickException: If *path* exists and *force* is False, or if
            *path* is one of *inputs*.
    """
    if path is None:
        return
    from metacooler.io._utils import check_output_path

    try:
        check_output_path(path, force, inputs=inputs, defer_removal=defer_removal)
    except (OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc


def option_given(name: str) -> bool:
    """Whether the option whose parameter is *name* was typed on the command line.

    For an option with a default, the value alone cannot tell "given" from
    "left at its default", and some options only mean something next to
    another one -- ``--contigs-format`` without ``--contigs``.
    """
    from click.core import ParameterSource

    ctx = click.get_current_context()
    assert ctx is not None  # only ever called from inside a running command
    return ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE


def refuse_contigs_format_without_contigs(contigs: object) -> None:
    """Refuse a ``--contigs-format`` given without the table it describes."""
    if contigs is None and option_given("contigs_format"):
        raise click.UsageError(
            "--contigs-format says how to read the --contigs table, which was "
            "not given. Pass the table, or leave both out to build a map with "
            "no binning yet."
        )


def force_option(f):
    """Attach the standard ``-f/--force`` flag to a command."""
    return click.option(
        "-f",
        "--force",
        is_flag=True,
        help="Overwrite output file if it already exists.",
    )(f)


def nproc_option(help_text: str):
    """Attach the standard ``--nproc`` option (default 1, at least 1)."""

    def decorate(f):
        return click.option(
            "--nproc",
            type=click.IntRange(min=1),
            default=1,
            show_default=True,
            help=help_text,
        )(f)

    return decorate


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


def _column_spec_callback(
    keys: tuple[str, ...],
    *,
    required: int | None = None,
    aliases: dict[str, str] | None = None,
):
    """Build a callback checking a column spec while the options are parsed.

    The spec is passed through unchanged -- the io layer parses it again when
    it reads the table -- but a typo is reported before any file is opened,
    and against the option that caused it.
    """

    def callback(ctx, param, value):
        if value is None:
            return value
        try:
            parse_column_spec(
                value,
                keys=keys,
                required=required,
                aliases=aliases,
                option=param.opts[0],
            )
        except ValueError as exc:
            raise click.BadParameter(str(exc)) from exc
        return value

    return callback


def contigs_format_option(f):
    """Attach ``--contigs-format`` for a contig-to-genome *assignment* table.

    Used by the commands whose ``--contigs`` table says which genome each
    contig belongs to (``create pairs``, ``create cool``, ``reassign``), so the
    spec names two columns.
    """
    return click.option(
        "--contigs-format",
        type=str,
        default=CONTIGS_FORMAT_DEFAULT,
        show_default=False,
        metavar="SPEC",
        callback=_column_spec_callback(ASSIGNMENT_KEYS, aliases=CONTIGS_FORMAT_ALIASES),
        help=CONTIGS_FORMAT_HELP,
    )(f)


def contig_metadata_format_option(f):
    """Attach ``--contigs-format`` for a per-contig *metadata* table.

    Used by the commands whose ``--contigs`` table only describes contigs whose
    genome is already known (``collect``, ``annotate``): the spec names the
    contig column.  A second, genome column is accepted -- it is only checked
    -- as are the assignment aliases, but neither is documented here.
    """
    return click.option(
        "--contigs-format",
        type=str,
        default=CONTIG_METADATA_FORMAT_DEFAULT,
        show_default=False,
        metavar="SPEC",
        callback=_column_spec_callback(
            ASSIGNMENT_KEYS, required=1, aliases=CONTIGS_FORMAT_ALIASES
        ),
        help=CONTIG_METADATA_FORMAT_HELP,
    )(f)


def genomes_format_option(f):
    """Attach ``--genomes-format``: which column of ``--genomes`` holds names."""
    return click.option(
        "--genomes-format",
        type=str,
        default=GENOMES_FORMAT_DEFAULT,
        show_default=False,
        metavar="SPEC",
        callback=_column_spec_callback(GENOMES_KEYS),
        help=GENOMES_FORMAT_HELP,
    )(f)
