"""metacooler dump command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli


def _given(ctx, name: str) -> bool:
    """Return True when *name* was supplied on the command line.

    Used to reject option/table combinations that the dump layer would
    otherwise accept and silently ignore.  Checking the *source* rather than
    the value is what distinguishes an explicit ``--weight count`` from the
    identical default.
    """
    from click.core import ParameterSource

    return ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE


@cli.command(name="dump", help="Extract tables from a metacool file.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-t",
    "--table",
    type=click.Choice(
        ["chroms", "bins", "pixels", "contigs", "genomes", "indexes", "all"]
    ),
    default="pixels",
    show_default=True,
    help=(
        "Which table to dump ('chroms', 'bins', 'contigs', 'genomes', 'pixels', "
        "'indexes'). Use 'all' to dump every table."
    ),
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help="Filter by genome name(s) (comma-separated or repeated).",
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help="Filter by contig name(s) (comma-separated or repeated).",
)
@click.option(
    "-H",
    "--header",
    is_flag=True,
    default=False,
    help="Include header row.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file (default: stdout).",
)
@click.option(
    "--prefix",
    type=str,
    help="File prefix (required when --table all).",
)
@click.option(
    "--pixels-format",
    type=click.Choice(["long", "wide"]),
    default="long",
    show_default=True,
    help="Pixel output format: long (default) or wide (2D matrix).",
)
@click.option(
    "--weight",
    type=str,
    default="count",
    show_default=True,
    help="Column name for wide format values.",
)
@click.option(
    "--join",
    is_flag=True,
    default=False,
    help=(
        "Annotate pixels with contig names and coordinates. Long format replaces "
        "bin1_id/bin2_id with contig1/start1/end1 and contig2/start2/end2; wide "
        "format labels the matrix rows/columns with contig:start-end. "
        "Applied only when --table pixels."
    ),
)
@click.option(
    "-r",
    "--range",
    "range1",
    type=str,
    default=None,
    help=(
        "Row-dimension range: genome name, contig name, or contig:start-end. "
        "Applied only when --table pixels."
    ),
)
@click.option(
    "-r2",
    "--range2",
    type=str,
    default=None,
    help=(
        "Column-dimension range (defaults to --range when omitted). "
        "Applied only when --table pixels."
    ),
)
@force_option
@verbose_option
@quiet_option
def dump(
    metacool_path: str,
    table: str,
    genomes: list[str],
    contigs: list[str],
    header: bool,
    output: Path | None,
    prefix: str | None,
    pixels_format: str,
    weight: str,
    join: bool,
    range1: str | None,
    range2: str | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Dump a table from a metacool file to TSV."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.dump import dump_table

    setup_logging(verbose, command_name="metacooler dump", quiet=quiet)

    if table == "all" and not prefix:
        raise click.UsageError("--prefix is required when --table all is specified.")

    # Reject option/table combinations that would otherwise be accepted and
    # silently ignored.  ``all`` is exempt from the pixel-option checks because
    # it does dump a pixels table (and forwards these options to it).
    ctx = click.get_current_context()
    if table not in ("pixels", "all"):
        for opt, name in (("--weight", "weight"), ("--join", "join")):
            if _given(ctx, name):
                raise click.UsageError(
                    f"{opt} only applies to the pixels table; "
                    f"it is meaningless with --table {table}."
                )
    if table == "indexes":
        for opt, values in (("--genomes", genomes), ("--contigs", contigs)):
            if values:
                raise click.UsageError(
                    f"{opt} cannot be used with --table indexes: the indexes "
                    "table has no genome or contig dimension."
                )
    if table == "all" and output is not None:
        raise click.UsageError(
            "-o/--output cannot be used with --table all, which writes one file "
            "per table. Use --prefix instead."
        )
    if table != "all" and prefix is not None:
        raise click.UsageError(
            "--prefix only applies to --table all; use -o/--output instead."
        )

    check_output(output, force)

    dump_table(
        uri=metacool_path,
        table=cast(
            Literal["bins", "chroms", "contigs", "genomes", "pixels", "indexes", "all"],
            table,
        ),
        genomes=genomes or None,
        contigs=contigs or None,
        range1=range1,
        range2=range2,
        header=header,
        out=output,
        prefix=prefix,
        pixels_format=cast(Literal["long", "wide"], pixels_format),
        weight=weight,
        join=join,
        force=force,
    )
