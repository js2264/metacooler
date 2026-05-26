"""metacooler dump command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


def _parse_comma_list(ctx, param, value):
    """Parse a comma-separated or multi-value option into a flat list."""
    if not value:
        return []
    result = []
    for v in value:
        result.extend(v.split(","))
    return [x.strip() for x in result if x.strip()]


@cli.command(name="dump", help="Extract tables from a metacool file.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-t",
    "--table",
    type=click.Choice(
        ["chroms", "bins", "pixels", "contigs", "genomes", "indexes", "all"]
    ),
    default="pixels",
    show_default=True,
    help="Which table to dump. Use 'all' to dump all tables.",
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=_parse_comma_list,
    help="Filter by genome name(s) (comma-separated or repeated).",
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=_parse_comma_list,
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
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
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
    range1: str | None,
    range2: str | None,
    verbose: bool,
) -> None:
    """Dump a table from a metacool file to TSV."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.dump import dump_table

    setup_logging(verbose, command_name="metacooler dump")

    if table == "all" and not prefix:
        raise click.UsageError("--prefix is required when --table all is specified.")

    dump_table(
        uri=metacool_path,
        table=table,  # type: ignore
        genomes=genomes or None,
        contigs=contigs or None,
        range1=range1,
        range2=range2,
        header=header,
        out=output,
        prefix=prefix,
        pixels_format=pixels_format,  # type: ignore
        weight=weight,
    )
