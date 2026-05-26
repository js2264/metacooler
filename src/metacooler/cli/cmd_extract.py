"""metacooler extract command."""

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


@cli.command(
    name="extract",
    help="Extract a standard cool file for specific genome(s)/contig(s).",
)
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=_parse_comma_list,
    help="Genome name(s) to extract (comma-separated or repeated).",
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=_parse_comma_list,
    help="Contig name(s) to extract (comma-separated or repeated).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .cool file path.",
)
@click.option(
    "--split",
    is_flag=True,
    help="Write one .cool per genome (output is used as prefix).",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite output file if it already exists (ignored when --split is set).",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_extract(
    metacool_path: Path,
    genomes: list[str],
    contigs: list[str],
    output: Path,
    split: bool,
    force: bool,
    verbose: bool,
) -> None:
    """Extract a standard cool file from a metacool file.

    \b
    Supports extracting one or more genomes and/or contigs.
    The output is a valid cooler file compatible with cooltools, HiGlass, etc.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.extract import extract

    setup_logging(verbose, command_name="metacooler extract")

    if not split and output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if not split and force and output.exists():
        output.unlink()

    genome_list = list(genomes)
    if not genome_list and not contigs:
        raise click.UsageError("At least one of --genomes or --contigs is required.")

    extract(
        uri=metacool_path,
        genomes=genome_list or None,
        contigs=list(contigs) or None,
        output=output,
        split=split,
    )
