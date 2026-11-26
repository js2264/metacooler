"""metacooler extract command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    filter_option,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(
    name="extract",
    help="Extract a standard cool file for specific genome(s)/contig(s).",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .cool file path.",
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help="Genome name(s) to extract (comma-separated or repeated).",
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help=(
        "Contig name(s) to extract, on top of those of --genomes "
        "(comma-separated or repeated)."
    ),
)
@filter_option
@click.option(
    "--split",
    is_flag=True,
    help="Write one .cool per genome (output is used as prefix).",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite existing output file(s), including per-genome files with --split.",
)
@verbose_option
@quiet_option
def cmd_extract(
    metacool_path: Path,
    genomes: list[str],
    contigs: list[str],
    filters: list[str],
    output: Path,
    split: bool,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Extract a standard cool file from a metacool file.

    \b
    Supports extracting one or more genomes and/or contigs.
    The output is a valid cooler file compatible with cooltools, HiGlass, etc.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.extract import extract

    setup_logging(verbose, command_name="metacooler extract", quiet=quiet)

    # With --split, *output* is a prefix rather than a file; the per-genome
    # destinations are guarded inside extract() instead.
    if not split:
        check_output(output, force, inputs=[metacool_path])

    genome_list = list(genomes)
    if not genome_list and not contigs and not filters:
        raise click.UsageError(
            "At least one of --genomes, --contigs or --filter is required."
        )
    if split and contigs:
        raise click.UsageError(
            "--contigs cannot be combined with --split (split is genome-only)."
        )

    extract(
        uri=metacool_path,
        genomes=genome_list or None,
        contigs=list(contigs) or None,
        filters=list(filters) or None,
        output=output,
        split=split,
        force=force,
    )
