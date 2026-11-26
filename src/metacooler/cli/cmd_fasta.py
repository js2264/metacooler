"""metacooler fasta command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    filter_option,
    option_given,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(
    name="fasta",
    help="Write the contigs of each genome or cluster as FASTA files.",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--fasta",
    "assembly",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help=(
        "The assembly the map was made from: FASTA, plain or gzip-compressed, "
        "or a pipe ([cyan]<(zcat ...)[/cyan]). Every contig written must be in "
        "it, once, with the map's length. Nothing is written next to it."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help=(
        "Directory to write one [cyan]<name>.fa[/cyan] per genome or cluster "
        "into; with [cyan]--single-file[/cyan], the FASTA file to write. "
        "Created if missing."
    ),
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help=(
        "Genome name(s) to write (comma-separated or repeated); with "
        "[cyan]--round[/cyan], cluster name(s) of that round. Default: every "
        "genome but __unbinned__, or every cluster of the round."
    ),
)
@filter_option(tables=("genomes", "contigs", "clusters"))
@click.option(
    "--round",
    "round_name",
    type=str,
    default=None,
    metavar="NAME",
    help=(
        "Write the clusters of this clustering round of the map "
        "([cyan]dump -t clusters[/cyan] lists them) instead of its genomes."
    ),
)
@click.option(
    "--unbinned",
    is_flag=True,
    help=(
        "Also write the contigs in no genome (in no cluster of "
        "[cyan]--round[/cyan]) as [cyan]__unbinned__[/cyan], as they are, even "
        "with [cyan]--stitch[/cyan]."
    ),
)
@click.option(
    "--stitch",
    is_flag=True,
    help=(
        "Write each genome or cluster as one record named after it: its "
        "contigs longest first (ties in map order), joined by "
        "[cyan]--spacer[/cyan] N, wrapped at 80 columns."
    ),
)
@click.option(
    "--spacer",
    type=click.IntRange(min=0),
    default=100,
    show_default=True,
    help="Number of N between two contigs of a stitched record.",
)
@click.option(
    "--single-file",
    is_flag=True,
    help=(
        "Write every record to one FASTA file, [cyan]-o[/cyan]: a contig "
        "record's header becomes [cyan]<genome>|<header>[/cyan]; a stitched "
        "record is named after its genome or cluster."
    ),
)
@click.option(
    "--temp-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help=(
        "Directory for the plain copy of a gzip-compressed or piped assembly "
        "(default: the system temporary directory), created if missing."
    ),
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite the output file(s) that already exist.",
)
@verbose_option
@quiet_option
def cmd_fasta(
    metacool_path: str,
    assembly: Path,
    output: Path,
    genomes: list[str],
    filters: list[str],
    round_name: str | None,
    unbinned: bool,
    stitch: bool,
    spacer: int,
    single_file: bool,
    temp_dir: Path | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Write the contig sequences of each genome, or of each cluster of a
    clustering round, from the assembly the map was made from.

    \b
    By default, one file per genome, its contigs' records as they are in the
    assembly (header and line breaks included), in map order.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.fasta import write_fasta

    setup_logging(verbose, command_name="metacooler fasta", quiet=quiet)
    if option_given("spacer") and not stitch:
        raise click.UsageError("--spacer separates stitched contigs: add --stitch.")
    if single_file:
        check_output(
            output, force, inputs=[metacool_path, assembly], defer_removal=True
        )
    write_fasta(
        metacool_path,
        assembly,
        output,
        genomes=genomes or None,
        filters=filters or None,
        round=round_name,
        unbinned=unbinned,
        stitch=stitch,
        spacer=spacer if stitch else None,
        single_file=single_file,
        force=force,
        temp_dir=temp_dir,
    )
