"""metacooler create command group with shared option decorators."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    filter_option,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import cli
from metacooler.qc import DEFAULT_DISTANCE_THRESHOLD

_CONTIGS_FORMAT_HELP = """
    Format of the contig assignment table. Four formats are supported: [dim]\\[default: contigs][/dim]\n
    - [yellow]contigs[/yellow]: TSV [white][bold]with a header row[/white][/bold]. Column 1 must be named `contig_name` and column 2 `genome_name`. Additional columns are stored verbatim.\n
    - [yellow]contig_genome[/yellow]: Two-column TSV (`contig<TAB>genome`) [white][bold]with no header[/white][/bold]. First value is contig name, second is genome assignment.\n
    - [yellow]metator[/yellow]: MetaTOR output table (ie `contig_data_final.txt`). Expects a header with columns `Name` (contig) and `Final_bin` (genome). Other columns ignored.\n
    - [yellow]I,J(,noheader)[/yellow]: Column-index spec (1-based). Use column `I` for contigs and `J` for genomes. Header expected by default. Append `,noheader` (e.g., `1,5,noheader`) to treat the first row as data.
    """


def shared_options(f):
    """Apply options shared by all ``create`` subcommands.

    Adds: ``-o/--output``, ``--contigs``, ``--contigs-format``, ``--genomes``,
    ``--filter``, ``--metadata``, ``--qc/--no-qc``, ``--qc-threshold``,
    ``-f/--force``, ``-v/--verbose``, ``-q/--quiet``.

    The filter/force/verbose/quiet options come from
    :mod:`metacooler.cli._options` so they stay identical to every other
    command in the CLI.
    """
    decorators = [
        click.option(
            "-o",
            "--output",
            type=click.Path(path_type=Path),
            required=True,
            help="Output .metacool file path.",
        ),
        click.option(
            "--contigs",
            type=click.Path(exists=True, path_type=Path),
            required=True,
            help="Contig-to-genome assignment table. The expected layout depends on --contigs-format.",
        ),
        click.option(
            "--contigs-format",
            type=str,
            default="contigs",
            show_default=False,
            help=_CONTIGS_FORMAT_HELP,
        ),
        click.option(
            "--genomes",
            "genomes_table",
            type=click.Path(exists=True, path_type=Path),
            help=(
                "Optional genomes metadata table with header. "
                "First column must be 'genome_name'; any further columns are stored in the /genomes table."
            ),
        ),
        filter_option,
        click.option(
            "--metadata",
            type=click.Path(exists=True, path_type=Path),
            help="Path to a JSON or YAML file containing metadata to ingest into the file.",
        ),
        click.option(
            "--qc/--no-qc",
            default=True,
            show_default=True,
            help=(
                "Tally contact-quality statistics during creation and store them "
                "in [cyan]metadata.qc[/cyan] and the [cyan]/genomes[/cyan] "
                "[cyan]qc_*[/cyan] columns. Read them back with "
                "[cyan]metacooler qc[/cyan]."
            ),
        ),
        click.option(
            "--qc-threshold",
            type=str,
            default=str(DEFAULT_DISTANCE_THRESHOLD),
            show_default=True,
            help=(
                "Distance splitting short- from long-range intra-contig signal. "
                "Accepts [cyan]bp/kb/Mb/Gb[/cyan] suffixes."
            ),
        ),
        force_option,
        verbose_option,
        quiet_option,
    ]
    for dec in reversed(decorators):
        f = dec(f)
    return f


@cli.group(name="create", help="Create a metacool file from various sources.")
def create() -> None:
    """Create metacool files from pairs or cool files."""
    pass
