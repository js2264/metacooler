"""metacooler create command group with shared option decorators."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import cli

_CONTIGS_FORMAT_HELP = """
    Format of the contig assignment table. Four formats are supported: [dim]\\[default: contigs][/dim]\n
    - [yellow]contigs[/yellow]: TSV [white][bold]with a header row[/white][/bold]. Column 1 must be named `contig_name` and column 2 `genome_name`. Additional columns are stored verbatim.\n
    - [yellow]contig_genome[/yellow]: Two-column TSV (`contig<TAB>genome`) [white][bold]with no header[/white][/bold]. First value is contig name, second is genome assignment.\n
    - [yellow]metator[/yellow]: MetaTOR output table (ie `contig_data_final.txt`). Expects a header with columns `Name` (contig) and `Final_bin` (genome). Other columns ignored.\n
    - [yellow]I,J(,noheader)[/yellow]: Column-index spec (1-based). Use column `I` for contigs and `J` for genomes. Header expected by default. Append `,noheader` (e.g., `1,5,noheader`) to treat the first row as data.
    """


def _reconstruct_command() -> str:
    """Reconstruct the exact command line invocation from the current click context."""
    ctx = click.get_current_context()
    if ctx is None:
        raise RuntimeError(
            "_reconstruct_command() called outside of a click command context."
        )
    parts = [ctx.command_path]
    for param in ctx.command.params:
        if isinstance(param, click.Option):
            value = ctx.params.get(param.name)
            if value is not None and value != param.default:
                opt_str = param.opts[0] if param.opts else f"--{param.name}"
                if param.is_flag:
                    if value:
                        parts.append(opt_str)
                elif param.multiple:
                    for v in value:
                        parts.extend([opt_str, str(v)])
                else:
                    parts.extend([opt_str, str(value)])
        elif isinstance(param, click.Argument):
            value = ctx.params.get(param.name)
            if value is not None:
                if isinstance(value, tuple):
                    for v in value:
                        parts.append(str(v))
                else:
                    parts.append(str(value))
    return " ".join(parts)


def shared_options(f):
    """Apply options shared by all ``create`` subcommands.

    Adds: ``-o/--output``, ``--contigs``, ``--contigs-format``, ``--genomes``,
    ``--metadata``, ``-f/--force``, ``-v/--verbose``.
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
        click.option(
            "--metadata",
            type=click.Path(exists=True, path_type=Path),
            help="Path to a JSON or YAML file containing metadata to ingest into the file.",
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Overwrite output file if it already exists.",
        ),
        click.option(
            "-v",
            "--verbose",
            is_flag=True,
            help="Increase verbosity (debug logging).",
        ),
    ]
    for dec in reversed(decorators):
        f = dec(f)
    return f


@cli.group(name="create", help="Create a metacool file from various sources.")
def create() -> None:
    """Create metacool files from pairs or cool files."""
    pass
