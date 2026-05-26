"""metacooler flatten command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="flatten", help="Flatten a metacool file into a standard cool file.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .cool file path.",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite output file if it already exists.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_flatten(
    metacool_path: Path,
    output: Path,
    force: bool,
    verbose: bool,
) -> None:
    """Flatten a metacool file by removing the genome hierarchy.

    \b
    All contigs become "chromosomes" in the output cool file.
    The resulting file is a valid cooler file compatible with
    cooltools, HiGlass, etc.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.flatten import flatten

    setup_logging(verbose, command_name="metacooler flatten")

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    flatten(uri=metacool_path, output=output)
