"""metacooler zoomify command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="zoomify", help="Create a multi-resolution mmetacool file.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .mmetacool file path.",
)
@click.option(
    "--levels",
    type=int,
    default=5,
    show_default=True,
    help="Maximum number of resolution levels.",
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
def cmd_zoomify(
    metacool_path: str,
    output: Path,
    levels: int,
    force: bool,
    verbose: bool,
) -> None:
    """Create a multi-resolution mmetacool file.

    \b
    Input must be a fixed-bin metacooler. Each level doubles the bin size.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.ops.zoomify import zoomify

    setup_logging(verbose, command_name="metacooler zoomify")

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    zoomify(uri=metacool_path, output=output, levels=levels)
    click.echo(f"Zoomified {metacool_path} → {output}")
