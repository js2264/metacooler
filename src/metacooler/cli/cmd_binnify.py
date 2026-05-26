"""metacooler binnify command."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli

if TYPE_CHECKING:
    from typing import Literal


@cli.command(name="binnify", help="Re-bin a metacool file to a new bin size.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .metacool file path.",
)
@click.option(
    "--binsize",
    type=int,
    required=True,
    help="New bin size in base pairs.",
)
@click.option(
    "--level",
    type=click.Choice(["contig", "genome"]),
    default="contig",
    show_default=True,
    help="Tiling level for re-binning.",
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
def cmd_binnify(
    metacool_path: str,
    output: Path,
    binsize: int,
    level: Literal["contig", "genome"],
    force: bool,
    verbose: bool,
) -> None:
    """Re-bin a metacool file to a larger bin size.

    \b
    The new bin size must be larger than the current bin size.
    Uses midpoint assignment to map old bins to new bins.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.ops.binnify import binnify

    setup_logging(verbose, command_name="metacooler binnify")

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    binnify(uri=metacool_path, output=output, binsize=binsize, level=level)
