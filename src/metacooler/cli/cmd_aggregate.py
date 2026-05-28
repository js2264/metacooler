"""metacooler aggregate command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="aggregate", help="Aggregate contacts to contig or genome level.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .cool file path.",
)
@click.option(
    "--level",
    type=click.Choice(["contig", "genome"]),
    default="contig",
    show_default=True,
    help=(
        "Aggregation level for the output chromosomes. "
        "[cyan]contig[/cyan]: one chromosome per contig (default). "
        "[cyan]genome[/cyan]: sum all contacts between genome pairs, "
        "producing one chromosome per genome."
    ),
)
@click.option(
    "--count-col",
    type=str,
    default="count",
    show_default=True,
    help="Pixel value column to aggregate (default: [cyan]count[/cyan]).",
)
@click.option(
    "--chunksize",
    type=int,
    default=1_000_000,
    show_default=True,
    help="Number of pixels read per chunk (controls peak memory usage).",
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
def cmd_aggregate(
    metacool_path: str,
    output: Path,
    level: str,
    count_col: str,
    chunksize: int,
    force: bool = False,
    verbose: bool = False,
) -> None:
    """Aggregate a metacool contact map to a standard ``.cool`` file.

    Writes a cooler-compatible output where each "chromosome" is either a
    contig (``--level contig``) or a genome (``--level genome``).  Pixels
    are processed in chunks so the operation scales to large files.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.ops.aggregate import aggregate

    setup_logging(verbose, command_name="metacooler aggregate")

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    aggregate(
        uri=metacool_path,
        output=output,
        level=cast(Literal["contig", "genome"], level),
        count_col=count_col,
        chunksize=chunksize,
    )
