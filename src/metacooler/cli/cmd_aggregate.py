"""metacooler aggregate command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.schema import DEFAULT_CHUNKSIZE


@cli.command(name="aggregate", help="Aggregate contacts to contig or genome level.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
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
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of pixels read and aggregated per chunk. Each chunk is written "
    "as a temporary cooler at a fixed per-call cost that scales with contig "
    "count, so smaller chunks cost disproportionately more; the trade-off is "
    "~24 bytes of peak memory per buffered pixel.",
)
@force_option
@verbose_option
@quiet_option
def cmd_aggregate(
    metacool_path: str,
    output: Path,
    level: str,
    count_col: str,
    chunksize: int,
    force: bool = False,
    verbose: bool = False,
    quiet: bool = False,
) -> None:
    """Aggregate a metacool contact map to a standard ``.cool`` file.

    Writes a cooler-compatible output where each "chromosome" is either a
    contig (``--level contig``) or a genome (``--level genome``).  Pixels
    are processed in chunks so the operation scales to large files.
    """
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.ops.aggregate import aggregate

    setup_logging(verbose, command_name="metacooler aggregate", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    with run_log(output, "Cool file written"):
        aggregate(
            uri=metacool_path,
            output=output,
            level=cast(Literal["contig", "genome"], level),
            count_col=count_col,
            chunksize=chunksize,
        )
