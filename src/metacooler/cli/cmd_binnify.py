"""metacooler binnify command."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from typing import Literal


@cli.command(name="binnify", help="Re-bin a metacool file to a new bin size.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
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
    help=(
        "Tiling strategy for re-binning. "
        "[cyan]contig[/cyan]: tile each contig independently (default). "
        "[cyan]genome[/cyan]: concatenate all contigs of each genome before tiling."
    ),
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of source pixels read and aggregated per chunk; peak memory "
    "is ~40 bytes per pixel of one chunk. At --level contig whole output rows "
    "are processed per chunk and written in order; at --level genome each chunk "
    "is written as a temporary cooler at a fixed per-call cost that scales with "
    "contig count, so smaller chunks cost disproportionately more there.",
)
@force_option
@verbose_option
@quiet_option
def cmd_binnify(
    metacool_path: str,
    output: Path,
    binsize: int,
    level: Literal["contig", "genome"],
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Re-bin a metacool file to a larger bin size.

    \b
    The new bin size must be an exact multiple of the current bin size.
    Uses midpoint assignment to map old bins to new bins.
    """
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.ops.binnify import binnify

    setup_logging(verbose, command_name="metacooler binnify", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    with run_log(output, "Map written"):
        binnify(
            uri=metacool_path,
            output=output,
            binsize=binsize,
            level=level,
            chunksize=chunksize,
        )
