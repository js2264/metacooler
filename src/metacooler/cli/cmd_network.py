"""metacooler network command."""

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


@cli.command(name="network", help="Export a network edge list from a metacool file.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--level",
    type=click.Choice(["bin", "contig", "genome"]),
    default="contig",
    show_default=True,
    help="Aggregation level for the network.",
)
@click.option(
    "--weight",
    type=str,
    default="count",
    show_default=True,
    help="Column name in the pixels table to use as edge weight.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file (default: stdout).",
)
@click.option(
    "--header",
    is_flag=True,
    help="Write a header line (node1, node2, value) to the output.",
)
@click.option(
    "--balance",
    type=str,
    default=None,
    metavar="NAME",
    help=(
        "Normalize the edges with the [cyan]/bins/NAME[/cyan] weights "
        "([cyan]metacooler normalize[/cyan], or another tool's): each pixel "
        "times the weights of its two bins, masked bins left out."
    ),
)
@click.option(
    "--drop-spurious",
    is_flag=True,
    help=(
        "Drop the contig pairs the [cyan]--balance[/cyan] normalization deems "
        "spurious: the cutoff [cyan]normalize normcc[/cyan] and "
        "[cyan]hiczin[/cyan] store, or another tool's. Contig or genome level."
    ),
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of pixels to process per chunk (memory control).",
)
@force_option
@verbose_option
@quiet_option
def cmd_network(
    metacool_path: str,
    level: str,
    weight: str,
    output: Path | None,
    header: bool,
    balance: str | None,
    drop_spurious: bool,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Export a network edge list (TSV, no header by default).

    \b
    Columns: node1, node2, value
    At --level bin: raw pixel pairs.
    At --level contig: aggregated by contig pairs.
    At --level genome: aggregated by genome pairs.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.network import network

    setup_logging(verbose, command_name="metacooler network", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    network(
        uri=metacool_path,
        level=cast(Literal["bin", "contig", "genome"], level),
        weight=weight,
        out=output,
        header=header,
        chunksize=chunksize,
        force=force,
        balance=balance,
        drop_spurious=drop_spurious,
    )
