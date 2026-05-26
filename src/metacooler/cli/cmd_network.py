"""metacooler network command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="network", help="Export a network edge list from a metacool file.")
@click.argument("metacool_path", type=MetaCoolURI())
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
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_network(
    metacool_path: str,
    level: str,
    weight: str,
    output: Path | None,
    header: bool,
    verbose: bool,
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

    setup_logging(verbose, command_name="metacooler network")

    network(uri=metacool_path, level=level, weight=weight, out=output, header=header)  # type: ignore
