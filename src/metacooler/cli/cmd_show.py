"""metacooler show command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="show", help="Visualize contact maps.")
@click.option(
    "-f",
    "--file",
    type=MetaCoolURI(),
    help="Path to a .metacool file or to a mmetacool entry",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Save to file instead of displaying.",
    required=True,
)
@click.option(
    "--level",
    type=click.Choice(["bin", "contig", "genome"]),
    default="genome",
    show_default=True,
    help="Heatmap resolution.",
)
@click.option(
    "--genome",
    type=str,
    help="Show only this genome (intra-genome view).",
)
@click.option(
    "--genomes",
    nargs=2,
    type=str,
    help="Show a pair of genomes (inter-genome view).",
)
@click.option(
    "--transform",
    type=click.Choice(["log1p", "log10", "none"]),
    default="log1p",
    show_default=True,
    help="Value transform.",
)
@click.option(
    "--cmap",
    type=str,
    default="afmhot_r",
    show_default=True,
    help="Matplotlib colormap.",
)
@click.option(
    "--vmax",
    type=float,
    default=None,
    show_default=True,
    help="Upper colormap bound. Defaults to the 99th percentile of non-zero values.",
)
@click.option(
    "--separators/--no-separators",
    default=True,
    show_default=True,
    help="Draw lines at genome/contig boundaries.",
)
@click.option(
    "--dpi",
    type=int,
    default=500,
    show_default=True,
    help="DPI for output image.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
@click.option(
    "-r",
    "--range",
    "range1",
    type=str,
    default=None,
    help=(
        "Row-dimension range: genome name, contig name, or contig:start-end. "
        "Implies --level bin."
    ),
)
@click.option(
    "-r2",
    "--range2",
    type=str,
    default=None,
    help="Column-dimension range (defaults to --range when omitted). Implies --level bin.",
)
def cmd_show(
    file: str,
    level: str,
    genome: str | None,
    genomes: tuple[str, str] | None,
    transform: str,
    cmap: str,
    vmax: float | None,
    separators: bool,
    output: Path | None,
    dpi: int,
    range1: str | None,
    range2: str | None,
    verbose: bool = False,
) -> None:
    """Render a contact map heatmap from a metacool file."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.ops.show import show

    setup_logging(verbose, command_name="metacooler show")

    # --range implies bin-level rendering
    if range1 is not None and level == "genome":
        level = "bin"

    show(
        uri=file,
        level=cast(Literal["bin", "contig", "genome"], level),
        genome=genome,
        genomes=genomes,
        range1=range1,
        range2=range2,
        transform=cast(Literal["log1p", "log10", "none"], transform),
        cmap=cmap,
        vmax=vmax,
        separators=separators,
        out=output,
        dpi=dpi,
    )
