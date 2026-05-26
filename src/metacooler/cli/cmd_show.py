"""metacooler show command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="show", help="Visualize contact maps.")
@click.argument("metacool_path", type=MetaCoolURI())
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
    type=click.Choice(["log1p", "none"]),
    default="log1p",
    show_default=True,
    help="Value transform.",
)
@click.option(
    "--cmap",
    type=str,
    default="YlOrRd",
    show_default=True,
    help="Matplotlib colormap.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Save to file instead of displaying.",
)
@click.option(
    "--dpi",
    type=int,
    default=150,
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
    metacool_path: str,
    level: str,
    genome: str | None,
    genomes: tuple[str, str] | None,
    transform: str,
    cmap: str,
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
        uri=metacool_path,
        level=level,  # type: ignore
        genome=genome,
        genomes=genomes,
        range1=range1,
        range2=range2,
        transform=transform,  # type: ignore
        cmap=cmap,
        out=output,
        dpi=dpi,
    )
