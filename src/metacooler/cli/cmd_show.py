"""metacooler show command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    filter_option,
    force_option,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="show", help="Visualize contact maps.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output image file path.",
    required=True,
)
@click.option(
    "--level",
    type=click.Choice(["bin", "contig", "genome"]),
    default="genome",
    show_default=True,
    help=(
        "Resolution of the contact heatmap. "
        "[cyan]bin[/cyan]: one pixel per bin (finest). "
        "[cyan]contig[/cyan]: aggregate contacts per contig. "
        "[cyan]genome[/cyan]: aggregate contacts per genome (coarsest, default)."
    ),
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help=(
        "Genome name(s) to show (comma-separated or repeated). Every name goes "
        "on [b]both[/b] axes, so the figure holds all the intra- and "
        "inter-genome contacts between them, blocked in the order given. "
        "Use [cyan]--range/--range2[/cyan] for a different selection per axis."
    ),
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help=(
        "Contig name(s) to show (comma-separated or repeated), on both axes "
        "like [cyan]--genomes[/cyan]. Bare names only — for coordinates use "
        "[cyan]--range[/cyan]. Upgrades [cyan]--level genome[/cyan] to "
        "[cyan]contig[/cyan]."
    ),
)
@filter_option
@click.option(
    "--transform",
    type=click.Choice(["log1p", "log10", "none"]),
    default="log1p",
    show_default=True,
    help=(
        "Pixel value transform applied before colouring. "
        "[cyan]log1p[/cyan]: log(1 + count) (default, good for sparse data). "
        "[cyan]log10[/cyan]: log10(1 + count). "
        "[cyan]none[/cyan]: raw counts."
    ),
)
@click.option(
    "--cmap",
    type=str,
    default="afmhot_r",
    show_default=True,
    help="Matplotlib colormap name (e.g. [cyan]afmhot_r[/cyan], [cyan]Reds[/cyan]).",
)
@click.option(
    "--vmax",
    type=float,
    default=None,
    show_default=True,
    help="Upper colormap bound. If not specificed, falls back to the 99th percentile of non-zero values.",
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
    "--balance",
    type=str,
    default=None,
    metavar="NAME",
    help=(
        "Show the map normalized with the [cyan]/bins/NAME[/cyan] weights "
        "([cyan]metacooler normalize[/cyan]). Normalized values are small: "
        "[cyan]--transform none[/cyan] draws them on a linear scale."
    ),
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
@force_option
@verbose_option
@quiet_option
def cmd_show(
    metacool_path: str,
    level: str,
    genomes: list[str],
    contigs: list[str],
    filters: list[str],
    transform: str,
    cmap: str,
    vmax: float | None,
    separators: bool,
    output: Path | None,
    dpi: int,
    balance: str | None,
    range1: str | None,
    range2: str | None,
    force: bool = False,
    verbose: bool = False,
    quiet: bool = False,
) -> None:
    """Render a contact map heatmap from a metacool file."""
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.ops.show import show

    setup_logging(verbose, command_name="metacooler show", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    # ``show`` upgrades the level for --range/--contigs itself, so the Python
    # API and the CLI cannot drift apart.
    with run_log(output, "Image written"):
        show(
            uri=metacool_path,
            level=cast(Literal["bin", "contig", "genome"], level),
            genomes=genomes or None,
            contigs=contigs or None,
            filters=filters or None,
            range1=range1,
            range2=range2,
            transform=cast(Literal["log1p", "log10", "none"], transform),
            cmap=cmap,
            vmax=vmax,
            separators=separators,
            out=output,
            dpi=dpi,
            force=force,
            balance=balance,
        )
