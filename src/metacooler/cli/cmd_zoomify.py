"""metacooler zoomify command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="zoomify", help="Create a multi-resolution mmetacool file.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
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
    help="Number of successive bin-size doublings (used when --resolutions is not set).",
)
@click.option(
    "-r",
    "--resolutions",
    type=str,
    default=None,
    help=(
        "Resolution ladder. Accepts an explicit list "
        "('--resolutions 1000,2000,5000'), or a single generator token: "
        "'1000B' (binary x2 progression from 1000). Every requested resolution must be a "
        "multiple of the base bin size. Takes precedence over --levels."
    ),
)
@click.option(
    "--temp-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help=(
        "Directory for the intermediate per-level files (default: the system "
        "temporary directory). Mirrors cooler's own [cyan]zoomify --temp-dir[/cyan]; "
        "point it at fast local scratch with room for one level."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_zoomify(
    metacool_path: str,
    output: Path,
    levels: int,
    resolutions: str | None,
    temp_dir: Path | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Create a multi-resolution mmetacool file.

    \b
    Input must be a fixed-bin metacooler. By default each level doubles the bin
    size; pass --resolutions for an explicit cooler-style resolution list.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.ops.zoomify import zoomify

    setup_logging(verbose, command_name="metacooler zoomify", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    zoomify(
        uri=metacool_path,
        output=output,
        levels=levels,
        resolutions=resolutions,
        temp_dir=temp_dir,
        force=force,
    )
    click.echo(f"Zoomified {metacool_path} → {output}")
