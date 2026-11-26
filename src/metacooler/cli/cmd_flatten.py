"""metacooler flatten command."""

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


@cli.command(name="flatten", help="Flatten a metacool file into a standard cool file.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .cool file path.",
)
@force_option
@verbose_option
@quiet_option
def cmd_flatten(
    metacool_path: Path,
    output: Path,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Flatten a metacool file by removing the genome hierarchy.

    \b
    All contigs become "chromosomes" in the output cool file.
    The resulting file is a valid cooler file compatible with
    cooltools, HiGlass, etc.
    """
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.io.flatten import flatten

    setup_logging(verbose, command_name="metacooler flatten", quiet=quiet)

    check_output(output, force, inputs=[metacool_path])

    with run_log(output, "Cool file written"):
        flatten(uri=metacool_path, output=output)
