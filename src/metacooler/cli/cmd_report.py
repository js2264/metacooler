"""metacooler report command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import force_option, quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(
    name="report",
    help="Write the report of a map: its numbers, its MultiQC row and its page.",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--outdir",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("."),
    show_default=True,
    help=(
        "Directory of the report, created if missing: "
        "[cyan]<sample>.report.txt[/cyan], [cyan]<sample>.report_mqc.yaml[/cyan], "
        "[cyan]<sample>.report.html[/cyan] and the log, [cyan]<sample>[/cyan] "
        "being the map's name without [cyan].metacool[/cyan], or its entry's."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_report(
    metacool_path: str,
    outdir: Path,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Write the report of a map.

    \b
    It reads the map, never writes it: what it holds (genomes, contigs, bins,
    contacts), where its contacts land (the QC block it stores, or one
    computed from its pixels), each genome, and its P(s). MultiQC finds the
    report_mqc.yaml files of several maps and merges their rows.
    """
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.io._utils import check_output_path
    from metacooler.io.report import report_paths, write_report
    from metacooler.uri import sample_name

    setup_logging(verbose, command_name="metacooler report", quiet=quiet)

    sample = sample_name(metacool_path)
    # An earlier report is refused before the log starts, as an existing
    # output of every other command is.
    for path in report_paths(outdir, sample):
        check_output_path(path, force, inputs=[metacool_path], defer_removal=True)
    with run_log(outdir, "Report written", directory=outdir, sample=sample, mkdir=True):
        write_report(metacool_path, outdir, force=force)
