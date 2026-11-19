"""metacooler ps command — contact probability vs genomic distance, P(s)."""

from __future__ import annotations

import sys
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
from metacooler.ops.ps import DEFAULT_BINS_PER_DECADE
from metacooler.schema import DEFAULT_CHUNKSIZE


@cli.command(name="ps", help="Export/plot P(s), contact probability vs distance.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="TSV output path. Defaults to stdout.",
)
@click.option(
    "--plot",
    type=click.Path(path_type=Path),
    default=None,
    help=(
        "Also render the curve to this image path "
        "([cyan].png[/cyan], [cyan].pdf[/cyan], [cyan].svg[/cyan])."
    ),
)
@click.option(
    "--pairs",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Measure from this pairs file instead of the binned pixel table "
        "(plain text or any compression [cyan]create pairs[/cyan] reads, "
        "e.g. gzip or bgzip — read exactly as [cyan]create pairs[/cyan] "
        "reads it). Distances then come from read "
        "positions, so the curve starts at 1 bp instead of at one bin. The "
        "metacool is still needed: it supplies the contig\u2192genome "
        "assignment, the contig lengths, and the selection."
    ),
)
@click.option(
    "--by",
    type=click.Choice(["all", "genome"]),
    default="all",
    show_default=True,
    help=(
        "Curve granularity. [cyan]all[/cyan]: one curve pooling every selected "
        "contig. [cyan]genome[/cyan]: one curve per genome (MAG)."
    ),
)
@click.option(
    "--genomes",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help="Restrict to these genome(s) (comma-separated or repeated).",
)
@click.option(
    "--contigs",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help="Restrict to these contig(s) (comma-separated or repeated).",
)
@filter_option
@click.option(
    "--bins-per-decade",
    type=int,
    default=DEFAULT_BINS_PER_DECADE,
    show_default=True,
    help="Number of log-spaced distance bins per decade.",
)
@click.option(
    "--min-distance",
    type=str,
    default=None,
    help=(
        "Upper edge of the leading [cyan][0, min)[/cyan] distance bin, which "
        "holds the diagonal. Defaults to 1 bp with [cyan]--pairs[/cyan], and "
        "to the file's bin size without. Accepts [cyan]bp/kb/Mb/Gb[/cyan] "
        "suffixes."
    ),
)
@click.option(
    "--max-distance",
    type=str,
    default=None,
    help=(
        "Largest separation to report. Defaults to the longest separation the "
        "selected contigs can offer. Accepts [cyan]bp/kb/Mb/Gb[/cyan] suffixes."
    ),
)
@click.option(
    "--value",
    type=click.Choice(["prob", "avg_contacts"]),
    default="prob",
    show_default=True,
    help=(
        "Column plotted on the y axis. [cyan]prob[/cyan]: unit-area normalised "
        "P(s). [cyan]avg_contacts[/cyan]: mean contacts per opportunity."
    ),
)
@click.option(
    "--header/--no-header",
    default=True,
    show_default=True,
    help="Write a header row in the TSV.",
)
@click.option(
    "--dpi",
    type=int,
    default=300,
    show_default=True,
    help="DPI for the output image.",
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of pixels (or pairs, under --pairs) to process per chunk.",
)
@force_option
@verbose_option
@quiet_option
def cmd_ps(
    metacool_path: str,
    output: Path | None,
    plot: Path | None,
    pairs: Path | None,
    by: str,
    genomes: list[str],
    contigs: list[str],
    filters: list[str],
    bins_per_decade: int,
    min_distance: str | None,
    max_distance: str | None,
    value: str,
    header: bool,
    dpi: int,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Export (and optionally plot) P(s) for a metacool file.

    Only contacts whose two ends sit on the same contig take part — those are
    the ones a genomic distance is defined for.  Each distance bin reports the
    contacts observed, the number of opportunities that could have produced one
    (bin pairs, or position pairs under ``--pairs``), and their ratio, so the
    decay is not confounded by the shrinking number of long-range pairs a
    contig offers.

    Read from the binned map by default, which cannot resolve anything below
    one bin; pass ``--pairs`` to measure the same curve from read positions.
    ``create pairs`` already writes that finer table beside its output as
    ``<output>.ps.tsv``.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.filters import parse_distance
    from metacooler.ops.ps import ps

    setup_logging(verbose, command_name="metacooler ps", quiet=quiet)

    check_output(output, force, inputs=[metacool_path, pairs])
    check_output(plot, force, inputs=[metacool_path, pairs])

    def _distance(raw: str | None, hint: str) -> int | None:
        if raw is None:
            return None
        try:
            return parse_distance(raw)
        except ValueError as exc:
            raise click.BadParameter(str(exc), param_hint=hint) from exc

    table = ps(
        uri=metacool_path,
        out=output,
        plot=plot,
        pairs=pairs,
        by=cast(Literal["all", "genome"], by),
        genomes=genomes or None,
        contigs=contigs or None,
        filters=filters or None,
        bins_per_decade=bins_per_decade,
        min_distance=_distance(min_distance, "--min-distance"),
        max_distance=_distance(max_distance, "--max-distance"),
        value=value,
        header=header,
        dpi=dpi,
        chunksize=chunksize,
        force=force,
    )
    if output is None:
        table.to_csv(sys.stdout, sep="\t", index=False, header=header)
