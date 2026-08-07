"""metacooler qc command — contact-quality statistics of a metacool file."""

from __future__ import annotations

import json
from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.qc import DEFAULT_DISTANCE_THRESHOLD
from metacooler.schema import DEFAULT_CHUNKSIZE

# Rows of the summary table, as (label, key) pairs, in reading order.
_TABLE_ROWS = (
    ("Source", "source"),
    ("Distance threshold (bp)", "distance_threshold"),
    ("Total contacts", "n_contacts"),
    ("Intra-contig, short range", "intra_contig_short"),
    ("Intra-contig, long range", "intra_contig_long"),
    ("Inter-contig, same genome", "inter_contig_intra_genome"),
    ("Inter-genome", "inter_genome"),
    ("Trans-contig ratio (3D ratio)", "trans_contig_ratio"),
    ("Intra-genome fraction", "intra_genome_fraction"),
    ("Inter-genome fraction", "inter_genome_fraction"),
    ("Pairs read", "n_pairs_read"),
    ("Pairs discarded", "n_pairs_discarded"),
)


@cli.command(name="qc", help="Report contact-quality statistics.")
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the JSON summary here instead of to stdout.",
)
@click.option(
    "--table",
    "table_out",
    is_flag=True,
    help="Print a rich table instead of JSON.",
)
@click.option(
    "--per-genome",
    "per_genome_out",
    type=click.Path(path_type=Path),
    default=None,
    help="Write the per-genome breakdown as a TSV to this path.",
)
@click.option(
    "--threshold",
    type=str,
    default=str(DEFAULT_DISTANCE_THRESHOLD),
    show_default=True,
    help=(
        "Distance splitting short- from long-range intra-contig signal. "
        "Accepts [cyan]bp/kb/Mb/Gb[/cyan] suffixes (e.g. [cyan]10kb[/cyan])."
    ),
)
@click.option(
    "--store",
    is_flag=True,
    help=(
        "Write the result back into the file: [cyan]metadata.qc[/cyan] plus the "
        "per-genome [cyan]qc_*[/cyan] columns of [cyan]/genomes[/cyan]."
    ),
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of pixels to process per chunk.",
)
@force_option
@verbose_option
@quiet_option
def cmd_qc(
    metacool_path: str,
    output: Path | None,
    table_out: bool,
    per_genome_out: Path | None,
    threshold: str,
    store: bool,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Report where a library's contacts land: the 3D ratio and its breakdown.

    Statistics are recomputed from the pixel table, so they describe the file
    as it is now.  A file created with ``create pairs`` already carries a more
    precise block (measured on read positions, before binning) under
    ``metadata.qc`` — read it with ``metacooler metadata get -k qc``.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.filters import parse_size
    from metacooler.ops.qc import qc as run_qc

    setup_logging(verbose, command_name="metacooler qc", quiet=quiet)

    check_output(output, force)
    check_output(per_genome_out, force)
    try:
        distance_threshold = int(parse_size(threshold))
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--threshold") from exc

    result = run_qc(
        metacool_path,
        distance_threshold=distance_threshold,
        chunksize=chunksize,
        store=store,
    )

    if per_genome_out is not None:
        result.per_genome.to_csv(per_genome_out, sep="\t", index=False)

    if table_out:
        from rich.console import Console
        from rich.table import Table

        table = Table(
            title=f"Contact QC: {Path(metacool_path.split('::')[0]).name}",
        )
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="magenta", justify="right")
        for label, key in _TABLE_ROWS:
            if key not in result.summary:
                continue
            value = result.summary[key]
            table.add_row(label, f"{value:,}" if isinstance(value, int) else str(value))
        Console().print(table)
        return

    text = json.dumps(result.summary, indent=2, default=str)
    if output is not None:
        output.write_text(text + "\n")
    else:
        click.echo(text)
