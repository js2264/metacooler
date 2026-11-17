"""metacooler collect command — one metacool file from many single-genome cools."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    contig_metadata_format_option,
    filter_option,
    force_option,
    genomes_format_option,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli._provenance import load_metadata_file, reconstruct_command
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.qc import DEFAULT_DISTANCE_THRESHOLD
from metacooler.schema import DEFAULT_CHUNKSIZE

# Spelled out rather than imported from metacooler.io.collect (which a test pins
# it to): that import would load cooler on every CLI start-up.
_ORDERINGS = ["cool", "genome-then-length", "genome-then-alpha"]


@cli.command(
    name="collect",
    help="Collect single-genome cool files into one metacool file, one genome each.",
)
@click.argument(
    "cool_paths",
    nargs=-1,
    required=True,
    type=MetaCoolURI(),
    metavar="COOL_PATH...",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .metacool file path.",
)
@click.option(
    "--names",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    help=(
        "Genome name for each input, in order (comma-separated or repeated). "
        "Defaults to each input's file stem "
        "([cyan]E_coli.cool[/cyan] → [cyan]E_coli[/cyan])."
    ),
)
@click.option(
    "--genomes",
    "genomes_table",
    type=click.Path(exists=True, path_type=Path),
    help=(
        "Optional genomes metadata table. Its genome-name column (see "
        "[cyan]--genomes-format[/cyan]) is matched against the genome names; "
        "every other column is stored in the /genomes table."
    ),
)
@genomes_format_option
@click.option(
    "--contigs",
    "contigs_table",
    type=click.Path(exists=True, path_type=Path),
    help=(
        "Optional per-contig metadata table. Its contig-name column (see "
        "[cyan]--contigs-format[/cyan]) names a contig as in its input cool, or "
        "as collected ([cyan]<genome>__<contig>[/cyan]). An optional "
        "'genome_name' column says which genome a row describes. Every other "
        "column is stored in the /contigs table."
    ),
)
@contig_metadata_format_option
@filter_option
@click.option(
    "--metadata",
    type=click.Path(exists=True, path_type=Path),
    help="Path to a JSON or YAML file containing metadata to ingest into the file.",
)
@click.option(
    "--contig-order",
    type=click.Choice(_ORDERINGS),
    default="cool",
    show_default=True,
    help=(
        """
        Order of the contigs within each genome; genomes are always laid out
        alphabetically.\n
        - [yellow]cool[/yellow]: keep each input's own contig order. Needs no
        pixel sorting.\n
        - [yellow]genome-then-length[/yellow]: longest contig first.\n
        - [yellow]genome-then-alpha[/yellow]: alphabetical.
        """
    ),
)
@click.option(
    "--assembly-name",
    type=str,
    help="Name for this collection, stored in the genome-assembly attribute.",
)
@click.option(
    "--qc/--no-qc",
    default=True,
    show_default=True,
    help=(
        "Tally contact-quality statistics while the pixels are copied and store "
        "them in [cyan]metadata.qc[/cyan] and the [cyan]qc_*[/cyan] columns of "
        "[cyan]/genomes[/cyan] and [cyan]/contigs[/cyan]. Inter-genome contacts "
        "are 0 by construction."
    ),
)
@click.option(
    "--qc-threshold",
    type=str,
    default=str(DEFAULT_DISTANCE_THRESHOLD),
    show_default=True,
    help=(
        "Distance splitting short- from long-range intra-contig signal. "
        "Accepts [cyan]bp/kb/Mb/Gb[/cyan] suffixes."
    ),
)
@click.option(
    "--sort-strategy",
    "sort_strategy",
    type=click.Choice(["auto", "memory", "external"]),
    default="auto",
    show_default=True,
    help=(
        "How the pixels of an input are re-sorted when --contig-order reorders "
        "its contigs (never needed under 'cool'). "
        "'auto': in memory, unless one such input exceeds 100 M pixels. "
        "'memory': in RAM, one input at a time. "
        "'external': cooler's external sort."
    ),
)
@click.option(
    "--chunksize",
    type=click.IntRange(min=1),
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help=(
        "Pixels read per chunk (peak-memory control), and the smallest chunk "
        "handed to cooler: under 'external' sorting each chunk is written as a "
        "temporary cooler at a fixed per-call cost, so small chunks cost "
        "disproportionately more."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_collect(
    cool_paths: tuple[str, ...],
    output: Path,
    names: list[str],
    genomes_table: Path | None,
    genomes_format: str,
    contigs_table: Path | None,
    contigs_format: str,
    filters: list[str],
    metadata: Path | None,
    contig_order: str,
    assembly_name: str | None,
    qc: bool,
    qc_threshold: str,
    sort_strategy: str,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Collect single-genome cool files into one metacool file.

    Every input becomes one genome — a ``.cool`` file, or one resolution of an
    ``.mcool`` (``file.mcool::resolutions/1000``) — and every contig is renamed
    ``<genome>__<contig>``.  Contacts between two inputs were never observed,
    so the inter-genome blocks of the output are empty.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.filters import parse_size
    from metacooler.io._utils import check_output_not_an_input
    from metacooler.io.collect import collect

    setup_logging(verbose, command_name="metacooler collect", quiet=quiet)

    # Before check_output: under --force that deletes the output up front,
    # which must never be one of the inputs.
    try:
        check_output_not_an_input(output, cool_paths)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="'-o' / '--output'") from exc
    check_output(output, force)
    try:
        qc_distance_threshold = int(parse_size(qc_threshold))
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--qc-threshold") from exc

    # ``commands`` is provenance and must survive the user's metadata file;
    # load_metadata_file() rejects a file that tries to set it.
    meta: dict = {"commands": [reconstruct_command()]}
    meta.update(load_metadata_file(metadata))

    collect(
        inputs=cool_paths,
        output=output,
        names=names or None,
        genomes_table=genomes_table,
        genomes_table_format=genomes_format,
        contigs_table=contigs_table,
        contigs_table_format=contigs_format,
        contig_ordering=cast(
            Literal["cool", "genome-then-length", "genome-then-alpha"], contig_order
        ),
        assembly_name=assembly_name,
        filters=list(filters) or None,
        metadata=meta,
        sort_in_memory={"auto": None, "memory": True, "external": False}[sort_strategy],
        chunksize=chunksize,
        qc=qc,
        qc_distance_threshold=qc_distance_threshold,
    )
