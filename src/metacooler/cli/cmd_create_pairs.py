"""metacooler create pairs command — create a metacool file from pairs files."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    nproc_option,
    refuse_contigs_format_without_contigs,
)
from metacooler.cli._provenance import load_metadata_file, reconstruct_command
from metacooler.cli.cmd_create import create, shared_options
from metacooler.ps import DEFAULT_BINS_PER_DECADE
from metacooler.schema import DEFAULT_CHUNKSIZE


@create.command(
    name="pairs", help="Create a metacool file from pairs and contig assignments."
)
@click.option(
    "--pairs",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    multiple=True,
    help=(
        "Pairs file: plain text, or compressed with gzip, bgzip, bzip2, zstd or "
        "lz4, read from its content. Repeat the option to read several "
        "files — one per Hi-C library of the same assembly — in one pass, as if "
        "concatenated: their contacts are summed, and the QC block, the P(s) "
        "table and the dropped-pairs report cover them all. Not comma-split: a "
        "path may contain a comma."
    ),
)
@click.option(
    "--bins",
    type=str,
    required=True,
    help="BED file of bins (without column names), or an integer bin size in bp (e.g. 1000). "
    "If an integer, contig lengths must be available from the pairs header or a --fasta should be provided.",
)
@shared_options
@click.option(
    "--contig-order",
    type=click.Choice(["genome-then-length", "genome-then-alpha", "fasta", "alpha"]),
    default="genome-then-length",
    show_default=True,
    help=(
        """
        Strategy for ordering contigs in the output metacool file.\n
        - [yellow]genome-then-length[/yellow]: group by genome, longest-first within each genome.\n
        - [yellow]genome-then-alpha[/yellow]: group by genome, alphabetical within each genome.\n
        - [yellow]fasta[/yellow]: preserve order from the input FASTA (or pairs header).\n
        - [yellow]alpha[/yellow]: sort all contigs alphabetically.
        """
    ),
)
@click.option(
    "--fasta",
    type=click.Path(exists=True, path_type=Path),
    help="Assembly FASTA. Used to extract contig lengths when --bins is an integer "
    "and contig sizes are absent from the pairs header.",
)
@click.option(
    "--assembly-name",
    type=str,
    help="Name for this assembly / metagenome.",
)
@click.option(
    "--drop-unassigned",
    is_flag=True,
    help="Drop contigs not assigned to any genome.",
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help="Number of pairs to process per chunk. This is chunk granularity, not "
    "read memory (the pyarrow read block is capped separately). Under "
    "'external' sorting each chunk is written as a temporary cooler at a fixed "
    "per-call cost that scales with contig count, so smaller chunks cost "
    "disproportionately more.",
)
@nproc_option(
    "Threads binning pairs while the file is read and decompressed (which "
    "stays on one thread). The output is identical whatever the value; each "
    "thread holds one more batch of pairs in memory."
)
@click.option(
    "--ps/--no-ps",
    default=True,
    show_default=True,
    help=(
        "Accumulate a per-MAG [b]P(s)[/b] curve from the same pass over the "
        "pairs and write it beside the output as [cyan]<output>.ps.tsv[/cyan]. "
        "Measured on read positions, so it resolves distances the binned map "
        "cannot. No plot is drawn — pass the table to "
        "[cyan]metacooler ps[/cyan] or your own plotting."
    ),
)
@click.option(
    "--ps-bins-per-decade",
    type=int,
    default=DEFAULT_BINS_PER_DECADE,
    show_default=True,
    help="Number of log-spaced distance bins per decade in the P(s) table.",
)
@click.option(
    "--sort-strategy",
    "sort_strategy",
    type=click.Choice(["auto", "memory", "external"]),
    default="auto",
    show_default=True,
    help=(
        "Pixel sorting strategy. "
        "'auto': use 'memory' unless the binned pixels exceed 100 M rows. "
        "'memory': collect all pixels in RAM and sort before writing (fast for small files). "
        "'external': stream unsorted pixels to cooler's external sort (scalable for large files)."
    ),
)
def cmd_create_pairs(
    output: Path,
    bins: str,
    pairs: tuple[Path, ...],
    contigs: Path | None,
    contigs_format: str,
    genomes_table: Path | None,
    genomes_format: str,
    fasta: Path | None,
    assembly_name: str | None,
    contig_order: str,
    drop_unassigned: bool,
    filters: list[str],
    chunksize: int,
    sort_strategy: str,
    force: bool,
    verbose: bool,
    quiet: bool,
    metadata: Path | None,
    qc: bool,
    qc_threshold: str,
    ps: bool,
    ps_bins_per_decade: int,
    nproc: int,
) -> None:
    """
    Create a metacool file from pairs file(s), bins, and contig assignments.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.filters import parse_distance
    from metacooler.io.create import from_pairs

    setup_logging(verbose, command_name="metacooler create pairs", quiet=quiet)
    refuse_contigs_format_without_contigs(contigs)

    check_output(
        output,
        force,
        inputs=[*pairs, bins, contigs, genomes_table, fasta, metadata],
        defer_removal=True,
    )
    try:
        qc_distance_threshold = parse_distance(qc_threshold)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--qc-threshold") from exc

    # Map --sort-strategy to sort_in_memory bool
    sort_in_memory: bool | None = {"auto": None, "memory": True, "external": False}[
        sort_strategy
    ]

    # Determine if bins is a path or an integer
    try:
        binsize = int(bins)
        bins_path = None
    except ValueError:
        binsize = None
        bins_path = Path(bins)
        if not bins_path.exists():
            raise click.BadParameter(f"Bins file not found: {bins}") from None

    # ``commands`` is provenance and must survive the user's metadata file;
    # load_metadata_file() rejects a file that tries to set it.
    meta: dict = {"commands": [reconstruct_command()]}
    meta.update(load_metadata_file(metadata))

    from_pairs(
        uri=output,
        bins=bins_path,
        pairs=list(pairs),
        contig_genome_map=contigs,
        contig_table_format=contigs_format,
        fasta=fasta,
        assembly_name=assembly_name,
        contig_ordering=cast(
            Literal["genome-then-length", "genome-then-alpha", "fasta", "alpha"],
            contig_order,
        ),
        drop_unassigned=drop_unassigned,
        filters=list(filters) or None,
        binsize=binsize,
        genomes_table=genomes_table,
        genomes_table_format=genomes_format,
        chunksize=chunksize,
        sort_in_memory=sort_in_memory,
        metadata=meta,
        qc=qc,
        qc_distance_threshold=qc_distance_threshold,
        ps=ps,
        ps_bins_per_decade=ps_bins_per_decade,
        force=force,
        nproc=nproc,
    )
