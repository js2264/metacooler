"""metacooler create cool command — create a metacool file from a cool file."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import rich_click as click

from metacooler.cli._options import check_output
from metacooler.cli._provenance import load_metadata_file, reconstruct_command
from metacooler.cli.cmd_create import create, shared_options


@create.command(name="cool", help="Create a metacool file from an existing cool file.")
@click.option(
    "--cool",
    "cool_path",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input cool file to convert.",
)
@shared_options
@click.option(
    "--contig-order",
    type=click.Choice(["genome-then-length", "genome-then-alpha", "cool", "alpha"]),
    default="genome-then-length",
    show_default=True,
    help=(
        """
        Strategy for ordering contigs in the output metacool file.\n
        - [yellow]genome-then-length[/yellow]: group by genome, longest-first within each genome.\n
        - [yellow]genome-then-alpha[/yellow]: group by genome, alphabetical within each genome.\n
        - [yellow]cool[/yellow]: preserve the exact contig order from the input cool file.\n
        - [yellow]alpha[/yellow]: sort all contigs alphabetically.
        """
    ),
)
@click.option(
    "--assembly-name",
    type=str,
    help="Name for this assembly / metagenome, stored in file attrs.",
)
@click.option(
    "--drop-unassigned",
    is_flag=True,
    help="Drop contigs not present in the contig table instead of collecting them under __unbinned__.",
)
@click.option(
    "--nproc",
    type=int,
    default=1,
    show_default=True,
    help="Number of parallel reader threads for pixel remapping.",
)
@click.option(
    "--sort-strategy",
    "sort_strategy",
    type=click.Choice(["auto", "memory", "external"]),
    default="auto",
    show_default=True,
    help=(
        "Pixel sorting strategy. "
        "'auto': detect based on pixel count (threshold: 100 M pixels). "
        "'memory': collect all pixels in RAM and sort before writing (fast for small files). "
        "'external': stream unsorted pixels to cooler's external sort (scalable for large files)."
    ),
)
def cmd_create_cool(
    cool_path: Path,
    output: Path,
    contigs: Path,
    contigs_format: str,
    genomes_table: Path | None,
    contig_order: str,
    assembly_name: str | None,
    drop_unassigned: bool,
    filters: list[str],
    nproc: int,
    sort_strategy: str,
    metadata: Path | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Create a metacool file from an existing cool file.

    Converts a standard cooler file to metacool format by adding the genome
    hierarchy. Bins and pixels from the input cool file are preserved exactly.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.create import from_cool

    setup_logging(verbose, command_name="metacooler create cool", quiet=quiet)

    check_output(output, force)

    # Map --sort-strategy to sort_in_memory bool
    sort_in_memory: bool | None = {"auto": None, "memory": True, "external": False}[
        sort_strategy
    ]

    # ``commands`` is provenance and must survive the user's metadata file;
    # load_metadata_file() rejects a file that tries to set it.
    meta: dict = {"commands": [reconstruct_command()]}
    meta.update(load_metadata_file(metadata))

    from_cool(
        cool_uri=cool_path,
        output=output,
        contig_genome_map=contigs,
        contig_table_format=cast(
            Literal["contigs", "contig_genome", "metator", "I,J"],
            contigs_format,
        ),
        genomes_table=genomes_table,
        contig_ordering=cast(
            Literal["genome-then-length", "genome-then-alpha", "cool", "alpha"],
            contig_order,
        ),
        assembly_name=assembly_name,
        drop_unassigned=drop_unassigned,
        filters=list(filters) or None,
        metadata=meta,
        nproc=nproc,
        sort_in_memory=sort_in_memory,
    )
