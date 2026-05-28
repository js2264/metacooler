"""metacooler create cool command — create a metacool file from a cool file."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal, cast

import rich_click as click
import yaml

from metacooler.cli.cmd_create import _reconstruct_command, create, shared_options
from metacooler.cli.cmd_metadata import _stringify_values


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
    type=click.Choice(
        ["cool", "genome-then-length", "genome-then-alpha", "alpha", "custom"]
    ),
    default="cool",
    show_default=True,
    help=(
        """
        Strategy for ordering contigs in the output metacool file.\n
        - [yellow]cool[/yellow]: preserve the exact contig order from the input cool file.\n
        - [yellow]genome-then-length[/yellow]: group by genome, longest-first within each genome.\n
        - [yellow]genome-then-alpha[/yellow]: group by genome, alphabetical within each genome.\n
        - [yellow]alpha[/yellow]: sort all contigs alphabetically.\n
        - [yellow]custom[/yellow]: explicit order supplied via the API only.
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
    "--sort-in-memory/--sort-external",
    "sort_in_memory",
    is_flag=True,
    default=None,
    show_default="auto-detect",
    help=(
        "Control pixel sorting strategy. "
        "--sort-in-memory: collect all pixels in RAM and sort before writing (fast for small files). "
        "--sort-external: stream unsorted pixels to cooler's external sort (scalable for large files). "
        "Default: auto-detect based on pixel count (threshold: 100 M pixels)."
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
    nproc: int,
    sort_in_memory: bool | None,
    metadata: Path | None,
    force: bool,
    verbose: bool,
) -> None:
    """Create a metacool file from an existing cool file.

    Converts a standard cooler file to metacool format by adding the genome
    hierarchy. Bins and pixels from the input cool file are preserved exactly.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.create import from_cool

    setup_logging(verbose)

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    meta: dict = {"commands": [_reconstruct_command()]}
    if metadata is not None:
        data = metadata.read_text()
        if metadata.suffix in (".yaml", ".yml"):
            loaded = yaml.safe_load(data)
        else:
            loaded = json.loads(data)
        if not isinstance(loaded, dict):
            raise click.BadParameter("Metadata file must contain a JSON/YAML object")
        meta.update(_stringify_values(loaded))

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
            Literal[
                "genome-then-length", "genome-then-alpha", "cool", "alpha", "custom"
            ],
            contig_order,
        ),
        assembly_name=assembly_name,
        drop_unassigned=drop_unassigned,
        metadata=meta,
        nproc=nproc,
        sort_in_memory=sort_in_memory,
    )
