"""metacooler create pairs command — create a metacool file from a pairs file."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal, cast

import rich_click as click
import yaml

from metacooler.cli.cmd_create import _reconstruct_command, create, shared_options
from metacooler.cli.cmd_metadata import _stringify_values


@create.command(
    name="pairs", help="Create a metacool file from pairs and contig assignments."
)
@click.option(
    "--pairs",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Pairs file (plain text or gzipped).",
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
    type=click.Choice(
        ["genome-then-length", "genome-then-alpha", "fasta", "alpha", "custom"]
    ),
    default="genome-then-length",
    show_default=True,
    help=(
        """
        Strategy for ordering contigs in the output metacool file.\n
        - [yellow]genome-then-length[/yellow]: group by genome, longest-first within each genome.\n
        - [yellow]genome-then-alpha[/yellow]: group by genome, alphabetical within each genome.\n
        - [yellow]fasta[/yellow]: preserve order from the input FASTA (or pairs header).\n
        - [yellow]alpha[/yellow]: sort all contigs alphabetically.\n
        - [yellow]custom[/yellow]: explicit order supplied via the API only.
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
    default=1_000_000,
    show_default=True,
    help="Number of pairs to process per chunk.",
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
        "Default: auto-detect based on pairs file size (threshold: 2 GiB)."
    ),
)
def cmd_create_pairs(
    output: Path,
    bins: str,
    pairs: Path,
    contigs: Path,
    contigs_format: str,
    genomes_table: Path | None,
    fasta: Path | None,
    assembly_name: str | None,
    contig_order: str,
    drop_unassigned: bool,
    chunksize: int,
    sort_in_memory: bool | None,
    force: bool,
    verbose: bool,
    metadata: Path | None,
) -> None:
    """
    Create a metacool file from a pairs file, bins, and contig assignments.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.create import from_pairs

    setup_logging(verbose)

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    # Determine if bins is a path or an integer
    try:
        binsize = int(bins)
        bins_path = None
    except ValueError:
        binsize = None
        bins_path = Path(bins)
        if not bins_path.exists():
            raise click.BadParameter(f"Bins file not found: {bins}") from None

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

    from_pairs(
        uri=output,
        bins=bins_path,
        pairs=pairs,
        contig_genome_map=contigs,
        contig_table_format=contigs_format,  # type: ignore
        fasta=fasta,
        assembly_name=assembly_name,
        contig_ordering=cast(
            Literal[
                "genome-then-length", "genome-then-alpha", "fasta", "alpha", "custom"
            ],
            contig_order,
        ),
        drop_unassigned=drop_unassigned,
        binsize=binsize,
        genomes_table=genomes_table,
        chunksize=chunksize,
        sort_in_memory=sort_in_memory,
        metadata=meta,
    )
