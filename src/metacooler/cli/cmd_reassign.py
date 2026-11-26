"""metacooler reassign command — re-attribute contigs to genomes (MAGs)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, cast

import rich_click as click

from metacooler.cli._options import (
    check_output,
    contigs_format_option,
    force_option,
    genomes_format_option,
    option_given,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.schema import DEFAULT_CHUNKSIZE, LAYOUT_CHOICES, Layout

if TYPE_CHECKING:
    from typing import Literal


@cli.command(
    name="reassign",
    help="Re-attribute contigs to genomes (MAGs) without re-binning contacts.",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--contigs",
    "contigs_table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Contig-to-genome assignment table. It need not list every contig — "
        "see [cyan]--missing[/cyan]. Only its two name columns are read; extra "
        "columns are ignored. Give it or [cyan]--round[/cyan]."
    ),
)
@contigs_format_option
@click.option(
    "--round",
    "round_name",
    default=None,
    metavar="NAME",
    help=(
        "Take the assignment from a clustering round of the map "
        "([cyan]dump -t clusters[/cyan] lists them) instead of a table: each "
        "contig of the round goes to the genome named after its cluster, and "
        "the clusters' own columns become the genomes' metadata. A round in "
        "which a contig is in two clusters is refused."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output .metacool file path. Required unless [cyan]--inplace[/cyan].",
)
@click.option(
    "--inplace",
    is_flag=True,
    help=(
        "Rewrite [cyan]METACOOL_PATH[/cyan] itself instead of writing a new "
        "file: through a copy beside it, which replaces it once it validates, "
        "under the map's lock. Not for an entry of a collection whose contigs "
        "have to move (see [cyan]--layout[/cyan])."
    ),
)
@click.option(
    "--genomes",
    "genomes_table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Genomes metadata table for the new genomes (genome names in the "
        "column [cyan]--genomes-format[/cyan] names). Replaces — does not "
        "merge with — the metadata columns carried over from the source "
        "[cyan]/genomes[/cyan] table."
    ),
)
@genomes_format_option
@click.option(
    "--missing",
    type=click.Choice(["keep", "unbinned"]),
    default=None,
    help=(
        "What to do with a contig the table (or the round) does not mention. "
        "[cyan]keep[/cyan]: leave it in its current genome, so a short table "
        "is an incremental edit. [cyan]unbinned[/cyan]: move it to "
        "[cyan]__unbinned__[/cyan], so the table is authoritative — use it for "
        "the output of a fresh binning run. Neither drops a contig; use "
        "[cyan]extract[/cyan] for that. [dim]\\[default: keep with "
        "--contigs, unbinned with --round][/dim]"
    ),
)
@click.option(
    "--layout",
    type=click.Choice(list(LAYOUT_CHOICES)),
    default="auto",
    show_default=True,
    help=(
        "How to lay the output's contigs out. [cyan]auto[/cyan]: re-group a "
        "file whose contig ordering is genome-grouped, leave any other file's "
        "order alone. [cyan]keep[/cyan]: never re-order (fast, but a genome's "
        "bins may stop being contiguous). "
        "[cyan]genome-then-alpha[/cyan]/[cyan]genome-then-length[/cyan]: force "
        "a re-grouping with that strategy. Re-ordering rewrites the whole "
        "pixel table."
    ),
)
@click.option(
    "--metadata",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to a JSON or YAML file with metadata to merge into the output.",
)
@click.option(
    "--qc/--no-qc",
    default=True,
    show_default=True,
    help=(
        "Keep the contact-quality block when the source carries one, for the "
        "new genomes: the per-contig intra-contig counts are carried (a "
        "[cyan]create pairs[/cyan] split measured on read positions survives) "
        "and the inter-contig ones re-measured from the pixels. A source "
        "without per-contig counts has its block recomputed from the pixels. "
        "With [cyan]--no-qc[/cyan] the block is dropped and not replaced."
    ),
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help=(
        "Pixels read per chunk: remapped on the re-ordering path, and "
        "classified when the QC block is stored. On the re-ordering path each "
        "chunk costs one [cyan]cooler.create[/cyan] call whose price "
        "scales with contig count, so smaller chunks cost disproportionately "
        "more; the trade-off is ~24 bytes of peak memory per buffered pixel."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_reassign(
    metacool_path: str,
    contigs_table: Path | None,
    contigs_format: str,
    round_name: str | None,
    output: Path | None,
    inplace: bool,
    genomes_table: Path | None,
    genomes_format: str,
    missing: Literal["keep", "unbinned"] | None,
    layout: str,
    metadata: Path | None,
    qc: bool,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Re-attribute the contigs of an existing metacool file to other genomes.

    \b
    The contact matrix is never recomputed: the contigs, their lengths and
    their bins are those of the source file, and every pixel keeps its count.
    Only the genome each contig belongs to changes, along with everything
    derived from it — the /genomes table, the genome_id column, and (when the
    layout has to be restored) the order the contigs are stored in. The new
    assignment is a --contigs table, or a --round of the map's own clusters.
    """
    from metacooler.cli._provenance import load_metadata_file, reconstruct_command
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.reassign import reassign

    setup_logging(verbose, command_name="metacooler reassign", quiet=quiet)

    if (contigs_table is None) == (round_name is None):
        raise click.UsageError(
            "Pass the new assignment as a --contigs table or as a --round of "
            "the map's clusters, and not both."
        )
    if round_name is not None and option_given("contigs_format"):
        raise click.UsageError(
            "--contigs-format says how to read a --contigs table; --round takes "
            "the assignment from the map itself."
        )
    if inplace and output is not None:
        raise click.UsageError(
            "--inplace rewrites the input file, so -o/--output is not expected. "
            "Pass one or the other."
        )
    if not inplace and output is None:
        raise click.UsageError(
            "Missing option '-o' / '--output'. Pass an output path, or "
            "--inplace to rewrite the input file."
        )

    check_output(
        output,
        force,
        inputs=[metacool_path, contigs_table, genomes_table, metadata],
        defer_removal=True,
    )

    meta: dict = {"commands": [reconstruct_command()]}
    meta.update(load_metadata_file(metadata))

    reassign(
        uri=metacool_path,
        contig_genome_map=contigs_table,
        output=output,
        round=round_name,
        contig_table_format=contigs_format,
        genomes_table=genomes_table,
        genomes_table_format=genomes_format,
        missing=missing,
        layout=cast("Layout", layout),
        inplace=inplace,
        metadata=meta,
        qc=qc,
        chunksize=chunksize,
        force=force,
    )
