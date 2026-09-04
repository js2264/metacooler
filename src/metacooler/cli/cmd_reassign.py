"""metacooler reassign command — re-attribute contigs to genomes (MAGs)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, cast

import rich_click as click

from metacooler.cli._options import (
    CONTIGS_FORMAT_HELP,
    check_output,
    force_option,
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
    required=True,
    help=(
        "Contig-to-genome assignment table. It need not list every contig — "
        "see [cyan]--missing[/cyan]. Only its two name columns are read; extra "
        "columns are ignored."
    ),
)
@click.option(
    "--contigs-format",
    type=str,
    default="contigs",
    show_default=False,
    help=CONTIGS_FORMAT_HELP,
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
        "file. Only possible when no contig has to move (see "
        "[cyan]--layout[/cyan]); keep a backup, since the source is modified."
    ),
)
@click.option(
    "--genomes",
    "genomes_table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Genomes metadata table for the new genomes (first column "
        "[cyan]genome_name[/cyan]). Replaces — does not merge with — the "
        "metadata columns carried over from the source [cyan]/genomes[/cyan] "
        "table."
    ),
)
@click.option(
    "--missing",
    type=click.Choice(["keep", "unbinned"]),
    default="keep",
    show_default=True,
    help=(
        "What to do with a contig the table does not mention. "
        "[cyan]keep[/cyan]: leave it in its current genome, so a short table "
        "is an incremental edit. [cyan]unbinned[/cyan]: move it to "
        "[cyan]__unbinned__[/cyan], so the table is authoritative — use it for "
        "the output of a fresh binning run. Neither drops a contig; use "
        "[cyan]extract[/cyan] for that."
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
        "Recompute the contact-quality block when the source carries one — the "
        "reassignment invalidates it. With [cyan]--no-qc[/cyan] the stale block "
        "is dropped and not replaced."
    ),
)
@click.option(
    "--chunksize",
    type=int,
    default=DEFAULT_CHUNKSIZE,
    show_default=True,
    help=(
        "Source pixels read and remapped per chunk, on the re-ordering path "
        "only. Each chunk costs one [cyan]cooler.create[/cyan] call whose price "
        "scales with contig count, so smaller chunks cost disproportionately "
        "more; the trade-off is ~24 bytes of peak memory per buffered pixel."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_reassign(
    metacool_path: str,
    contigs_table: Path,
    contigs_format: str,
    output: Path | None,
    inplace: bool,
    genomes_table: Path | None,
    missing: Literal["keep", "unbinned"],
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
    layout has to be restored) the order the contigs are stored in.
    """
    from metacooler.cli._provenance import load_metadata_file, reconstruct_command
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.reassign import reassign

    setup_logging(verbose, command_name="metacooler reassign", quiet=quiet)

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

    check_output(output, force)

    meta: dict = {"commands": [reconstruct_command()]}
    meta.update(load_metadata_file(metadata))

    reassign(
        uri=metacool_path,
        contig_genome_map=contigs_table,
        output=output,
        contig_table_format=contigs_format,
        genomes_table=genomes_table,
        missing=missing,
        layout=cast("Layout", layout),
        inplace=inplace,
        metadata=meta,
        qc=qc,
        chunksize=chunksize,
        force=force,
    )
