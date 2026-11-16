"""metacooler annotate command — add, replace or drop metadata columns."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import rich_click as click

from metacooler.cli._options import (
    check_output,
    contig_metadata_format_option,
    force_option,
    genomes_format_option,
    parse_comma_list,
    quiet_option,
    verbose_option,
)
from metacooler.cli.main import MetaCoolURI, cli

if TYPE_CHECKING:
    from typing import Literal


@cli.command(
    name="annotate",
    help="Add, replace or drop metadata columns of /genomes and /contigs.",
)
# Not require_entry: a bare .mmetacool annotates every entry of the collection.
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "--genomes",
    "genomes_table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Genomes metadata table. Its genome-name column (see "
        "[cyan]--genomes-format[/cyan]) is matched against the map's genomes; "
        "every other column is added to [cyan]/genomes[/cyan]."
    ),
)
@genomes_format_option
@click.option(
    "--contigs",
    "contigs_table",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Per-contig metadata table. Its contig-name column (see "
        "[cyan]--contigs-format[/cyan]) is matched against the map's contigs; "
        "an optional 'genome_name' column is checked against the genome "
        "assignment, never applied; every other column is added to "
        "[cyan]/contigs[/cyan] (ASCII only: cooler reads that table)."
    ),
)
@contig_metadata_format_option
@click.option(
    "--columns",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    metavar="TABLE.COLUMN",
    help=(
        "Take only these columns of the tables, as [cyan]genomes.<name>[/cyan] "
        "or [cyan]contigs.<name>[/cyan] (comma-separated or repeated). "
        "Default: every metadata column."
    ),
)
@click.option(
    "--drop",
    type=str,
    multiple=True,
    callback=parse_comma_list,
    metavar="TABLE.COLUMN",
    help=(
        "Remove these metadata columns, as [cyan]genomes.<name>[/cyan] or "
        "[cyan]contigs.<name>[/cyan] (comma-separated or repeated). Columns "
        "metacooler derives cannot be dropped."
    ),
)
@click.option(
    "--replace",
    is_flag=True,
    help=(
        "Allow changing the values of a column the map already has. A column "
        "given the values it already holds is never a change, so re-running "
        "is safe without it."
    ),
)
@click.option(
    "--missing",
    type=click.Choice(["keep", "na"]),
    default="keep",
    show_default=True,
    help=(
        "For a replaced column, what the rows the table does not list get. "
        "[cyan]keep[/cyan]: their current value, so the table is a patch. "
        "[cyan]na[/cyan]: an empty value, so the table is authoritative. New "
        "columns always leave unlisted rows empty."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help=(
        "Write an annotated copy here. A collection stays a collection; an "
        "entry URI gives a standalone .metacool. Required unless "
        "[cyan]--inplace[/cyan]."
    ),
)
@click.option(
    "--inplace",
    is_flag=True,
    help=(
        "Edit [cyan]METACOOL_PATH[/cyan] itself. Only /genomes and /contigs "
        "are rewritten, so this costs seconds even on a large map."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_annotate(
    metacool_path: str,
    genomes_table: Path | None,
    genomes_format: str,
    contigs_table: Path | None,
    contigs_format: str,
    columns: list[str],
    drop: list[str],
    replace: bool,
    missing: Literal["keep", "na"],
    output: Path | None,
    inplace: bool,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Add, replace or drop metadata columns of an existing map.

    \b
    Rows are matched by name and never added, removed or reordered; the
    columns metacooler derives, the bins and the pixels are left alone. A bare
    .mmetacool path annotates every entry of the collection.
    """
    from metacooler.cli._provenance import reconstruct_command
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.annotate import annotate

    setup_logging(verbose, command_name="metacooler annotate", quiet=quiet)

    if inplace and output is not None:
        raise click.UsageError(
            "--inplace edits the input file, so -o/--output is not expected. "
            "Pass one or the other."
        )
    if not inplace and output is None:
        raise click.UsageError(
            "Missing option '-o' / '--output'. Pass an output path, or "
            "--inplace to edit the input file."
        )
    if genomes_table is None and contigs_table is None and not drop:
        raise click.UsageError("Nothing to do: pass --genomes, --contigs, or --drop.")

    check_output(
        output,
        force,
        inputs=[metacool_path, genomes_table, contigs_table],
        defer_removal=True,
    )

    annotate(
        metacool_path,
        output,
        genomes_table=genomes_table,
        contigs_table=contigs_table,
        genomes_table_format=genomes_format,
        contigs_table_format=contigs_format,
        columns=columns or None,
        drop=drop or None,
        replace=replace,
        missing=missing,
        inplace=inplace,
        metadata={"commands": [reconstruct_command()]},
        force=force,
    )
