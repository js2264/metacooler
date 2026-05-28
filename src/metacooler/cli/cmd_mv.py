"""metacooler mv command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli.main import cli


@cli.command(name="mv", help="Move (rename) an entry within an mmetacool file.")
@click.argument(
    "src",
    type=str,
    metavar="SRC",
    help="Source URI: [cyan]file.mmetacool::entries/entry_name[/cyan].",
)
@click.argument(
    "dst",
    type=str,
    metavar="DST",
    help="Destination URI (same file): [cyan]file.mmetacool::entries/new_name[/cyan].",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite destination without confirmation.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_mv(src: str, dst: str, force: bool, verbose: bool) -> None:
    """Move (rename) an entry within an mmetacool file.

    \b
    SRC and DST use URI syntax: file.mmetacool::entries/entry_name
    Both must be in the same file.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import list_entries, move_entry
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler mv")

    src_path, src_group = parse_uri(src)
    dst_path, dst_group = parse_uri(dst)

    if src_path != dst_path:
        raise click.ClickException(
            "Move is only supported within the same file. Use 'cp' + 'rm' for cross-file moves."
        )

    src_name = src_group.rstrip("/").split("/")[-1]
    dst_name = dst_group.rstrip("/").split("/")[-1]

    if not src_path.exists():
        raise click.ClickException(f"File not found: {src_path}")

    entries = list_entries(src_path)
    if src_name not in entries:
        raise click.ClickException(
            f"Source entry {src_name!r} not found. Available: {entries}"
        )
    if (
        dst_name in entries
        and not force
        and not click.confirm(
            f"Destination entry {dst_name!r} already exists. Overwrite?"
        )
    ):
        click.echo("Aborted.")
        return

    move_entry(src_path, src_name, dst_name, overwrite=force)
