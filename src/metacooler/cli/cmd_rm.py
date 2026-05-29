"""metacooler rm command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli.main import cli


@cli.command(
    name="rm",
    help="""Remove an entry from an mmetacool file.

    \b
    TARGET   Entry URI: file.mmetacool::entries/entry_name
    """,
)
@click.argument(
    "target",
    type=str,
    metavar="TARGET",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Remove without confirmation.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_rm(target: str, force: bool, verbose: bool) -> None:
    """Remove an entry from an mmetacool file."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import list_entries, remove_entry
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler rm")

    file_path, group_path = parse_uri(target)
    entry_name = group_path.rstrip("/").split("/")[-1]

    if not file_path.exists():
        raise click.ClickException(f"File not found: {file_path}")

    entries = list_entries(file_path)
    if entry_name not in entries:
        raise click.ClickException(
            f"Entry {entry_name!r} not found. Available: {entries}"
        )

    if not force and not click.confirm(
        f"Remove entry {entry_name!r} from {file_path}?"
    ):
        click.echo("Aborted.")
        return

    remove_entry(file_path, entry_name)
