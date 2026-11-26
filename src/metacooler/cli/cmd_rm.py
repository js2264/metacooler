"""metacooler rm command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.cmd_cp import parse_entry_uri
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(
    name="rm",
    help="""Remove an entry from an mmetacool file.

    \b
    TARGET   Entry URI: file.mmetacool::entries/entry_name
    """,
)
@click.argument(
    "target",
    type=MetaCoolURI(require_entry=True),
    metavar="TARGET",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Remove without confirmation.",
)
@verbose_option
@quiet_option
def cmd_rm(target: str, force: bool, verbose: bool, quiet: bool) -> None:
    """Remove an entry from an mmetacool file."""
    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.io.mmetacool import list_entries, remove_entry

    setup_logging(verbose, command_name="metacooler rm", quiet=quiet)

    # File existence and "this is not a bare collection path" are enforced by
    # MetaCoolURI; parse_entry_uri rejects the remaining empty-name cases.
    file_path, entry_name = parse_entry_uri(target, "Target")

    entries = list_entries(file_path)
    if entry_name not in entries:
        raise click.ClickException(
            f"Entry {entry_name!r} not found. Available: {entries}"
        )

    # abort=True so a declined delete exits non-zero — `rm … && echo ok` used to
    # print "ok" after the user said no.
    if not force:
        click.confirm(f"Remove entry {entry_name!r} from {file_path}?", abort=True)

    with run_log(file_path, "Collection written"):
        remove_entry(file_path, entry_name)
