"""metacooler mv command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.cmd_cp import parse_entry_uri
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(
    name="mv",
    help="""Move (rename) an entry within an mmetacool file.

    \b
    SRC   Source URI: file.mmetacool::entries/entry_name
    DST   Destination URI (same file): file.mmetacool::entries/new_name
    Both must be in the same file.
""",
)
@click.argument(
    "src",
    type=MetaCoolURI(require_entry=True),
    metavar="SRC",
)
@click.argument(
    "dst",
    type=str,
    metavar="DST",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite destination without confirmation.",
)
@verbose_option
@quiet_option
def cmd_mv(src: str, dst: str, force: bool, verbose: bool, quiet: bool) -> None:
    """Move (rename) an entry within an mmetacool file."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import list_entries, move_entry

    setup_logging(verbose, command_name="metacooler mv", quiet=quiet)

    src_path, src_name = parse_entry_uri(src, "Source")
    dst_path, dst_name = parse_entry_uri(dst, "Destination")

    if src_path.resolve() != dst_path.resolve():
        raise click.ClickException(
            "Move is only supported within the same file. "
            "Use 'cp' + 'rm' for cross-file moves."
        )

    # Source file existence is guaranteed by MetaCoolURI.
    entries = list_entries(src_path)
    if src_name not in entries:
        raise click.ClickException(
            f"Source entry {src_name!r} not found. Available: {entries}"
        )

    # Answering "y" must actually overwrite; declining must exit non-zero.
    overwrite = force
    if dst_name in entries and not force:
        click.confirm(
            f"Destination entry {dst_name!r} already exists. Overwrite?",
            abort=True,
        )
        overwrite = True

    move_entry(src_path, src_name, dst_name, overwrite=overwrite)
