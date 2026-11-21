"""metacooler cp command."""

from __future__ import annotations

import os
from pathlib import Path

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli


def parse_entry_uri(uri: str, role: str) -> tuple[Path, str]:
    """Split an ``file::entries/name`` URI into ``(file path, entry name)``.

    A bare path (no ``::group``) used to fall through the old
    ``group_path.rstrip("/").split("/")[-1]`` idiom as an empty entry name,
    producing the useless ``Entry '' not found``.  This is the single guard
    that turns it into an actionable message.  Destination URIs cannot be
    validated by :class:`~metacooler.cli.main.MetaCoolURI` because the group
    they name does not exist yet, so they go through here instead.

    Args:
        uri: Raw CLI URI string.
        role: Human-readable role of the argument (``"Source"`` /
            ``"Destination"``), used in the error message.

    Returns:
        A ``(file_path, entry_name)`` tuple.

    Raises:
        click.ClickException: If the URI does not name an entry.
    """
    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(uri)
    entry_name = group_path.rstrip("/").rsplit("/", 1)[-1]
    if group_path == "/" or not entry_name:
        raise click.ClickException(
            f"{role} must name an entry with URI syntax, e.g. "
            f"'{file_path}::entries/<name>' (got {uri!r})."
        )
    return file_path, entry_name


@cli.command(
    name="cp",
    help="""Copy an entry within or between mmetacool files.

    \b
    SRC   Source URI: file.mmetacool::entries/entry_name
    DST   Destination URI: file.mmetacool::entries/entry_name
    Cross-file copy is supported.
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
def cmd_cp(src: str, dst: str, force: bool, verbose: bool, quiet: bool) -> None:
    """Copy an entry within or between mmetacool files."""
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import (
        copy_entry,
        copy_entry_cross_file,
        init_mmetacool_attrs,
        list_entries,
    )

    setup_logging(verbose, command_name="metacooler cp", quiet=quiet)

    src_path, src_name = parse_entry_uri(src, "Source")
    dst_path, dst_name = parse_entry_uri(dst, "Destination")

    # Source existence is guaranteed by MetaCoolURI; the entry inside it is not.
    src_entries = list_entries(src_path)
    if src_name not in src_entries:
        raise click.ClickException(
            f"Source entry {src_name!r} not found. Available: {src_entries}"
        )

    same_file = src_path.resolve() == dst_path.resolve()

    if same_file:
        dst_entries = src_entries
    else:
        # Only now that the source is known good may we touch the destination —
        # creating it first left a stray empty .mmetacool behind on a typo.
        if not os.path.lexists(dst_path):
            import h5py

            from metacooler.io.staging import writing_output

            # Created exclusively: a collection another writer created
            # meanwhile -- another cp, an entry copied into it -- is copied
            # into, never truncated.
            try:
                with (
                    writing_output(dst_path, force=False),
                    h5py.File(dst_path, "w-") as f,
                ):
                    init_mmetacool_attrs(f)
            except FileExistsError:
                pass
            else:
                click.echo(f"Initialised new mmetacool file: {dst_path}")
        dst_entries = list_entries(dst_path)

    # Answering "y" must actually overwrite: the old code fell through with
    # overwrite=force (still False) and the io layer raised anyway.
    overwrite = force
    if dst_name in dst_entries and not force:
        where = "" if same_file else f" in {dst_path}"
        click.confirm(
            f"Destination entry {dst_name!r} already exists{where}. Overwrite?",
            abort=True,
        )
        overwrite = True

    if same_file:
        copy_entry(src_path, src_name, dst_name, overwrite=overwrite)
    else:
        copy_entry_cross_file(
            src_path, src_name, dst_path, dst_name, overwrite=overwrite
        )
