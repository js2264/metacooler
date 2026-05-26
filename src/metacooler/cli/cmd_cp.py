"""metacooler cp command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli.main import cli


@cli.command(name="cp", help="Copy an entry within or between mmetacool files.")
@click.argument("src", type=str)
@click.argument("dst", type=str)
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
def cmd_cp(src: str, dst: str, force: bool, verbose: bool) -> None:
    """Copy an entry within or between mmetacool files.

    \b
    SRC and DST use URI syntax: file.mmetacool::entries/entry_name
    Cross-file copy is supported.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import (
        copy_entry,
        copy_entry_cross_file,
        list_entries,
    )
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler cp")

    src_path, src_group = parse_uri(src)
    dst_path, dst_group = parse_uri(dst)

    # Extract entry name from group path (last component)
    src_name = src_group.rstrip("/").split("/")[-1]
    dst_name = dst_group.rstrip("/").split("/")[-1]

    if not src_path.exists():
        raise click.ClickException(f"Source file not found: {src_path}")

    if src_path == dst_path:
        # Same-file copy — dst must already exist
        if not dst_path.exists():
            raise click.ClickException(f"Destination file not found: {dst_path}")
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
        copy_entry(src_path, src_name, dst_name, overwrite=force)
    else:
        # Cross-file copy — initialise a new mmetacool if dst does not exist
        if not dst_path.exists():
            import h5py as _h5py

            with _h5py.File(dst_path, "w") as _f:
                from metacooler.io.mmetacool import _init_mmetacool_attrs

                _init_mmetacool_attrs(_f)
            click.echo(f"Initialised new mmetacool file: {dst_path}")
        entries = list_entries(dst_path)
        if (
            dst_name in entries
            and not force
            and not click.confirm(
                f"Destination entry {dst_name!r} already exists in {dst_path}. Overwrite?"
            )
        ):
            click.echo("Aborted.")
            return
        copy_entry_cross_file(src_path, src_name, dst_path, dst_name, overwrite=force)
