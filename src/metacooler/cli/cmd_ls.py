"""metacooler ls command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli


def _is_mcool(file_path: Path) -> bool:
    """Return ``True`` if *file_path* is a cooler multi-resolution (``.mcool``) file.

    A file is considered an mcool if it contains a ``/resolutions`` group
    and at least one child group inside it with the cooler ``format``
    attribute set to ``HDF5::Cooler``.

    Args:
        file_path: Path to the HDF5 file to inspect.

    Returns:
        ``True`` if the file is a valid mcool collection, ``False``
        otherwise (including if the file cannot be opened).
    """
    import h5py

    from metacooler.schema import COOLER_FORMAT

    try:
        with h5py.File(file_path, "r") as f:
            if "resolutions" not in f:
                return False
            # At least one child must look like a cooler
            for name in f["resolutions"]:
                child = f["resolutions"][name]
                if (
                    isinstance(child, h5py.Group)
                    and child.attrs.get("format") == COOLER_FORMAT
                ):
                    return True
            return False
    except Exception:
        return False


@cli.command(name="ls", help="List entries in a metacool/mmetacool/mcool file.")
@click.argument("metacool_path", type=MetaCoolURI(), metavar="METACOOL_PATH")
@verbose_option
@quiet_option
def cmd_ls(metacool_path: str, verbose: bool, quiet: bool) -> None:
    """List entries in a metacool, mmetacool, or mcool file.

    \b
    For .mmetacool: lists all entries with full URI path.
    For .mcool: lists resolutions under /resolutions/.
    For .metacool: returns the root path.

    A URI that already addresses a single entry (``file::entries/name``) is
    echoed back unchanged — there is nothing below it to enumerate.
    """
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import list_entries
    from metacooler.schema import is_mmetacooler
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler ls", quiet=quiet)

    file_path, group_path = parse_uri(metacool_path)

    if group_path != "/":
        # The URI already names a single map — echo it back so `ls` output is
        # always a list of URIs this CLI can consume.
        click.echo(f"{file_path}::{group_path.lstrip('/')}")
        return

    if is_mmetacooler(file_path):
        entries = list_entries(file_path)
        if not entries:
            click.echo(f"No entries in {file_path}")
            return
        for entry in entries:
            click.echo(f"{file_path}::entries/{entry}")
    elif _is_mcool(file_path):
        with h5py.File(file_path, "r") as f:
            for name in f["resolutions"]:
                click.echo(f"{file_path}::resolutions/{name}")
    else:
        click.echo(f"{file_path}::/")
