"""metacooler ls command."""

from __future__ import annotations

from pathlib import Path

import h5py
import rich_click as click

from metacooler.cli.main import cli
from metacooler.schema import COOLER_FORMAT


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
@click.argument("file_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_ls(file_path: Path, verbose: bool) -> None:
    """List entries in a metacool, mmetacool, or mcool file.

    \b
    For .mmetacool: lists all entries with full URI path.
    For .mcool: lists resolutions under /resolutions/.
    For .metacool: returns the root path.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import list_entries
    from metacooler.schema import is_mmetacooler

    setup_logging(verbose, command_name="metacooler ls")

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
