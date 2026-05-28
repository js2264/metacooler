"""metacooler tree command - Display the data hierarchy."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import rich_click as click
from rich.console import Console
from rich.tree import Tree

if TYPE_CHECKING:
    import h5py

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="tree", help="Display the data hierarchy.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_tree(metacool_path: str, verbose: bool) -> None:
    """Display the HDF5 hierarchy as a tree."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.schema import is_metacooler, is_mmetacooler
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler tree")
    console = Console()

    file_path, group_path = parse_uri(metacool_path)

    def _build_tree(name: str, obj, parent: Tree) -> None:
        if isinstance(obj, h5py.Group):
            node = parent.add(f"[bold]{name}/[/bold]")
            for key, child in obj.items():
                _build_tree(key, child, node)
        elif isinstance(obj, h5py.Dataset):
            shape = " × ".join(str(s) for s in obj.shape)
            dtype = str(obj.dtype)
            parent.add(f"{name} [dim]({shape}, {dtype})[/dim]")

    with h5py.File(file_path, "r") as f:
        tree = Tree(f"[bold]{file_path.name}[/bold]")

        if group_path != "/":
            # URI points to a specific entry — show its tree
            entry = f[group_path]
            for key in ("genomes", "contigs", "bins", "pixels", "indexes"):
                if key in entry:
                    _build_tree(key, entry[key], tree)
        elif is_mmetacooler(file_path):
            # Show top-level entries, then tree of each entry
            if "entries" in f:
                entries_node = tree.add("[bold]entries/[/bold]")
                for entry_name in f["entries"]:
                    entry_node = entries_node.add(f"[bold]{entry_name}/[/bold]")
                    entry = f["entries"][entry_name]
                    for key in ("genomes", "contigs", "bins", "pixels", "indexes"):
                        if key in entry:
                            _build_tree(key, entry[key], entry_node)
        elif is_metacooler(file_path):
            # Metacool file: show core tables in fixed order.
            # chroms is a shim (soft link to contigs) and is omitted.
            for key in ("genomes", "contigs", "bins", "pixels", "indexes"):
                if key in f:
                    _build_tree(key, f[key], tree)
        else:
            # Plain cooler or other HDF5: show all top-level groups.
            for key in f:
                _build_tree(key, f[key], tree)

        console.print(tree)
