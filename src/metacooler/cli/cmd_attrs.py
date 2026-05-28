"""metacooler attrs command - Display HDF5 attributes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import rich_click as click
import yaml

if TYPE_CHECKING:
    import h5py

from metacooler.cli.main import MetaCoolURI, cli


def _convert_attr(val: object) -> object:
    """Convert an HDF5 attribute value to a YAML-friendly Python type.

    Handles numpy scalars (→ int/float), bytes (→ str), and strings that
    encode JSON objects, ISO datetimes, or the literal ``"null"`` / ``"{}"``.

    Args:
        val: Raw attribute value as returned by h5py.

    Returns:
        A plain Python int, float, str, dict, list, datetime, or ``None``.
    """
    import json
    from datetime import datetime

    import numpy as np

    if isinstance(val, np.integer):
        return int(val)
    if isinstance(val, np.floating):
        return float(val)
    if isinstance(val, bytes):
        val = val.decode()
    if isinstance(val, str):
        if val == "null":
            return None
        if val == "{}":
            return {}

        try:
            return json.loads(val)
        except ValueError:
            pass
        try:
            return datetime.fromisoformat(val)
        except ValueError:
            pass
    return val


def _collect_attrs_tree(obj: h5py.HLObject) -> dict:
    """Recursively build a nested dict of HDF5 attributes.

    Mirrors the output format of ``cooler attrs``: each node in the tree
    has an ``"@attrs"`` key containing its attribute dictionary, and child
    groups are nested under their HDF5 key names.

    Args:
        obj: An open ``h5py.Group`` or ``h5py.Dataset``.

    Returns:
        Nested dict suitable for YAML serialisation.  Inaccessible
        children (broken soft-links, permission errors) are represented
        as ``{"@error": "inaccessible"}``.
    """
    import h5py as _h5py

    result: dict = {}
    result["@attrs"] = {k: _convert_attr(v) for k, v in obj.attrs.items()}
    if isinstance(obj, _h5py.Group):
        for key in obj:
            try:
                child = obj[key]
                result[key] = _collect_attrs_tree(child)
            except Exception:
                # Skip broken soft links or inaccessible objects
                result[key] = {"@error": "inaccessible"}
    return result


@cli.command(
    name="attrs", help="Display HDF5 attributes (YAML, same format as cooler attrs)."
)
@click.argument("metacool_path", type=MetaCoolURI())
@click.argument(
    "group_path",
    type=str,
    default="/",
    metavar="GROUP_PATH",
    help="HDF5 group path to inspect (default: [cyan]/[/cyan] — the file root).",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_attrs(metacool_path: str, group_path: str, verbose: bool) -> None:
    """Dump HDF5 attributes as YAML, matching cooler attrs output format.

    \b
    GROUP_PATH    HDF5 group path (default: /)
    """
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.schema import is_mmetacooler
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler attrs")

    file_path, uri_group = parse_uri(metacool_path)
    # Merge URI group path with explicit group_path argument
    resolved_group = uri_group if group_path == "/" else group_path

    with h5py.File(file_path, "r") as f:
        if is_mmetacooler(file_path) and resolved_group == "/":
            # Show attrs for ALL contained metacoolers
            all_data = {}
            all_data["@attrs"] = {k: _convert_attr(v) for k, v in f.attrs.items()}
            if "entries" in f:
                entries_data = {}
                for entry_name in f["entries"]:
                    entry = f["entries"][entry_name]
                    entries_data[entry_name] = _collect_attrs_tree(entry)
                all_data["entries"] = entries_data
            click.echo(
                yaml.dump(all_data, default_flow_style=False, allow_unicode=True),
                nl=False,
            )
        else:
            grp = f[resolved_group]
            data = _collect_attrs_tree(grp)
            click.echo(
                yaml.dump(data, default_flow_style=False, allow_unicode=True),
                nl=False,
            )
