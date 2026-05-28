"""metacooler metadata command."""

from __future__ import annotations

import json
from pathlib import Path

import rich_click as click
import yaml

from metacooler.cli.main import MetaCoolURI, cli


@cli.group(name="metadata", help="Manipulate metadata in a metacool file.")
def metadata() -> None:
    """Add, remove, or list metadata in a metacool file."""
    pass


def _check_mmetacool(metacool_path: str) -> None:
    """Raise a ClickException if the URI points to an mmetacool without an entry.

    When a bare mmetacool file path is given (no ``::entries/…`` suffix),
    this helper prints the available entry URIs so the user knows what to
    specify.

    Args:
        metacool_path: Raw CLI argument string — plain path or URI.

    Raises:
        click.ClickException: If the file is an mmetacool and no entry
            path was specified in the URI.
    """
    from metacooler.schema import is_mmetacooler, list_metacoolers
    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(metacool_path)
    if group_path != "/":
        # URI already specifies an entry — no need to check
        return
    if is_mmetacooler(file_path):
        entries = list_metacoolers(file_path)
        entry_uris = [f"{file_path}::{e}" for e in entries]
        raise click.ClickException(
            f"This is an mmetacool file with {len(entries)} entries. "
            f"Specify an entry with URI syntax:\n"
            + "\n".join(f"  {u}" for u in entry_uris)
        )


@metadata.command(name="list")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def metadata_list(metacool_path: str, verbose: bool) -> None:
    """List only the user metadata entries in a metacool file."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata list")
    _check_mmetacool(metacool_path)

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        raw = grp.attrs.get("metadata", "{}")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        meta = json.loads(raw) if isinstance(raw, str) else dict(raw)
    click.echo(json.dumps(meta, indent=2, default=str))


def _read_metadata(f):
    """Read and parse the ``metadata`` attribute from an open HDF5 group.

    Args:
        f: An open ``h5py.File`` or ``h5py.Group`` object.

    Returns:
        Parsed metadata as a plain Python dict.  Returns ``{}`` if the
        attribute is absent or is an empty JSON object.
    """
    raw = f.attrs.get("metadata", "{}")
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return json.loads(raw) if isinstance(raw, str) else dict(raw)


@metadata.command(name="get")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option("-k", "--key", required=True, type=str, help="Metadata key to retrieve.")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def metadata_get(metacool_path: str, key: str, verbose: bool) -> None:
    """Get a single metadata value by key."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata get")

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = _read_metadata(grp)
        if key not in meta:
            raise click.BadParameter(f"Key not found: {key}")
        value = meta[key]
    click.echo(json.dumps(value, indent=2, default=str))


@metadata.command(name="add")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-k", "--key", required=True, type=str, help="Metadata key to add or update."
)
@click.option("--value", required=True, type=str, help="Value to store.")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def metadata_add(metacool_path: str, key: str, value: str, verbose: bool) -> None:
    """Add or update a metadata key-value pair."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata add")

    # Try to parse as JSON for structured values
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        parsed = value

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = _read_metadata(grp)
        meta[key] = parsed
        grp.attrs["metadata"] = json.dumps(meta)


@metadata.command(name="remove")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option("-k", "--key", required=True, type=str, help="Metadata key to remove.")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def metadata_remove(metacool_path: str, key: str, verbose: bool) -> None:
    """Remove a metadata key."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata remove")

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = _read_metadata(grp)
        if key not in meta:
            raise click.BadParameter(f"Key not found: {key}")
        del meta[key]
        grp.attrs["metadata"] = json.dumps(meta)


def _stringify_values(obj):
    """Recursively convert every leaf value in a nested structure to a string.

    Args:
        obj: Any Python object.  Dicts and lists are traversed recursively;
            all other values are converted with ``str()``.

    Returns:
        A new structure of the same shape with all leaf values replaced by
        their string representations.
    """
    if isinstance(obj, dict):
        return {k: _stringify_values(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_stringify_values(v) for v in obj]
    return str(obj)


@metadata.command(name="load")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def metadata_load(input_file: Path, metacool_path: str, verbose: bool) -> None:
    """Load metadata from a JSON or YAML file into the 'metadata' attribute."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata load")

    data = input_file.read_text()
    if input_file.suffix in (".yaml", ".yml"):
        meta = yaml.safe_load(data)
    else:
        meta = json.loads(data)

    if not isinstance(meta, dict):
        raise click.BadParameter("Input file must contain a JSON/YAML object")

    meta = _stringify_values(meta)

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        existing = _read_metadata(grp)
        existing.update(meta)
        grp.attrs["metadata"] = json.dumps(existing)
