"""metacooler metadata command."""

from __future__ import annotations

import json
from pathlib import Path

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli

# Every metadata subcommand reads or writes the attributes of exactly one map.
# ``require_entry=True`` is what stops ``metadata add coll.mmetacool -k x`` from
# writing to the collection root, where the value is invisible from every entry.
_ENTRY_URI = MetaCoolURI(require_entry=True)


@cli.group(name="metadata", help="Manipulate metadata in a metacool file.")
def metadata() -> None:
    """Add, remove, or list metadata in a metacool file."""
    pass


@metadata.command(name="list")
@click.argument("metacool_path", type=_ENTRY_URI)
@verbose_option
@quiet_option
def metadata_list(metacool_path: str, verbose: bool, quiet: bool) -> None:
    """List only the user metadata entries in a metacool file."""
    import h5py

    from metacooler.cli._provenance import read_metadata
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata list", quiet=quiet)

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = read_metadata(grp)
    click.echo(json.dumps(meta, indent=2, default=str))


@metadata.command(name="get")
@click.argument("metacool_path", type=_ENTRY_URI)
@click.option("-k", "--key", required=True, type=str, help="Metadata key to retrieve.")
@verbose_option
@quiet_option
def metadata_get(metacool_path: str, key: str, verbose: bool, quiet: bool) -> None:
    """Get a single metadata value by key."""
    import h5py

    from metacooler.cli._provenance import read_metadata
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata get", quiet=quiet)

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = read_metadata(grp)
        if key not in meta:
            raise click.BadParameter(
                f"Key not found: {key}. Available keys: {list(meta.keys())}"
            )
        value = meta[key]
    click.echo(json.dumps(value, indent=2, default=str))


@metadata.command(name="add")
@click.argument("metacool_path", type=_ENTRY_URI)
@click.option(
    "-k", "--key", required=True, type=str, help="Metadata key to add or update."
)
@click.option("--value", required=True, type=str, help="Value to store.")
@verbose_option
@quiet_option
def metadata_add(
    metacool_path: str, key: str, value: str, verbose: bool, quiet: bool
) -> None:
    """Add or update a metadata key-value pair."""
    import h5py

    from metacooler.cli._provenance import (
        read_metadata,
        refuse_reserved_metadata_key,
    )
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata add", quiet=quiet)
    refuse_reserved_metadata_key(key, "set")

    # Try to parse as JSON for structured values
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        parsed = value

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = read_metadata(grp)
        meta[key] = parsed
        grp.attrs["metadata"] = json.dumps(meta)


@metadata.command(name="remove")
@click.argument("metacool_path", type=_ENTRY_URI)
@click.option("-k", "--key", required=True, type=str, help="Metadata key to remove.")
@verbose_option
@quiet_option
def metadata_remove(metacool_path: str, key: str, verbose: bool, quiet: bool) -> None:
    """Remove a metadata key."""
    import h5py

    from metacooler.cli._provenance import (
        read_metadata,
        refuse_reserved_metadata_key,
    )
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata remove", quiet=quiet)
    refuse_reserved_metadata_key(key, "remove")

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        meta = read_metadata(grp)
        if key not in meta:
            raise click.BadParameter(
                f"Key not found: {key}. Available keys: {list(meta.keys())}"
            )
        del meta[key]
        grp.attrs["metadata"] = json.dumps(meta)


@metadata.command(name="load")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("metacool_path", type=_ENTRY_URI)
@verbose_option
@quiet_option
def metadata_load(
    input_file: Path, metacool_path: str, verbose: bool, quiet: bool
) -> None:
    """Load metadata from a JSON or YAML file into the 'metadata' attribute."""
    import h5py

    from metacooler.cli._provenance import load_metadata_file, read_metadata
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler metadata load", quiet=quiet)

    # Shared loader: parses JSON/YAML, stringifies leaves, and refuses the
    # reserved 'commands' provenance key.
    meta = load_metadata_file(input_file)

    file_path, group_path = parse_uri(metacool_path)
    with h5py.File(file_path, "a") as f:
        grp = f[group_path] if group_path != "/" else f
        existing = read_metadata(grp)
        existing.update(meta)
        grp.attrs["metadata"] = json.dumps(existing)
