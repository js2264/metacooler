"""Provenance and metadata helpers shared across metacooler CLI commands.

A few keys of a file's ``metadata`` attribute belong to metacooler, which
writes them itself (:data:`RESERVED_METADATA_KEYS`): ``commands`` records every
metacooler invocation that produced the file, ``qc`` the contact-quality block,
``reassign`` and ``annotate`` those commands' histories.  Every user path that
writes the attribute -- a ``--metadata`` file, ``metadata load``, ``metadata
add`` and ``metadata remove`` -- must refuse to set or remove them, and this
module is the single place that enforces it.
"""

from __future__ import annotations

import json
import shlex
from pathlib import Path
from typing import Any

import rich_click as click

#: Keys of the ``metadata`` attribute that metacooler writes itself, with what
#: each holds.  Anyone may read them (``metadata list``/``get``); no user path
#: may set or remove them.
RESERVED_METADATA_KEYS: dict[str, str] = {
    "commands": "reserved for provenance tracking",
    "qc": "reserved for the contact-quality statistics metacooler records",
    "reassign": "reserved for the history 'metacooler reassign' records",
    "annotate": "reserved for the history 'metacooler annotate' records",
}


def reconstruct_command() -> str:
    """Reconstruct the exact command-line invocation from the current Click context.

    Called from within a running Click command to produce a string
    representation of all non-default options and arguments, which is then
    stored in the output file's ``metadata.commands`` list for provenance.

    Returns:
        A shell-escaped string of the form ``metacooler <subcommand>
        [OPTIONS] [ARGS]``.

    Raises:
        RuntimeError: If called outside of an active Click command context.
    """
    ctx = click.get_current_context(silent=True)
    if ctx is None:
        raise RuntimeError(
            "reconstruct_command() called outside of a click command context."
        )
    parts = [ctx.command_path]
    for param in ctx.command.params:
        if isinstance(param, click.Option):
            value = ctx.params.get(param.name)
            if value is not None and value != param.default:
                opt_str = param.opts[0] if param.opts else f"--{param.name}"
                if param.is_flag:
                    if value:
                        parts.append(opt_str)
                elif param.multiple:
                    for v in value:
                        parts.extend([opt_str, str(v)])
                else:
                    parts.extend([opt_str, str(value)])
        elif isinstance(param, click.Argument):
            value = ctx.params.get(param.name)
            if value is not None:
                if isinstance(value, tuple):
                    for v in value:
                        parts.append(str(v))
                else:
                    parts.append(str(value))
    return shlex.join(parts)


def read_metadata(group) -> dict[str, Any]:
    """Read and parse the ``metadata`` attribute from an open HDF5 group.

    Args:
        group: An open ``h5py.File`` or ``h5py.Group`` object.

    Returns:
        Parsed metadata as a plain Python dict.  ``{}`` if the attribute is
        absent.
    """
    raw = group.attrs.get("metadata", "{}")
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return json.loads(raw) if isinstance(raw, str) else dict(raw)


def read_metadata_uri(uri: str) -> dict[str, Any]:
    """Read the ``metadata`` attribute of the group addressed by a metacool URI.

    Args:
        uri: ``file.metacool`` or ``file.mmetacool::entries/name`` URI.

    Returns:
        Parsed metadata dict, or ``{}`` if the file has none.
    """
    import h5py

    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(uri)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        return read_metadata(grp)


def stringify_values(obj: Any) -> Any:
    """Recursively convert every leaf value in a nested structure to a string.

    Args:
        obj: Any Python object.  Dicts and lists are traversed recursively;
            all other values are converted with ``str()``.

    Returns:
        A new structure of the same shape with all leaves stringified.
    """
    if isinstance(obj, dict):
        return {k: stringify_values(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [stringify_values(v) for v in obj]
    return str(obj)


def load_metadata_file(path: str | Path | None) -> dict[str, Any]:
    """Load a user-supplied JSON/YAML metadata file, rejecting reserved keys.

    Args:
        path: Path to a ``.json``, ``.yaml`` or ``.yml`` file.  ``None``
            returns an empty dict.

    Returns:
        The parsed mapping with every leaf stringified.

    Raises:
        click.BadParameter: If the file is not a JSON/YAML object, or contains
            a key of :data:`RESERVED_METADATA_KEYS`.
    """
    if path is None:
        return {}
    path = Path(path)
    text = path.read_text()
    if path.suffix in (".yaml", ".yml"):
        import yaml

        meta = yaml.safe_load(text)
    else:
        meta = json.loads(text)

    if not isinstance(meta, dict):
        raise click.BadParameter(
            f"Metadata file {str(path)!r} must contain a JSON/YAML object."
        )
    # 'qc' used to be overwritten by create's contact-quality block without a
    # word, and qc --store crashed on a non-object value.
    for key, why in RESERVED_METADATA_KEYS.items():
        if key in meta:
            raise click.BadParameter(
                f"Metadata file must not contain a {key!r} key — that key is {why}."
            )
    return stringify_values(meta)


def refuse_reserved_metadata_key(key: str, action: str) -> None:
    """Refuse to set or remove a key of :data:`RESERVED_METADATA_KEYS`.

    ``metadata add -k commands --value '[]'`` used to wipe the provenance
    record (and ``metadata remove`` to delete it), although a ``--metadata``
    file was already refused the same key.

    Args:
        key: The key the user asked for.
        action: What was asked, for the message: ``"set"`` or ``"remove"``.

    Raises:
        click.BadParameter: If *key* is reserved.
    """
    why = RESERVED_METADATA_KEYS.get(key)
    if why is not None:
        raise click.BadParameter(
            f"Cannot {action} {key!r}: that key is {why}. Read it with "
            f"'metacooler metadata get -k {key}'.",
            param_hint="'-k' / '--key'",
        )
