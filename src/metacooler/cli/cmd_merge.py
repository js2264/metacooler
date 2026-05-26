"""metacooler merge command."""

from __future__ import annotations

import json
from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


def _reconstruct_command() -> str:
    """Reconstruct the exact command line invocation from click context."""
    ctx = click.get_current_context()
    assert ctx is not None
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
    return " ".join(parts)


def _read_metadata_dict(uri: str) -> dict:
    """Read and parse the 'metadata' attribute from a metacool file."""
    import h5py

    from metacooler.uri import parse_uri

    file_path, group_path = parse_uri(uri)
    with h5py.File(file_path, "r") as f:
        grp = f[group_path] if group_path != "/" else f
        raw = grp.attrs.get("metadata", "{}")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        return json.loads(raw) if isinstance(raw, str) else dict(raw)


@cli.command(name="merge", help="Merge metacool files with identical axes.")
@click.argument("metacool_paths", nargs=-1, type=MetaCoolURI())
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .metacool file path.",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite output file if it already exists.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_merge(
    metacool_paths: tuple[str, ...],
    output: Path,
    force: bool,
    verbose: bool,
) -> None:
    """Merge multiple metacooler files."""
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.merge import merge

    setup_logging(verbose, command_name="metacooler merge")

    if len(metacool_paths) < 2:
        raise click.UsageError("At least two input files required.")

    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )
    if force and output.exists():
        output.unlink()

    merge(inputs=metacool_paths, output=output)

    # Merge metadata from inputs and append command
    combined_meta: dict = {}
    all_commands: list[str] = []
    for inp in metacool_paths:
        meta = _read_metadata_dict(inp)
        for key, value in meta.items():
            if key == "commands":
                if isinstance(value, list):
                    for cmd in value:
                        if cmd not in all_commands:
                            all_commands.append(cmd)
            else:
                combined_meta[key] = value
    all_commands.append(_reconstruct_command())
    combined_meta["commands"] = all_commands

    with h5py.File(output, "a") as f:
        f.attrs["metadata"] = json.dumps(combined_meta)
