"""metacooler stack command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="stack", help="Combine metacool files into a single mmetacool file.")
@click.argument(
    "input_files",
    nargs=-1,
    type=MetaCoolURI(),
    required=True,
    metavar="INPUT_FILES...",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .mmetacool file path.",
)
@click.option(
    "--names",
    type=str,
    default=None,
    help=(
        "Comma-separated entry names for the output collection "
        "(must match number of inputs). "
        "Defaults to the file stem of each input (e.g. [cyan]sample.metacool[/cyan] → [cyan]sample[/cyan])."
    ),
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Overwrite output file if it exists.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase verbosity (debug logging).",
)
def cmd_stack(
    input_files: tuple[str, ...],
    output: Path,
    names: str | None,
    force: bool,
    verbose: bool,
) -> None:
    """Combine multiple .metacool files into a single .mmetacool file.

    \b
    INPUT_FILES    One or more input .metacool file paths
    OUTPUT         Output .mmetacool file path
    """
    import logging

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.mmetacool import create_mmetacool
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler stack")
    logger = logging.getLogger("metacooler")

    # Validate: at least 2 inputs
    if len(input_files) < 2:
        raise click.UsageError("At least 2 input files required.")

    # Validate: output doesn't exist (unless --force)
    if output.exists() and not force:
        raise click.ClickException(
            f"Output file already exists: {output}. Use --force to overwrite."
        )

    # Determine entry names
    if names:
        entry_names = [n.strip() for n in names.split(",")]
        if len(entry_names) != len(input_files):
            raise click.UsageError(
                f"Number of --names ({len(entry_names)}) must match number of inputs ({len(input_files)})."
            )
    else:
        # Auto-generate names from filenames
        entry_names = []
        seen = set()
        for inp_uri in input_files:
            # Parse URI to get file path
            file_path, _ = parse_uri(inp_uri)
            # Get stem (filename without any extensions)
            basename = file_path.stem

            # Check for collisions
            if basename in seen:
                raise click.ClickException(
                    f"Entry name collision: {basename!r}. Use --names to specify unique names."
                )
            seen.add(basename)
            entry_names.append(basename)
            logger.debug(f"Auto-named entry: {inp_uri} -> {basename}")

    # Convert input URIs to file paths and build entries dict
    entries = {}
    for inp_uri, name in zip(input_files, entry_names, strict=True):
        file_path, _ = parse_uri(inp_uri)
        entries[name] = file_path
        logger.debug(f"Entry {name!r}: {file_path}")

    # Create mmetacool
    try:
        create_mmetacool(output, entries, force=force)
        click.echo(f"Created {output} with {len(entries)} entries")
    except FileNotFoundError as e:
        raise click.ClickException(str(e)) from e
    except Exception as e:
        raise click.ClickException(f"Failed to create stack: {e}") from e
