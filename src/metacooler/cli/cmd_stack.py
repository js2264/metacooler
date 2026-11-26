"""metacooler stack command."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
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
        "Defaults to the entry name for an entry URI, otherwise the file stem "
        "of each input (e.g. [cyan]sample.metacool[/cyan] → [cyan]sample[/cyan])."
    ),
)
@force_option
@verbose_option
@quiet_option
def cmd_stack(
    input_files: tuple[str, ...],
    output: Path,
    names: str | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Combine multiple .metacool files into a single .mmetacool file.

    \b
    INPUT_FILES    One or more input .metacool paths or entry URIs
    OUTPUT         Output .mmetacool file path

    An input may address a single entry of an existing collection
    (``file.mmetacool::entries/name``); only that entry is copied in.
    """
    import logging

    from metacooler.cli.logging_utils import run_log, setup_logging
    from metacooler.io.mmetacool import copy_entry_cross_file, create_mmetacool
    from metacooler.uri import parse_uri

    setup_logging(verbose, command_name="metacooler stack", quiet=quiet)
    logger = logging.getLogger("metacooler")

    # Validate: at least 2 inputs
    if len(input_files) < 2:
        raise click.UsageError("At least 2 input files required.")

    check_output(output, force, inputs=input_files, defer_removal=True)

    # (file path, group path) for every input. A group path other than "/" used
    # to be discarded, silently copying the whole collection in as one entry.
    parsed = [parse_uri(uri) for uri in input_files]

    def _source_entry(group_path: str) -> str | None:
        """Entry name addressed by a URI group path, or None for a whole file."""
        if group_path == "/":
            return None
        return group_path.rstrip("/").rsplit("/", 1)[-1]

    # Determine entry names
    if names:
        entry_names = [n.strip() for n in names.split(",")]
        if len(entry_names) != len(input_files):
            raise click.UsageError(
                f"Number of --names ({len(entry_names)}) must match number of "
                f"inputs ({len(input_files)})."
            )
        if any(not n for n in entry_names):
            raise click.UsageError("--names contains an empty entry name.")
        duplicates = sorted({n for n in entry_names if entry_names.count(n) > 1})
        if duplicates:
            raise click.UsageError(
                f"--names contains duplicate entry names: {duplicates}. "
                "Entry names must be unique."
            )
    else:
        # Auto-generate names from the entry name, or the filename stem
        entry_names = []
        seen: set[str] = set()
        for uri, (file_path, group_path) in zip(input_files, parsed, strict=True):
            basename = _source_entry(group_path) or file_path.stem

            # Check for collisions
            if basename in seen:
                raise click.ClickException(
                    f"Entry name collision: {basename!r}. "
                    "Use --names to specify unique names."
                )
            seen.add(basename)
            entry_names.append(basename)
            logger.debug(f"Auto-named entry: {uri} -> {basename}")

    # Whole-file inputs are copied wholesale; entry URIs are pulled in one by
    # one afterwards so only the addressed group lands in the output.
    whole_files: dict[str, Path] = {}
    sub_entries: list[tuple[Path, str, str]] = []
    for (file_path, group_path), name in zip(parsed, entry_names, strict=True):
        src_entry = _source_entry(group_path)
        if src_entry is None:
            whole_files[name] = file_path
            logger.debug(f"Entry {name!r}: {file_path}")
        else:
            sub_entries.append((file_path, src_entry, name))
            logger.debug(f"Entry {name!r}: {file_path}::{group_path.lstrip('/')}")

    with run_log(output, "Collection written"):
        from metacooler.io.create import _staged_output
        from metacooler.io.staging import _writing_directly
        from metacooler.schema import validate_mmetacooler

        # Built in one staged file, which only this run knows of: each entry
        # pulled in is written to it directly, and it is validated once, then
        # renamed onto the output.
        with (
            _staged_output(
                Path(output), keep_partial=False, inputs=input_files
            ) as staged,
            _writing_directly(staged),
        ):
            create_mmetacool(staged, whole_files, force=True)
            for src_path, src_entry, name in sub_entries:
                copy_entry_cross_file(src_path, src_entry, staged, name, overwrite=True)
            validate_mmetacooler(staged)

        n_entries = len(whole_files) + len(sub_entries)
        click.echo(f"Created {output} with {n_entries} entries")
