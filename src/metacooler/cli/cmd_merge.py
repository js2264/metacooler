"""metacooler merge command."""

from __future__ import annotations

import json
from pathlib import Path

import rich_click as click

from metacooler.cli._options import (
    check_output,
    force_option,
    quiet_option,
    verbose_option,
)
from metacooler.cli._provenance import (
    read_metadata,
    read_metadata_uri,
    reconstruct_command,
)
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="merge", help="Merge metacool files with identical axes.")
@click.argument("metacool_paths", nargs=-1, type=MetaCoolURI(require_entry=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output .metacool file path.",
)
@force_option
@verbose_option
@quiet_option
def cmd_merge(
    metacool_paths: tuple[str, ...],
    output: Path,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Merge multiple ``.metacool`` files with identical axes.

    All input files must have the same contigs, genomes, and
    ``contig-ordering`` attribute.  Pixel tables are concatenated and
    re-sorted externally, and the inputs' QC blocks are summed when they were
    counted alike.  Metadata from all inputs is merged and the command line
    is appended for provenance.
    """
    import h5py

    from metacooler.cli.logging_utils import setup_logging
    from metacooler.io.merge import merge
    from metacooler.qc import QC_METADATA_KEY, strip_qc_metadata

    setup_logging(verbose, command_name="metacooler merge", quiet=quiet)

    if len(metacool_paths) < 2:
        raise click.UsageError("At least two input files required.")

    check_output(output, force, inputs=metacool_paths)

    merge(inputs=metacool_paths, output=output)

    # Merge metadata from inputs and append command
    combined_meta: dict = {}
    all_commands: list[str] = []
    for inp in metacool_paths:
        # Each input's qc block describes that input alone, not the summed
        # pixel table: the block merge() wrote, if it could sum them, is kept.
        meta = strip_qc_metadata(read_metadata_uri(inp))
        for key, value in meta.items():
            if key == "commands":
                if isinstance(value, list):
                    for cmd in value:
                        if cmd not in all_commands:
                            all_commands.append(cmd)
            else:
                combined_meta[key] = value
    all_commands.append(reconstruct_command())
    combined_meta["commands"] = all_commands

    with h5py.File(output, "a") as f:
        merged_qc = read_metadata(f).get(QC_METADATA_KEY)
        if merged_qc is not None:
            combined_meta[QC_METADATA_KEY] = merged_qc
        f.attrs["metadata"] = json.dumps(combined_meta)
