"""Write the report of a map, and read its report.txt back.

The report of one map (CONVENTIONS.md § 12), named after its sample (the
map's name, or its entry's, as :func:`metacooler.uri.sample_name` gives it),
in an output directory:

- ``<sample>.report.txt``: every number of the report, for scripts, which
  :func:`read_report` reads back;
- ``<sample>.report_mqc.yaml``: the map's row of MultiQC's General
  Statistics table;
- ``<sample>.report.html``: a page to read the map with.

:func:`metacooler.ops.report.collect` gathers it; chromokit's report kit
(``metacooler._report``) writes it. ``docs/format.md`` documents every key.
"""

from __future__ import annotations

import logging
from pathlib import Path

from metacooler.schema import DEFAULT_CHUNKSIZE

logger = logging.getLogger("metacooler")


def report_paths(outdir: str | Path, sample: str) -> list[Path]:
    """The files of the report of *sample* in *outdir*: its report.txt, its
    page and its report_mqc.yaml, in that order."""
    from metacooler._report import report_files

    return [Path(outdir) / name for name in report_files(sample)]


def write_report(
    uri: str | Path,
    outdir: str | Path = ".",
    *,
    force: bool = False,
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> list[Path]:
    """Write the report of a map in *outdir*, which is made if missing.

    Args:
        uri: ``file.metacool`` path or ``file.mmetacool::entries/name`` URI.
            It is read only.
        outdir: The report's directory.
        force: Replace the report of the same sample in *outdir*; without it,
            an existing one is refused before the map is read.
        chunksize: Pixels read per iteration, where the QC block or the P(s)
            is computed.

    Returns:
        The report's files: its report.txt, its page and its report_mqc.yaml.

    Raises:
        FileExistsError: If a file of the report exists and *force* is not
            set.
        ValueError: If one of them is the map itself.
    """
    from metacooler._report import write_report as write_files
    from metacooler.io._utils import check_output_path
    from metacooler.ops.report import collect, tool
    from metacooler.uri import sample_name

    paths = report_paths(outdir, sample_name(uri))
    for path in paths:
        check_output_path(path, force, inputs=[uri], defer_removal=True)
    collected = collect(uri, chunksize=chunksize)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    written = write_files(
        collected.report,
        outdir,
        tool(),
        figures=collected.figures,
        inputs=collected.inputs,
    )
    logger.info("The report of %s: %s", uri, ", ".join(written))
    return [Path(path) for path in written]


def read_report(path: str | Path) -> tuple[dict, dict]:
    """Read a report.txt back: its sections, by name, each its keys and their
    values, as text; and its tables, by name, as ``pandas.DataFrame`` of text.
    ``NA`` reads as None.

    Args:
        path: A ``<sample>.report.txt``.

    Returns:
        ``(sections, tables)``.
    """
    from metacooler._report import read_txt

    return read_txt(path)
