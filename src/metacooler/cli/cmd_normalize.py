"""metacooler normalize — balance a map and list the weights it carries.

``normalize ice`` balances a map, by default genome by genome, and stores the
weights as a column of ``/bins``, as ``cooler balance`` does: the pixels are
never rewritten, and the column is read back with ``show``/``dump``/``network
--balance NAME`` -- or by cooler, cooltools and HiGlass, which apply
``/bins/NAME`` the same way.  ``normalize list`` lists every weight column of a
map, including those another tool stored with
:func:`metacooler.io.weights.store_weights`.
"""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.schema import DEFAULT_CHUNKSIZE


def _common_options(default_name: str):
    """Options every subcommand shares: the column, chunking, force, logging."""

    def decorate(f):
        f = quiet_option(f)
        f = verbose_option(f)
        f = click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Replace an existing [cyan]/bins[/cyan] column of that name.",
        )(f)
        f = click.option(
            "--chunksize",
            type=click.IntRange(min=1),
            default=DEFAULT_CHUNKSIZE,
            show_default=True,
            help="Pixels read per chunk.",
        )(f)
        f = click.option(
            "--name",
            type=str,
            default=default_name,
            show_default=True,
            help="Column of [cyan]/bins[/cyan] to store the weights as.",
        )(f)
        return f

    return decorate


def _start(verbose: bool, quiet: bool, command: str) -> list[str]:
    """Set up logging and return the command line to record as provenance."""
    from metacooler.cli._provenance import reconstruct_command
    from metacooler.cli.logging_utils import setup_logging

    setup_logging(verbose, command_name=f"metacooler normalize {command}", quiet=quiet)
    return [reconstruct_command()]


@cli.group(
    name="normalize",
    help="""
Balance a map and store the weights in [cyan]/bins[/cyan].

A weight column holds one weight per bin, as [cyan]cooler balance[/cyan] writes
it: the normalized contact of two bins is the raw count times both weights, NaN
for a masked bin. The pixels are not rewritten. Read the weights back with
[cyan]show[/cyan], [cyan]dump[/cyan] or [cyan]network --balance NAME[/cyan], or
with cooler, cooltools and HiGlass. Another tool's normalization is stored the
same way with [cyan]metacooler.store_weights[/cyan] in Python.
[cyan]normalise[/cyan] is an alias.
""",
)
def normalize() -> None:
    """Group for ``normalize ice`` and ``normalize list``."""


@normalize.command(
    name="ice",
    help="""
Iterative correction (matrix balancing), by default genome by genome.

[yellow]genome[/yellow] balances each MAG's block of the matrix on its own and
leaves the contacts between genomes out of the fit; each [cyan]__unbinned__[/cyan]
contig is a genome of its own. [yellow]contig[/yellow] and
[yellow]global[/yellow] are cooler's [cyan]--cis-only[/cyan] and whole-matrix
balancing.
""",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--scope",
    type=click.Choice(["genome", "contig", "global"]),
    default="genome",
    show_default=True,
    help="What is balanced on its own.",
)
@click.option(
    "--ignore-diags",
    type=click.IntRange(min=0),
    default=2,
    show_default=True,
    help="Diagonals of each contig left out of the fit.",
)
@click.option(
    "--mad-max",
    type=click.IntRange(min=0),
    default=5,
    show_default=True,
    help="Mask bins whose marginal is this many median absolute deviations "
    "below the median of their genome (0 disables).",
)
@click.option(
    "--min-nnz",
    type=click.IntRange(min=0),
    default=10,
    show_default=True,
    help="Mask bins with fewer non-zero pixels.",
)
@click.option(
    "--min-count",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="Mask bins with a smaller total count.",
)
@click.option(
    "--tol",
    type=click.FloatRange(min=0, min_open=True),
    default=1e-5,
    show_default=True,
    help="Variance of the marginals at which balancing has converged.",
)
@click.option(
    "--max-iters",
    type=click.IntRange(min=0),
    default=200,
    show_default=True,
    help="Iteration limit.",
)
@_common_options("weight")
def normalize_ice_cmd(
    metacool_path: str,
    scope: str,
    ignore_diags: int,
    mad_max: int,
    min_nnz: int,
    min_count: int,
    tol: float,
    max_iters: int,
    name: str,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Balance a map and store the weights."""
    commands = _start(verbose, quiet, "ice")
    from metacooler.ops.normalize import normalize_ice

    normalize_ice(
        metacool_path,
        scope=scope,
        name=name,
        ignore_diags=ignore_diags,
        mad_max=mad_max,
        min_nnz=min_nnz,
        min_count=min_count,
        tol=tol,
        max_iters=max_iters,
        chunksize=chunksize,
        force=force,
        commands=commands,
    )


#: Attributes shown by ``normalize list``, in order, when a column has them.
_LIST_DETAILS = ("scope", "converged")


@normalize.command(
    name="list",
    help="List the weight columns of a map and how they were computed.",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@verbose_option
@quiet_option
def normalize_list_cmd(metacool_path: str, verbose: bool, quiet: bool) -> None:
    """Print one row per weight column of ``/bins``."""
    from metacooler.cli.logging_utils import setup_logging

    setup_logging(verbose, command_name="metacooler normalize list", quiet=quiet)
    from rich.console import Console
    from rich.table import Table

    from metacooler.api import MetaCooler
    from metacooler.io.weights import describe_weights

    rows = describe_weights(MetaCooler(metacool_path))
    if not rows:
        click.echo(
            f"No weight columns in {metacool_path}; compute one with "
            "'metacooler normalize ice'."
        )
        return
    table = Table(title=f"Weights: {Path(metacool_path.split('::')[0]).name}")
    for column in ("Name", "Method", "Details", "Masked bins", "Depends on"):
        # Fold rather than truncate: a listing that cuts names is no listing.
        table.add_column(column, overflow="fold")
    for row in rows:
        attrs = row["attrs"]
        assert isinstance(attrs, dict)
        details = ", ".join(
            f"{key}={attrs[key]}" for key in _LIST_DETAILS if key in attrs
        )
        depends = attrs.get("depends_on")
        table.add_row(
            str(row["name"]),
            str(row["method"] or "-"),
            details or "-",
            f"{row['n_masked']:,}",
            ", ".join(depends)
            if isinstance(depends, list)
            else str(depends or "pixels"),
        )
    Console().print(table)


# British spelling: the same group under a second, hidden name.  The two share
# one ``commands`` dict, so every subcommand above answers to both.
cli.add_command(
    click.RichGroup(
        name="normalise",
        commands=normalize.commands,
        help=normalize.help,
        hidden=True,
    )
)
