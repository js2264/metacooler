"""metacooler normalize — compute normalization weights and store them in /bins.

One subcommand per family of methods.  Each stores a weight column in the map
itself, as ``cooler balance`` does: the pixels are never rewritten, and the
column is read back with ``show``/``dump``/``network --balance NAME`` -- or by
cooler, cooltools and HiGlass, which apply ``/bins/NAME`` the same way.
"""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from metacooler.cli._options import parse_comma_list, quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli
from metacooler.schema import DEFAULT_CHUNKSIZE

#: Defaults shown in ``--help``; the API's own live in metacooler.ops.normalize,
#: which is not imported here so that ``--help`` stays fast.
_MIN_LENGTH = "1000"
_MIN_SIGNAL = 2.0
_SPURIOUS = 0.05
_MAX_PAIRS = 1_000_000


def _common_options(default_name: str | None):
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
            show_default=default_name is not None,
            help=(
                "Column of [cyan]/bins[/cyan] to store the weights as."
                if default_name is not None
                else "Column of [cyan]/bins[/cyan] to store the weights as "
                "[default: depends on the method]."
            ),
        )(f)
        return f

    return decorate


def _sites_options(f):
    """``--sites`` / ``--enzyme`` / ``--fasta``: where restriction sites come from."""
    f = click.option(
        "--fasta",
        type=click.Path(exists=True, dir_okay=False, path_type=Path),
        default=None,
        help="Assembly FASTA (plain or [cyan].gz[/cyan]) to count the sites of "
        "[cyan]--enzyme[/cyan] in.",
    )(f)
    f = click.option(
        "--enzyme",
        "enzymes",
        type=str,
        multiple=True,
        callback=parse_comma_list,
        help="Restriction enzyme(s) of the library (comma-separated or repeated): "
        "names such as [cyan]DpnII[/cyan], [cyan]HinfI[/cyan], [cyan]Arima[/cyan], "
        "or IUPAC sites such as [cyan]GATC[/cyan]. Needs [cyan]--fasta[/cyan].",
    )(f)
    f = click.option(
        "--sites",
        type=str,
        default=None,
        metavar="COLUMN",
        help="[cyan]/contigs[/cyan] column holding each contig's restriction-site "
        "count, used as is.",
    )(f)
    return f


def _filter_options(f):
    """The contig filters of MetaCC and HiCBin, and the spurious-contact quantile."""
    f = click.option(
        "--spurious",
        type=click.FloatRange(min=0, max=1, max_open=True),
        default=_SPURIOUS,
        show_default=True,
        help="Quantile of the normalized contig-pair contacts stored as the "
        "spurious-contact cutoff (see [cyan]network --drop-spurious[/cyan]).",
    )(f)
    f = click.option(
        "--min-signal",
        type=click.FloatRange(min=0),
        default=_MIN_SIGNAL,
        show_default=True,
        help="Leave out contigs whose largest contact with another contig is smaller.",
    )(f)
    f = click.option(
        "--min-length",
        type=str,
        default=_MIN_LENGTH,
        show_default=True,
        help="Leave out shorter contigs. Accepts [cyan]bp/kb/Mb[/cyan] suffixes.",
    )(f)
    return f


def _distance(value: str, option: str) -> int:
    """Parse a distance option, reporting a bad value against the option."""
    from metacooler.filters import parse_distance

    try:
        return parse_distance(value)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint=option) from exc


def _metator_mode(ctx, param, value: str) -> str:
    """Accept MetaTOR's spellings (``theoritical_hit``, ``RS``) of a mode."""
    from metacooler.normalize import canonical_metator_mode

    try:
        return canonical_metator_mode(value)
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc


def _start(verbose: bool, quiet: bool, command: str) -> list[str]:
    """Set up logging and return the command line to record as provenance."""
    from metacooler.cli._provenance import reconstruct_command
    from metacooler.cli.logging_utils import setup_logging

    setup_logging(verbose, command_name=f"metacooler normalize {command}", quiet=quiet)
    return [reconstruct_command()]


@cli.group(
    name="normalize",
    help="""
Compute normalization weights and store them in [cyan]/bins[/cyan].

Every method stores one weight per bin, as [cyan]cooler balance[/cyan] does: the
normalized contact of two bins is the raw count times both weights, NaN for a
masked bin. The pixels are not rewritten. Read the weights back with
[cyan]show[/cyan], [cyan]dump[/cyan] or [cyan]network --balance NAME[/cyan], or
with cooler, cooltools and HiGlass. [cyan]normalise[/cyan] is an alias.
""",
)
def normalize() -> None:
    """Group for the normalization methods."""


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


@normalize.command(
    name="metator",
    help="""
MetaTOR's network normalizations, with MetaTOR's formulas and constants.

[yellow]empirical_hit[/yellow] (MetaTOR's default) divides the contacts of two
contigs by the geometric mean of their Hi-C hits; [yellow]abundance[/yellow] by
that of their shotgun coverage ([cyan]--coverage[/cyan]);
[yellow]length[/yellow] of their hits per bp; [yellow]rs[/yellow] of their hits
per restriction site; [yellow]theoretical_hit[/yellow] of 10 x coverage x
sqrt(length x sites). Sites come from [cyan]--sites[/cyan] (MetaTOR's
[cyan]Restriction_site[/cyan] column) or are counted from [cyan]--enzyme[/cyan]
in [cyan]--fasta[/cyan] (plus one, as MetaTOR counts them).
""",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--mode",
    type=str,
    default="empirical_hit",
    show_default=True,
    callback=_metator_mode,
    help="[yellow]empirical_hit[/yellow], [yellow]abundance[/yellow], "
    "[yellow]length[/yellow], [yellow]rs[/yellow] or "
    "[yellow]theoretical_hit[/yellow] (MetaTOR's spellings are accepted).",
)
@click.option(
    "--coverage",
    type=str,
    default=None,
    metavar="COLUMN",
    help="[cyan]/contigs[/cyan] column holding each contig's shotgun coverage "
    "(MetaTOR's [cyan]Shotgun_coverage[/cyan]). For [yellow]abundance[/yellow] "
    "and [yellow]theoretical_hit[/yellow].",
)
@_sites_options
@_common_options(None)
def normalize_metator_cmd(
    metacool_path: str,
    mode: str,
    coverage: str | None,
    sites: str | None,
    enzymes: list[str],
    fasta: Path | None,
    name: str | None,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Normalize a map as MetaTOR normalizes its contig network."""
    commands = _start(verbose, quiet, "metator")
    from metacooler.ops.normalize import normalize_metator

    normalize_metator(
        metacool_path,
        mode=mode,
        coverage=coverage,
        sites=sites,
        enzymes=enzymes or None,
        fasta=fasta,
        name=name,
        chunksize=chunksize,
        force=force,
        commands=commands,
    )


@normalize.command(
    name="mge",
    help="""
The MGE score of MetaTOR's mge module (Bignaud et al., Nat. Microbiol. 2025).

The contacts of two contigs over the geometric mean of their own intra-contig
contacts, pairs closer than [cyan]--min-distance[/cyan] left out:
x_ij / sqrt(x_ii x_jj). Exact when [cyan]create pairs[/cyan] measured the QC
block at that threshold ([cyan]--qc-threshold[/cyan]); otherwise measured
between bin midpoints. [yellow]genome[/yellow] takes the contacts inside each
genome instead.
""",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--min-distance",
    type=str,
    default="1000",
    show_default=True,
    help="Intra-contig pairs closer than this are left out. Accepts "
    "[cyan]bp/kb/Mb[/cyan] suffixes.",
)
@click.option(
    "--scope",
    type=click.Choice(["contig", "genome"]),
    default="contig",
    show_default=True,
    help="Whose intra contacts divide: each contig's or each genome's.",
)
@_common_options(None)
def normalize_mge_cmd(
    metacool_path: str,
    min_distance: str,
    scope: str,
    name: str | None,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Normalize a map by the intra contacts of each contig."""
    commands = _start(verbose, quiet, "mge")
    from metacooler.ops.normalize import normalize_mge

    normalize_mge(
        metacool_path,
        min_distance=_distance(min_distance, "--min-distance"),
        scope=scope,
        name=name,
        chunksize=chunksize,
        force=force,
        commands=commands,
    )


@normalize.command(
    name="normcc",
    help="""
NormCC, the normalization of MetaCC (Du & Sun, Nat. Commun. 2023).

A negative binomial regression of each contig's contacts with the other contigs
on its restriction sites, length and intra-contig contacts (which stand in for
coverage). Without [cyan]--sites[/cyan] or [cyan]--enzyme[/cyan], the site term
is left out, as in MetaCC.
""",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@_sites_options
@_filter_options
@_common_options("normcc")
def normalize_normcc_cmd(
    metacool_path: str,
    sites: str | None,
    enzymes: list[str],
    fasta: Path | None,
    min_length: str,
    min_signal: float,
    spurious: float,
    name: str,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Normalize a map with NormCC."""
    commands = _start(verbose, quiet, "normcc")
    from metacooler.ops.normalize import normalize_normcc

    normalize_normcc(
        metacool_path,
        sites=sites,
        enzymes=enzymes or None,
        fasta=fasta,
        min_length=_distance(min_length, "--min-length"),
        min_signal=min_signal,
        spurious=spurious,
        name=name,
        chunksize=chunksize,
        force=force,
        commands=commands,
    )


@normalize.command(
    name="hiczin",
    help="""
HiCzin, the normalization of HiCBin (Du et al., J. Comput. Biol. 2022).

A zero-inflated negative binomial regression of the contacts between contigs of
one genome on their restriction sites, lengths and coverages. HiCzin learns
those pairs from a taxonomic labelling; here they come from the genome
assignment. Without [cyan]--sites[/cyan] or [cyan]--enzyme[/cyan], the site term
is left out (HiCzin_LC).
""",
)
@click.argument("metacool_path", type=MetaCoolURI(require_entry=True))
@click.option(
    "--coverage",
    type=str,
    required=True,
    metavar="COLUMN",
    help="[cyan]/contigs[/cyan] column holding each contig's shotgun coverage.",
)
@_sites_options
@_filter_options
@click.option(
    "--max-pairs",
    type=click.IntRange(min=0),
    default=_MAX_PAIRS,
    show_default=True,
    help="Fit on a uniform sample of at most this many pairs of contigs "
    "(0: all of them).",
)
@click.option(
    "--seed",
    type=int,
    default=0,
    show_default=True,
    help="Seed of that sample.",
)
@_common_options("hiczin")
def normalize_hiczin_cmd(
    metacool_path: str,
    coverage: str,
    sites: str | None,
    enzymes: list[str],
    fasta: Path | None,
    min_length: str,
    min_signal: float,
    spurious: float,
    max_pairs: int,
    seed: int,
    name: str,
    chunksize: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Normalize a map with HiCzin."""
    commands = _start(verbose, quiet, "hiczin")
    from metacooler.ops.normalize import normalize_hiczin

    normalize_hiczin(
        metacool_path,
        coverage=coverage,
        sites=sites,
        enzymes=enzymes or None,
        fasta=fasta,
        min_length=_distance(min_length, "--min-length"),
        min_signal=min_signal,
        spurious=spurious,
        max_pairs=max_pairs,
        seed=seed,
        name=name,
        chunksize=chunksize,
        force=force,
        commands=commands,
    )


#: Attributes shown by ``normalize list``, in order, when a column has them.
_LIST_DETAILS = ("scope", "mode", "converged", "spurious_cutoff")


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
            "'metacooler normalize'."
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
