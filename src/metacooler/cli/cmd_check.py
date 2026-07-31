"""metacooler check command."""

from __future__ import annotations

import rich_click as click

from metacooler.cli._options import quiet_option, verbose_option
from metacooler.cli.main import MetaCoolURI, cli


@cli.command(name="check", help="Check the integrity of a metacool file.")
@click.argument("metacool_path", type=MetaCoolURI())
@click.option(
    "--deep",
    is_flag=True,
    default=False,
    help=(
        "Also scan the [b]whole[/b] pixel table instead of the three sampled "
        "windows the fast check uses, and reconcile [cyan]attrs['sum'][/cyan] "
        "with the stored counts. O(nnz) — slow on large files."
    ),
)
@verbose_option
@quiet_option
def cmd_check(metacool_path: str, deep: bool, verbose: bool, quiet: bool) -> None:
    """Report on the integrity of a metacool file or mmetacool collection.

    Prints one line per check and exits non-zero if any of them failed.  A bare
    ``.mmetacool`` path is accepted and checks every entry.
    """
    from metacooler.cli.logging_utils import setup_logging
    from metacooler.schema import check

    setup_logging(verbose, command_name="metacooler check", quiet=quiet)

    results = check(metacool_path, deep=deep)
    for result in results:
        mark, colour = ("✓", "green") if result.ok else ("✗", "red")
        click.secho(f"{mark} {result.name}", fg=colour, nl=False)
        click.echo(f": {result.detail}")

    failed = [r for r in results if not r.ok]
    if failed:
        raise click.ClickException(
            f"{metacool_path} failed {len(failed)} of {len(results)} checks: "
            + ", ".join(r.name for r in failed)
        )
    click.secho(
        f"{metacool_path}: OK ({len(results)} checks passed)", fg="green", bold=True
    )
