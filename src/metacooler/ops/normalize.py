"""Normalize a metacool file that already exists: compute and store bin weights.

Each function computes one normalization of the map and stores it as a weight
column of ``/bins`` (see :mod:`metacooler.io.weights`): the normalized value of
a pixel is ``raw * w[bin1] * w[bin2]``, the convention of cooler's balancing
weights, so cooler, cooltools and HiGlass read the result as they read
``cooler balance``'s.  The pixels are never rewritten.  The arithmetic lives in
:mod:`metacooler.normalize`; this module reads what it needs from the file.

- :func:`normalize_ice` -- iterative correction (matrix balancing), by default
  **per genome**: each MAG's block of the matrix balanced on its own, contacts
  between genomes left out of the fit.  ``scope="contig"`` and
  ``scope="global"`` are cooler's ``--cis-only`` and whole-matrix balancing.
- :func:`normalize_normcc` -- NormCC (MetaCC; Du & Sun, Nat. Commun. 2023).
- :func:`normalize_hiczin` -- HiCzin (Du et al., J. Comput. Biol. 2022, the
  normalization of HiCBin).

NormCC and HiCzin give every bin of a contig its contig's weight, so a
normalized map summed back to contigs (``network --level contig --balance``)
holds exactly their contig-level values.  They store their weights with
:func:`metacooler.io.weights.store_weights`, as any other tool does.
"""

from __future__ import annotations

import contextlib
import logging
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd

from metacooler.api import MetaCooler
from metacooler.io.weights import (
    CONTIG_COLUMN_PREFIX,
    DEPENDS_ON_ASSIGNMENT,
    DEPENDS_ON_PIXELS,
    SPURIOUS_CUTOFF_ATTR,
    SPURIOUS_RULE_ATTR,
    check_weight_name,
    store_weights,
    weight_attributes,
)
from metacooler.normalize import (
    compact_groups,
    fit_nb_glm,
    fit_zinb,
    hiczin_contig_weights,
    ice_iterate,
    ice_mad_drop,
    ice_marginals,
    pairs_within_groups,
)
from metacooler.qc import unbinned_genome_id
from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from scipy import sparse

logger = logging.getLogger("metacooler")

#: What :func:`normalize_ice` balances separately.
ICE_SCOPES = ("genome", "contig", "global")

#: Default name of an ICE weight column: cooler's, so tools that read
#: ``weight`` by default (``cooler show``, HiGlass, cooltools) pick it up.
DEFAULT_ICE_NAME = "weight"

#: Contig filters of MetaCC and HiCBin: contigs shorter than ``min_len`` are
#: left out, and so are those whose largest contact with another contig is
#: below ``min_signal``.
DEFAULT_MIN_LENGTH = 1_000
DEFAULT_MIN_SIGNAL = 2

#: Share of the normalized contig-pair contacts NormCC and HiCzin call
#: spurious (their ``--thres``).
DEFAULT_SPURIOUS_QUANTILE = 0.05

#: Contig pairs HiCzin's regression is fitted on, at most (a uniform sample when
#: the genomes hold more).
DEFAULT_MAX_PAIRS = 1_000_000

#: Filtered pixels genome-scoped balancing keeps in memory between iterations
#: (about 16 bytes each); a larger map is re-read at every iteration, as
#: cooler does.
_ICE_CACHE_PIXELS = 1 << 25


class Normalization(NamedTuple):
    """A computed normalization.

    Attributes:
        name: The ``/bins`` column it is (or would be) stored as.
        weights: One multiplicative weight per bin, ``NaN`` for a masked bin.
        attrs: The attributes stored with the column: method, parameters,
            fit statistics and dependencies.
    """

    name: str
    weights: np.ndarray
    attrs: dict[str, object]


# ----------------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------------


def _open_target(uri: str | Path, name: str, *, store: bool, force: bool) -> MetaCooler:
    """Open the map, refusing up front a column that could not be stored.

    A balancing run can take a while; finding out at the end that the name is
    taken, or not a valid column name, would waste it.
    """
    mc = MetaCooler(uri)
    if store:
        check_weight_name(name)
        if not force:
            with mc.open() as grp:
                if name in grp["bins"]:
                    raise ValueError(
                        f"/bins/{name} already exists in {mc.uri}. Pass --force to "
                        "replace it, or --name to store this normalization beside it."
                    )
    return mc


def _held_for_store(uri: str | Path, store: bool) -> contextlib.AbstractContextManager:
    """The map's lock when the weights are to be stored (see
    :func:`metacooler.lock`), else nothing."""
    if not store:
        return contextlib.nullcontext()
    from metacooler.io.staging import lock, map_file

    return lock(map_file(uri))


def _finish(
    mc: MetaCooler,
    name: str,
    weights: np.ndarray,
    *,
    method: str,
    depends_on: Sequence[str],
    parameters: Mapping[str, object],
    store: bool,
    force: bool,
    commands: Sequence[str],
) -> Normalization:
    """Log, store and return a computed normalization."""
    n_bins = len(weights)
    n_masked = int(np.count_nonzero(~np.isfinite(weights)))
    logger.info("%s: %s of %s bin(s) masked", method, f"{n_masked:,}", f"{n_bins:,}")
    if n_bins and n_masked == n_bins:
        logger.warning(
            "Every bin is masked: the map has no data this normalization can use."
        )
    if store:
        store_weights(
            mc.uri,
            name,
            weights,
            method=method,
            depends_on=depends_on,
            parameters=parameters,
            force=force,
            commands=commands,
            source=mc,
        )
        logger.info("Stored the weights as /bins/%s in %s", name, mc.uri)
    return Normalization(
        name, weights, weight_attributes(method, depends_on, parameters)
    )


def _contig_groups(mc: MetaCooler) -> np.ndarray:
    """Group of every contig: its genome, or itself when it is unbinned.

    ``__unbinned__`` is the absence of a genome, so -- as in the QC
    classification -- each of its contigs is a genome of its own.
    """
    genome = mc.contigs["genome_id"].to_numpy(dtype=np.int64)
    unbinned = unbinned_genome_id(mc.genomes["name"])
    key = genome.copy()
    if unbinned is not None:
        alone = genome == unbinned
        key[alone] = len(mc.genomes) + np.flatnonzero(alone)
    return compact_groups(key)


def _numeric_contig_column(
    contigs: pd.DataFrame, column: str, *, what: str
) -> np.ndarray:
    """Read a ``/contigs`` column as numbers, ``NaN`` where a value is unusable."""
    if column not in contigs.columns:
        listing = ", ".join(
            c for c in contigs.columns if c not in ("name", "genome_id")
        )
        raise ValueError(
            f"/contigs has no column {column!r} to read the {what} from. "
            f"Its columns: {listing}. Add one with 'metacooler annotate'."
        )
    values = np.array(pd.to_numeric(contigs[column], errors="coerce"), dtype=np.float64)
    bad = ~np.isfinite(values) | (values < 0)
    if bad.any():
        logger.warning(
            "%s contig(s) have no usable %s in /contigs/%s (missing, non-numeric "
            "or negative); their bins are masked.",
            f"{int(bad.sum()):,}",
            what,
            column,
        )
        values[bad] = np.nan
    return values


class _Sites(NamedTuple):
    """Restriction-site counts and where they came from."""

    counts: np.ndarray
    attrs: dict[str, object]
    depends: list[str]


def _contig_sites(
    mc: MetaCooler,
    *,
    sites: str | None,
    enzymes: Sequence[str] | None,
    fasta: str | Path | None,
) -> _Sites | None:
    """Restriction-site count of every contig, from a column or the FASTA.

    Args:
        mc: The open map.
        sites: ``/contigs`` column holding the counts, used as they are.
        enzymes: Enzymes whose sites are counted in *fasta*.
        fasta: The assembly.

    Returns:
        A :class:`_Sites`, or ``None`` when neither source is given.

    Raises:
        ValueError: If both sources are given, or a FASTA without enzymes or
            enzymes without a FASTA.
    """
    if sites is not None and enzymes:
        raise ValueError(
            "Give the restriction sites either as a /contigs column (--sites) or as "
            "enzymes to count in the assembly (--enzyme with --fasta), not both."
        )
    if fasta is not None and not enzymes:
        raise ValueError(
            "--fasta is only read to count the restriction sites of --enzyme; "
            "pass the enzyme too."
        )
    contigs = mc.contigs
    if sites is not None:
        values = _numeric_contig_column(contigs, sites, what="restriction-site count")
        return _Sites(
            values,
            {"sites": f"{CONTIG_COLUMN_PREFIX}{sites}"},
            [f"{CONTIG_COLUMN_PREFIX}{sites}"],
        )
    if enzymes:
        if fasta is None:
            raise ValueError(
                "--enzyme counts restriction sites in the assembly: pass the "
                "FASTA the map was built from with --fasta."
            )
        from metacooler.io.sites import count_sites

        counts = count_sites(
            fasta,
            contigs["name"].astype(str).tolist(),
            contigs["length"].to_numpy(),
            list(enzymes),
        ).astype(np.float64)
        return _Sites(counts, {"enzymes": list(enzymes), "fasta": str(fasta)}, [])
    return None


def _contig_matrix(mc: MetaCooler, chunksize: int) -> sparse.csr_matrix:
    """The symmetric contig-by-contig contact matrix, as ``float64``."""
    from metacooler.ops.aggregate import _aggregate_sparse

    return (
        _aggregate_sparse(mc, level="contig", chunksize=chunksize)
        .astype(np.float64)
        .tocsr()
    )


def _off_diagonal(matrix: sparse.csr_matrix) -> sparse.csr_matrix:
    """*matrix* with its diagonal (the intra-contig contacts) removed."""
    off = matrix.copy()
    off.setdiag(0)
    off.eliminate_zeros()
    return off


def _largest_contact(off: sparse.csr_matrix, among: np.ndarray) -> np.ndarray:
    """Each contig's largest contact with another contig of *among*.

    MetaCC and HiCBin drop the short contigs before they build the contact
    matrix, so a contig's signal is measured against the long ones only.
    """
    from scipy import sparse

    restricted = off @ sparse.diags(among.astype(np.float64))
    return np.asarray(restricted.max(axis=1).todense()).ravel()


def _spurious_cutoff(
    off: sparse.csr_matrix, weights: np.ndarray, keep: np.ndarray, quantile: float
) -> tuple[float, int]:
    """Quantile of the normalized contacts between two kept contigs.

    Returns:
        ``(cutoff, n_contacts)`` -- ``NaN`` when there is no such contact.
    """
    from scipy import sparse

    upper = sparse.triu(off, k=1).tocoo()
    sel = keep[upper.row] & keep[upper.col] & (upper.data > 0)
    values = upper.data[sel] * weights[upper.row[sel]] * weights[upper.col[sel]]
    if not values.size:
        return float("nan"), 0
    return float(np.quantile(values, quantile)), int(values.size)


def _check_quantile(value: float, option: str) -> None:
    if not 0.0 <= value < 1.0:
        raise ValueError(f"{option} must be a fraction in [0, 1), got {value}.")


# ----------------------------------------------------------------------------
# ICE
# ----------------------------------------------------------------------------


class _GroupPixels:
    """The pixels genome-scoped balancing sees, chunk by chunk.

    A pixel takes part when its two bins are in the same group and it is not on
    one of the first *ignore_diags* diagonals of a contig (cooler's
    ``ignore_diags``, which also stops at a contig boundary).  The filtered
    chunks are kept in memory after the first pass when they fit in
    :data:`_ICE_CACHE_PIXELS`, so the iterations do not re-read the file.
    """

    def __init__(
        self,
        mc: MetaCooler,
        group_of_bin: np.ndarray,
        contig_of_bin: np.ndarray,
        ignore_diags: int,
        chunksize: int,
    ) -> None:
        self._mc = mc
        self._group = group_of_bin
        self._contig = contig_of_bin
        self._ignore_diags = ignore_diags
        self._chunksize = chunksize
        self._n_bins = len(group_of_bin)
        self._cache: list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None = None

    def _chunks(self):
        if self._cache is not None:
            yield from self._cache
            return
        bin_dtype = np.int32 if self._n_bins < np.iinfo(np.int32).max else np.int64
        cache: list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None = []
        kept = 0
        for chunk in self._mc.iter_pixels(self._chunksize):
            b1 = chunk["bin1_id"].to_numpy(dtype=np.int64)
            b2 = chunk["bin2_id"].to_numpy(dtype=np.int64)
            counts = chunk["count"].to_numpy(dtype=np.float64)
            keep = self._group[b1] == self._group[b2]
            if self._ignore_diags:
                keep &= ~(
                    (np.abs(b2 - b1) < self._ignore_diags)
                    & (self._contig[b1] == self._contig[b2])
                )
            part = (
                b1[keep].astype(bin_dtype),
                b2[keep].astype(bin_dtype),
                counts[keep],
            )
            if cache is not None:
                kept += len(part[0])
                if kept > _ICE_CACHE_PIXELS:
                    cache = None
                else:
                    cache.append(part)
            yield part
        self._cache = cache

    def marginals(
        self, bias: np.ndarray | None = None, *, binarize: bool = False
    ) -> np.ndarray:
        """Marginals of the filtered matrix, optionally weighted by *bias*."""
        marg = np.zeros(self._n_bins)
        for b1, b2, counts in self._chunks():
            values = (counts != 0).astype(np.float64) if binarize else counts
            if bias is not None:
                values = bias[b1] * bias[b2] * values
            marg += ice_marginals(b1, b2, values, self._n_bins)
        return marg


def _ice_by_genome(
    mc: MetaCooler,
    *,
    ignore_diags: int,
    mad_max: int,
    min_nnz: int,
    min_count: int,
    tol: float,
    max_iters: int,
    chunksize: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """Balance each genome's block of the matrix on its own."""
    contig_of_bin = mc._bin_contig_ids
    group_of_bin = _contig_groups(mc)[contig_of_bin]
    n_bins = len(contig_of_bin)
    pixels = _GroupPixels(mc, group_of_bin, contig_of_bin, ignore_diags, chunksize)

    bias = np.ones(n_bins)
    if min_nnz > 0:
        bias[pixels.marginals(binarize=True) < min_nnz] = 0
    marg = pixels.marginals()
    if min_count:
        bias[marg < min_count] = 0
    if mad_max > 0 and n_bins:
        bias[ice_mad_drop(marg, group_of_bin, mad_max)] = 0

    stats = ice_iterate(
        bias,
        group_of_bin,
        lambda b: pixels.marginals(b),
        tol=tol,
        max_iters=max_iters,
    )
    n_groups = len(stats.scale)
    unconverged = int(np.count_nonzero(~stats.converged))
    if unconverged:
        logger.warning(
            "Iteration limit (%d) reached without convergence on %d of %d "
            "genome(s); raise --max-iters or --tol.",
            max_iters,
            unconverged,
            n_groups,
        )
    finite = np.isfinite(stats.var)
    attrs: dict[str, object] = {
        "cis_only": False,
        "trans_only": False,
        "converged": unconverged == 0,
        "var": float(np.max(stats.var[finite])) if finite.any() else float("nan"),
        "n_groups": n_groups,
        "n_unconverged": unconverged,
        "n_iters": stats.n_iters,
    }
    return bias, attrs


def _ice_with_cooler(
    mc: MetaCooler,
    *,
    cis_only: bool,
    ignore_diags: int,
    mad_max: int,
    min_nnz: int,
    min_count: int,
    tol: float,
    max_iters: int,
    chunksize: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """cooler's own balancing: per contig (*cis_only*) or of the whole matrix."""
    import cooler

    if not len(mc._clr.pixels()):
        return np.full(len(mc._bin_contig_ids), np.nan), {"converged": True}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        bias, stats = cooler.balance_cooler(
            mc._clr,
            cis_only=cis_only,
            ignore_diags=ignore_diags if ignore_diags else False,
            mad_max=mad_max,
            min_nnz=min_nnz,
            min_count=min_count,
            tol=tol,
            max_iters=max_iters,
            chunksize=chunksize,
            rescale_marginals=True,
            store=False,
        )
    for warning in caught:
        logger.warning("cooler: %s", warning.message)
    return np.asarray(bias, dtype=np.float64), _summarize_stats(stats)


def _summarize_stats(stats: dict[str, object]) -> dict[str, object]:
    """cooler's balancing statistics, with per-contig arrays summarised.

    Balanced per contig, cooler reports a scale, a variance and a convergence
    flag for every contig.  Two million values cannot be an HDF5 attribute
    (64 kB at most), so they are reduced to the file-wide figures the other
    scopes report.
    """
    out: dict[str, object] = {}
    for key, value in stats.items():
        if key == "divisive_weights":
            continue  # store_weights writes it, always False
        array = np.asarray(value)
        if array.ndim == 0:
            out[key] = value.item() if isinstance(value, np.generic) else value
            continue
        if key == "converged":
            out["converged"] = bool(array.all())
            out["n_unconverged"] = int(np.count_nonzero(~array.astype(bool)))
            out["n_groups"] = int(array.size)
        elif key == "var":
            finite = array[np.isfinite(array)]
            out["var"] = float(finite.max()) if finite.size else float("nan")
        # Per-contig scales mean nothing beside the weights; leave them out.
    return out


def normalize_ice(
    uri: str | Path,
    *,
    scope: str = "genome",
    name: str = DEFAULT_ICE_NAME,
    ignore_diags: int = 2,
    mad_max: int = 5,
    min_nnz: int = 10,
    min_count: int = 0,
    tol: float = 1e-5,
    max_iters: int = 200,
    chunksize: int = DEFAULT_CHUNKSIZE,
    store: bool = True,
    force: bool = False,
    commands: Sequence[str] = (),
) -> Normalization:
    """Balance a map by iterative correction (ICE) and store the weights.

    Matrix balancing gives every bin the same total of normalized contacts.
    Applied to a whole metagenome, that equalises genomes of very different
    abundance and inflates the stray contacts of rare ones; so by default each
    *genome* is balanced on its own, contacts between genomes left out of the
    fit (they are still weighted when the map is read).  An ``__unbinned__``
    contig is a genome of its own.  The bin filters are cooler's, except that
    the MAD filter compares a bin with the median of its genome, not of its
    contig.

    Args:
        uri: ``file.metacool`` or ``file.mmetacool::entries/<name>``.
        scope: ``genome`` (default), ``contig`` (cooler's ``--cis-only``) or
            ``global`` (cooler's default, the whole matrix at once).
        name: ``/bins`` column to store the weights as.
        ignore_diags: Diagonals of each contig left out (cooler's default 2).
        mad_max: Drop bins whose marginal is this many median absolute
            deviations below the median (0 disables).
        min_nnz: Drop bins with fewer non-zero pixels.
        min_count: Drop bins with a smaller total count.
        tol: Variance of the marginals at which balancing has converged.
        max_iters: Iteration limit.
        chunksize: Pixels read per chunk.
        store: Write the weights into the file.
        force: Replace an existing column of that name.
        commands: Command lines to record in ``metadata["commands"]``.

    Returns:
        The :class:`Normalization`.

    Raises:
        ValueError: On an unknown *scope*, a negative parameter, or a name that
            is taken (without *force*) or cannot hold a weight column.
    """
    if scope not in ICE_SCOPES:
        raise ValueError(f"scope must be one of {', '.join(ICE_SCOPES)}, got {scope!r}")
    for option, value in (
        ("ignore_diags", ignore_diags),
        ("mad_max", mad_max),
        ("min_nnz", min_nnz),
        ("min_count", min_count),
        ("max_iters", max_iters),
    ):
        if value < 0:
            raise ValueError(f"{option} must be >= 0, got {value}")
    if tol <= 0:
        raise ValueError(f"tol must be > 0, got {tol}")
    # Held from before the map is read to after the weights are stored:
    # no other writer changes it meanwhile.
    with _held_for_store(uri, store):
        mc = _open_target(uri, name, store=store, force=force)
        options = {
            "ignore_diags": int(ignore_diags),
            "mad_max": int(mad_max),
            "min_nnz": int(min_nnz),
            "min_count": int(min_count),
            "tol": tol,
            "max_iters": int(max_iters),
            "chunksize": chunksize,
        }
        if scope == "genome":
            weights, stats = _ice_by_genome(mc, **options)
            depends = [DEPENDS_ON_PIXELS, DEPENDS_ON_ASSIGNMENT]
        else:
            weights, stats = _ice_with_cooler(mc, cis_only=scope == "contig", **options)
            depends = [DEPENDS_ON_PIXELS]
        parameters: dict[str, object] = {
            "scope": scope,
            **{k: v for k, v in options.items() if k != "chunksize"},
            **stats,
        }
        return _finish(
            mc,
            name,
            weights,
            method="ice",
            depends_on=depends,
            parameters=parameters,
            store=store,
            force=force,
            commands=commands,
        )


# ----------------------------------------------------------------------------
# NormCC
# ----------------------------------------------------------------------------


def normalize_normcc(
    uri: str | Path,
    *,
    sites: str | None = None,
    enzymes: Sequence[str] | None = None,
    fasta: str | Path | None = None,
    min_length: int = DEFAULT_MIN_LENGTH,
    min_signal: float = DEFAULT_MIN_SIGNAL,
    spurious: float = DEFAULT_SPURIOUS_QUANTILE,
    name: str = "normcc",
    chunksize: int = DEFAULT_CHUNKSIZE,
    store: bool = True,
    force: bool = False,
    commands: Sequence[str] = (),
) -> Normalization:
    """Normalize a map with NormCC, MetaCC's regression normalization.

    NormCC (Du & Sun, Nat. Commun. 2023) fits one negative binomial regression
    (dispersion fixed at 1) of each contig's total contacts with the other
    contigs on the logs of its restriction-site count, its length and its
    intra-contig contacts -- which stand in for its coverage, so no coverage is
    needed.  A contact between contigs ``i`` and ``j`` is normalized as
    ``raw * max(mu) / sqrt(mu_i mu_j)``, ``mu`` being the fitted means: the
    weight of a contig is ``sqrt(max(mu) / mu_i)``.  Without sites (no
    *sites* nor *enzymes*), the site term is left out, as in MetaCC's mode for
    an unknown enzyme.

    As in MetaCC, contigs shorter than *min_length*, with no intra-contig
    contact, no site, or no contact of at least *min_signal* with another long
    contig are left out (masked).  The *spurious* quantile of the normalized
    contacts between two kept contigs is stored as the spurious-contact cutoff
    (``network --drop-spurious`` drops the contig pairs below it).

    Args:
        uri: ``file.metacool`` or ``file.mmetacool::entries/<name>``.
        sites: ``/contigs`` column of restriction-site counts, used as they
            are.  Or count them with *enzymes*.
        enzymes: Enzymes whose sites are counted in *fasta*.
        fasta: The assembly the map was built from.
        min_length: Shorter contigs are left out.
        min_signal: Contigs whose largest contact with another contig is
            smaller are left out.
        spurious: Quantile stored as the spurious-contact cutoff.
        name: ``/bins`` column.
        chunksize: Pixels read per chunk.
        store: Write the weights into the file.
        force: Replace an existing column of that name.
        commands: Command lines to record in ``metadata["commands"]``.

    Returns:
        The :class:`Normalization`.

    Raises:
        ValueError: If too few contigs pass the filters to fit the regression.
    """
    _check_quantile(spurious, "spurious")
    # Held from before the map is read to after the weights are stored:
    # no other writer changes it meanwhile.
    with _held_for_store(uri, store):
        mc = _open_target(uri, name, store=store, force=force)
        contigs = mc.contigs
        lengths = contigs["length"].to_numpy(dtype=np.float64)
        found = _contig_sites(mc, sites=sites, enzymes=enzymes, fasta=fasta)

        matrix = _contig_matrix(mc, chunksize)
        intra = matrix.diagonal()
        off = _off_diagonal(matrix)
        long_enough = lengths >= min_length
        keep = (
            long_enough
            & (_largest_contact(off, long_enough) >= min_signal)
            & (intra > 0)
        )
        covariates = ["intercept"]
        columns = [np.ones(len(lengths))]
        if found is not None:
            keep &= np.isfinite(found.counts) & (found.counts > 0)
            covariates.append("log_sites")
            columns.append(np.log(np.where(keep, found.counts, 1.0)))
        covariates += ["log_length", "log_intra_contacts"]
        columns += [
            np.log(np.maximum(lengths, 1.0)),
            np.log(np.where(keep, intra, 1.0)),
        ]
        n_kept = int(keep.sum())
        if n_kept <= len(covariates):
            raise ValueError(
                f"Only {n_kept} contig(s) pass NormCC's filters (length >= {min_length} "
                f"bp, an intra-contig contact, a contact of at least {min_signal} with "
                "another contig"
                + (", a restriction site" if found is not None else "")
                + f"), too few to fit {len(covariates)} coefficients."
            )
        signal = off @ keep.astype(np.float64)
        X = np.column_stack(columns)[keep]
        beta = fit_nb_glm(signal[keep], X, alpha=1.0)
        mu = np.exp(X @ beta)
        contig_weights = np.full(len(lengths), np.nan)
        contig_weights[keep] = np.sqrt(mu.max() / mu)
        cutoff, n_contacts = _spurious_cutoff(off, contig_weights, keep, spurious)
        logger.info(
            "NormCC fitted on %s contig(s): %s; spurious-contact cutoff %.6g "
            "(the %g quantile of %s normalized contig-pair contacts)",
            f"{n_kept:,}",
            ", ".join(f"{c}={b:.4g}" for c, b in zip(covariates, beta, strict=True)),
            cutoff,
            spurious,
            f"{n_contacts:,}",
        )
        parameters: dict[str, object] = {
            "covariates": covariates,
            "coefficients": beta,
            "dispersion": 1.0,
            "min_length": int(min_length),
            "min_signal": min_signal,
            "n_contigs_fitted": n_kept,
            "spurious_quantile": spurious,
            SPURIOUS_CUTOFF_ATTR: cutoff,
            SPURIOUS_RULE_ATTR: "<",
            **(found.attrs if found is not None else {}),
        }
        return _finish(
            mc,
            name,
            contig_weights[mc._bin_contig_ids],
            method="normcc",
            depends_on=[
                DEPENDS_ON_PIXELS,
                *(found.depends if found is not None else []),
            ],
            parameters=parameters,
            store=store,
            force=force,
            commands=commands,
        )


# ----------------------------------------------------------------------------
# HiCzin
# ----------------------------------------------------------------------------


def normalize_hiczin(
    uri: str | Path,
    *,
    coverage: str,
    sites: str | None = None,
    enzymes: Sequence[str] | None = None,
    fasta: str | Path | None = None,
    min_length: int = DEFAULT_MIN_LENGTH,
    min_signal: float = DEFAULT_MIN_SIGNAL,
    spurious: float = DEFAULT_SPURIOUS_QUANTILE,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    seed: int = 0,
    name: str = "hiczin",
    chunksize: int = DEFAULT_CHUNKSIZE,
    store: bool = True,
    force: bool = False,
    commands: Sequence[str] = (),
) -> Normalization:
    """Normalize a map with HiCzin, the normalization of HiCBin.

    HiCzin (Du et al., J. Comput. Biol. 2022) fits a zero-inflated negative
    binomial regression of the contacts between two contigs of the *same
    species* -- the contacts that are genuine -- on the logs of the products of
    their restriction-site counts, lengths and coverages, each standardised.  It
    learns those pairs from a taxonomic labelling; a metacool file has a better
    one, its genome assignment, so the pairs here are the pairs of contigs of
    one genome (``__unbinned__`` contigs excluded), zero contacts included.  A
    contact is normalized by its fitted mean (the count model), which is one
    factor per contig.  Without sites, the site term is left out (HiCzin_LC).

    As in HiCBin, contigs shorter than *min_length*, or whose largest contact
    with another long contig is below *min_signal*, are left out; a zero
    coverage is raised to the smallest positive one and a zero site count to
    1.  The *spurious* quantile of the positive normalized contacts of the
    fitted pairs is stored as the spurious-contact cutoff (a contact at or
    below it is spurious).

    Args:
        uri: ``file.metacool`` or ``file.mmetacool::entries/<name>``.
        coverage: ``/contigs`` column of shotgun coverages.
        sites: ``/contigs`` column of restriction-site counts, or count them
            with *enzymes*.
        enzymes: Enzymes whose sites are counted in *fasta*.
        fasta: The assembly the map was built from.
        min_length: Shorter contigs are left out.
        min_signal: Contigs whose largest contact with another contig is
            smaller are left out.
        spurious: Quantile stored as the spurious-contact cutoff.
        max_pairs: Fit on at most this many pairs, a uniform sample of them
            when the genomes hold more (0: all of them).
        seed: Seed of that sample.
        name: ``/bins`` column.
        chunksize: Pixels read per chunk.
        store: Write the weights into the file.
        force: Replace an existing column of that name.
        commands: Command lines to record in ``metadata["commands"]``.

    Returns:
        The :class:`Normalization`.

    Raises:
        ValueError: If no genome has two contigs passing the filters, or a
            covariate is constant across the fitted pairs.
    """
    _check_quantile(spurious, "spurious")
    if max_pairs < 0:
        raise ValueError(f"max_pairs must be >= 0, got {max_pairs}")
    # Held from before the map is read to after the weights are stored:
    # no other writer changes it meanwhile.
    with _held_for_store(uri, store):
        mc = _open_target(uri, name, store=store, force=force)
        contigs = mc.contigs
        lengths = contigs["length"].to_numpy(dtype=np.float64)
        cov = _numeric_contig_column(contigs, coverage, what="coverage")
        found = _contig_sites(mc, sites=sites, enzymes=enzymes, fasta=fasta)

        matrix = _contig_matrix(mc, chunksize)
        off = _off_diagonal(matrix)
        long_enough = lengths >= min_length
        keep = (
            long_enough
            & (_largest_contact(off, long_enough) >= min_signal)
            & np.isfinite(cov)
        )
        if found is not None:
            keep &= np.isfinite(found.counts)
        positive = keep & (cov > 0)
        if not positive.any():
            raise ValueError(
                "No contig passing HiCzin's filters has a positive coverage in "
                f"/contigs/{coverage}."
            )
        cov = np.where(keep & (cov == 0), cov[positive].min(), cov)

        covariates = []
        logs = []
        if found is not None:
            covariates.append("log_sites")
            logs.append(np.log(np.where(found.counts > 0, found.counts, 1.0)))
        covariates += ["log_length", "log_coverage"]
        logs += [np.log(np.maximum(lengths, 1.0)), np.log(np.where(cov > 0, cov, 1.0))]
        log_cov = np.column_stack(logs)

        # The pairs of kept contigs of one genome, enumerated genome by genome.
        genome = contigs["genome_id"].to_numpy(dtype=np.int64)
        unbinned = unbinned_genome_id(mc.genomes["name"])
        member = keep if unbinned is None else keep & (genome != unbinned)
        candidates = np.flatnonzero(member)
        order = candidates[np.argsort(genome[candidates], kind="stable")]
        _, starts, sizes = np.unique(
            genome[order], return_index=True, return_counts=True
        )
        total = int(np.sum(sizes * (sizes - 1) // 2))
        if total == 0:
            raise ValueError(
                "HiCzin fits the contacts between contigs of one genome, but no genome "
                f"has two contigs passing its filters (length >= {min_length} bp, a "
                f"contact of at least {min_signal} with another contig, a coverage)."
            )
        sampled = bool(max_pairs) and total > max_pairs
        if sampled:
            rng = np.random.default_rng(seed)
            index = np.sort(rng.choice(total, size=max_pairs, replace=False))
        else:
            index = np.arange(total, dtype=np.int64)
        group, a, b = pairs_within_groups(sizes, index)
        ci = order[starts[group] + a]
        cj = order[starts[group] + b]
        y = np.asarray(matrix[ci, cj]).ravel()

        pair_logs = log_cov[ci] + log_cov[cj]
        means = pair_logs.mean(axis=0)
        sds = pair_logs.std(axis=0, ddof=1) if len(y) > 1 else np.zeros(len(covariates))
        # Identical values can leave a standard deviation of rounding error.
        constant = [
            c
            for c, s, m in zip(covariates, sds, means, strict=True)
            if not s > 1e-9 * max(1.0, abs(m))
        ]
        if constant:
            raise ValueError(
                f"HiCzin standardises its covariates, but {', '.join(constant)} is "
                "constant across the fitted pairs."
            )
        X = np.column_stack([np.ones(len(y)), (pair_logs - means) / sds])
        fit = fit_zinb(y, X)

        contig_weights = np.full(len(lengths), np.nan)
        contig_weights[keep] = hiczin_contig_weights(
            fit.beta, means, sds, log_cov[keep]
        )
        normalized = y * contig_weights[ci] * contig_weights[cj]
        positive_norm = normalized[normalized > 0]
        cutoff = (
            float(np.quantile(positive_norm, spurious))
            if positive_norm.size
            else float("nan")
        )
        logger.info(
            "HiCzin fitted on %s pair(s) of contigs of one genome%s (%s with a "
            "contact): count model %s, theta=%.4g; spurious-contact cutoff %.6g",
            f"{len(y):,}",
            f" sampled from {total:,}" if sampled else "",
            f"{int(np.count_nonzero(y)):,}",
            ", ".join(
                f"{c}={v:.4g}"
                for c, v in zip(["intercept", *covariates], fit.beta, strict=True)
            ),
            fit.theta,
            cutoff,
        )
        parameters: dict[str, object] = {
            "covariates": ["intercept", *covariates],
            "coefficients": fit.beta,
            "zero_inflation_coefficients": fit.gamma,
            "theta": fit.theta,
            "covariate_means": means,
            "covariate_sds": sds,
            "converged": fit.converged,
            "min_length": int(min_length),
            "min_signal": min_signal,
            "coverage": f"{CONTIG_COLUMN_PREFIX}{coverage}",
            "n_pairs_fitted": len(y),
            "n_pairs_total": total,
            "seed": seed if sampled else None,
            "spurious_quantile": spurious,
            SPURIOUS_CUTOFF_ATTR: cutoff,
            SPURIOUS_RULE_ATTR: "<=",
            **(found.attrs if found is not None else {}),
        }
        return _finish(
            mc,
            name,
            contig_weights[mc._bin_contig_ids],
            method="hiczin",
            # Fitted on the pairs of contigs of one genome: a reassignment
            # changes the pairs, hence the fit.
            depends_on=[
                DEPENDS_ON_PIXELS,
                DEPENDS_ON_ASSIGNMENT,
                f"{CONTIG_COLUMN_PREFIX}{coverage}",
                *(found.depends if found is not None else []),
            ],
            parameters=parameters,
            store=store,
            force=force,
            commands=commands,
        )
