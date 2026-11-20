"""ICE arithmetic for any grouping of the bins, pure numpy.

``metacooler normalize ice`` (:mod:`metacooler.ops.normalize`) balances a map
by iterative correction, and this module holds the arithmetic that produces
its weights.  It never opens a file: the file-facing half reads the pixels,
hands this module plain arrays, and stores what it returns (see
:mod:`metacooler.io.weights`).

cooler balances either the whole matrix or each *chromosome* -- a contig, in
a metacool file -- on its own; a metagenome wants each *genome* balanced on
its own, which cooler cannot express.  :func:`ice_iterate`,
:func:`ice_marginals` and :func:`ice_mad_drop` are cooler's iterative
correction written for an arbitrary grouping of the bins.  The iteration
follows cooler 0.12 step for step, so a grouping that coincides with cooler's
gives cooler's weights.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger("metacooler")


class IceStats(NamedTuple):
    """How an iterative correction ended, per group of bins.

    Attributes:
        scale: Mean marginal of each group at its last iteration (``NaN`` for a
            group with nothing to balance).
        var: Variance of the marginals of each group at its last iteration.
        converged: Whether each group reached ``var < tol``.  A group with
            nothing left to balance counts as converged unless every group is
            empty, as in cooler.
        n_iters: Iterations run.
    """

    scale: np.ndarray
    var: np.ndarray
    converged: np.ndarray
    n_iters: int


def ice_marginals(
    bin1: np.ndarray, bin2: np.ndarray, values: np.ndarray, n_bins: int
) -> np.ndarray:
    """Row sums of a symmetric matrix stored as its upper triangle.

    A pixel adds its value to the marginal of both of its bins, except a pixel
    on the main diagonal, which adds it once (cooler 0.12's convention).
    """
    offdiag = bin1 != bin2
    return np.bincount(
        bin1[offdiag], weights=values[offdiag], minlength=n_bins
    ) + np.bincount(bin2, weights=values, minlength=n_bins)


def ice_mad_drop(marg: np.ndarray, groups: np.ndarray, mad_max: float) -> np.ndarray:
    """Bins whose marginal is an outlier on the low side: cooler's MAD filter.

    Each marginal is first divided by the median positive marginal of its
    group, then bins more than *mad_max* median absolute deviations below the
    median of the log marginals are dropped, as are empty bins.  cooler takes
    the medians per chromosome -- a contig here, often a single bin, which
    always scores exactly 1 and so is never filtered -- where a metagenome
    wants them per genome.

    Returns:
        A boolean mask of the bins to drop.
    """
    import pandas as pd

    marg = np.array(marg, dtype=np.float64)
    positive = marg > 0
    if positive.any():
        medians = pd.Series(marg[positive]).groupby(groups[positive]).median()
        per_group = np.full(int(groups.max()) + 1, np.nan)
        per_group[medians.index.to_numpy(dtype=np.int64)] = medians.to_numpy()
        marg[positive] /= per_group[groups[positive]]
    log_nz = np.log(marg[marg > 0])
    if not len(log_nz):
        return np.zeros(len(marg), dtype=bool)
    median = np.median(log_nz)
    deviation = np.median(np.abs(log_nz - median))
    cutoff = np.exp(median - mad_max * deviation)
    return marg < cutoff


def ice_iterate(
    bias: np.ndarray,
    groups: np.ndarray,
    get_marginals: Callable[[np.ndarray], np.ndarray],
    *,
    tol: float,
    max_iters: int,
    rescale_marginals: bool = True,
) -> IceStats:
    """Balance *bias* in place, separately for each group of bins.

    cooler 0.12's ``_iterate``: every iteration computes all the marginals in
    one pass over the pixels (so pixels between two groups must already be
    filtered out by *get_marginals*), divides each bin's bias by its marginal
    over its group's mean, and freezes a group once the variance of its
    marginals drops below *tol*.  Bins left with no data are ``NaN``; with
    *rescale_marginals*, each group's balanced rows then sum to 1.

    Args:
        bias: Initial weights, one per bin; 0 for a bin filtered out.
            Modified in place.
        groups: Group of every bin, numbered from 0.
        get_marginals: Marginals of the matrix weighted by a bias vector.
        tol: Variance of the marginals at which a group has converged.
        max_iters: Iteration limit.
        rescale_marginals: Scale each group so its balanced rows sum to 1.

    Returns:
        An :class:`IceStats`.
    """
    n_groups = int(groups.max()) + 1 if len(groups) else 0
    scale = np.full(n_groups, np.nan)
    var = np.full(n_groups, np.nan)
    active = np.ones(n_groups, dtype=bool)
    last_marg = np.full(len(bias), np.nan)
    n_iters = 0
    for i in range(max_iters):
        n_iters = i + 1
        marg = get_marginals(bias)
        sel = active[groups]
        last_marg[sel] = marg[sel]

        nz = sel & (marg != 0)
        g = groups[nz]
        n = np.bincount(g, minlength=n_groups)
        active &= n > 0
        if not active.any():
            break

        with np.errstate(divide="ignore", invalid="ignore"):
            mean = np.bincount(g, weights=marg[nz], minlength=n_groups) / n
            dev = (marg[nz] - mean[g]) ** 2
            v = np.bincount(g, weights=dev, minlength=n_groups) / n
        scale[active] = mean[active]
        var[active] = v[active]

        sel = active[groups]
        factor = marg[sel] / mean[groups[sel]]
        factor[factor == 0] = 1
        bias[sel] /= factor

        logger.debug(
            "Iteration %d: %d of %d group(s) left to converge, maximum variance "
            "of marginals %.6g",
            i + 1,
            int(np.count_nonzero(v[active] >= tol)),
            n_groups,
            float(np.max(v[active])),
        )
        active &= ~(v < tol)
        if not active.any():
            break

    bias[last_marg == 0] = 0
    bias[bias == 0] = np.nan
    if rescale_marginals:
        s = scale[groups]
        ok = np.isfinite(s)
        bias[ok] /= np.sqrt(s[ok])

    empty = np.isnan(scale)
    if max_iters > 0 and not empty.all():
        var[empty] = 0.0
    with np.errstate(invalid="ignore"):
        converged = var < tol
    return IceStats(scale, var, converged, n_iters)


def compact_groups(keys: np.ndarray) -> np.ndarray:
    """Renumber arbitrary integer group keys as ``0 .. n_groups - 1``."""
    return np.unique(np.asarray(keys), return_inverse=True)[1].astype(np.int64)
