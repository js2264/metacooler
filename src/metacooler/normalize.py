"""Normalization engine: per-bin weights from contact sums, pure numpy.

Every normalization of :mod:`metacooler.ops.normalize` factorises -- the
normalized value of a pixel is ``raw[i, j] * w[i] * w[j]`` -- and this module
holds the arithmetic that produces ``w``.  It never opens a file: the
file-facing half reads the pixels, hands this module plain arrays, and stores
what it returns (see :mod:`metacooler.io.weights`).

Four families live here:

- **Matrix balancing** (:func:`ice_iterate`, :func:`ice_marginals`,
  :func:`ice_mad_drop`): cooler's iterative correction, written for an
  arbitrary grouping of the bins.  cooler balances either the whole matrix or
  each *chromosome* -- a contig, in a metacool file -- on its own; a metagenome
  wants each *genome* balanced on its own, which cooler cannot express.  The
  iteration follows cooler 0.12 step for step, so a grouping that coincides
  with cooler's gives cooler's weights.
- **Geometric-mean normalizations** (:func:`metator_contig_weights`,
  :func:`inverse_sqrt`): MetaTOR's five network normalizations and the MGE
  score of Bignaud et al. (2025), each a division by the square root of the
  product of a per-contig quantity.
- **NormCC** (:func:`fit_nb_glm`): MetaCC's negative binomial regression of a
  contig's inter-contig signal on its site count, length and intra-contig
  contacts (Du & Sun 2023).
- **HiCzin** (:func:`fit_zinb`, :func:`hiczin_contig_weights`): the
  zero-inflated negative binomial regression of the contacts between contigs
  of one genome (Du et al. 2022), with :func:`pairs_within_groups` to list or
  sample those pairs without materialising them.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger("metacooler")


# ----------------------------------------------------------------------------
# Geometric-mean normalizations
# ----------------------------------------------------------------------------

#: MetaTOR's network normalizations, by the name ``metator network`` gives them
#: (its ``theoritical_hit`` spelled correctly; :func:`canonical_metator_mode`
#: accepts both).
METATOR_MODES = ("empirical_hit", "abundance", "length", "rs", "theoretical_hit")

#: Modes that need each contig's shotgun coverage.
METATOR_COVERAGE_MODES = frozenset({"abundance", "theoretical_hit"})

#: Modes that need each contig's restriction-site count.
METATOR_SITES_MODES = frozenset({"rs", "theoretical_hit"})

#: Modes that need each contig's Hi-C hits (counted from the pixels).
METATOR_HIT_MODES = frozenset({"empirical_hit", "length", "rs"})


def canonical_metator_mode(mode: str) -> str:
    """Return the canonical spelling of a MetaTOR mode.

    MetaTOR names one mode ``theoritical_hit`` and spells another ``RS``;
    either spelling, in any case, is accepted.

    Raises:
        ValueError: If *mode* is not a MetaTOR normalization.
    """
    key = mode.strip().lower()
    if key == "theoritical_hit":
        key = "theoretical_hit"
    if key not in METATOR_MODES:
        raise ValueError(
            f"Unknown MetaTOR normalization {mode!r}; choose one of "
            f"{', '.join(METATOR_MODES)}."
        )
    return key


def inverse_sqrt(values: np.ndarray) -> np.ndarray:
    """Return ``1 / sqrt(values)``, ``NaN`` wherever that is not a finite number.

    A zero, negative, missing or infinite value has no meaningful factor, so
    its bin is masked rather than given an infinite weight.
    """
    v = np.asarray(values, dtype=np.float64)
    out = np.full(v.shape, np.nan)
    ok = np.isfinite(v) & (v > 0)
    out[ok] = 1.0 / np.sqrt(v[ok])
    return out


def metator_contig_weights(
    mode: str,
    *,
    hits: np.ndarray | None = None,
    lengths: np.ndarray | None = None,
    coverage: np.ndarray | None = None,
    sites: np.ndarray | None = None,
) -> np.ndarray:
    """Per-contig weights of one of MetaTOR's network normalizations.

    MetaTOR divides the contacts of two contigs by a geometric mean
    (``metator.network.normalize_pair``); each such division is the product of
    one factor per contig.  The constants MetaTOR uses (``0.01`` in ``RS``,
    ``10`` in ``theoritical_hit``) are kept, so the normalized values are
    MetaTOR's own:

    ============== ======================================== =====================
    mode           MetaTOR divides the contacts by          weight of a contig
    ============== ======================================== =====================
    empirical_hit  sqrt(hit_i hit_j)                        1 / sqrt(hit)
    abundance      sqrt(cov_i cov_j)                        1 / sqrt(cov)
    length         sqrt(hit_i/len_i * hit_j/len_j)          sqrt(len / hit)
    rs             0.01 sqrt(hit_i/RS_i * hit_j/RS_j)       10 sqrt(RS / hit)
    theoretical_hit sqrt(10 cov_i sqrt(len_i RS_i) * ...)   1 / sqrt(10 cov sqrt(len RS))
    ============== ======================================== =====================

    Where MetaTOR's factor is zero it reports a normalized contact of 0; here
    the contig's weight is ``NaN`` (masked), cooler's convention.

    Args:
        mode: One of :data:`METATOR_MODES`.
        hits: Hi-C hits of each contig: the read ends that fall on it, an
            intra-contig pair counting twice.
        lengths: Contig lengths, in bp.
        coverage: Shotgun coverage of each contig.
        sites: Restriction-site count of each contig, as MetaTOR counts it
            (sites + 1).

    Raises:
        ValueError: If an input the mode needs is missing.
    """
    mode = canonical_metator_mode(mode)

    def need(value: np.ndarray | None, what: str) -> np.ndarray:
        if value is None:
            raise ValueError(f"MetaTOR's {mode!r} normalization needs {what}.")
        return np.asarray(value, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        match mode:
            case "empirical_hit":
                return inverse_sqrt(need(hits, "Hi-C hits"))
            case "abundance":
                return inverse_sqrt(need(coverage, "coverage"))
            case "length":
                return inverse_sqrt(need(hits, "Hi-C hits") / need(lengths, "lengths"))
            case "rs":
                return 10.0 * inverse_sqrt(
                    need(hits, "Hi-C hits") / need(sites, "restriction sites")
                )
            case _:  # theoretical_hit
                return inverse_sqrt(
                    10.0
                    * need(coverage, "coverage")
                    * np.sqrt(
                        need(lengths, "lengths") * need(sites, "restriction sites")
                    )
                )


def contig_hits(
    contig1: np.ndarray, contig2: np.ndarray, counts: np.ndarray, n_contigs: int
) -> np.ndarray:
    """Hi-C hits per contig: every contact counts once for each of its two ends.

    A contact inside one contig therefore counts twice for it -- as each read
    pair of ``metator network`` adds one hit to the contig of each of its ends.
    """
    w = np.asarray(counts, dtype=np.float64)
    return np.bincount(contig1, weights=w, minlength=n_contigs) + np.bincount(
        contig2, weights=w, minlength=n_contigs
    )


# ----------------------------------------------------------------------------
# Matrix balancing
# ----------------------------------------------------------------------------


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


# ----------------------------------------------------------------------------
# NormCC: negative binomial regression with a fixed dispersion
# ----------------------------------------------------------------------------


def _nb_deviance(y: np.ndarray, mu: np.ndarray, alpha: float) -> float:
    """Deviance of a negative binomial fit with variance ``mu + alpha mu**2``."""
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(y > 0, y * np.log(y / mu), 0.0)
    return float(
        2.0
        * np.sum(
            term - (y + 1.0 / alpha) * np.log((1.0 + alpha * y) / (1.0 + alpha * mu))
        )
    )


def fit_nb_glm(
    y: np.ndarray,
    X: np.ndarray,
    *,
    alpha: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 200,
) -> np.ndarray:
    """Fit a log-link negative binomial GLM with a known dispersion *alpha*.

    Iteratively reweighted least squares on the model of NormCC (MetaCC's
    ``statsmodels`` ``GLM(family=NegativeBinomial(alpha=1))``): the variance of
    an observation is ``mu + alpha * mu**2``.  IRLS and the L-BFGS fit MetaCC
    runs maximise the same likelihood, so they reach the same coefficients.

    Args:
        y: Non-negative responses.
        X: Design matrix, intercept column included.
        alpha: The fixed dispersion.
        tol: Relative change in deviance at which the fit has converged.
        max_iter: Iteration limit.

    Returns:
        The coefficients, one per column of *X*.

    Raises:
        ValueError: If *y* holds no positive value -- there is nothing to fit.
    """
    y = np.asarray(y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    if not (y > 0).any():
        raise ValueError("Every response is zero; there is no signal to fit.")
    mu = (y + y.mean()) / 2.0
    eta = np.log(mu)
    beta = np.zeros(X.shape[1])
    deviance = np.inf
    for _ in range(max_iter):
        weights = mu / (1.0 + alpha * mu)
        z = eta + (y - mu) / mu
        xtw = X.T * weights
        try:
            beta = np.linalg.solve(xtw @ X, xtw @ z)
        except np.linalg.LinAlgError:
            beta = np.linalg.lstsq(xtw @ X, xtw @ z, rcond=None)[0]
        eta = X @ beta
        mu = np.exp(eta)
        new = _nb_deviance(y, mu, alpha)
        if abs(new - deviance) <= tol * (abs(new) + 0.1):
            return beta
        deviance = new
    logger.warning(
        "The negative binomial fit did not converge in %d iterations.", max_iter
    )
    return beta


# ----------------------------------------------------------------------------
# HiCzin: zero-inflated negative binomial regression
# ----------------------------------------------------------------------------


#: Bounds of the zero-inflation logits and of ``log(theta)`` during a fit: a
#: probability of 1e-13 is as good as zero, and a dispersion beyond these
#: limits is a Poisson or a degenerate model for any practical count.
_ZINB_LOGIT_BOUND = 30.0
_ZINB_LOG_THETA_BOUND = 12.0

#: Largest per-observation gradient component of a converged fit.
_ZINB_GRADIENT_TOL = 1e-5


class ZinbFit(NamedTuple):
    """A zero-inflated negative binomial fit.

    Attributes:
        beta: Coefficients of the count model (log link), intercept first.
        gamma: Coefficients of the zero-inflation model (logit link).
        theta: Dispersion of the count model (``nbinom2``: variance
            ``mu + mu**2 / theta``).
        loglik: Log-likelihood at the fit.
        converged: Whether the optimiser reported convergence.
    """

    beta: np.ndarray
    gamma: np.ndarray
    theta: float
    loglik: float
    converged: bool


def zinb_loglik(
    params: np.ndarray, y: np.ndarray, X: np.ndarray, Z: np.ndarray
) -> tuple[float, np.ndarray]:
    """Log-likelihood of a zero-inflated negative binomial model and its gradient.

    The model of HiCzin (glmmTMB's ``nbinom2`` family with a ``ziformula``):
    a contact count is 0 with probability ``pi`` and otherwise negative
    binomial with mean ``mu`` and dispersion ``theta``, where
    ``log(mu) = X beta`` and ``logit(pi) = Z gamma``.

    Args:
        params: ``beta`` (``X.shape[1]`` values), ``gamma`` (``Z.shape[1]``),
            then ``log(theta)``.
        y: Counts.
        X: Count-model design matrix.
        Z: Zero-inflation design matrix.

    Returns:
        ``(loglik, gradient)``, the gradient with respect to *params*.
    """
    from scipy.special import digamma, gammaln

    kx, kz = X.shape[1], Z.shape[1]
    beta, gamma, tau = params[:kx], params[kx : kx + kz], float(params[-1])
    eta = np.clip(X @ beta, -700.0, 700.0)
    zeta = Z @ gamma
    theta = np.exp(tau)
    log_theta_mu = np.logaddexp(tau, eta)  # log(theta + mu)
    log_pi = -np.logaddexp(0.0, -zeta)
    log_1mpi = -np.logaddexp(0.0, zeta)
    pi = np.exp(log_pi)

    zero = y == 0
    pos = ~zero
    ll = np.empty_like(eta)
    d_eta = np.empty_like(eta)
    d_zeta = np.empty_like(eta)
    d_tau = np.empty_like(eta)

    yp, ep, lp = y[pos], eta[pos], log_theta_mu[pos]
    ll[pos] = (
        log_1mpi[pos]
        + gammaln(yp + theta)
        - gammaln(theta)
        - gammaln(yp + 1.0)
        + theta * (tau - lp)
        + yp * (ep - lp)
    )
    d_eta[pos] = yp - (yp + theta) * np.exp(ep - lp)
    d_zeta[pos] = -pi[pos]
    d_tau[pos] = theta * (
        digamma(yp + theta)
        - digamma(theta)
        + tau
        + 1.0
        - lp
        - (yp + theta) * np.exp(-lp)
    )

    l0 = log_theta_mu[zero]
    log_f0 = theta * (tau - l0)  # log P(NB = 0)
    a = log_pi[zero]
    b = log_1mpi[zero] + log_f0
    ll0 = np.logaddexp(a, b)
    nb_share = np.exp(b - ll0)
    ll[zero] = ll0
    d_eta[zero] = -nb_share * theta * np.exp(eta[zero] - l0)
    d_zeta[zero] = np.exp(a - ll0) - pi[zero]
    d_tau[zero] = theta * nb_share * (tau - l0 + 1.0 - np.exp(tau - l0))

    grad = np.concatenate([X.T @ d_eta, Z.T @ d_zeta, [d_tau.sum()]])
    return float(ll.sum()), grad


def fit_zinb(
    y: np.ndarray,
    X: np.ndarray,
    Z: np.ndarray | None = None,
    *,
    max_iter: int = 5000,
) -> ZinbFit:
    """Fit a zero-inflated negative binomial model by maximum likelihood.

    The likelihood of :func:`zinb_loglik` -- HiCzin's glmmTMB model -- is
    maximised with L-BFGS-B and its analytic gradient, starting from a
    negative binomial fit of the count model.

    Args:
        y: Counts.
        X: Count-model design matrix, intercept column included.
        Z: Zero-inflation design matrix (defaults to *X*, as in HiCzin).
        max_iter: Optimiser iteration limit.

    Returns:
        A :class:`ZinbFit`.
    """
    from scipy.optimize import minimize

    y = np.asarray(y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    Z = X if Z is None else np.asarray(Z, dtype=np.float64)
    n = max(len(y), 1)
    kx, kz = X.shape[1], Z.shape[1]
    beta0 = fit_nb_glm(y, X, alpha=1.0)
    gamma0 = np.zeros(kz)
    gamma0[0] = -2.0
    p0 = np.concatenate([beta0, gamma0, [0.0]])
    # With few or no zeros the zero-inflation is not identified and its
    # likelihood keeps rising as pi -> 0; bounding the logit (and the log
    # dispersion) lets the optimiser settle there instead of running off.
    bounds = (
        [(None, None)] * kx
        + [(-_ZINB_LOGIT_BOUND, _ZINB_LOGIT_BOUND)] * kz
        + [(-_ZINB_LOG_THETA_BOUND, _ZINB_LOG_THETA_BOUND)]
    )

    def objective(params: np.ndarray) -> tuple[float, np.ndarray]:
        ll, grad = zinb_loglik(params, y, X, Z)
        return -ll / n, -grad / n

    result = minimize(
        objective,
        p0,
        jac=True,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": max_iter, "ftol": 1e-13, "gtol": 1e-9},
    )
    params = result.x
    ll, grad = zinb_loglik(params, y, X, Z)
    # A line search that can no longer improve at the optimum ends L-BFGS-B
    # "abnormally"; what matters is whether the projected gradient is flat.
    lower = np.array([-np.inf if lo is None else lo for lo, _ in bounds])
    upper = np.array([np.inf if hi is None else hi for _, hi in bounds])
    step = -grad / n  # gradient of the objective
    free = ~(((params <= lower) & (step > 0)) | ((params >= upper) & (step < 0)))
    converged = bool(result.success) or bool(
        np.max(np.abs(step[free]), initial=0.0) < _ZINB_GRADIENT_TOL
    )
    if not converged:
        logger.warning(
            "The zero-inflated negative binomial fit stopped before converging: %s",
            result.message,
        )
    return ZinbFit(
        beta=params[:kx],
        gamma=params[kx : kx + kz],
        theta=float(np.exp(params[-1])),
        loglik=ll,
        converged=converged,
    )


def hiczin_contig_weights(
    beta: np.ndarray,
    means: np.ndarray,
    sds: np.ndarray,
    log_covariates: np.ndarray,
) -> np.ndarray:
    """Per-contig weights equivalent to HiCzin's normalization.

    HiCzin divides the contacts of contigs ``i`` and ``j`` by
    ``exp(b0 + sum_k b_k (log(x_ik x_jk) - m_k) / s_k)``, where each covariate
    ``x_k`` (sites, length, coverage) enters as the log of a product and is
    standardised with the mean ``m_k`` and standard deviation ``s_k`` of that
    log over the fitted pairs.  The log of a product is a sum, so the divisor is
    one factor per contig; its inverse is the weight returned here.

    Args:
        beta: Count-model coefficients, intercept first.
        means: Mean of each standardised covariate.
        sds: Standard deviation of each standardised covariate.
        log_covariates: ``(n_contigs, k)`` array of ``log x_ik``.

    Returns:
        One weight per contig.
    """
    slopes = np.asarray(beta[1:], dtype=np.float64) / np.asarray(sds, dtype=np.float64)
    intercept = float(beta[0]) - float(np.sum(slopes * np.asarray(means)))
    return np.exp(-intercept / 2.0 - np.asarray(log_covariates) @ slopes)


def pairs_within_groups(
    sizes: np.ndarray, index: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decode pair indices into ``(group, a, b)`` with ``a < b`` inside the group.

    The pairs of members of every group are numbered group after group, and
    within a group of ``n`` members in row-major order: ``(0, 1), (0, 2), ...,
    (0, n-1), (1, 2), ...``.  Decoding an index rather than listing the pairs
    is what lets a sample of them be drawn without ever materialising all
    ``sum(n * (n - 1) / 2)`` of them.

    Args:
        sizes: Number of members of each group.
        index: Pair indices, each in ``[0, total pairs)``.

    Returns:
        ``(group, a, b)`` arrays: the group of each pair and the positions of
        its two members inside it.
    """
    sizes = np.asarray(sizes, dtype=np.int64)
    index = np.asarray(index, dtype=np.int64)
    n_pairs = sizes * (sizes - 1) // 2
    offsets = np.concatenate([[0], np.cumsum(n_pairs)])
    group = np.searchsorted(offsets, index, side="right") - 1
    k = index - offsets[group]
    n = sizes[group]
    disc = np.maximum((2 * n - 1) ** 2 - 8 * k, 0).astype(np.float64)
    a = np.floor(((2 * n - 1) - np.sqrt(disc)) / 2).astype(np.int64)
    a = np.clip(a, 0, np.maximum(n - 2, 0))

    def row_start(r: np.ndarray) -> np.ndarray:
        return r * n - r * (r + 1) // 2

    # Floating-point rounding can land one row off either way.
    a = np.where(row_start(a) > k, a - 1, a)
    a = np.where(row_start(a + 1) <= k, a + 1, a)
    b = k - row_start(a) + a + 1
    return group, a, b
