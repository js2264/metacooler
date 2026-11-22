"""Normalize a metacool file that already exists: ICE weights in ``/bins``.

:func:`normalize_ice` balances a map by iterative correction and stores the
weights as a column of ``/bins`` (see :mod:`metacooler.io.weights`): the
normalized value of a pixel is ``raw * w[bin1] * w[bin2]``, the convention of
cooler's balancing weights, so cooler, cooltools and HiGlass read the result
as they read ``cooler balance``'s.  The pixels are never rewritten.

By default each **genome** is balanced on its own: each MAG's block of the
matrix, with the contacts between genomes left out of the fit.
``scope="contig"`` and ``scope="global"`` are cooler's ``--cis-only`` and
whole-matrix balancing.  The arithmetic lives in :mod:`metacooler.normalize`;
this module reads what it needs from the file.

Other tools store their own normalizations with
:func:`metacooler.io.weights.store_weights`.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

import numpy as np

from metacooler.api import MetaCooler
from metacooler.io.weights import (
    DEPENDS_ON_ASSIGNMENT,
    DEPENDS_ON_PIXELS,
    check_weight_name,
    store_weights,
    weight_attributes,
)
from metacooler.normalize import (
    compact_groups,
    ice_iterate,
    ice_mad_drop,
    ice_marginals,
)
from metacooler.qc import unbinned_genome_id
from metacooler.schema import DEFAULT_CHUNKSIZE

logger = logging.getLogger("metacooler")

#: What :func:`normalize_ice` balances separately.
ICE_SCOPES = ("genome", "contig", "global")

#: Default name of an ICE weight column: cooler's, so tools that read
#: ``weight`` by default (``cooler show``, HiGlass, cooltools) pick it up.
DEFAULT_ICE_NAME = "weight"

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
