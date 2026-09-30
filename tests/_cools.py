"""Tiny synthetic single-genome cool files for the ``collect`` tests.

Two "species" share contig names, as real assemblies do (every other one has a
``chromosome``), and E. coli's contigs are deliberately *not* longest-first, so
``--contig-order genome-then-length`` has to reorder — and re-sort — one of
them.
"""

from __future__ import annotations

from pathlib import Path

import cooler
import numpy as np
import pandas as pd

ECOLI = {"chromosome": 1000, "plasmid1": 250, "plasmid2": 400}
BSUB = {"chromosome": 800, "plasmid1": 300}


def make_cool(
    path: Path,
    chroms: dict[str, int],
    binsize: int = 100,
    *,
    seed: int = 0,
    n_contacts: int = 400,
    count_dtype: type = np.int32,
    count_scale: int = 1,
    bins: pd.DataFrame | None = None,
    symmetric_upper: bool = True,
) -> Path:
    """Write a small random-contact cool file and return its path.

    Args:
        path: Where to write it; its stem becomes the default genome name.
        chroms: Contig name to length, in the order the cool should list them.
        binsize: Fixed bin width, used unless *bins* is given.
        seed: Seed of the random contacts.
        n_contacts: Contacts drawn before aggregation into pixels; 0 writes an
            empty map.
        count_dtype: dtype of the ``count`` column.
        count_scale: Factor applied to every count, to reach wide values.
        bins: Explicit (e.g. variable-width) bins instead of a fixed grid.
        symmetric_upper: ``False`` writes a ``square`` storage-mode cool.

    Returns:
        *path*.
    """
    if bins is None:
        bins = cooler.binnify(pd.Series(chroms), binsize)
    rng = np.random.default_rng(seed)
    n = len(bins)
    b1 = rng.integers(0, n, n_contacts)
    b2 = rng.integers(0, n, n_contacts)
    pixels = (
        pd.DataFrame({"bin1_id": np.minimum(b1, b2), "bin2_id": np.maximum(b1, b2)})
        .value_counts()
        .rename("count")
        .reset_index()
        .sort_values(["bin1_id", "bin2_id"], ignore_index=True)
    )
    pixels["count"] = (pixels["count"] * count_scale).astype(count_dtype)
    cooler.create_cooler(
        str(path),
        bins,
        pixels,
        ordered=True,
        symmetric_upper=symmetric_upper,
        dtypes={"count": count_dtype},
    )
    return path
