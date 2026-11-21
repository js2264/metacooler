"""Comparing the edge tables of :func:`metacooler.io.network.network`."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal


def assert_ids_name(by_id: pd.DataFrame, by_name: pd.DataFrame, names) -> None:
    """*by_id* (``network(ids=True)``) is *by_name* row for row, in its order,
    once each id is named through *names* (``/contigs`` or ``/genomes``)."""
    names = np.asarray(names, dtype=object)
    assert list(by_id.columns) == ["node1", "node2", "value"]
    assert by_id["node1"].dtype == np.int64
    assert by_id["node2"].dtype == np.int64
    pd.testing.assert_index_equal(by_id.index, by_name.index)
    for node in ("node1", "node2"):
        named = names[by_id[node].to_numpy()].tolist()
        assert named == by_name[node].tolist()
    assert by_id["value"].dtype == by_name["value"].dtype
    assert_array_equal(by_id["value"].to_numpy(), by_name["value"].to_numpy())
