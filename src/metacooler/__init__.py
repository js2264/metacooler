"""metacooler: Hi-C contact maps for metagenomes.

Public names are loaded lazily on first access via ``__getattr__`` to keep
import time fast and avoid circular dependencies.  All symbols listed in
``__all__`` are importable directly from ``metacooler``:

.. code-block:: python

    import metacooler
    mc = metacooler.MetaCooler("file.metacool")
    metacooler.from_pairs(...)

Each ``if name == "..."`` branch inside :func:`__getattr__` imports the real
symbol from its submodule and returns it; Python then caches the resolved
attribute on the module object so subsequent accesses skip the lookup.
"""

from metacooler._version import __version__

__all__ = [
    "__version__",
    "MetaCooler",
    # creation
    "from_pairs",
    "from_cool",
    "collect",
    # transforms and readers
    "reassign",
    "extract",
    "flatten",
    "merge",
    "binnify",
    "zoomify",
    "aggregate",
    "network",
    "dump_table",
    # quality control
    # Spelled ``contact_qc``, not ``qc``: ``metacooler.qc`` is the tally engine
    # *module*, and importing it anywhere (``io.create`` does, at import time)
    # binds that name on the package — which would shadow a lazily-exported
    # function of the same name, making ``metacooler.qc`` mean different things
    # depending on what had been imported first.
    "contact_qc",
    "contact_probability",
    "plot_ps",
    # format
    "is_metacooler",
    "is_mmetacooler",
    "validate",
    "validate_mmetacooler",
    "check",
    "CheckResult",
    "parse_uri",
    # mmetacool collections
    "create_mmetacool",
    "add_entry",
    "remove_entry",
    "copy_entry",
    "move_entry",
    "list_entries",
    "copy_entry_cross_file",
]


def __getattr__(name: str) -> object:
    """Lazily import public symbols from their submodules on first access.

    Args:
        name: Attribute name requested from the ``metacooler`` package.

    Returns:
        The resolved symbol from its source submodule.

    Raises:
        AttributeError: If *name* is not a recognised public attribute.
    """
    if name == "MetaCooler":
        from metacooler.api import MetaCooler as _MetaCooler

        return _MetaCooler
    if name == "from_pairs":
        from metacooler.io.create import from_pairs as _from_pairs

        return _from_pairs
    if name == "reassign":
        from metacooler.io.reassign import reassign as _reassign

        return _reassign
    if name == "from_cool":
        from metacooler.io.create import from_cool as _from_cool

        return _from_cool
    if name == "collect":
        from metacooler.io.collect import collect as _collect

        return _collect
    if name == "extract":
        from metacooler.io.extract import extract as _extract

        return _extract
    if name == "flatten":
        from metacooler.io.flatten import flatten as _flatten

        return _flatten
    if name == "merge":
        from metacooler.io.merge import merge as _merge

        return _merge
    if name == "network":
        from metacooler.io.network import network as _network

        return _network
    if name == "dump_table":
        from metacooler.io.dump import dump_table as _dump_table

        return _dump_table
    if name == "binnify":
        from metacooler.ops.binnify import binnify as _binnify

        return _binnify
    if name == "zoomify":
        from metacooler.ops.zoomify import zoomify as _zoomify

        return _zoomify
    if name == "aggregate":
        from metacooler.ops.aggregate import aggregate as _aggregate

        return _aggregate
    if name == "contact_qc":
        from metacooler.ops.qc import qc as _qc

        return _qc
    if name == "contact_probability":
        from metacooler.ops.ps import contact_probability as _contact_probability

        return _contact_probability
    if name == "plot_ps":
        from metacooler.ops.ps import plot_ps as _plot_ps

        return _plot_ps
    if name == "validate":
        from metacooler.schema import validate as _validate

        return _validate
    if name == "validate_mmetacooler":
        from metacooler.schema import validate_mmetacooler as _validate_mmetacooler

        return _validate_mmetacooler
    if name == "check":
        from metacooler.schema import check as _check

        return _check
    if name == "CheckResult":
        from metacooler.schema import CheckResult as _CheckResult

        return _CheckResult
    if name == "is_metacooler":
        from metacooler.schema import is_metacooler as _is_metacooler

        return _is_metacooler
    if name == "is_mmetacooler":
        from metacooler.schema import is_mmetacooler as _is_mmetacooler

        return _is_mmetacooler
    if name == "parse_uri":
        from metacooler.uri import parse_uri as _parse_uri

        return _parse_uri
    if name == "create_mmetacool":
        from metacooler.io.mmetacool import create_mmetacool as _create_mmetacool

        return _create_mmetacool
    if name == "add_entry":
        from metacooler.io.mmetacool import add_entry as _add_entry

        return _add_entry
    if name == "remove_entry":
        from metacooler.io.mmetacool import remove_entry as _remove_entry

        return _remove_entry
    if name == "copy_entry":
        from metacooler.io.mmetacool import copy_entry as _copy_entry

        return _copy_entry
    if name == "move_entry":
        from metacooler.io.mmetacool import move_entry as _move_entry

        return _move_entry
    if name == "list_entries":
        from metacooler.io.mmetacool import list_entries as _list_entries

        return _list_entries
    if name == "copy_entry_cross_file":
        from metacooler.io.mmetacool import (
            copy_entry_cross_file as _copy_entry_cross_file,
        )

        return _copy_entry_cross_file
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
