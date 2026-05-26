"""metacooler: Hi-C contact maps for metagenomes."""

from __future__ import annotations

from metacooler._version import __version__

__all__ = [
    "__version__",
    "MetaCooler",
    "from_pairs",
    "from_cool",
    "is_metacooler",
    "is_mmetacooler",
    "parse_uri",
]


def __getattr__(name: str) -> object:
    if name == "MetaCooler":
        from metacooler.api import MetaCooler as _MetaCooler

        return _MetaCooler
    if name == "from_pairs":
        from metacooler.io.create import from_pairs as _from_pairs

        return _from_pairs
    if name == "from_cool":
        from metacooler.io.create import from_cool as _from_cool

        return _from_cool
    if name == "is_metacooler":
        from metacooler.schema import is_metacooler as _is_metacooler

        return _is_metacooler
    if name == "is_mmetacooler":
        from metacooler.schema import is_mmetacooler as _is_mmetacooler

        return _is_mmetacooler
    if name == "parse_uri":
        from metacooler.uri import parse_uri as _parse_uri

        return _parse_uri
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
