"""Tests for metacooler.__init__ (lazy attribute exports)."""

from __future__ import annotations

import pytest


class TestInit:
    def test_version_string(self):
        import metacooler

        assert isinstance(metacooler.__version__, str)
        assert len(metacooler.__version__) > 0

    def test_metacooler_class_accessible(self):
        from metacooler import MetaCooler

        assert MetaCooler is not None

    def test_metacooler_class_correct_type(self):
        from metacooler import MetaCooler
        from metacooler.api import MetaCooler as _MetaCooler

        assert MetaCooler is _MetaCooler

    def test_from_pairs_callable(self):
        from metacooler import from_pairs

        assert callable(from_pairs)

    def test_from_cool_callable(self):
        from metacooler import from_cool

        assert callable(from_cool)

    def test_is_metacooler_callable(self):
        from metacooler import is_metacooler

        assert callable(is_metacooler)

    def test_is_mmetacooler_callable(self):
        from metacooler import is_mmetacooler

        assert callable(is_mmetacooler)

    def test_parse_uri_callable(self):
        from metacooler import parse_uri

        assert callable(parse_uri)

    def test_unknown_attribute_raises_attribute_error(self):
        import metacooler

        with pytest.raises(AttributeError, match="no attribute"):
            _ = metacooler.nonexistent_attribute_xyz

    def test_all_exports_in_dunder_all(self):
        import metacooler

        expected = {
            "__version__",
            "MetaCooler",
            # creation
            "from_pairs",
            "from_cool",
            "collect",
            # transforms and readers
            "reassign",
            "annotate",
            "write_clusters",
            "drop_clusters",
            "extract",
            "write_fasta",
            "flatten",
            "merge",
            "binnify",
            "zoomify",
            "aggregate",
            "network",
            "dump_table",
            # normalization
            "normalize_ice",
            "normalize_normcc",
            "normalize_hiczin",
            "store_weights",
            "read_weights",
            # quality control
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
            # in-place writes
            "lock",
            "staged_edit",
            "MapLockError",
            "MapLockedError",
            "MapChangedError",
            # mmetacool collections
            "create_mmetacool",
            "add_entry",
            "remove_entry",
            "copy_entry",
            "move_entry",
            "list_entries",
            "copy_entry_cross_file",
        }
        assert expected == set(metacooler.__all__)

    def test_all_exports_have_a_getattr_branch(self):
        """Every name in ``__all__`` must resolve through the lazy dispatch.

        Guards the convention in CLAUDE.md: public API is added by extending
        ``__all__`` *and* adding an ``if name == ...`` branch.  Forgetting the
        branch is otherwise only caught at first use.
        """
        import metacooler

        for name in metacooler.__all__:
            assert hasattr(metacooler, name), f"{name} has no __getattr__ branch"

    def test_type_checkers_see_every_export_at_its_source(self):
        """The ``if TYPE_CHECKING:`` block names each export as __getattr__ does.

        Without it a type checker types ``metacooler.from_pairs`` and every
        other lazy export as what ``__getattr__`` returns, ``object``, and
        flags each call in client code as not callable.  This keeps the block
        in step with ``__all__`` and with the dispatch.
        """
        import ast
        import importlib
        import inspect

        import metacooler

        block = next(
            node
            for node in ast.parse(inspect.getsource(metacooler)).body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Name)
            and node.test.id == "TYPE_CHECKING"
        )
        seen = {}
        for node in block.body:
            assert isinstance(node, ast.ImportFrom) and node.module
            for alias in node.names:
                seen[alias.asname or alias.name] = (node.module, alias.name)
        assert set(seen) == set(metacooler.__all__) - {"__version__"}
        for public, (module, name) in seen.items():
            source = getattr(importlib.import_module(module), name)
            assert getattr(metacooler, public) is source, public

    def test_the_package_says_it_is_typed(self):
        """PEP 561: without py.typed, mypy skips an installed package entirely."""
        from pathlib import Path

        import metacooler

        assert (Path(metacooler.__file__).parent / "py.typed").is_file()

    def test_the_type_checking_block_never_runs(self):
        """The block is for type checkers: at run time it binds no export."""
        import metacooler

        assert metacooler.TYPE_CHECKING is False
        assert "from_pairs" not in vars(metacooler)

    def test_all_exports_accessible(self):
        import metacooler

        for name in metacooler.__all__:
            obj = getattr(metacooler, name)
            assert obj is not None

    def test_from_pairs_is_create_from_pairs(self):
        from metacooler import from_pairs
        from metacooler.io.create import from_pairs as _from_pairs

        assert from_pairs is _from_pairs

    def test_from_cool_is_create_from_cool(self):
        from metacooler import from_cool
        from metacooler.io.create import from_cool as _from_cool

        assert from_cool is _from_cool

    def test_collect_is_io_collect(self):
        from metacooler import collect
        from metacooler.io.collect import collect as _collect

        assert collect is _collect

    def test_is_metacooler_is_schema_is_metacooler(self):
        from metacooler import is_metacooler
        from metacooler.schema import is_metacooler as _is_metacooler

        assert is_metacooler is _is_metacooler

    def test_is_mmetacooler_is_schema_is_mmetacooler(self):
        from metacooler import is_mmetacooler
        from metacooler.schema import is_mmetacooler as _is_mmetacooler

        assert is_mmetacooler is _is_mmetacooler

    def test_parse_uri_is_uri_parse_uri(self):
        from metacooler import parse_uri
        from metacooler.uri import parse_uri as _parse_uri

        assert parse_uri is _parse_uri

    def test_is_metacooler_works_with_real_file(self):
        from pathlib import Path

        import metacooler

        data_dir = Path(__file__).parent / "data"
        real_file = data_dir / "assembly_1_rep1_1kb.metacool"
        if real_file.exists():
            assert metacooler.is_metacooler(real_file)

    def test_is_mmetacooler_works_with_real_file(self):
        from pathlib import Path

        import metacooler

        data_dir = Path(__file__).parent / "data"
        real_file = data_dir / "assembly_1_1kb.mmetacool"
        if real_file.exists():
            assert metacooler.is_mmetacooler(real_file)

    def test_parse_uri_basic(self):
        from metacooler import parse_uri

        path, group = parse_uri("file.mmetacool::groups/entry1")
        assert str(path) == "file.mmetacool"
        assert "groups/entry1" in group

    def test_metacooler_class_instantiation_with_real_file(self):
        from pathlib import Path

        from metacooler import MetaCooler

        data_dir = Path(__file__).parent / "data"
        real_file = data_dir / "assembly_1_rep1_1kb.metacool"
        if real_file.exists():
            mc = MetaCooler(real_file)
            assert mc is not None
