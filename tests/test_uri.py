"""Tests for URI parsing utility."""

from __future__ import annotations

from pathlib import Path

from metacooler.uri import parse_uri


class TestParseUri:
    """Tests for parse_uri()."""

    def test_plain_metacool(self):
        path, group = parse_uri("file.metacool")
        assert path == Path("file.metacool")
        assert group == "/"

    def test_mmetacool_with_group(self):
        path, group = parse_uri("file.mmetacool::entries/entry1")
        assert path == Path("file.mmetacool")
        assert group == "/entries/entry1"

    def test_mmetacool_empty_group(self):
        path, group = parse_uri("file.mmetacool::")
        assert path == Path("file.mmetacool")
        assert group == "/"

    def test_mmetacool_root_group(self):
        path, group = parse_uri("file.mmetacool::/")
        assert path == Path("file.mmetacool")
        assert group == "/"

    def test_mmetacool_group_with_leading_slash(self):
        path, group = parse_uri("file.mmetacool::/entries/entry1")
        assert path == Path("file.mmetacool")
        assert group == "/entries/entry1"

    def test_absolute_path(self):
        path, group = parse_uri("/tmp/data/file.metacool")
        assert path == Path("/tmp/data/file.metacool")
        assert group == "/"

    def test_absolute_path_with_group(self):
        path, group = parse_uri("/tmp/data/file.mmetacool::entries/sample1")
        assert path == Path("/tmp/data/file.mmetacool")
        assert group == "/entries/sample1"

    def test_nested_group_path(self):
        path, group = parse_uri("file.mmetacool::resolutions/1000")
        assert path == Path("file.mmetacool")
        assert group == "/resolutions/1000"

    def test_whitespace_in_group_stripped(self):
        path, group = parse_uri("file.mmetacool::  entries/entry1  ")
        assert path == Path("file.mmetacool")
        assert group == "/entries/entry1"
