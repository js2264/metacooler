"""Tests for URI parsing utility."""

from __future__ import annotations

from pathlib import Path

import pytest

from metacooler.uri import parse_uri, sample_name


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


class TestSampleName:
    """The name of a map's sample, which names its report: metaTOR's rule."""

    @pytest.mark.parametrize(
        ("uri", "name"),
        [
            ("out/gut_T0.metacool", "gut_T0"),
            ("gut_T0.mmetacool", "gut_T0"),
            # An entry of a collection names its sample, as metaTOR names it.
            ("reps.mmetacool::entries/rep1", "rep1"),
            ("reps.mmetacool::/entries/rep1/", "rep1"),
            ("reps.mmetacool::", "reps"),
            # Only a map's extensions go, as in metaTOR.
            ("gut_T0.cool", "gut_T0.cool"),
            ("gut_T0.metacool.bak", "gut_T0.metacool.bak"),
            # A name that is an extension and nothing else stays as it is.
            (".metacool", ".metacool"),
        ],
    )
    def test_names(self, uri: str, name: str):
        assert sample_name(uri) == name

    def test_a_path(self, tmp_path):
        assert sample_name(tmp_path / "x.metacool") == "x"

    def test_an_empty_file_part_is_refused(self):
        with pytest.raises(ValueError, match="empty file path"):
            sample_name("::entries/x")
