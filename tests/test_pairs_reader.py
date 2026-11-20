"""Pairs files are read as their content says, header lines wherever they are.

The compression of a pairs file comes from its first bytes, not from its name,
and every line that starts with ``#`` is a header line: pairs files
concatenated with ``cat`` read as the files given one after the other.
"""

from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path
from typing import cast

import cooler
import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

from metacooler.io.create import (
    _HeaderLinesBelow,
    _HeaderLinesSkipped,
    _scan_pairs_header,
    _skip_header_lines,
    from_pairs,
)

LENGTHS = {"c1": 5_000, "c2": 3_000, "c3": 2_000}

HEADER = (
    "## pairs format v1.0\n"
    + "".join(f"#chromsize: {contig} {length}\n" for contig, length in LENGTHS.items())
    + "#samheader: @SQ\tSN:c1\tLN:5000\n"
    + "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"
)
#: A header line with as many tab-separated fields as a pair (7).
RG_LINE = "#samheader: @RG\tID:L1\tSM:S\tLB:L\tPL:ILLUMINA\tPU:U\tCN:C\n"


def _pairs(seed: int, n: int = 300) -> str:
    rng = np.random.default_rng(seed)
    names = list(LENGTHS)
    rows = []
    for i in range(n):
        a, b = names[rng.integers(3)], names[rng.integers(3)]
        pos_a = rng.integers(1, LENGTHS[a])
        pos_b = rng.integers(1, LENGTHS[b])
        rows.append(f"r{seed}_{i}\t{a}\t{pos_a}\t{b}\t{pos_b}\t+\t-\n")
    return "".join(rows)


def _write(path: Path, text: str, codec: str | None = None) -> Path:
    with pa.output_stream(str(path), compression=codec) as out:
        out.write(text.encode())
    return path


def _map(path: Path, pairs: Path | list[Path]) -> Path:
    from_pairs(path, binsize=1_000, pairs=pairs, ps=False)
    return path


def _pixels(path: Path) -> pd.DataFrame:
    return cast(pd.DataFrame, cooler.Cooler(str(path)).pixels()[:])


@pytest.mark.parametrize(
    ("codec", "name"),
    [
        ("gzip", "lib.pairs"),
        ("bz2", "lib.txt"),
        ("zstd", "lib.pairs"),
        ("lz4", "lib.pairs"),
        (None, "lib.pairs.gz"),
        (None, "lib.pairs.zst"),
    ],
)
def test_the_compression_is_read_from_the_content(
    tmp_path: Path, codec: str | None, name: str
):
    """A file named otherwise than its compression reads as the plain file.

    By its suffix, a gzip file named ``lib.pairs`` was parsed as raw
    compressed bytes, and a plain file named ``lib.pairs.gz`` failed to
    decompress.
    """
    text = HEADER + _pairs(1)
    expected = _map(tmp_path / "plain.metacool", _write(tmp_path / "plain.pairs", text))
    other = tmp_path / "other"
    other.mkdir()
    result = _map(tmp_path / "result.metacool", _write(other / name, text, codec))

    pd.testing.assert_frame_equal(_pixels(result), _pixels(expected))


def test_a_zip_archive_is_refused_by_name(tmp_path: Path):
    """pyarrow cannot read one: the error says so, and what to do."""
    plain = _write(tmp_path / "lib.pairs", HEADER + _pairs(1))
    archive = tmp_path / "lib.pairs.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.write(plain, "lib.pairs")

    with pytest.raises(ValueError, match="is a zip archive"):
        _map(tmp_path / "out.metacool", archive)


@pytest.mark.parametrize(
    "header", [HEADER, HEADER + RG_LINE], ids=["header", "header-like-a-pair"]
)
@pytest.mark.parametrize("codec", [None, "gzip"])
def test_concatenated_pairs_files_read_as_the_files_one_after_the_other(
    tmp_path: Path, header: str, codec: str | None
):
    """``cat a.pairs b.pairs``, whose second header sits below the top.

    pyarrow parsed that header as pairs and failed on its first line. A header
    line with as many fields as a pair (the ``@RG`` line) is the case pyarrow's
    row handler never sees: it is left out before parsing. ``cat`` of two gzip
    files is a valid gzip file, of two members.
    """
    first = _write(tmp_path / "a.pairs", header + _pairs(1), codec)
    second = _write(tmp_path / "b.pairs", header + _pairs(2), codec)
    together = tmp_path / "ab.pairs"
    together.write_bytes(first.read_bytes() + second.read_bytes())

    expected = _map(tmp_path / "separate.metacool", [first, second])
    result = _map(tmp_path / "concatenated.metacool", together)

    pd.testing.assert_frame_equal(_pixels(result), _pixels(expected))
    assert cooler.Cooler(str(result)).info["sum"] == 600


@pytest.mark.parametrize(
    "header", [HEADER, HEADER + RG_LINE], ids=["header", "header-like-a-pair"]
)
def test_a_malformed_pair_is_still_an_error(tmp_path: Path, header: str):
    """Only ``#`` lines are skipped: a row with too few fields still stops the run.

    A guard: this held before header lines below the top were skipped.
    """
    pairs = _write(tmp_path / "lib.pairs", header + _pairs(1) + "r_bad\tc1\t10\n")

    with pytest.raises(pa.ArrowInvalid, match="Expected 7 columns, got 3"):
        _map(tmp_path / "out.metacool", pairs)


@pytest.mark.parametrize(
    ("head", "what"),
    [
        (b"PK\x05\x06" + bytes(18), "an empty zip archive"),
        (b"\xfd7zXZ\x00\x00", "xz-compressed"),
        (b"\x02\x21\x4c\x18\x00\x00", "a legacy LZ4 file"),
        (b"7z\xbc\xaf\x27\x1c", "a 7z archive"),
    ],
)
def test_other_formats_pyarrow_cannot_read_are_refused_by_name(
    tmp_path: Path, head: bytes, what: str
):
    """They were parsed as text, and failed with an error about columns."""
    pairs = tmp_path / "lib.pairs"
    pairs.write_bytes(head + b"rest of the file\n")
    with pytest.raises(ValueError, match=f"is {what}, which metacooler cannot read"):
        _map(tmp_path / "out.metacool", pairs)


def _concatenated(tmp_path: Path, first: str, second: str) -> tuple[Path, Path, Path]:
    """Two pairs files, and their concatenation."""
    a = _write(tmp_path / "a.pairs", first)
    b = _write(tmp_path / "b.pairs", second)
    together = tmp_path / "ab.pairs"
    together.write_bytes(a.read_bytes() + b.read_bytes())
    return a, b, together


def test_a_quote_in_a_header_line_below_the_top_swallows_no_pair(tmp_path: Path):
    """A field starting with a quote made pyarrow read every line up to the
    next quote as one field: the pairs after it were lost, without a word."""
    a, b, together = _concatenated(
        tmp_path,
        HEADER + _pairs(1),
        HEADER + '#samheader: @CO\t"unclosed\n' + _pairs(2),
    )
    result = _map(tmp_path / "concatenated.metacool", together)
    expected = _map(tmp_path / "separate.metacool", [a, b])
    pd.testing.assert_frame_equal(_pixels(result), _pixels(expected))
    assert cooler.Cooler(str(result)).info["sum"] == 600


def test_a_read_name_starting_with_a_quote_is_read_as_it_is(tmp_path: Path):
    """Pairs are read without quoting."""
    body = _pairs(1).replace("r1_0\t", '"r1_0\t', 1)
    pairs = _write(tmp_path / "a.pairs", HEADER + body)
    assert cooler.Cooler(str(_map(tmp_path / "out.metacool", pairs))).info["sum"] == 300


@pytest.mark.parametrize(
    "header", [HEADER, HEADER + RG_LINE], ids=["header", "header-like-a-pair"]
)
def test_lines_of_white_space_are_blank_lines(tmp_path: Path, header: str):
    """At the top, as in a header below it, a line of spaces or tabs only is
    a blank line: skipped. With a header line like a pair, the line at the top
    became the first row once the filter had left the header out."""
    _, _, together = _concatenated(
        tmp_path, "  \n" + header + _pairs(1), header + " \t \n" + _pairs(2)
    )
    result = _map(tmp_path / "concatenated.metacool", together)
    assert cooler.Cooler(str(result)).info["sum"] == 600


def test_a_header_line_like_a_pair_below_the_top_alone_is_named(tmp_path: Path):
    """When only a header below the top holds a line with a pair's fields, it
    is parsed as a pair: the error says what it may be, and what to do."""
    _, _, together = _concatenated(
        tmp_path, HEADER + _pairs(1), HEADER + RG_LINE + _pairs(2)
    )
    with pytest.raises(pa.ArrowInvalid, match="give the pairs files one by one"):
        _map(tmp_path / "out.metacool", together)


@pytest.mark.parametrize(
    ("chromsize", "shown"),
    [
        ("#chromsize: c1 4000\n", "c1: 5000 at the top, 4000 below"),
        ("#chromsize: c9 900\n", "c9: none at the top, 900 below"),
    ],
)
def test_headers_of_other_assemblies_are_refused(
    tmp_path: Path, chromsize: str, shown: str
):
    """The contig lengths come from the top header: a header below it that
    gives another length, or another contig, is refused rather than read with
    the top header's lengths."""
    _, _, together = _concatenated(
        tmp_path, HEADER + _pairs(1), "## pairs format v1.0\n" + chromsize + _pairs(2)
    )
    with pytest.raises(ValueError, match=re.escape(f"disagree ({shown})")):
        _map(tmp_path / "out.metacool", together)


@pytest.mark.parametrize(
    ("length", "refused"), [("05000", False), ("5000.0", True), ("five", True)]
)
def test_lengths_below_the_top_are_compared_as_numbers(
    tmp_path: Path, length: str, refused: bool
):
    """``05000`` is the length the top header gives as ``5000``, as the
    header's own parser reads it; a length that is not an integer is not."""
    _, _, together = _concatenated(
        tmp_path,
        HEADER + _pairs(1),
        f"## pairs format v1.0\n#chromsize: c1 {length}\n" + _pairs(2),
    )
    if refused:
        with pytest.raises(ValueError, match=f"c1: 5000 at the top, {length} below"):
            _map(tmp_path / "out.metacool", together)
    else:
        result = _map(tmp_path / "out.metacool", together)
        assert cooler.Cooler(str(result)).info["sum"] == 600


# pyarrow reports the row it could not decode as an exception ignored in the
# handler, then fails the read, which metacooler names.
@pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning")
def test_a_header_line_below_the_top_that_is_not_text_is_named(tmp_path: Path):
    """pyarrow decodes a row as UTF-8 for the row handler: a header line below
    the top that is not (a Latin-1 byte) never reaches it, and the read
    stopped with pyarrow's "Expected 7 columns, got 2" alone."""
    a = _write(tmp_path / "a.pairs", HEADER + _pairs(1))
    b = tmp_path / "b.pairs"
    b.write_bytes(
        HEADER.encode()
        + "#samheader: @CO\tcaf\xe9\n".encode("latin-1")
        + _pairs(2).encode()
    )
    together = tmp_path / "ab.pairs"
    together.write_bytes(a.read_bytes() + b.read_bytes())
    with pytest.raises(pa.ArrowInvalid, match="that is not UTF-8 text, which cannot"):
        _map(tmp_path / "out.metacool", together)


def test_chromsize_lines_below_a_top_header_without_any_are_not_checked(
    tmp_path: Path,
):
    """With no lengths at the top (they come from elsewhere), there is nothing
    to contradict. A guard: the check of the headers below the top must not
    refuse this."""
    plain = "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n"
    _, _, together = _concatenated(tmp_path, plain + _pairs(1), HEADER + _pairs(2))
    bins = pd.DataFrame(
        {"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())}
    )
    out = tmp_path / "out.metacool"
    from_pairs(out, bins=bins, pairs=together, ps=False)
    assert cooler.Cooler(str(out)).info["sum"] == 600


def test_the_header_scan_says_whether_a_header_line_looks_like_a_pair(
    tmp_path: Path,
):
    plain = _write(tmp_path / "plain.pairs", HEADER + _pairs(1))
    like = _write(tmp_path / "like.pairs", HEADER + RG_LINE + _pairs(1))
    empty = _write(tmp_path / "empty.pairs", HEADER + RG_LINE)

    assert _scan_pairs_header(plain) == (6, True, False)
    assert _scan_pairs_header(like) == (7, True, True)
    assert _scan_pairs_header(empty) == (7, False, False)


@pytest.mark.parametrize("block_size", [1, 3, 7, 1 << 20])
@pytest.mark.parametrize("nbytes", [-1, 1, 5, 1 << 20])
def test_header_lines_are_left_out_at_any_block_size(
    tmp_path: Path, block_size: int, nbytes: int
):
    """Lines that straddle blocks, a block of header lines only, a last line
    with no newline, and a header line as the last line of the file. The top
    header is left out as it is, a line of spaces included; the header lines
    below it are noted."""
    text = (
        b"#h1\n  \n#h2\tx\nr1\tc1\t1\n"
        + b"#mid\n" * 3
        + b"r2\tc2\t"
        + b"9" * 20
        + b"\n\n#end-header\nr3\tc3\t3"
    )
    below = _HeaderLinesBelow(tmp_path / "unread.pairs")
    stream = _HeaderLinesSkipped(
        io.BytesIO(text + b"\n#trailing"), 3, below, block_size
    )

    out = b""
    while chunk := stream.read(nbytes):
        out += chunk
    stream.close()

    assert out == b"r1\tc1\t1\nr2\tc2\t" + b"9" * 20 + b"\n\nr3\tc3\t3\n"
    assert below.count == 5
    assert stream.closed
    assert stream.readable()
    assert not stream.seekable()


def test_a_top_header_that_ends_the_file_without_a_newline_is_left_out(
    tmp_path: Path,
):
    """Its last line is the top header's, not a header line below it: it was
    left out of the pairs, but counted as one."""
    below = _HeaderLinesBelow(tmp_path / "unread.pairs")
    stream = _HeaderLinesSkipped(io.BytesIO(b"#h1\n#h2"), 2, below)
    assert stream.read() == b""
    assert below.count == 0


def test_the_row_handler_skips_header_and_blank_lines_and_refuses_the_rest(
    tmp_path: Path,
):
    """pyarrow calls it from its own threads, out of coverage's sight."""
    from types import SimpleNamespace

    below = _HeaderLinesBelow(tmp_path / "unread.pairs")
    handler = _skip_header_lines(below)
    assert handler(SimpleNamespace(text="## pairs format v1.0")) == "skip"
    assert handler(SimpleNamespace(text=" \t ")) == "skip"
    assert handler(SimpleNamespace(text="r1\tc1\t10")) == "error"
    assert below.count == 1


def test_a_chromsize_line_with_no_length_is_not_checked(tmp_path: Path):
    below = _HeaderLinesBelow(tmp_path / "unread.pairs")
    below.note("#chromsize: c1")
    below.check()
    assert below.count == 1


def test_a_length_that_is_not_an_integer_is_a_conflict(tmp_path: Path):
    """The check itself, called here as pyarrow calls the row handler, from
    its own threads, which coverage does not trace. A guard: compared as
    text, such a length was a conflict too."""
    top = _write(tmp_path / "top.pairs", HEADER + _pairs(1))
    below = _HeaderLinesBelow(top)
    below.note("#chromsize: c1 five")
    with pytest.raises(ValueError, match="c1: 5000 at the top, five below"):
        below.check()


def test_an_unreadable_top_header_leaves_nothing_to_check(tmp_path: Path):
    """A top header whose own #chromsize lines cannot be read fails where it is
    read for the lengths, not here."""
    pairs = tmp_path / "bad.pairs"
    pairs.write_text("#chromsize: c1 many\nr1\tc1\t1\tc1\t2\n")
    below = _HeaderLinesBelow(pairs)
    below.note("#chromsize: c1 100")
    below.check()


def test_many_conflicts_are_shortened(tmp_path: Path):
    pairs = _write(tmp_path / "a.pairs", HEADER + _pairs(1))
    below = _HeaderLinesBelow(pairs)
    for i in range(5):
        below.note(f"#chromsize: x{i} 10")
    with pytest.raises(ValueError, match=r"x2: none at the top, 10 below; \.\.\.\)"):
        below.check()
