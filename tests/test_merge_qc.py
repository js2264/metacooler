"""``merge`` sums the inputs' QC instead of dropping it.

The guarantee: ``from_pairs`` on library A and on library B, over the same
bins and contig-to-genome map, then ``merge``, stores the very map that
``from_pairs`` on A and B together stores -- the same pixels and bins, the
same ``/contigs`` and ``/genomes`` with their ``qc_*`` columns, and the same QC
block; only provenance differs.  That holds because every count of the QC adds
up over libraries -- the four categories, in total, per genome and per contig,
and the pairs read and discarded -- and every fraction and ratio is recomputed
from the sums by the code that computed each library's.

The libraries below have all of it: intra-contig pairs below and above the
distance threshold, inter-contig pairs within and across genomes, unbinned
contigs, and pairs ``create pairs`` discards.

When the inputs' blocks were not counted alike, or one has none, the merged
map has no QC block, as before, and a warning says why.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from metacooler.api import MetaCooler
from metacooler.cli.main import cli
from metacooler.io.create import from_pairs
from metacooler.io.merge import merge
from metacooler.ops.qc import qc
from metacooler.qc import QC_CONTIG_COLUMNS, QC_GENOME_COLUMNS, QC_METADATA_KEY

# Distinct lengths, so the genome-then-length layout is fully determined.
LENGTHS = {
    "c1": 9_000,
    "c2": 8_000,
    "c3": 7_000,
    "c4": 6_000,
    "c5": 5_000,
    "c6": 4_000,
    "c7": 3_000,
    "c8": 2_000,
}
THRESHOLD = 2_500
# c7 and c8 are left unbinned.
ASSIGNMENT = {
    "c1": "MAG_1",
    "c2": "MAG_1",
    "c3": "MAG_2",
    "c4": "MAG_2",
    "c5": "MAG_2",
    "c6": "MAG_3",
}
# Pairs per library that create pairs discards: an unknown contig, and a
# position past the end of its contig.
DISCARDED = 3

# One bin per contig, as a map built to be binned by MetaTOR; and 1 kb bins,
# finer than the threshold.
BINNINGS = ("one-bin-per-contig", "1kb")
ONE_BIN_PER_CONTIG = pd.DataFrame(
    {"chrom": list(LENGTHS), "start": 0, "end": list(LENGTHS.values())}
)


def _library(path: Path, seed: int) -> Path:
    """One library: random pairs over every QC category, plus the discarded."""
    rng = np.random.default_rng(seed)
    names = list(LENGTHS)
    rows = []
    for _ in range(1_500):
        c1 = names[rng.integers(len(names))]
        c2 = c1 if rng.random() < 0.6 else names[rng.integers(len(names))]
        rows.append(
            (
                c1,
                int(rng.integers(1, LENGTHS[c1])),
                c2,
                int(rng.integers(1, LENGTHS[c2])),
            )
        )
    rows += [("cX", 10, "c1", 20)] * 2
    rows += [("c1", LENGTHS["c1"] + 50, "c2", 10)]
    path.write_text(
        "## pairs format v1.0\n"
        + "".join(f"#chromsize: {c} {n}\n" for c, n in LENGTHS.items())
        + "#columns: readID chrom1 pos1 chrom2 pos2\n"
        + "".join(f"r{i}\t{a}\t{b}\t{c}\t{d}\n" for i, (a, b, c, d) in enumerate(rows))
    )
    return path


@pytest.fixture(scope="module")
def libraries(tmp_path_factory: pytest.TempPathFactory) -> list[Path]:
    """Three libraries of one assembly."""
    directory = tmp_path_factory.mktemp("libraries")
    return [
        _library(directory / f"lib{i}.pairs", seed) for i, seed in enumerate((1, 2, 3))
    ]


def _create(
    path: Path,
    pairs: Path | list[Path],
    assignment: dict[str, str] | None = ASSIGNMENT,
    binning: str = "1kb",
    **kwargs,
) -> Path:
    one_bin = binning == "one-bin-per-contig"
    kwargs.setdefault("qc_distance_threshold", THRESHOLD)
    from_pairs(
        path,
        bins=ONE_BIN_PER_CONTIG if one_bin else None,
        binsize=None if one_bin else 1_000,
        pairs=pairs,
        contig_genome_map=assignment,
        ps=False,
        **kwargs,
    )
    return path


def _contents(path: Path) -> dict:
    """Every dataset, attribute and metadata key of a map, bar its provenance.

    Provenance is the ``creation-date`` attribute and the ``commands`` list.
    ``pixels/count`` is compared by value only: ``merge`` sums counts in
    ``int64`` where ``create pairs`` stores ``int32``.
    """
    out: dict = {}
    with h5py.File(path, "r") as f:

        def visit(name: str, obj: object) -> None:
            if isinstance(obj, h5py.Dataset):
                values = obj[()].tolist()
                out[name] = values if name == "pixels/count" else (obj.dtype, values)

        f.visititems(visit)
        out["attrs"] = {
            key: value.tolist() if isinstance(value, np.generic) else value
            for key, value in f.attrs.items()
            if key not in ("creation-date", "metadata")
        }
        metadata = json.loads(f.attrs["metadata"])
    metadata.pop("commands", None)
    out["metadata"] = metadata
    return out


def _assert_same_map(result: Path, expected: Path) -> None:
    got, want = _contents(result), _contents(expected)
    assert sorted(got) == sorted(want)
    for key in want:
        assert got[key] == want[key], key


def _assert_no_qc(path: Path) -> None:
    mc = MetaCooler(path)
    assert QC_METADATA_KEY not in mc.metadata
    assert not [c for c in mc.genomes.columns if c.startswith("qc_")]
    assert not set(QC_CONTIG_COLUMNS) & set(mc.contigs.columns)


# ---------------------------------------------------------------------------
# merge(from_pairs(A), from_pairs(B)) == from_pairs([A, B])
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("binning", BINNINGS)
@pytest.mark.parametrize(
    "assignment", [ASSIGNMENT, None], ids=["binned", "no-binning-yet"]
)
def test_merged_libraries_are_the_map_of_both_at_once(
    tmp_path: Path,
    libraries: list[Path],
    assignment: dict[str, str] | None,
    binning: str,
):
    """The proof: with a genome map, and with none (every contig in
    ``__unbinned__``, each counted as a genome of its own); at one bin per
    contig, and at 1 kb."""
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a, assignment, binning)
    second = _create(tmp_path / "b.metacool", b, assignment, binning)
    both = _create(tmp_path / "both.metacool", [a, b], assignment, binning)
    out = tmp_path / "merged.metacool"

    merge([first, second], out)

    _assert_same_map(out, both)
    block = MetaCooler(out).metadata[QC_METADATA_KEY]
    assert block["source"] == "pairs"
    assert block["n_pairs_discarded"] == 2 * DISCARDED
    for key in ("intra_contig_short", "intra_contig_long", "inter_genome"):
        assert block[key] > 0, key
    assert (block["inter_contig_intra_genome"] > 0) == (assignment is not None)
    with h5py.File(out, "r") as f:
        assert f["pixels/count"].dtype == np.int64


def test_three_libraries(tmp_path: Path, libraries: list[Path]):
    maps = [_create(tmp_path / f"{p.stem}.metacool", p) for p in libraries]
    every = _create(tmp_path / "every.metacool", list(libraries))
    out = tmp_path / "merged.metacool"

    merge(maps, out)

    _assert_same_map(out, every)


@pytest.mark.parametrize("binning", BINNINGS)
def test_blocks_measured_on_the_pixels_add_up_too(
    tmp_path: Path, libraries: list[Path], binning: str
):
    """``qc --store`` on each library, then merge: ``qc --store`` on both."""
    a, b = libraries[:2]
    maps = [
        _create(tmp_path / "a.metacool", a, binning=binning),
        _create(tmp_path / "b.metacool", b, binning=binning),
        _create(tmp_path / "both.metacool", [a, b], binning=binning),
    ]
    for path in maps:
        qc(path, distance_threshold=THRESHOLD, store=True)
    out = tmp_path / "merged.metacool"

    merge(maps[:2], out)

    _assert_same_map(out, maps[2])
    assert MetaCooler(out).metadata[QC_METADATA_KEY]["source"] == "pixels"


def _scale_counts_to_float(path: Path, factor: float) -> None:
    """Rewrite ``pixels/count`` as float64 times *factor* (a rescaled map)."""
    with h5py.File(path, "a") as f:
        counts = f["pixels"]["count"][:].astype(np.float64) * factor
        del f["pixels"]["count"]
        f["pixels"].create_dataset("count", data=counts)
        f.attrs["sum"] = float(counts.sum())


def test_float_counts_add_up_to_what_qc_measures_on_the_merged_map(
    tmp_path: Path, libraries: list[Path]
):
    """Up to rounding: float addition does not associate, so a sum of sums
    need not have the bits of one sum."""
    maps = [_create(tmp_path / f"{p.stem}.metacool", p) for p in libraries[:2]]
    for path, factor in zip(maps, (0.3, 0.7), strict=True):
        _scale_counts_to_float(path, factor)
        qc(path, distance_threshold=THRESHOLD, store=True)
    out = tmp_path / "merged.metacool"

    merge(maps, out)

    mc = MetaCooler(out)
    block = mc.metadata[QC_METADATA_KEY]
    assert isinstance(block["n_contacts"], float)
    measured = qc(out, distance_threshold=THRESHOLD)
    floats = [k for k, v in measured.summary.items() if isinstance(v, float)]
    assert {k: v for k, v in block.items() if k not in floats} == {
        k: v for k, v in measured.summary.items() if k not in floats
    }
    assert [block[k] for k in floats] == pytest.approx(
        [measured.summary[k] for k in floats], rel=1e-9, abs=2e-6
    )
    np.testing.assert_allclose(
        mc.genomes[list(QC_GENOME_COLUMNS)].to_numpy(),
        measured.per_genome[list(QC_GENOME_COLUMNS)].to_numpy(),
        rtol=1e-9,
    )
    assert mc.contigs[list(QC_CONTIG_COLUMNS)].dtypes.tolist() == [np.float64] * 2


def test_float_maps_reassigned_then_merged_or_merged_then_reassigned(
    tmp_path: Path, libraries: list[Path]
):
    """Both orders keep the QC, and give the same block and columns.

    ``reassign`` keeps a block measured on pairs and re-measures its two
    inter-contig categories from the pixels, which a rescaled map holds as
    floats: the block then counts in floats, some of them fractional, and
    ``merge`` sums it as it is.  Counts scaled by halves stay exact in
    binary, so the two orders give the same bits.
    """
    from metacooler.io.reassign import reassign

    a, b = libraries[:2]
    maps = [_create(tmp_path / "a.metacool", a), _create(tmp_path / "b.metacool", b)]
    for path, factor in zip(maps, (0.5, 1.5), strict=True):
        _scale_counts_to_float(path, factor)
    regrouped = {**ASSIGNMENT, "c2": "MAG_2", "c7": "MAG_3"}
    reassigned = [path.with_name(f"{path.stem}.reassigned.metacool") for path in maps]
    for path, out in zip(maps, reassigned, strict=True):
        reassign(path, regrouped, out)
    then_merged = tmp_path / "reassigned-then-merged.metacool"
    merged = tmp_path / "merged.metacool"
    merged_then = tmp_path / "merged-then-reassigned.metacool"

    merge(reassigned, then_merged)
    merge(maps, merged)
    reassign(merged, regrouped, merged_then)

    blocks = [MetaCooler(path).metadata[QC_METADATA_KEY] for path in reassigned]
    assert all(block["source"] == "pairs" for block in blocks)
    assert any(block["inter_genome"] % 1 for block in blocks)
    _assert_same_qc(then_merged, merged_then)
    block = MetaCooler(then_merged).metadata[QC_METADATA_KEY]
    assert block["source"] == "pairs"
    assert block["n_pairs_read"] == 2 * 1_503


def test_entries_of_a_collection(tmp_path: Path, libraries: list[Path]):
    from metacooler.io.mmetacool import create_mmetacool

    a, b = libraries[:2]
    collection = tmp_path / "libraries.mmetacool"
    create_mmetacool(
        collection,
        {
            "a": _create(tmp_path / "a.metacool", a),
            "b": _create(tmp_path / "b.metacool", b),
        },
    )
    both = _create(tmp_path / "both.metacool", [a, b])
    out = tmp_path / "merged.metacool"

    merge([f"{collection}::entries/a", f"{collection}::entries/b"], out)

    _assert_same_map(out, both)


def test_the_cli_sums_the_blocks(tmp_path: Path, libraries: list[Path]):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    both = _create(tmp_path / "both.metacool", [a, b])
    out = tmp_path / "merged.metacool"

    result = CliRunner().invoke(cli, ["merge", str(first), str(second), "-o", str(out)])

    assert result.exit_code == 0, result.output
    _assert_same_map(out, both)
    assert "merge" in MetaCooler(out).metadata["commands"][-1]


# ---------------------------------------------------------------------------
# Blocks that do not add up are dropped, as before, and the warning says why
# ---------------------------------------------------------------------------


def _assert_dropped(
    caplog: pytest.LogCaptureFixture,
    out: Path,
    expected: Path,
    reason: str,
    *,
    on_pairs: bool = True,
) -> None:
    """*out* is *expected* without QC, and one warning gave *reason*.

    The warning names ``qc --store``, and ``create pairs`` first when a block
    was measured on pairs, since ``qc --store`` cannot give back what only
    read positions do.
    """
    record = next(r for r in caplog.records if "Dropping the inputs' QC" in r.message)
    assert record.levelname == "WARNING"
    message = record.getMessage()
    assert reason in message
    assert f"'metacooler qc --store {out}'" in message
    assert ("give 'metacooler create pairs' every library's" in message) == on_pairs
    _assert_no_qc(out)
    _assert_same_map(out, expected)


@pytest.mark.parametrize(
    ("second_kwargs", "reason"),
    [
        ({"qc": False}, "{b} has none"),
        (
            {"qc_distance_threshold": 1_000},
            "not counted alike: distance_threshold is 2500 in {a} but 1000 in {b}",
        ),
    ],
    ids=["no-block", "another-threshold"],
)
def test_blocks_that_do_not_add_up_are_dropped(
    tmp_path: Path,
    libraries: list[Path],
    caplog: pytest.LogCaptureFixture,
    second_kwargs: dict,
    reason: str,
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b, **second_kwargs)
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    _assert_dropped(caplog, out, expected, reason.format(a=first, b=second))


def test_blocks_measured_on_pixels_point_to_qc_store_only(
    tmp_path: Path, libraries: list[Path], caplog: pytest.LogCaptureFixture
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a, qc=False)
    second = _create(tmp_path / "b.metacool", b, qc=False)
    qc(first, distance_threshold=THRESHOLD, store=True)
    qc(second, distance_threshold=1_000, store=True)
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    _assert_dropped(
        caplog,
        out,
        expected,
        f"distance_threshold is 2500 in {first} but 1000 in {second}",
        on_pairs=False,
    )


def test_blocks_from_another_source_are_dropped(
    tmp_path: Path, libraries: list[Path], caplog: pytest.LogCaptureFixture
):
    """A block measured on read positions and one between bin midpoints."""
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    qc(second, distance_threshold=THRESHOLD, store=True)
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    _assert_dropped(
        caplog, out, expected, f"source is 'pairs' in {first} but 'pixels' in {second}"
    )
    # So are the pairs counts that only a block from pairs records.
    assert f"n_pairs_read is 1503 in {first} but absent in {second}" in caplog.text


def _qc_of(path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """A map's QC block and its ``qc_*`` columns."""
    mc = MetaCooler(path)
    return (
        mc.metadata[QC_METADATA_KEY],
        mc.genomes[["name", *QC_GENOME_COLUMNS]],
        mc.contigs[["name", *QC_CONTIG_COLUMNS]],
    )


def _assert_same_qc(result: Path, expected: Path) -> None:
    got, want = _qc_of(result), _qc_of(expected)
    assert got[0] == want[0]
    pd.testing.assert_frame_equal(got[1], want[1])
    pd.testing.assert_frame_equal(got[2], want[2])


def test_bins_given_by_size_or_as_a_bed_add_up(
    tmp_path: Path, libraries: list[Path], caplog: pytest.LogCaptureFixture
):
    """The same 1 kb bins, once tiled (``binsize`` recorded) and once given.

    ``binsize`` changes no count -- merge has already checked the bins are the
    same -- so the blocks add up, and the merged one, without ``binsize``, is
    that of a map of both libraries made with the bins given.
    """
    a, b = libraries[:2]
    starts = [(c, s) for c, n in LENGTHS.items() for s in range(0, n, 1_000)]
    grid = pd.DataFrame(
        {
            "chrom": [c for c, _ in starts],
            "start": [s for _, s in starts],
            "end": [min(s + 1_000, LENGTHS[c]) for c, s in starts],
        }
    )

    def from_grid(path: Path, pairs: Path | list[Path]) -> Path:
        from_pairs(
            path,
            bins=grid,
            pairs=pairs,
            contig_genome_map=ASSIGNMENT,
            qc_distance_threshold=THRESHOLD,
            ps=False,
        )
        return path

    first = _create(tmp_path / "a.metacool", a)
    second = from_grid(tmp_path / "b.metacool", b)
    direct = from_grid(tmp_path / "both.metacool", [a, b])
    assert MetaCooler(first).metadata[QC_METADATA_KEY]["binsize"] == 1_000
    assert "binsize" not in MetaCooler(direct).metadata[QC_METADATA_KEY]
    out = tmp_path / "merged.metacool"
    caplog.clear()

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    assert not caplog.records
    _assert_same_qc(out, direct)


_DELETE = object()


@pytest.mark.parametrize(
    ("change", "binsize"),
    [({"binsize": 2_000}, _DELETE), ({"binsize": _DELETE}, _DELETE), ({}, 1_000)],
    ids=["another-value", "absent-in-one", "the-same"],
)
def test_binsize_is_kept_only_when_every_input_records_the_same(
    tmp_path: Path, libraries: list[Path], change: dict, binsize: object
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    with h5py.File(second, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        for key, value in change.items():
            if value is _DELETE:
                del metadata[QC_METADATA_KEY][key]
            else:
                metadata[QC_METADATA_KEY][key] = value
        f.attrs["metadata"] = json.dumps(metadata)
    direct = _create(tmp_path / "both.metacool", [a, b])
    out = tmp_path / "merged.metacool"

    merge([first, second], out)

    expected = dict(MetaCooler(direct).metadata[QC_METADATA_KEY])
    if binsize is _DELETE:
        del expected["binsize"]
    assert MetaCooler(out).metadata[QC_METADATA_KEY] == expected


def test_a_key_merge_does_not_know_is_carried_once(
    tmp_path: Path, libraries: list[Path]
):
    """Equal in every input, as a later release might record: kept, not summed.

    A key that adds up over libraries comes with a new QC_VERSION, which this
    release's merge does not sum (see QC_VERSION).
    """
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    for path in (first, second):
        _edit_block(path, x_later=7)
    direct = _create(tmp_path / "both.metacool", [a, b])
    out = tmp_path / "merged.metacool"

    merge([first, second], out)

    block = MetaCooler(out).metadata[QC_METADATA_KEY]
    assert block == {**MetaCooler(direct).metadata[QC_METADATA_KEY], "x_later": 7}


@pytest.mark.parametrize(
    ("first_value", "reason"),
    [
        (7, "x_later is 7 in {a} but 8 in {b}"),
        (_DELETE, "x_later is absent in {a} but 8 in {b}"),
    ],
    ids=["differs", "absent-in-one"],
)
def test_a_key_merge_does_not_know_drops_the_qc_unless_the_inputs_agree(
    tmp_path: Path,
    libraries: list[Path],
    caplog: pytest.LogCaptureFixture,
    first_value: object,
    reason: str,
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    if first_value is not _DELETE:
        _edit_block(first, x_later=first_value)
    _edit_block(second, x_later=8)
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    _assert_dropped(caplog, out, expected, reason.format(a=first, b=second))


def test_maps_without_blocks_merge_without_a_word(
    tmp_path: Path, libraries: list[Path], caplog: pytest.LogCaptureFixture
):
    """Guard, which passes without the change too: with no block anywhere,
    nothing is dropped and nothing is said about it."""
    a, b = libraries[:2]
    maps = [
        _create(tmp_path / "a.metacool", a, qc=False),
        _create(tmp_path / "b.metacool", b, qc=False),
    ]
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.INFO, logger="metacooler"):
        merge(maps, out)

    assert "QC" not in caplog.text
    _assert_same_map(out, expected)


# ---------------------------------------------------------------------------
# Blocks that cannot be read back
# ---------------------------------------------------------------------------


def _edit_block(path: Path, **changes: object) -> None:
    with h5py.File(path, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        metadata[QC_METADATA_KEY].update(changes)
        f.attrs["metadata"] = json.dumps(metadata)


def _set_not_an_object(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        metadata = json.loads(f.attrs["metadata"])
        metadata[QC_METADATA_KEY] = "not a block"
        f.attrs["metadata"] = json.dumps(metadata)


def _drop_genome_columns(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        del f["genomes"]["qc_inter_genome"]


def _set_version(path: Path) -> None:
    _edit_block(path, version=1)


def _set_pairs_count_text(path: Path) -> None:
    _edit_block(path, n_pairs_read="many")


def _set_pairs_count_fractional(path: Path) -> None:
    _edit_block(path, n_pairs_read=1503.5)


def _set_pairs_count_negative(path: Path) -> None:
    _edit_block(path, n_pairs_discarded=-5)


def _set_category_negative(path: Path) -> None:
    _edit_block(path, inter_genome=-100)


def _set_metadata_not_json(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        f.attrs["metadata"] = "{not json"


def _drop_metadata(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        del f.attrs["metadata"]


def _rewrite_column(path: Path, column: str, change: Callable) -> None:
    """Replace ``/genomes`` or ``/contigs`` column *column* with change(values)."""
    with h5py.File(path, "r+") as f:
        table, name = column.split("/")
        values = change(f[table][name][:])
        del f[table][name]
        f[table].create_dataset(name, data=values)


def _shift(values: np.ndarray, by: int) -> np.ndarray:
    values[0] += by
    return values


def _make_negative(values: np.ndarray) -> np.ndarray:
    """Row 0 below zero, the column's total kept."""
    shift = values[0] + 1
    values[0] -= shift
    values[1] += shift
    return values


def _as_text(values: np.ndarray) -> np.ndarray:
    return np.array([str(v).encode() for v in values])


_GENOME_COUNT = "genomes/qc_inter_genome"
_CONTIG_COUNT = "contigs/qc_intra_contig_short"


@pytest.mark.parametrize(
    ("damage", "both", "reason"),
    [
        (_set_not_an_object, False, "the QC block of {b} is not a JSON object"),
        (
            _set_metadata_not_json,
            False,
            "the metadata attribute of {b} is not a JSON object",
        ),
        (_drop_metadata, False, "{b} has none"),
        (
            _drop_genome_columns,
            False,
            "{b} lacks the qc_* count columns of /genomes that go with it",
        ),
        (_set_version, False, "version is 2 in {a} but 1 in {b}"),
        (
            _set_version,
            True,
            "the QC of {a} cannot be read back (a version 1 QC block, not 2)",
        ),
        (
            _set_pairs_count_text,
            False,
            "the QC of {b} cannot be read back (not a finite number: n_pairs_read)",
        ),
        (
            _set_pairs_count_fractional,
            False,
            "the QC of {b} cannot be read back (not an integer, as a count of "
            "pairs is: n_pairs_read)",
        ),
        (
            _set_pairs_count_negative,
            False,
            "the QC of {b} cannot be read back (negative: n_pairs_discarded)",
        ),
        (
            _set_category_negative,
            False,
            "the QC of {b} cannot be read back (negative: inter_genome)",
        ),
        (
            lambda path: _rewrite_column(
                path, _GENOME_COUNT, lambda v: _shift(v, 1000)
            ),
            False,
            "the QC of {b} cannot be read back (per_genome does not add up to the "
            "block: ",
        ),
        (
            lambda path: _rewrite_column(path, _GENOME_COUNT, _make_negative),
            False,
            "the QC of {b} cannot be read back (per_genome holds a negative count)",
        ),
        (
            lambda path: _rewrite_column(
                path, _GENOME_COUNT, lambda v: v.astype(np.float64)
            ),
            False,
            "the QC of {b} cannot be read back (per_genome holds floats, where the "
            "block counts in integers)",
        ),
        (
            lambda path: _rewrite_column(path, _GENOME_COUNT, _as_text),
            False,
            "the qc_* columns of /genomes of {b} are not numbers",
        ),
    ],
    ids=[
        "not-an-object",
        "metadata-not-json",
        "no-metadata",
        "no-genome-columns",
        "versions",
        "old-version",
        "pairs-count-text",
        "pairs-count-fractional",
        "pairs-count-negative",
        "category-negative",
        "genome-column-edited",
        "genome-column-negative",
        "genome-column-float",
        "genome-column-text",
    ],
)
def test_blocks_that_cannot_be_read_back_are_dropped(
    tmp_path: Path,
    libraries: list[Path],
    caplog: pytest.LogCaptureFixture,
    damage: Callable[[Path], None],
    both: bool,
    reason: str,
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    for path in (first, second) if both else (second,):
        damage(path)
    expected = _create(tmp_path / "both.metacool", [a, b], qc=False)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.WARNING, logger="metacooler"):
        merge([first, second], out)

    _assert_dropped(caplog, out, expected, reason.format(a=first, b=second))


@pytest.mark.parametrize("outcome", ["summed", "dropped"])
def test_the_qc_is_read_before_any_pixel_is_merged(
    tmp_path: Path,
    libraries: list[Path],
    caplog: pytest.LogCaptureFixture,
    monkeypatch: pytest.MonkeyPatch,
    outcome: str,
):
    """Whether the blocks add up is known, and said, before the long step."""
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    if outcome == "dropped":
        _set_version(second)
    said = "Summed the QC blocks" if outcome == "summed" else "Dropping the inputs' QC"

    def pixel_merge(*args, **kwargs):
        assert said in caplog.text
        raise RuntimeError("stopped at the pixel merge")

    monkeypatch.setattr(cooler, "merge_coolers", pixel_merge)
    with (
        caplog.at_level(logging.INFO, logger="metacooler"),
        pytest.raises(RuntimeError, match="stopped at the pixel merge"),
    ):
        merge([first, second], tmp_path / "merged.metacool")


# ---------------------------------------------------------------------------
# Without per-contig counts, the block and the /genomes columns still add up
# ---------------------------------------------------------------------------


def _drop_contig_columns(path: Path) -> None:
    """Make *path* look like a map written by metacooler 0.2.0."""
    with h5py.File(path, "r+") as f:
        for column in QC_CONTIG_COLUMNS:
            del f["contigs"][column]


def _shift_one_count(path: Path) -> None:
    """Leave per-contig counts that no longer add up to the block."""
    with h5py.File(path, "r+") as f:
        column = f["contigs"][QC_CONTIG_COLUMNS[0]]
        values = column[:]
        values[0] += 1
        column[...] = values


def _drop_one_contig_column(path: Path) -> None:
    with h5py.File(path, "r+") as f:
        del f["contigs"][QC_CONTIG_COLUMNS[1]]


def _a_user_text_column(path: Path) -> None:
    """A 0.2.0 map with a user's column named like a QC one, stored as given.

    0.2.0 reserved no QC name in ``/contigs``, and wrote no QC column there.
    """
    _drop_contig_columns(path)
    with h5py.File(path, "r+") as f:
        n = len(f["contigs"]["name"])
        f["contigs"].create_dataset(
            QC_CONTIG_COLUMNS[0], data=np.array([b"pass", b"fail"] * n)[:n]
        )


@pytest.mark.parametrize(
    ("damage", "both", "level", "message"),
    [
        (_drop_contig_columns, True, "INFO", "No input has per-contig QC columns"),
        (
            _drop_contig_columns,
            False,
            "WARNING",
            "since {b} has none (maps written by metacooler 0.2.0 and earlier",
        ),
        (
            _shift_one_count,
            False,
            "WARNING",
            "since those of {b} do not add up to the block they go with",
        ),
        (
            _drop_one_contig_column,
            False,
            "WARNING",
            "since {b} lacks qc_intra_contig_long. Without them",
        ),
        (
            lambda path: _rewrite_column(path, _CONTIG_COUNT, _as_text),
            False,
            "WARNING",
            "since the qc_intra_contig_short of {b} is not numbers (a user's "
            "column of that name, as metacooler 0.2.0 stored one)",
        ),
        (
            _a_user_text_column,
            False,
            "WARNING",
            "since {b} lacks qc_intra_contig_long; the qc_intra_contig_short of "
            "{b} is not numbers",
        ),
    ],
    ids=[
        "none-has-them",
        "one-lacks-them",
        "one-stale",
        "one-column-missing",
        "one-column-text",
        "a-user-text-column-of-0.2.0",
    ],
)
def test_without_per_contig_counts_the_rest_still_adds_up(
    tmp_path: Path,
    libraries: list[Path],
    caplog: pytest.LogCaptureFixture,
    damage: Callable[[Path], None],
    both: bool,
    level: str,
    message: str,
):
    a, b = libraries[:2]
    first = _create(tmp_path / "a.metacool", a)
    second = _create(tmp_path / "b.metacool", b)
    for path in (first, second) if both else (second,):
        damage(path)
    direct = _create(tmp_path / "both.metacool", [a, b])
    _drop_contig_columns(direct)
    out = tmp_path / "merged.metacool"

    with caplog.at_level(logging.INFO, logger="metacooler"):
        merge([first, second], out)

    record = next(r for r in caplog.records if "per-contig" in r.message)
    assert record.levelname == level
    assert message.format(b=second) in record.getMessage()
    if level == "WARNING":
        assert "recomputes its whole QC block from the pixels" in record.getMessage()
    _assert_same_map(out, direct)
