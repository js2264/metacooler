"""Test fixtures for metacooler."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from metacooler.io.create import from_pairs


@pytest.fixture(autouse=True)
def _isolate_metacooler_logger():
    """Undo ``cli.setup_logging``'s global side effects around every test.

    ``setup_logging`` reconfigures the shared ``"metacooler"`` logger: it
    replaces the handler list and sets ``propagate = False`` so a host
    application does not receive every record a second time in its own format.
    That mutation outlives the CLI invocation that caused it, so without this
    fixture any test that runs a CLI command would silently stop ``caplog`` —
    which listens on the *root* logger — from seeing records in every test that
    happens to run afterwards, making log assertions pass or fail on ordering
    alone.
    """
    logger = logging.getLogger("metacooler")
    saved = (logger.handlers[:], logger.propagate, logger.level)
    logger.handlers = []
    logger.propagate = True
    logger.setLevel(logging.NOTSET)
    yield
    logger.handlers, logger.propagate, logger.level = saved


@pytest.fixture
def tiny_fasta(tmp_path: Path) -> Path:
    """A tiny 4-contig fasta."""
    path = tmp_path / "tiny.fa"
    with open(path, "w") as f:
        f.write(">contig_A\n" + "ACGT" * 250 + "\n")  # 1000 bp
        f.write(">contig_B\n" + "TGCA" * 150 + "\n")  # 600 bp
        f.write(">contig_C\n" + "AATT" * 100 + "\n")  # 400 bp
        f.write(">contig_D\n" + "GGCC" * 50 + "\n")  # 200 bp
    return path


@pytest.fixture
def tiny_contigs_table(tmp_path: Path) -> Path:
    """A 2-MAG assignment table (contig_D unassigned). Has header."""
    path = tmp_path / "mag_table.tsv"
    with open(path, "w") as f:
        f.write("contig_name\tgenome_name\n")
        f.write("contig_A\tMAG_1\n")
        f.write("contig_B\tMAG_1\n")
        f.write("contig_C\tMAG_2\n")
        # contig_D is unassigned -> __unbinned__
    return path


@pytest.fixture
def tiny_bins_bed(tmp_path: Path) -> Path:
    """Fixed 100-bp bins for the 4 contigs."""
    path = tmp_path / "bins.bed"
    rows = []
    for contig, length in [
        ("contig_A", 1000),
        ("contig_B", 600),
        ("contig_C", 400),
        ("contig_D", 200),
    ]:
        for start in range(0, length, 100):
            end = min(start + 100, length)
            rows.append(f"{contig}\t{start}\t{end}\n")
    with open(path, "w") as f:
        f.writelines(rows)
    return path


@pytest.fixture
def tiny_pairs(tmp_path: Path) -> Path:
    """Synthetic pairs file (~200 contacts)."""
    path = tmp_path / "pairs.pairs"
    np.random.seed(42)
    contigs = ["contig_A", "contig_B", "contig_C", "contig_D"]
    lengths = {"contig_A": 1000, "contig_B": 600, "contig_C": 400, "contig_D": 200}

    rows = ["#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n"]
    for i in range(200):
        c1 = np.random.choice(contigs)
        c2 = np.random.choice(contigs)
        p1 = np.random.randint(1, lengths[c1])
        p2 = np.random.randint(1, lengths[c2])
        rows.append(f"read_{i}\t{c1}\t{p1}\t{c2}\t{p2}\t+\t+\n")

    with open(path, "w") as f:
        f.writelines(rows)
    return path


@pytest.fixture
def tiny_genomes_table(tmp_path: Path) -> Path:
    """A genomes metadata table with extra columns (with header)."""
    path = tmp_path / "genomes_table.tsv"
    with open(path, "w") as f:
        f.write("genome_name\tquality\tphylum\n")
        f.write("MAG_1\thigh_quality\tBacteroidetes\n")
        f.write("MAG_2\tmedium_quality\tFirmicutes\n")
    return path


@pytest.fixture
def tiny_metacooler(
    tmp_path: Path,
    tiny_fasta: Path,
    tiny_contigs_table: Path,
    tiny_bins_bed: Path,
    tiny_pairs: Path,
) -> Path:
    """A tiny metacool file built from fixtures."""
    out = tmp_path / "tiny.metacool"
    from_pairs(
        uri=out,
        bins=tiny_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=tiny_contigs_table,
        contig_table_format="contigs",
        fasta=tiny_fasta,
        assembly_name="tiny_test",
        contig_ordering="genome-then-alpha",
        chunksize=1000,
    )
    return out


@pytest.fixture
def tiny_cool(tmp_path: Path, tiny_pairs: Path, tiny_bins_bed: Path) -> Path:
    """A tiny .cool file created directly from the pairs and bins fixtures."""
    import cooler
    import pandas as pd

    bins = pd.read_csv(
        tiny_bins_bed, sep="\t", header=None, names=["chrom", "start", "end"]
    )
    pairs = pd.read_csv(
        tiny_pairs,
        sep="\t",
        comment="#",
        header=None,
        names=["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"],
    )

    def _assign_bin_ids(chrom_col: pd.Series, pos_col: pd.Series) -> np.ndarray:
        result = np.full(len(chrom_col), -1, dtype=np.int64)
        for chrom, grp in bins.groupby("chrom"):
            mask = (chrom_col == chrom).values
            if not mask.any():
                continue
            starts = grp["start"].values
            ends = grp["end"].values
            gids = grp.index.values
            positions = pos_col.values[mask]
            idx = np.searchsorted(starts, positions, side="right") - 1
            clip_idx = np.clip(idx, 0, len(ends) - 1)
            valid = (idx >= 0) & (positions < ends[clip_idx])
            result[mask] = np.where(valid, gids[np.where(valid, idx, 0)], -1)
        return result

    pairs["bin1_id"] = _assign_bin_ids(pairs["chrom1"], pairs["pos1"])
    pairs["bin2_id"] = _assign_bin_ids(pairs["chrom2"], pairs["pos2"])
    pairs = pairs[(pairs["bin1_id"] >= 0) & (pairs["bin2_id"] >= 0)].copy()

    # Enforce upper triangle
    swap = pairs["bin1_id"] > pairs["bin2_id"]
    pairs.loc[swap, ["bin1_id", "bin2_id"]] = pairs.loc[
        swap, ["bin2_id", "bin1_id"]
    ].values

    pixels = (
        pairs.groupby(["bin1_id", "bin2_id"], sort=True)
        .size()
        .reset_index(name="count")
    )

    out = tmp_path / "tiny.cool"
    cooler.create_cooler(str(out), bins=bins, pixels=pixels, ordered=True)
    return out


@pytest.fixture
def tiny_mmetacool(tmp_path: Path, tiny_metacooler: Path) -> Path:
    """A 2-entry .mmetacool from the tiny_metacooler fixture."""
    from metacooler.io.mmetacool import create_mmetacool

    out = tmp_path / "tiny.mmetacool"
    create_mmetacool(out, {"entry1": tiny_metacooler, "entry2": tiny_metacooler})
    return out


@pytest.fixture
def tiny_cool_reversed_order(tmp_path: Path) -> tuple[Path, Path]:
    """A tiny cool + contigs table where cool contig order differs from genome-then-length order.

    Cool order:    contig_D (200bp), contig_C (400bp), contig_B (600bp), contig_A (1000bp)
    MAG ordering:  MAG_1 = [contig_A, contig_B], MAG_2 = [contig_C], unbinned = [contig_D]
    Genome-then-length reorder: contig_A (1000), contig_B (600), contig_C (400), contig_D (200)

    Using reversed bins exercises the chrom-ID remapping in ``io.create.from_cool``.
    """
    import cooler
    import pandas as pd

    contig_lengths = {
        "contig_D": 200,
        "contig_C": 400,
        "contig_B": 600,
        "contig_A": 1000,
    }

    # Bins in REVERSED order (contig_D first, contig_A last)
    rows = []
    for contig in ["contig_D", "contig_C", "contig_B", "contig_A"]:
        length = contig_lengths[contig]
        for start in range(0, length, 100):
            end = min(start + 100, length)
            rows.append({"chrom": contig, "start": start, "end": end})
    bins = pd.DataFrame(rows)

    # Synthetic pixels (a few contacts per contig pair)
    np.random.seed(0)
    bin_ids = bins.index.tolist()
    pixel_rows = []
    for _ in range(300):
        b1 = int(np.random.choice(bin_ids))
        b2 = int(np.random.choice(bin_ids))
        if b1 > b2:
            b1, b2 = b2, b1
        pixel_rows.append({"bin1_id": b1, "bin2_id": b2})
    pixels = (
        pd.DataFrame(pixel_rows)
        .groupby(["bin1_id", "bin2_id"])
        .size()
        .reset_index(name="count")  # ty:ignore[no-matching-overload]
    )

    cool_path = tmp_path / "reversed.cool"
    cooler.create_cooler(str(cool_path), bins=bins, pixels=pixels, ordered=True)

    # Contigs table: contig_A and contig_B -> MAG_1, contig_C -> MAG_2
    contigs_table = tmp_path / "reversed_contigs.tsv"
    with open(contigs_table, "w") as f:
        f.write("contig_name\tgenome_name\n")
        f.write("contig_A\tMAG_1\n")
        f.write("contig_B\tMAG_1\n")
        f.write("contig_C\tMAG_2\n")

    return cool_path, contigs_table


@pytest.fixture
def tiny_contigs_contig_genome(tmp_path: Path) -> Path:
    """Contig-to-genome table in ``contig_genome`` format (no header, 2 cols)."""
    path = tmp_path / "contigs_contig_genome.tsv"
    with open(path, "w") as f:
        f.write("contig_A\tMAG_1\n")
        f.write("contig_B\tMAG_1\n")
        f.write("contig_C\tMAG_2\n")
    return path


@pytest.fixture
def tiny_contigs_metator(tmp_path: Path) -> Path:
    """Contig-to-genome table in ``metator`` format (header: Name, Final_bin)."""
    path = tmp_path / "contigs_metator.tsv"
    with open(path, "w") as f:
        f.write("Name\tFinal_bin\tsize\n")
        f.write("contig_A\tMAG_1\t1000\n")
        f.write("contig_B\tMAG_1\t600\n")
        f.write("contig_C\tMAG_2\t400\n")
    return path


@pytest.fixture
def tiny_contigs_4col(tmp_path: Path) -> Path:
    """4-column contig table (header: id, contig_name, extra, genome_name).

    Suitable for ``--contigs-format 2,4`` (1-based column indices).
    """
    path = tmp_path / "contigs_4col.tsv"
    with open(path, "w") as f:
        f.write("id\tcontig_name\textra\tgenome_name\n")
        f.write("1\tcontig_A\tabc\tMAG_1\n")
        f.write("2\tcontig_B\tdef\tMAG_1\n")
        f.write("3\tcontig_C\tghi\tMAG_2\n")
    return path


@pytest.fixture
def tiny_variable_bin_metacooler(
    tmp_path: Path,
    tiny_fasta: Path,
    tiny_contigs_table: Path,
    tiny_pairs: Path,
) -> Path:
    """A metacooler with non-uniform (variable) bin sizes.

    Bins are built with irregular fragment-like sizes so that
    ``cooler.Cooler.binsize`` returns ``None``.
    """
    import pandas as pd

    # Build irregular bins: two large + one small fragment per contig
    rows = []
    for contig, length in [
        ("contig_A", 1000),
        ("contig_B", 600),
        ("contig_C", 400),
        ("contig_D", 200),
    ]:
        rows += [
            {"chrom": contig, "start": 0, "end": length // 3},
            {"chrom": contig, "start": length // 3, "end": 2 * length // 3},
            {"chrom": contig, "start": 2 * length // 3, "end": length},
        ]
    bins_df = pd.DataFrame(rows)

    var_bins_bed = tmp_path / "var_bins.bed"
    bins_df.to_csv(var_bins_bed, sep="\t", header=False, index=False)

    out = tmp_path / "tiny_varbin.metacool"
    from_pairs(
        uri=out,
        bins=var_bins_bed,
        pairs=tiny_pairs,
        contig_genome_map=tiny_contigs_table,
        contig_table_format="contigs",
        fasta=tiny_fasta,
        assembly_name="tiny_varbin_test",
        contig_ordering="genome-then-alpha",
        chunksize=1000,
    )
    return out


@pytest.fixture
def tiny_genomes_wrong_col_name(tmp_path: Path) -> Path:
    """Genomes table where first column is named 'name' instead of 'genome_name'."""
    path = tmp_path / "genomes_wrong_col.tsv"
    with open(path, "w") as f:
        f.write("name\tquality\n")
        f.write("MAG_1\thigh\n")
        f.write("MAG_2\tmedium\n")
    return path


@pytest.fixture
def tiny_genomes_no_header(tmp_path: Path) -> Path:
    """Genomes data file with no header row (first row is data, not column names)."""
    path = tmp_path / "genomes_no_header.tsv"
    with open(path, "w") as f:
        f.write("MAG_1\thigh\n")
        f.write("MAG_2\tmedium\n")
    return path


# NOTE: a ``tiny_fasta_valid`` fixture used to live here.  It was byte-identical
# to ``tiny_fasta`` and its docstring described a pyfaidx "Duplicate key" bug
# that does not exist — ``tiny_fasta`` already writes one single-line record per
# contig.  Use ``tiny_fasta``.
