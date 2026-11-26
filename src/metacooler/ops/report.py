"""The report of a map: what it holds, and where its contacts land.

:func:`collect` gathers the report of one map, a ``.metacool`` file or an
entry of a collection, as chromokit's report kit takes it (CONVENTIONS.md
§ 12): sections, each a key and its value per line, and tables, with the
plots of its page. :mod:`metacooler.io.report` writes it.

- **run**: how the map was made: what wrote it, when, its assembly, its bins
  and its contig ordering; **commands**: the command lines that wrote it, in
  order (its ``metadata["commands"]``).
- **genomes**, **contigs**, **bins**, **contacts**: what it holds.
- **qc**: where its contacts land, in the four categories of
  :mod:`metacooler.qc`: the block the map stores (``metacooler qc --store``,
  or ``create pairs``), or one computed from its pixels when it stores none
  of this version.
- tables **genomes** (each genome, with its contacts by category) and **ps**
  (the contact probability by distance, the map's contigs pooled).

The report's format is :data:`FORMAT_VERSION`: a key, a column, a unit or a
MultiQC column that changes or goes raises it; one that comes does not.
``docs/format.md`` documents every key.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd

from metacooler.schema import DEFAULT_CHUNKSIZE

if TYPE_CHECKING:
    from metacooler._report import Report, Tool
    from metacooler.api import MetaCooler

logger = logging.getLogger("metacooler")

#: The version of the format of report.txt and report_mqc.yaml.
FORMAT_VERSION = 1

#: Where the report's keys are documented.
DOCS = "https://js2264.github.io/metacooler/format/#the-report-of-a-map"

#: The page's template, which extends the kit's base.
TEMPLATE = Path(__file__).resolve().parents[1] / "templates" / "report.html.j2"

#: metacooler's columns in MultiQC's General Statistics table, which are also
#: the cards of the page's overview: (key, section, pconfig).
MQC_COLUMNS = (
    (
        "genomes",
        "genomes",
        {
            "title": "Genomes",
            "description": "Genomes of the map, its unbinned contigs aside",
            "format": "{:,.0f}",
            "min": 0,
        },
    ),
    (
        "contigs",
        "contigs",
        {
            "title": "Contigs",
            "description": "Contigs of the map",
            "format": "{:,.0f}",
            "min": 0,
        },
    ),
    (
        "contacts",
        "contacts",
        {
            "title": "Contacts",
            "description": "Contacts of the map: the sum of its pixels' counts",
            "format": "{:,.0f}",
            "min": 0,
        },
    ),
    (
        "inter_genome_pct",
        "qc",
        {
            "title": "Inter-genome",
            "description": "Share of the contacts between two genomes",
            "suffix": "%",
            "min": 0,
            "max": 100,
            "format": "{:,.2f}",
        },
    ),
)

#: The columns of the genomes table, from the map's /genomes and its QC.
GENOME_COLUMNS = (
    "name",
    "length_bp",
    "contigs",
    "bins",
    "intra_contig_short",
    "intra_contig_long",
    "inter_contig_intra_genome",
    "inter_genome",
    "intra_genome_ratio",
)

#: The columns of the P(s) table.
PS_COLUMNS = ("distance_bp", "contacts", "bin_pairs", "probability")


class Collected(NamedTuple):
    """The report of a map, and what its page takes besides.

    Attributes:
        report: The report: its sections and tables.
        figures: The page's plots, by the id of their place in it.
        inputs: What the page's header names: the map.
    """

    report: Report
    figures: dict
    inputs: dict


def tool() -> Tool:
    """metacooler, as its reports name it: its version, the report's format,
    its MultiQC columns and its page's template."""
    from metacooler._report import Tool
    from metacooler._version import __version__

    return Tool(
        name="metacooler",
        version=__version__,
        format_version=FORMAT_VERSION,
        columns=MQC_COLUMNS,
        description="metacooler maps: their genomes, contigs and contacts",
        docs=DOCS,
        template=TEMPLATE,
    )


def n50(lengths: np.ndarray) -> tuple[int, int]:
    """The N50 and the L50 of *lengths*: the length of the contig that, the
    contigs taken longest first, brings their sum to half the total, and how
    many contigs that takes. ``(0, 0)`` for none."""
    ordered = np.sort(np.asarray(lengths, dtype=np.int64))[::-1]
    if not len(ordered) or not ordered.sum():
        return 0, 0
    position = int(np.argmax(np.cumsum(ordered) * 2 >= ordered.sum()))
    return int(ordered[position]), position + 1


def _number(value: object) -> int | float | None:
    """A value of the map as a Python number: None for a missing one (None,
    NaN), an integer as an int."""
    if isinstance(value, (int, float, np.integer, np.floating, str, bytes)):
        try:
            number = float(value)
        except ValueError:
            return None
    else:
        return None
    if math.isnan(number):
        return None
    return int(number) if number.is_integer() else number


def _pct(part: float | None, whole: float | None) -> float | None:
    """*part* as a percent of *whole*; None when either is missing or whole is 0."""
    if part is None or not whole:
        return None
    return 100 * part / whole


def _qc(mc: MetaCooler, chunksize: int) -> tuple[dict, pd.DataFrame, bool]:
    """The map's QC block, its per-genome table, and whether the map stored
    them: the block of this version the map stores with its ``/genomes``
    columns, or else one computed from its pixels, which is not stored."""
    from metacooler.qc import QC_GENOME_COLUMNS, QC_METADATA_KEY, QC_VERSION

    block = mc.metadata.get(QC_METADATA_KEY)
    genomes = mc.genomes
    if (
        isinstance(block, dict)
        and block.get("version") == QC_VERSION
        and all(column in genomes.columns for column in QC_GENOME_COLUMNS)
    ):
        return block, genomes[["name", *QC_GENOME_COLUMNS]], True
    from metacooler.ops.qc import qc

    logger.info("The map stores no QC block of version %s: computing one", QC_VERSION)
    result = qc(mc.uri, chunksize=chunksize)
    return result.summary, result.per_genome, False


def _ps(mc: MetaCooler, chunksize: int) -> pd.DataFrame:
    """The contact probability of the map's contigs, pooled, by distance."""
    from metacooler.ops.ps import contact_probability

    try:
        table = contact_probability(mc.uri, chunksize=chunksize)
    except ValueError as error:
        logger.warning("No P(s) for this map: %s", error)
        return pd.DataFrame(columns=pd.Index(PS_COLUMNS))
    # Distances to the base pair, probabilities to 6 significant digits.
    probability = table["prob"].to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "distance_bp": np.round(table["dist"].to_numpy(dtype=float)).astype(
                np.int64
            ),
            "contacts": table["n_contacts"].to_numpy(),
            "bin_pairs": table["n_bin_pairs"].to_numpy(),
            "probability": [float(f"{value:.6g}") for value in probability],
        }
    )


def _total_contacts(mc: MetaCooler, chunksize: int) -> int | float:
    """The sum of the counts of the map's pixels: its ``sum`` attribute, as
    cooler writes it, or else the sum of its pixels, read chunk by chunk."""
    total = mc.info.get("sum")
    if total is not None:
        return _number(total) or 0
    selector = mc.pixels()
    count = 0
    for start in range(0, int(mc.info["nnz"]), chunksize):
        count += selector[start : start + chunksize]["count"].sum()
    return _number(count) or 0


def collect(uri: str | Path, *, chunksize: int = DEFAULT_CHUNKSIZE) -> Collected:
    """Gather the report of a map.

    Args:
        uri: ``file.metacool`` path or ``file.mmetacool::entries/name`` URI.
        chunksize: Pixels read per iteration, where the QC block or the P(s)
            is computed.

    Returns:
        A :class:`Collected`: the report (see the module's docstring for its
        sections and tables), the plots of its page, and its header's inputs.
    """
    from metacooler._report import Report, report_section
    from metacooler.api import MetaCooler
    from metacooler.qc import QC_CATEGORIES, QC_COLUMN_PREFIX, UNBINNED_GENOME
    from metacooler.uri import sample_name

    mc = MetaCooler(uri)
    info = mc.info
    genomes = mc.genomes
    contigs = mc.contigs
    names = genomes["name"].astype(str).to_numpy()
    unbinned = names == UNBINNED_GENOME
    lengths = contigs["length"].to_numpy(dtype=np.int64)
    in_unbinned = np.isin(contigs["genome_id"].to_numpy(), np.flatnonzero(unbinned))
    contigs_bp = int(lengths.sum())
    binned_bp = int(lengths[~in_unbinned].sum())
    contig_n50, contig_l50 = n50(lengths)

    block, per_genome, stored = _qc(mc, chunksize)
    total = _number(block.get("n_contacts"))
    commands = mc.metadata.get("commands")

    sections = {
        "report": report_section(tool(), sample_name(uri), {"map": str(uri)}),
        "run": {
            "generated_by": info.get("generated-by"),
            "created": info.get("creation-date"),
            "assembly": info.get("assembly"),
            "bin_type": info.get("bin-type"),
            "bin_size_bp": _number(info.get("bin-size")),
            "contig_ordering": info.get("ordering"),
            "storage_mode": info.get("storage-mode"),
        },
    }
    if isinstance(commands, list) and commands:
        sections["commands"] = {
            f"command_{number}": str(line) for number, line in enumerate(commands, 1)
        }
    sections["genomes"] = {
        "genomes": int((~unbinned).sum()),
        "genomes_bp": int(genomes.loc[~unbinned, "length"].sum()),
        "unbinned_contigs": int(in_unbinned.sum()),
        "unbinned_bp": contigs_bp - binned_bp,
        "binned_pct": _pct(binned_bp, contigs_bp),
    }
    sections["contigs"] = {
        "contigs": len(contigs),
        "contigs_bp": contigs_bp,
        "contig_n50_bp": contig_n50,
        "contig_l50": contig_l50,
        "longest_contig_bp": int(lengths.max()) if len(lengths) else None,
    }
    sections["bins"] = {"bins": _number(info.get("nbins"))}
    sections["contacts"] = {
        "contacts": _total_contacts(mc, chunksize),
        "pixels": _number(info.get("nnz")),
    }
    sections["qc"] = {
        "stored": stored,
        "source": block.get("source"),
        "distance_threshold_bp": _number(block.get("distance_threshold")),
        "classified_contacts": total,
        **{category: _number(block.get(category)) for category in QC_CATEGORIES},
        "trans_contig_ratio": _number(block.get("trans_contig_ratio")),
        "intra_genome_pct": _pct(_number(block.get("intra_genome")), total),
        "inter_genome_pct": _pct(_number(block.get("inter_genome")), total),
    }

    table = pd.DataFrame(
        {
            "name": names,
            "length_bp": genomes["length"].to_numpy(),
            "contigs": genomes["n_contigs"].to_numpy(),
            "bins": genomes["n_bins"].to_numpy(),
        }
    )
    qc_columns = per_genome.set_index(per_genome["name"].astype(str))
    for category in (*QC_CATEGORIES, "intra_genome_ratio"):
        column = qc_columns[f"{QC_COLUMN_PREFIX}{category}"]
        table[category] = column.reindex(names).to_numpy()
    ps = _ps(mc, chunksize)
    report = Report(
        sections=sections, tables={"genomes": table[list(GENOME_COLUMNS)], "ps": ps}
    )
    return Collected(
        report=report,
        figures=figures(report, lengths),
        inputs={"Map": str(uri)},
    )


#: The four categories of the QC, as the page names them.
QC_LABELS = {
    "intra_contig_short": "Within a contig, short range",
    "intra_contig_long": "Within a contig, long range",
    "inter_contig_intra_genome": "Between contigs of a genome",
    "inter_genome": "Between genomes",
}


#: What every plot's layout starts from: its margins, its height, and the
#: closest point under the pointer (metaTOR's).
_LAYOUT = {"margin": {"l": 60, "r": 20, "t": 30, "b": 50}, "hovermode": "closest"}
#: A logarithmic axis, ticked at each power of ten.
_LOG = {"type": "log", "dtick": 1, "exponentformat": "power"}


def _layout(height: int, **axes: object) -> dict:
    """A plot's layout: :data:`_LAYOUT`, *height*, and its axes."""
    return {**_LAYOUT, "height": height, **axes}


def _length_steps(lengths: np.ndarray) -> tuple[list, list]:
    """The contigs' lengths as a histogram over 24 bins, evenly spaced on a
    logarithmic scale from the shortest to the longest: the bins' edges, and
    the number of contigs from each edge on (the last repeated, for a step
    line that ends on the last edge)."""
    if not len(lengths):
        return [], []
    low, high = max(float(lengths.min()), 1.0), float(lengths.max())
    if high <= low:
        low, high = low * 0.9, low * 1.1
    edges = np.logspace(math.log10(low), math.log10(high), 25)
    counts, _ = np.histogram(lengths, bins=edges)
    return edges.tolist(), [*counts.tolist(), int(counts[-1])]


def figures(report: Report, lengths: np.ndarray) -> dict:
    """The plots of the page, as Plotly figures (data and layout), by the id
    of their place: where the contacts land (``plot-qc``), the genomes'
    lengths against their intra-genome ratio (``plot-genomes``), the contigs'
    lengths (``plot-contigs``) and the P(s) (``plot-ps``). A colour left out
    follows the page's palette."""
    qc = report.sections["qc"]
    total = qc.get("classified_contacts") or 0
    shares = [100 * (qc.get(key) or 0) / total if total else 0 for key in QC_LABELS]
    genomes = report.tables["genomes"]
    genomes = genomes[genomes["name"] != "__unbinned__"]
    ps = report.tables["ps"]
    edges, steps = _length_steps(lengths)
    return {
        "plot-qc": {
            "data": [
                {
                    "type": "bar",
                    "orientation": "h",
                    "x": shares,
                    "y": list(QC_LABELS.values()),
                    "text": [f"{share:.1f}%" for share in shares],
                    "textposition": "auto",
                    "hovertemplate": "%{y}: %{x:.2f}% of the contacts<extra></extra>",
                }
            ],
            "layout": _layout(
                300,
                margin={"l": 200, "r": 20, "t": 30, "b": 50},
                xaxis={
                    "title": {"text": "Share of the contacts (%)"},
                    "range": [0, 100],
                },
                yaxis={"autorange": "reversed"},
            ),
        },
        "plot-genomes": {
            "data": [
                {
                    "type": "scatter",
                    "mode": "markers",
                    "x": genomes["length_bp"].tolist(),
                    "y": genomes["intra_genome_ratio"].tolist(),
                    "text": genomes["name"].tolist(),
                    "customdata": genomes["contigs"].tolist(),
                    "hovertemplate": (
                        "<b>%{text}</b><br>%{x:,} bp, %{customdata} contigs"
                        "<br>intra-genome ratio %{y:.3f}<extra></extra>"
                    ),
                }
            ],
            "layout": _layout(
                300,
                xaxis={"title": {"text": "Length (bp)"}, "type": "log"},
                yaxis={"title": {"text": "Intra-genome ratio"}, "range": [0, 1.05]},
            ),
        },
        "plot-contigs": {
            "data": [
                {
                    "type": "scatter",
                    "mode": "lines",
                    "line": {"shape": "hv"},
                    "fill": "tozeroy",
                    "x": edges,
                    "y": steps,
                    "hovertemplate": "%{y} contigs from %{x:,.0f} bp<extra></extra>",
                }
            ],
            "layout": _layout(
                300,
                xaxis={"title": {"text": "Length (bp)"}, **_LOG},
                yaxis={"title": {"text": "Contigs"}, "rangemode": "tozero"},
            ),
        },
        "plot-ps": {
            "data": [
                {
                    "type": "scatter",
                    "mode": "lines+markers",
                    "x": ps["distance_bp"].tolist(),
                    "y": ps["probability"].tolist(),
                    "hovertemplate": "%{x:,.0f} bp: %{y:.3g}<extra></extra>",
                }
            ],
            "layout": _layout(
                380,
                xaxis={"title": {"text": "Distance (bp)"}, **_LOG},
                yaxis={"title": {"text": "Contact probability"}, **_LOG},
            ),
        },
    }
