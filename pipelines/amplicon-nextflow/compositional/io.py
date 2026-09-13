"""Reading QIIME2 exports and the sample metadata, and shaping them for analysis.

A `.qza` is a zip file with a UUID-named directory inside; `qiime tools
export` unpacks the payload. That payload is BIOM for a feature table, which
`biom convert` turns into a TSV with two header lines. These readers exist
so the rest of the package never has to know any of that.

Orientation note: QIIME2 writes features as rows and samples as columns.
Every function here returns the transpose — samples as rows — because that
is what the transforms and tests expect, and because getting it backwards
produces a CLR that is silently meaningless rather than an error.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .utils import get_logger

logger = get_logger(__name__)

# Greengenes/SILVA rank prefixes, in order.
RANK_PREFIXES = {
    "kingdom": "k__", "phylum": "p__", "class": "c__",
    "order": "o__", "family": "f__", "genus": "g__", "species": "s__",
}


def read_feature_table(path: str | Path) -> pd.DataFrame:
    """Read a `biom convert`-produced TSV into a samples x features frame.

    The file starts with a `# Constructed from biom file` comment line, then
    a header row whose first column is literally `#OTU ID` — a name QIIME2
    keeps for backwards compatibility even though the rows are ASVs, not
    OTUs, in any DADA2-based workflow.
    """
    table = pd.read_csv(path, sep="\t", skiprows=1, index_col=0)
    table.index.name = "feature_id"
    table = table.T
    table.index.name = "sample_id"
    logger.info("Feature table: %d samples x %d features, %d total reads",
                table.shape[0], table.shape[1], int(table.values.sum()))
    return table


def read_taxonomy(path: str | Path) -> pd.DataFrame:
    """Read the classifier's output: Feature ID -> Taxon string, Confidence."""
    tax = pd.read_csv(path, sep="\t", index_col=0)
    tax.index.name = "feature_id"
    logger.info("Taxonomy: %d features classified, median confidence %.3f",
                tax.shape[0], tax["Confidence"].median() if "Confidence" in tax else float("nan"))
    return tax


def read_metadata(path: str | Path) -> pd.DataFrame:
    """Read a QIIME2 sample metadata TSV, dropping the `#q2:types` directive row.

    That second row declares each column as categorical or numeric to QIIME2.
    Pandas would otherwise read it as a data row and coerce every numeric
    column to object dtype — a failure that shows up much later as a
    confusing comparison error, so it is worth handling at the door.
    """
    meta = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    meta = meta[~meta.index.astype(str).str.startswith("#")]
    meta.index.name = "sample_id"
    for col in meta.columns:
        converted = pd.to_numeric(meta[col], errors="coerce")
        if converted.notna().sum() > 0.8 * len(meta):
            meta[col] = converted
    logger.info("Metadata: %d samples x %d columns", meta.shape[0], meta.shape[1])
    return meta


def align(table: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict both frames to the samples they share, in a common order.

    Demultiplexing routinely drops samples entirely — a barcode that picked
    up too few reads produces no per-sample file at all — so the feature
    table is usually a subset of the metadata sheet. Intersecting explicitly
    (rather than letting pandas broadcast NaNs) means the sample count in
    the report is the number of samples actually analysed.
    """
    shared = table.index.intersection(metadata.index)
    dropped_table = set(table.index) - set(shared)
    dropped_meta = set(metadata.index) - set(shared)
    if dropped_table:
        logger.warning("%d sample(s) in the table have no metadata: %s",
                       len(dropped_table), sorted(dropped_table)[:5])
    if dropped_meta:
        logger.info("%d metadata sample(s) absent from the table (dropped in demux/denoise)",
                    len(dropped_meta))
    if len(shared) == 0:
        raise ValueError("No overlap between feature table and metadata sample IDs")
    return table.loc[shared], metadata.loc[shared]


def collapse_to_rank(table: pd.DataFrame, taxonomy: pd.DataFrame, rank: str = "genus") -> pd.DataFrame:
    """Sum ASV counts into taxonomic bins at the requested rank.

    Collapsing is a real analytical choice, not just tidying. ASVs are exact
    sequences, so two ASVs differing by one base are separate features even
    when they are the same organism for every practical purpose; that splits
    a real signal across columns and costs power. Collapsing to genus pools
    them. The cost is resolution — strain-level differences vanish — so this
    pipeline reports ASV-level results too rather than only the collapsed view.

    Unclassified features at the requested rank are binned into a single
    `Unassigned_<rank>` column rather than dropped, so read totals (and
    therefore every ratio downstream) stay intact.
    """
    if rank not in RANK_PREFIXES:
        raise ValueError(f"Unknown rank {rank!r}; expected one of {list(RANK_PREFIXES)}")
    prefix = RANK_PREFIXES[rank]

    def extract(taxon: str) -> str:
        for part in str(taxon).split(";"):
            part = part.strip()
            if part.startswith(prefix):
                name = part[len(prefix):].strip()
                return name if name else f"Unassigned_{rank}"
        return f"Unassigned_{rank}"

    labels = taxonomy.reindex(table.columns)["Taxon"].fillna("").map(extract)
    collapsed = table.T.groupby(labels.values).sum().T
    collapsed.index.name = "sample_id"
    n_unassigned = int((labels == f"Unassigned_{rank}").sum())
    logger.info("Collapsed %d ASVs -> %d %s-level bins (%d ASVs unassigned at this rank)",
                table.shape[1], collapsed.shape[1], rank, n_unassigned)
    return collapsed


def filter_features(table: pd.DataFrame, min_prevalence: float = 0.1,
                    min_total_count: int = 10) -> pd.DataFrame:
    """Drop features seen in too few samples, or with too few reads overall.

    Rare features are mostly zeros, and every zero has to be imputed before
    the CLR. Keep thousands of near-empty columns and the multiplicative
    replacement — not the data — starts driving the geometric mean, which
    shifts every CLR value in the sample.

    This does change the composition, and under Aitchison geometry that is
    permitted: the transform is subcompositionally coherent, so distances
    computed from the retained parts remain valid statements about those
    parts. What it is *not* is free — the results describe the filtered
    subcomposition, which is why the thresholds are parameters and are
    recorded in the run summary.
    """
    n_samples = table.shape[0]
    prevalence = (table > 0).sum(axis=0) / n_samples
    totals = table.sum(axis=0)
    keep = (prevalence >= min_prevalence) & (totals >= min_total_count)
    filtered = table.loc[:, keep]
    logger.info(
        "Feature filter: %d -> %d features (prevalence >= %.0f%% of %d samples, "
        "total count >= %d); retained %.1f%% of reads",
        table.shape[1], filtered.shape[1], 100 * min_prevalence, n_samples,
        min_total_count, 100 * filtered.values.sum() / table.values.sum(),
    )
    if filtered.shape[1] < 2:
        raise ValueError("Fewer than 2 features survived filtering — thresholds are too strict")
    return filtered


def drop_empty_samples(table: pd.DataFrame, min_depth: int = 1) -> pd.DataFrame:
    """Remove samples whose remaining depth is too low to form a composition."""
    depths = table.sum(axis=1)
    keep = depths >= min_depth
    if (~keep).any():
        logger.warning("Dropping %d sample(s) below depth %d: %s",
                       int((~keep).sum()), min_depth, list(table.index[~keep]))
    return table.loc[keep]
