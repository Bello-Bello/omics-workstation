"""Differential expression and dose-response summarization.

The chr19 notebook's marker-gene section asked "which genes define this
cluster vs every other cell" (Wilcoxon, cluster vs rest). A perturbation
screen asks the related but distinct question "which genes changed in THIS
drug's treated cells vs THIS drug's own control" — same underlying test,
different groupby/reference semantics, which is why this isn't just reusing
the notebook's marker-gene cell directly.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd
import scanpy as sc

from .utils import get_logger

logger = get_logger(__name__)


def differential_expression(
    adata: ad.AnnData,
    reference: str = "control",
    groupby: str = "condition",
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """Run condition-vs-control DE for every non-control condition at once.

    `sc.tl.rank_genes_groups` with an explicit `reference` runs each group
    against that one reference rather than against "all other cells" (the
    marker-gene default) — the correct comparison for a perturbation screen,
    where "all other cells" would incorrectly include cells from every other
    drug too.
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, reference=reference, method=method)

    groups = [g for g in adata.obs[groupby].unique() if g != reference]
    frames = []
    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group)
        df["condition"] = group
        frames.append(df)

    result = pd.concat(frames, ignore_index=True)
    logger.info("DE computed for %d conditions vs %r reference", len(groups), reference)
    return result


def n_significant_by_dose(
    de_by_dose: dict[tuple[str, float], pd.DataFrame],
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """Count significant DE genes per (drug, dose), for a dose-response curve.

    Takes a pre-computed dict rather than raw AnnData + a loop internally:
    dose-stratified DE (re-running rank_genes_groups once per dose bin) is
    expensive enough that callers should be able to cache/parallelize it
    (e.g. one Nextflow task per dose) rather than this function forcing a
    specific, non-parallel loop structure.
    """
    rows = []
    for (drug, dose), df in de_by_dose.items():
        n_sig = int(((df["pvals_adj"] < padj_threshold) & (df["logfoldchanges"].abs() > lfc_threshold)).sum())
        rows.append({"drug": drug, "dose": dose, "n_significant_genes": n_sig})
    return pd.DataFrame(rows).sort_values(["drug", "dose"]).reset_index(drop=True)
