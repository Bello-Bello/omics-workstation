"""Per-cell QC metrics and filtering for perturbation screens.

Same core logic as the chr19 mechanics notebook (../analysis/scanpy_analysis.ipynb)
— total_counts, n_genes_by_counts, pct_counts_mt, then filter — but exposed as
functions so the same code runs identically from Nextflow, Snakemake, and the
notebook rather than being copy-pasted into each.
"""

from __future__ import annotations

import anndata as ad
import scanpy as sc

from .utils import get_logger

logger = get_logger(__name__)


def compute_qc_metrics(adata: ad.AnnData, mt_prefix: str = "MT-") -> ad.AnnData:
    """Annotate mitochondrial genes and compute standard per-cell QC metrics.

    `mt_prefix` defaults to the human convention ("MT-", uppercase) — sciPlex
    is A549 (human), the opposite of the chr19 notebook's mouse data ("mt-",
    lowercase). Getting this wrong silently zeroes out pct_counts_mt instead
    of erroring, which is exactly the failure mode worth a named parameter
    rather than a hardcoded string.
    """
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    n_mt = int(adata.var["mt"].sum())
    logger.info("Mitochondrial genes flagged (prefix=%r): %d", mt_prefix, n_mt)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    return adata


def filter_cells_and_genes(
    adata: ad.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mt: float = 20.0,
) -> ad.AnnData:
    """Apply standard gene- and cell-level QC filters.

    Defaults (min_genes=200, max_pct_mt=20%) are the conventional
    whole-transcriptome thresholds — safe to use as-is here since sciPlex is
    a full human transcriptome, unlike the chr19-only reference in the
    mechanics notebook where those defaults had to be abandoned.
    """
    n_before = adata.n_obs
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)

    if adata.var["mt"].sum() > 0:
        adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    else:
        logger.warning("No mitochondrial genes found — skipping pct_counts_mt filter")

    logger.info(
        "Filtering: %d -> %d cells (min_genes=%d, max_pct_mt=%.1f)",
        n_before, adata.n_obs, min_genes, max_pct_mt,
    )
    return adata
