"""Normalization, dimensionality reduction, and clustering.

Same steps as the chr19 mechanics notebook's Sections 4-5, factored out so
the perturbation case study and any future dataset can reuse them without
copy-pasting notebook cells into pipeline code.
"""

from __future__ import annotations

import anndata as ad
import scanpy as sc

from .utils import get_logger

logger = get_logger(__name__)


def normalize_and_reduce(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    n_comps: int = 30,
) -> ad.AnnData:
    """Library-size normalize, log-transform, select HVGs, and run PCA."""
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_top_genes = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    logger.info("Highly variable genes: %d of %d", int(adata.var["highly_variable"].sum()), adata.n_vars)

    n_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
    return adata


def cluster(
    adata: ad.AnnData,
    n_pcs: int = 15,
    n_neighbors: int = 30,
    resolution: float = 0.5,
) -> ad.AnnData:
    """Build the neighbor graph, run Leiden clustering, and compute UMAP.

    Unlike the chr19 notebook, this data has enough real per-cell signal
    (full transcriptome, thousands of counts/cell) that fragmentation into
    dozens of disconnected components isn't expected — worth checking
    anyway rather than assuming, since it's a two-line diagnostic.
    """
    from scipy.sparse.csgraph import connected_components

    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    n_components, _ = connected_components(adata.obsp["connectivities"], directed=False)
    logger.info("Neighbor graph connected components: %d", n_components)

    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)
    logger.info("Leiden clusters: %d", adata.obs["leiden"].nunique())

    sc.tl.umap(adata)
    return adata
