"""Dimensionality reduction and clustering on CLR-transformed data.

Both of the unsupervised methods the role description names, applied to
compositional data the way compositional data requires.

The important idea: once CLR has moved the samples into real space,
ordinary Euclidean methods become correct rather than merely convenient.
A PCA on CLR values is a *principal component analysis of log-ratios* —
what the literature calls an Aitchison biplot — and it is equivalent to
PCoA on the Aitchison distance matrix, but with the loadings kept, so you
can read which taxa drive each axis.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

from .utils import get_logger

logger = get_logger(__name__)


@dataclass
class PCAResult:
    coordinates: np.ndarray          # samples x n_components
    explained_variance_ratio: np.ndarray
    loadings: np.ndarray             # features x n_components
    feature_ids: list[str]

    def top_loadings(self, component: int = 0, n: int = 10):
        """Features contributing most to a component, by absolute loading."""
        order = np.argsort(-np.abs(self.loadings[:, component]))[:n]
        return [(self.feature_ids[i], float(self.loadings[i, component])) for i in order]


def pca(clr_values: np.ndarray, feature_ids: list[str], n_components: int = 5) -> PCAResult:
    """PCA via SVD on column-centered CLR values.

    SVD rather than eigendecomposition of the covariance matrix, for a
    reason specific to this data: CLR rows sum to zero by construction, so
    the covariance matrix is rank-deficient by exactly one and its smallest
    eigenvalue should be 0. Eigensolvers can return a small negative number
    there and then a square root produces NaN. SVD on the centered matrix
    sidesteps the issue — it never forms the covariance matrix at all.

    The last component is therefore structurally empty; `n_components` is
    capped below the rank so it is never returned as if it meant something.
    """
    values = np.asarray(clr_values, dtype=float)
    n_samples, n_features = values.shape
    centered = values - values.mean(axis=0, keepdims=True)

    max_components = min(n_samples - 1, n_features - 1)
    n_components = int(min(n_components, max_components))
    if n_components < 1:
        raise ValueError(f"Not enough samples/features for a PCA ({n_samples}x{n_features})")

    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    explained = (s ** 2) / (s ** 2).sum()

    logger.info("PCA: PC1=%.1f%%, PC2=%.1f%% of variance (%d components kept)",
                100 * explained[0], 100 * explained[1] if explained.size > 1 else float("nan"),
                n_components)
    return PCAResult(
        coordinates=u[:, :n_components] * s[:n_components],
        explained_variance_ratio=explained[:n_components],
        loadings=vt[:n_components].T,
        feature_ids=list(feature_ids),
    )


def hierarchical_clustering(distances: np.ndarray, method: str = "average"):
    """Hierarchical clustering of samples on the Aitchison distance matrix.

    Returns (linkage_matrix, leaf_order).

    `average` (UPGMA) by default rather than `ward`: Ward's criterion
    minimises within-cluster variance and is defined for Euclidean distances
    specifically. Aitchison distance *is* Euclidean in CLR space, so Ward is
    not wrong here — but the linkage is computed from the condensed distance
    matrix, where SciPy cannot verify that, and average linkage makes no
    such assumption. Choosing the method that stays valid if someone later
    swaps in Bray-Curtis costs nothing today.
    """
    distances = np.asarray(distances, dtype=float)
    # squareform requires exact symmetry; floating-point error in the
    # distance computation can leave d[i,j] and d[j,i] differing in the
    # last bit, which SciPy rejects outright.
    symmetric = (distances + distances.T) / 2.0
    np.fill_diagonal(symmetric, 0.0)
    condensed = squareform(symmetric, checks=False)

    link = linkage(condensed, method=method)
    order = leaves_list(link)
    logger.info("Hierarchical clustering: %d samples, %s linkage", distances.shape[0], method)
    return link, order


def cut_clusters(link, n_clusters: int) -> np.ndarray:
    """Cut the dendrogram into `n_clusters` flat groups."""
    return fcluster(link, t=n_clusters, criterion="maxclust")


def cluster_group_agreement(cluster_labels: np.ndarray, group_labels: np.ndarray) -> float:
    """Adjusted Rand index between unsupervised clusters and a known grouping.

    Implemented directly (it is four counts and a formula) to give the report
    a single number for "did the clustering rediscover the study design?".
    Chance-corrected: 0 is what random labels score, 1 is a perfect match,
    and negative values mean worse than chance.
    """
    from itertools import combinations

    cluster_labels = np.asarray(cluster_labels)
    group_labels = np.asarray(group_labels)
    n = cluster_labels.size

    same_cluster_same_group = same_cluster = same_group = 0
    for i, j in combinations(range(n), 2):
        sc = cluster_labels[i] == cluster_labels[j]
        sg = group_labels[i] == group_labels[j]
        same_cluster += sc
        same_group += sg
        same_cluster_same_group += sc and sg

    n_pairs = n * (n - 1) / 2
    expected = same_cluster * same_group / n_pairs
    maximum = (same_cluster + same_group) / 2
    if maximum == expected:
        return 0.0
    return float((same_cluster_same_group - expected) / (maximum - expected))
