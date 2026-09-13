#!/usr/bin/env python3
"""Cross-check the from-scratch compositional methods against reference implementations.

The point of implementing CLR, PERMANOVA and Benjamini-Hochberg by hand is
to understand them. The point of this script is to make sure understanding
them did not come at the cost of getting them wrong.

Three kinds of check, in increasing order of what they can catch:

1. **Invariants** — properties the maths guarantees, tested on data where
   the answer is known by construction (CLR rows sum to zero; CLR is
   invariant to rescaling a sample; Aitchison distance is a metric).
2. **Reference agreement** — the same statistic computed by scikit-bio and
   statsmodels, which are independently written and widely used.
3. **Behaviour under a known truth** — planted differential taxa that the
   pipeline should recover, and a null dataset where it should find nothing.

Check 3 is the one that matters most. Agreeing with scikit-bio only proves
both implementations do the same thing; recovering planted signal and
staying quiet on null data proves the thing is the right thing.

Run:  python bin/validate_compositional.py [--project-root .]
Exits non-zero on the first failure, so it works as a CI gate.
"""

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", default=str(Path(__file__).resolve().parent.parent))
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--monte-carlo", type=int, default=64)
    return p.parse_args()


class Checks:
    """Minimal assertion harness that reports every check, not just failures."""

    def __init__(self) -> None:
        self.failures: list[str] = []

    def ok(self, name: str, condition: bool, detail: str = "") -> None:
        mark = "PASS" if condition else "FAIL"
        print(f"  [{mark}] {name}" + (f" — {detail}" if detail else ""))
        if not condition:
            self.failures.append(name)

    def close(self, name: str, a: float, b: float, tol: float = 1e-9) -> None:
        self.ok(name, abs(a - b) <= tol, f"{a:.10g} vs {b:.10g} (tol {tol:g})")


def main() -> int:
    args = parse_args()
    sys.path.insert(0, args.project_root)
    from compositional import ordination, stats, transforms

    checks = Checks()
    rng = np.random.default_rng(20260913)

    # ---- 1. Invariants ------------------------------------------------
    print("\nInvariants")
    counts = rng.gamma(2.0, 40.0, size=(30, 80))

    clr_values = transforms.clr(counts)
    checks.ok("CLR rows sum to zero",
              np.abs(clr_values.sum(axis=1)).max() < 1e-8,
              f"max |rowsum| = {np.abs(clr_values.sum(axis=1)).max():.2e}")

    rescaled = counts.copy()
    rescaled[0] *= 1_000_000
    checks.ok("CLR is invariant to sequencing depth",
              np.allclose(transforms.clr(counts)[0], transforms.clr(rescaled)[0]),
              "a sample scaled 1e6x gives identical CLR values")

    sparse = counts.copy()
    sparse[rng.random(sparse.shape) < 0.6] = 0.0
    replaced = transforms.multiplicative_replacement(sparse)
    checks.ok("Zero replacement keeps rows closed and positive",
              np.allclose(replaced.sum(axis=1), 1.0) and replaced.min() > 0,
              f"min part = {replaced.min():.2e}")

    # The regression test for the bug this pipeline actually hit: a sample
    # so shallow that n_zeros * delta would exceed the whole composition.
    shallow = np.zeros((1, 12)); shallow[0, 0] = 3.0; shallow[0, -1] = 1.0
    shallow_replaced = transforms.multiplicative_replacement(shallow)
    checks.ok("Shallow sample does not produce negative parts",
              shallow_replaced.min() > 0 and abs(shallow_replaced.sum() - 1.0) < 1e-12,
              "4 reads across 12 features")

    dist = transforms.aitchison_distance(counts)
    triangle = all(
        dist[i, j] <= dist[i, k] + dist[k, j] + 1e-9
        for i, j, k in rng.integers(0, 30, size=(300, 3))
    )
    checks.ok("Aitchison distance satisfies the triangle inequality", triangle,
              "300 random triples")
    checks.ok("Aitchison distance is symmetric with zero diagonal",
              np.allclose(dist, dist.T) and np.allclose(np.diag(dist), 0))

    # ---- 2. Reference agreement ---------------------------------------
    print("\nAgreement with reference implementations")
    from statsmodels.stats.multitest import multipletests
    pvals = rng.uniform(0, 1, 2000)
    pvals[:50] /= 5000
    mine = stats.benjamini_hochberg(pvals)
    reference = multipletests(pvals, method="fdr_bh")[1]
    checks.ok("Benjamini-Hochberg matches statsmodels",
              np.abs(mine - reference).max() < 1e-12,
              f"max |diff| = {np.abs(mine - reference).max():.2e}")
    checks.ok("Adjusted p-values are monotone in raw p-values",
              bool(np.all(np.diff(mine[np.argsort(pvals)]) >= -1e-12)))

    import pandas as pd
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.distance import permanova as skbio_permanova

    groups = np.array(["veg"] * 15 + ["bare"] * 15)
    for label, matrix in [("null", counts), ("with planted signal", None)]:
        if matrix is None:
            matrix = counts.copy()
            matrix[:15, :12] *= 6.0
        distances = transforms.aitchison_distance(matrix)
        mine_result = stats.permanova(distances, groups, permutations=args.permutations)
        ids = [f"s{i}" for i in range(matrix.shape[0])]
        reference_result = skbio_permanova(
            DistanceMatrix(distances, ids=ids),
            pd.Series(groups, index=ids, name="group"),
            permutations=args.permutations,
        )
        checks.close(f"PERMANOVA pseudo-F matches scikit-bio ({label})",
                     mine_result.pseudo_f, float(reference_result["test statistic"]),
                     tol=1e-10)

    # ---- 3. Behaviour under a known truth ------------------------------
    print("\nRecovery of planted signal")
    n_planted = 12
    truth = counts.copy()
    truth[:15, :n_planted] *= 6.0
    truth = np.round(truth)

    null_result = stats.permanova(transforms.aitchison_distance(counts), groups,
                                  permutations=args.permutations)
    signal_result = stats.permanova(transforms.aitchison_distance(truth), groups,
                                    permutations=args.permutations)
    checks.ok("PERMANOVA finds nothing in null data", null_result.p_value > 0.05,
              f"p = {null_result.p_value:.3f}")
    checks.ok("PERMANOVA detects planted community shift", signal_result.p_value < 0.05,
              f"p = {signal_result.p_value:.3f}, R^2 = {signal_result.r_squared:.3f}")

    feature_ids = [f"ASV{i:03d}" for i in range(truth.shape[1])]
    records = stats.differential_abundance(truth, groups, "veg", "bare", feature_ids,
                                           n_monte_carlo=args.monte_carlo)
    hits = [r for r in records if r["p_adj_monte_carlo"] < 0.05]
    recovered = sum(1 for r in hits if int(r["feature_id"][3:]) < n_planted)
    checks.ok("All planted taxa recovered", recovered == n_planted,
              f"{recovered}/{n_planted} recovered, {len(hits)} total hits")
    checks.ok("Planted taxa are called as enriched in the right group",
              all(r["enriched_in"] == "veg" for r in hits
                  if int(r["feature_id"][3:]) < n_planted))

    n_naive = sum(1 for r in records if r["p_adj"] < 0.05)
    checks.ok("Monte-Carlo pass is no less conservative than the naive test",
              len(hits) <= n_naive, f"{len(hits)} MC hits vs {n_naive} naive hits")

    print("\nUnsupervised structure")
    link, order = ordination.hierarchical_clustering(transforms.aitchison_distance(truth))
    cluster_labels = ordination.cut_clusters(link, n_clusters=2)
    ari = ordination.cluster_group_agreement(cluster_labels, groups)
    checks.ok("Clustering rediscovers the planted grouping", ari > 0.5,
              f"adjusted Rand index = {ari:.2f}")
    checks.ok("Leaf order is a permutation of the samples",
              sorted(order.tolist()) == list(range(truth.shape[0])))

    pca_result = ordination.pca(transforms.clr(truth), feature_ids, n_components=5)
    checks.ok("PCA variance ratios are ordered and sum to <= 1",
              bool(np.all(np.diff(pca_result.explained_variance_ratio) <= 1e-12))
              and pca_result.explained_variance_ratio.sum() <= 1 + 1e-9)
    top_pc1 = {fid for fid, _ in pca_result.top_loadings(0, n_planted)}
    planted = {f"ASV{i:03d}" for i in range(n_planted)}
    overlap = len(top_pc1 & planted)
    checks.ok("PC1 is driven by the planted taxa", overlap >= n_planted // 2,
              f"{overlap}/{n_planted} of the top PC1 loadings are planted taxa")

    print()
    if checks.failures:
        print(f"FAILED: {len(checks.failures)} check(s): {', '.join(checks.failures)}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
