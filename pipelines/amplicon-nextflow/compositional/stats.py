"""Association testing on compositional data: PERMANOVA and differential abundance.

Two questions, two tests:

* "Does the whole community differ between groups?" -> PERMANOVA on the
  Aitchison distance matrix. A permutation test, because the null
  distribution of the pseudo-F statistic has no tractable closed form.
* "Which taxa differ?" -> a per-feature test on CLR values, with
  Benjamini-Hochberg control across the thousands of features tested.

The per-feature test is run twice: once directly on the point-estimate CLR,
and once as a Monte-Carlo average over Dirichlet draws from each sample's
count posterior (the ALDEx2 idea). The second is slower and almost always
more conservative — see `differential_abundance` for why that difference is
the interesting part, not an inconvenience.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats as sp_stats

from .transforms import clr
from .utils import get_logger

logger = get_logger(__name__)


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR correction, returning adjusted p-values.

    Same procedure implemented in ../../stats/rnaseq-stats-notebook for
    RNA-seq; reimplemented here rather than imported across pipelines
    because the two are deliberately independent checks on each other.

    The cumulative-minimum pass from the largest p-value downwards is the
    step-up part, and it is the bit that is easy to get wrong: without it
    the adjusted values are not monotone, and you can end up rejecting a
    larger p-value while accepting a smaller one.
    """
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out


@dataclass
class PermanovaResult:
    pseudo_f: float
    p_value: float
    permutations: int
    r_squared: float
    group_sizes: dict[str, int]


def permanova(
    distances: np.ndarray,
    groups: np.ndarray,
    permutations: int = 999,
    random_state: int = 0,
) -> PermanovaResult:
    """Permutational MANOVA on a distance matrix (Anderson 2001).

    Partitions total variance into between- and within-group components
    using only the pairwise distances — no need to place samples in a
    coordinate space first, which is what lets it work on any metric.

        SS_total   = (1/N) * sum_{i<j} d_ij^2
        SS_within  = sum_g (1/n_g) * sum_{i<j in g} d_ij^2
        pseudo-F   = (SS_between / (a-1)) / (SS_within / (N-a))

    The statistic has the shape of an F, but not its null distribution —
    distances between the same N points are not independent. So the p-value
    comes from shuffling the group labels `permutations` times and asking
    how often chance produces a pseudo-F at least this large. The observed
    value is counted as one of the draws, which is why the denominator is
    `permutations + 1`: it keeps the test valid (and means the smallest
    achievable p-value with 999 permutations is 0.001, never 0).

    Caveat worth stating in any report that uses this: PERMANOVA is
    sensitive to differences in within-group *dispersion*, not only in
    location. A significant result means "these groups differ", not
    necessarily "these groups are centered in different places".
    """
    distances = np.asarray(distances, dtype=float)
    groups = np.asarray(groups)
    n = distances.shape[0]
    if distances.shape != (n, n):
        raise ValueError(f"Distance matrix must be square, got {distances.shape}")
    if groups.shape[0] != n:
        raise ValueError(f"Got {groups.shape[0]} labels for {n} samples")

    labels, counts = np.unique(groups, return_counts=True)
    n_groups = labels.size
    if n_groups < 2:
        raise ValueError(f"PERMANOVA needs at least 2 groups, got {n_groups}")
    if n - n_groups <= 0:
        raise ValueError("Not enough samples for the within-group degrees of freedom")

    squared = distances ** 2

    def sum_squares(labels_vec: np.ndarray) -> tuple[float, float]:
        total = squared.sum() / (2.0 * n)
        within = 0.0
        for label in labels:
            idx = np.flatnonzero(labels_vec == label)
            if idx.size:
                block = squared[np.ix_(idx, idx)]
                within += block.sum() / (2.0 * idx.size)
        return total, within

    total_ss, within_ss = sum_squares(groups)
    between_ss = total_ss - within_ss
    observed_f = (between_ss / (n_groups - 1)) / (within_ss / (n - n_groups))

    rng = np.random.default_rng(random_state)
    shuffled = groups.copy()
    n_at_least = 1  # the observed value counts as one draw
    for _ in range(permutations):
        rng.shuffle(shuffled)
        _, perm_within = sum_squares(shuffled)
        perm_between = total_ss - perm_within
        perm_f = (perm_between / (n_groups - 1)) / (perm_within / (n - n_groups))
        if perm_f >= observed_f:
            n_at_least += 1

    p_value = n_at_least / (permutations + 1)
    logger.info(
        "PERMANOVA: pseudo-F=%.3f, R^2=%.3f, p=%.4f (%d permutations)",
        observed_f, between_ss / total_ss, p_value, permutations,
    )
    return PermanovaResult(
        pseudo_f=float(observed_f),
        p_value=float(p_value),
        permutations=permutations,
        r_squared=float(between_ss / total_ss),
        group_sizes={str(l): int(c) for l, c in zip(labels, counts)},
    )


def _welch_by_feature(values: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray):
    """Welch t-test per column, vectorised. Returns (t, p, mean_diff)."""
    a = values[mask_a]
    b = values[mask_b]
    t_stat, p_val = sp_stats.ttest_ind(a, b, equal_var=False, axis=0)
    return t_stat, p_val, a.mean(axis=0) - b.mean(axis=0)


def differential_abundance(
    counts: np.ndarray,
    groups: np.ndarray,
    group_a: str,
    group_b: str,
    feature_ids: list[str],
    n_monte_carlo: int = 128,
    random_state: int = 0,
):
    """Per-feature differential abundance between two groups, on CLR values.

    Returns a list of dicts (one per feature), sorted by adjusted p-value.

    Why a t-test on CLR rather than a count model
    ---------------------------------------------
    The negative-binomial GLMs used for RNA-seq (see
    ../../stats/rnaseq-stats-notebook) model counts as draws whose mean
    scales with library size. That framing is a poor fit here: 16S counts
    are not independent across taxa within a sample, because the simplex
    constraint couples them. Once CLR has mapped each sample into real
    space, a location test on the transformed values is both simpler and
    better matched to what the data actually are.

    Why the Monte-Carlo pass
    ------------------------
    A point-estimate CLR pretends the observed counts *are* the composition.
    For a taxon seen 3 times in a sample of 8,000 reads, that is a strong
    claim to make from very little evidence. Drawing the composition from
    its Dirichlet posterior instead — Dirichlet(counts + 0.5), the Jeffreys
    prior — and repeating the whole test across draws propagates that
    sequencing-depth uncertainty into the p-value. The returned
    `p_adj_monte_carlo` is the expected adjusted p-value over draws.

    Expect the Monte-Carlo column to be *less* significant than the naive
    one, and expect the gap to be widest for the rarest taxa. That gap is
    the honest answer to "how much of this hit is real signal and how much
    is a small number that happened to land well" — which is why both
    columns are reported rather than only the better-looking one.
    """
    counts = np.asarray(counts, dtype=float)
    groups = np.asarray(groups)
    mask_a = groups == group_a
    mask_b = groups == group_b
    n_a, n_b = int(mask_a.sum()), int(mask_b.sum())
    if n_a < 2 or n_b < 2:
        raise ValueError(f"Need >=2 samples per group, got {group_a}={n_a}, {group_b}={n_b}")
    if counts.shape[1] != len(feature_ids):
        raise ValueError("feature_ids length does not match the count matrix")

    logger.info("Differential abundance: %s (n=%d) vs %s (n=%d), %d features",
                group_a, n_a, group_b, n_b, counts.shape[1])

    # --- Pass 1: point-estimate CLR -------------------------------------
    clr_values = clr(counts)
    t_stat, p_val, mean_diff = _welch_by_feature(clr_values, mask_a, mask_b)
    p_val = np.nan_to_num(p_val, nan=1.0)
    p_adj = benjamini_hochberg(p_val)

    # Standardised effect size: CLR difference over pooled within-group SD.
    pooled_sd = np.sqrt(
        (clr_values[mask_a].var(axis=0, ddof=1) + clr_values[mask_b].var(axis=0, ddof=1)) / 2.0
    )
    effect = np.divide(mean_diff, pooled_sd, out=np.zeros_like(mean_diff), where=pooled_sd > 0)

    # --- Pass 2: Monte-Carlo over the Dirichlet posterior ----------------
    rng = np.random.default_rng(random_state)
    p_adj_accumulator = np.zeros(counts.shape[1], dtype=float)
    diff_accumulator = np.zeros(counts.shape[1], dtype=float)
    for _ in range(n_monte_carlo):
        # Jeffreys prior; gamma-normalise is the standard way to draw
        # Dirichlet rows without a python-level loop over samples.
        gamma = rng.gamma(shape=counts + 0.5, scale=1.0)
        draw = gamma / gamma.sum(axis=1, keepdims=True)
        draw_clr = clr(draw, already_replaced=True)
        _, draw_p, draw_diff = _welch_by_feature(draw_clr, mask_a, mask_b)
        p_adj_accumulator += benjamini_hochberg(np.nan_to_num(draw_p, nan=1.0))
        diff_accumulator += draw_diff
    p_adj_mc = p_adj_accumulator / n_monte_carlo
    mean_diff_mc = diff_accumulator / n_monte_carlo

    prevalence_a = (counts[mask_a] > 0).mean(axis=0)
    prevalence_b = (counts[mask_b] > 0).mean(axis=0)

    records = [
        {
            "feature_id": feature_ids[i],
            "mean_clr_diff": float(mean_diff[i]),
            "effect_size": float(effect[i]),
            "t_statistic": float(t_stat[i]),
            "p_value": float(p_val[i]),
            "p_adj": float(p_adj[i]),
            "mean_clr_diff_monte_carlo": float(mean_diff_mc[i]),
            "p_adj_monte_carlo": float(p_adj_mc[i]),
            f"prevalence_{group_a}": float(prevalence_a[i]),
            f"prevalence_{group_b}": float(prevalence_b[i]),
            "total_reads": float(counts[:, i].sum()),
            "enriched_in": group_a if mean_diff[i] > 0 else group_b,
        }
        for i in range(counts.shape[1])
    ]
    records.sort(key=lambda r: (r["p_adj_monte_carlo"], r["p_adj"]))
    return records
