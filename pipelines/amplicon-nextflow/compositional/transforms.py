"""Compositional transforms: closure, zero replacement, CLR, Aitchison distance.

Why any of this is needed
-------------------------
A 16S feature table does not measure abundance. Sequencing returns a fixed
number of reads per sample — whatever the library prep and the flow cell
happened to give you — so a column of counts carries information only about
*ratios between taxa within that sample*, not about how many cells were in
the soil. The data live on a simplex, not in real space.

Two consequences that bite in practice:

1. **Spurious negative correlation.** Because the parts must sum to a
   constant, one taxon going up forces every other proportion down. Run
   Pearson correlation on raw proportions and you manufacture negative
   associations that do not exist biologically.

2. **Differences are not distances.** "Taxon A went from 1% to 2%" and
   "taxon B went from 50% to 51%" are both +1 percentage point, but the
   first is a doubling and the second is a 2% relative change. Euclidean
   distance on proportions treats them as equal. They are not.

Aitchison's answer is to work with log-ratios, which are invariant to the
arbitrary total. The centered log-ratio (CLR) divides each part by the
geometric mean of its own sample before taking logs, which removes the
sample-total dependence and lands the data back in ordinary real space
where Euclidean geometry, PCA, t-tests and clustering all behave.

Implemented here from first principles rather than imported, for the same
reason as ../../stats/rnaseq-stats-notebook: knowing the formula is not the
same as knowing what it does to your data. `bin/validate_compositional.py` cross-checks
these against scikit-bio and statsmodels where an equivalent exists.

Reference: Aitken, J. (1982) is the usual citation via Aitchison, J.,
"The Statistical Analysis of Compositional Data", JRSS-B 44(2):139-177.
"""

from __future__ import annotations

import numpy as np

from .utils import get_logger

logger = get_logger(__name__)


def closure(mat: np.ndarray) -> np.ndarray:
    """Rescale each row (sample) to sum to 1.

    This is the formal definition of "make it a composition". It is also the
    step that makes explicit what sequencing already did to you implicitly:
    after closure, the total carries no information, so anything you conclude
    afterwards must be a statement about ratios.
    """
    mat = np.asarray(mat, dtype=float)
    if mat.ndim != 2:
        raise ValueError(f"Expected a 2-D samples x features matrix, got shape {mat.shape}")
    totals = mat.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        bad = int(np.sum(totals <= 0))
        raise ValueError(f"{bad} sample(s) have zero total counts — filter them before closure")
    return mat / totals


def multiplicative_replacement(
    mat: np.ndarray,
    delta: float | None = None,
    max_imputed_mass: float = 0.5,
) -> np.ndarray:
    """Replace zeros with a small value, then re-close the composition.

    The CLR takes a logarithm, and log(0) is undefined — so every
    compositional workflow has to decide what a zero means. In amplicon data
    a zero is almost always a *rounded* zero: the taxon was below the
    detection limit of this library's sequencing depth, not truly absent.

    The multiplicative replacement substitutes `delta` for each zero and
    shrinks the non-zero parts proportionally so the row still sums to 1.
    The "multiplicative" part matters: naively adding a pseudocount to
    every cell (`counts + 1`) distorts the ratios between the non-zero
    parts, and does so hardest in the shallowest samples — which is exactly
    backwards, since those are the samples you already trust least.

    Default `delta` is 65% of the smallest detectable proportion (1 / depth),
    the conventional choice from Martin-Fernandez et al. (2003).

    The `max_imputed_mass` cap, and why it is not optional
    -----------------------------------------------------
    That conventional delta is derived for samples with a sensible number of
    reads, and it silently breaks when they do not have one. The imputed
    mass is `n_zeros * delta`; with a shallow sample and a sparse table,
    that product can exceed 1, at which point re-closure assigns the
    observed parts a *negative* proportion and the CLR fails downstream with
    an error that points nowhere near the real cause.

    Found on real data, not in theory: the Atacama soil table at 1%
    subsampling gives some samples ~4 reads across 12 retained features —
    10 zeros x 0.1625 = 1.625 of imputed mass for a composition that must
    sum to 1. The cap bounds the total imputed mass at `max_imputed_mass`
    (shrinking delta for that row only), which keeps re-closure well-posed.

    A row that hits the cap is telling you something real: more than half of
    that sample's composition is now imputation rather than measurement. The
    warning below says so, and the right response is usually to raise the
    minimum depth and drop the sample, not to trust the transformed values.
    """
    mat = np.asarray(mat, dtype=float)
    if np.any(mat < 0):
        raise ValueError("Negative values in count matrix")
    if not 0 < max_imputed_mass < 1:
        raise ValueError(f"max_imputed_mass must be in (0, 1), got {max_imputed_mass}")

    closed = closure(mat)
    n_zero = int(np.sum(closed == 0))
    if n_zero == 0:
        return closed

    zeros = closed == 0
    n_zeros_per_row = zeros.sum(axis=1, keepdims=True)

    if delta is None:
        # Smallest proportion this table could have resolved, per sample.
        depth = mat.sum(axis=1, keepdims=True)
        delta_vec = 0.65 / depth
    else:
        delta_vec = np.full((mat.shape[0], 1), float(delta))

    cap = max_imputed_mass / np.maximum(n_zeros_per_row, 1)
    capped_rows = delta_vec > cap
    if np.any(capped_rows & (n_zeros_per_row > 0)):
        logger.warning(
            "%d sample(s) too shallow for the standard zero replacement: imputed mass "
            "would have exceeded %.0f%% of the composition, so delta was reduced for "
            "those rows. Their CLR values are mostly imputation — consider raising the "
            "minimum sample depth.",
            int(np.sum(capped_rows & (n_zeros_per_row > 0))), 100 * max_imputed_mass,
        )
    delta_vec = np.minimum(delta_vec, cap)

    replaced = np.where(zeros, np.broadcast_to(delta_vec, closed.shape), closed)

    # Shrink the observed (non-zero) parts so each row still sums to 1.
    zero_mass = (zeros * delta_vec).sum(axis=1, keepdims=True)
    observed_mass = np.where(zeros, 0.0, closed).sum(axis=1, keepdims=True)
    if np.any(observed_mass <= 0):
        raise ValueError("Sample(s) with no non-zero features — raise the minimum depth")
    scale = (1.0 - zero_mass) / observed_mass
    replaced = np.where(zeros, replaced, replaced * scale)

    if np.any(replaced <= 0):
        raise AssertionError("Zero replacement produced a non-positive part — this is a bug")

    logger.info(
        "Multiplicative replacement: %d zeros of %d cells (%.1f%%) replaced",
        n_zero, closed.size, 100 * n_zero / closed.size,
    )
    return replaced


def clr(mat: np.ndarray, already_replaced: bool = False) -> np.ndarray:
    """Centered log-ratio transform.

        clr(x)_i = log( x_i / g(x) ),  where g(x) is the geometric mean of x

    Dividing by the sample's own geometric mean is what removes the
    dependence on sequencing depth: scale every count in a sample by k and
    the geometric mean scales by k too, so the ratio — and therefore the CLR
    value — is unchanged. That invariance is the whole point.

    Interpretation: a CLR value is "how much more of this taxon than the
    typical taxon in this sample", on a log scale. Positive means enriched
    relative to the sample's own average; zero means average. Because the
    values are centered per sample, each row sums to zero by construction —
    which is also why CLR data is singular and a covariance matrix built
    from it is rank-deficient by one. PCA handles that fine; anything that
    needs to invert the covariance matrix does not.
    """
    mat = np.asarray(mat, dtype=float)
    parts = mat if already_replaced else multiplicative_replacement(mat)
    if np.any(parts <= 0):
        raise ValueError("Non-positive values remain — run multiplicative_replacement first")

    log_parts = np.log(parts)
    geometric_mean_log = log_parts.mean(axis=1, keepdims=True)
    out = log_parts - geometric_mean_log

    row_sums = np.abs(out.sum(axis=1))
    if np.max(row_sums) > 1e-8:
        raise AssertionError(
            f"CLR rows should sum to 0 by construction; max |rowsum| = {np.max(row_sums):.2e}"
        )
    return out


def aitchison_distance(mat: np.ndarray) -> np.ndarray:
    """Pairwise Aitchison distance = Euclidean distance in CLR space.

    Returned as a full square matrix (n_samples x n_samples).

    This is the compositional analogue of Euclidean distance, and the reason
    it is worth the trouble over Bray-Curtis: it is subcompositionally
    coherent. Drop half the taxa from the table and the distances between
    the samples you kept are computed from the remaining ratios only — they
    do not lurch around because the total changed. Bray-Curtis, computed on
    proportions, has no such guarantee, so filtering rare taxa can move your
    ordination.
    """
    coords = clr(mat)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    np.fill_diagonal(dist, 0.0)
    return dist
