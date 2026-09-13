"""Report figures: scatter, box, bar, heatmap and dendrogram.

Deliberately matplotlib-only and Agg-backed so the same code produces the
same PNGs on a laptop, in a container, and on a CI runner with no display.

Every figure here is meant to be readable by someone who did not run the
analysis: axis labels carry units and variance explained, group colours are
consistent across panels, and n is stated on the figure rather than left in
the log.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from .utils import get_logger

logger = get_logger(__name__)

# Colour-blind-safe pair (Okabe-Ito), used for the two-group comparisons
# throughout so a colour means the same thing in every panel.
GROUP_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]


def _color_map(groups: np.ndarray) -> dict[str, str]:
    labels = sorted(set(map(str, groups)))
    return {label: GROUP_COLORS[i % len(GROUP_COLORS)] for i, label in enumerate(labels)}


def pca_scatter(pca_result, groups: np.ndarray, outpath: Path, title: str = "") -> Path:
    """Sample ordination: PC1 vs PC2 of the CLR values, coloured by group."""
    colors = _color_map(groups)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for label, color in colors.items():
        mask = np.asarray(list(map(str, groups))) == label
        ax.scatter(pca_result.coordinates[mask, 0], pca_result.coordinates[mask, 1],
                   c=color, s=55, alpha=0.85, edgecolor="white", linewidth=0.8,
                   label=f"{label} (n={int(mask.sum())})")
    ev = pca_result.explained_variance_ratio
    ax.set_xlabel(f"PC1 — {100 * ev[0]:.1f}% of variance")
    ax.set_ylabel(f"PC2 — {100 * ev[1]:.1f}% of variance")
    ax.set_title(title or "Aitchison PCA of CLR-transformed abundances")
    ax.axhline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath


def scree_bar(pca_result, outpath: Path, n: int = 10) -> Path:
    """Bar chart of variance explained per component."""
    ev = pca_result.explained_variance_ratio[:n]
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.bar(np.arange(1, ev.size + 1), 100 * ev, color="#0072B2", alpha=0.9)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title("Scree plot")
    ax.set_xticks(np.arange(1, ev.size + 1))
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath


def clustered_heatmap(clr_values: np.ndarray, sample_ids: list[str], feature_ids: list[str],
                      link, leaf_order: np.ndarray, groups: np.ndarray,
                      outpath: Path, n_features: int = 30) -> Path:
    """Dendrogram of samples beside a heatmap of the most variable taxa.

    Rows are ordered by the dendrogram, not by group, so the figure shows
    what the clustering found rather than what the study design asserts —
    the group colour strip on the left is the only place the design appears,
    which makes agreement (or disagreement) visible at a glance.

    The colour scale is diverging and centered at 0 because these are CLR
    values: 0 is "average for this sample", not "absent".
    """
    values = np.asarray(clr_values, dtype=float)
    variance = values.var(axis=0)
    top = np.argsort(-variance)[:n_features]
    sub = values[np.ix_(leaf_order, top)]
    colors = _color_map(groups)

    fig = plt.figure(figsize=(12, max(5.0, 0.22 * len(sample_ids))))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.1, 0.09, 5.0], wspace=0.02)

    ax_dendro = fig.add_subplot(grid[0])
    dendrogram(link, orientation="left", ax=ax_dendro, no_labels=True,
               color_threshold=0, above_threshold_color="#555555")
    ax_dendro.invert_yaxis()
    ax_dendro.set_xticks([])
    for spine in ax_dendro.spines.values():
        spine.set_visible(False)
    ax_dendro.set_title("Aitchison\nlinkage", fontsize=9)

    ax_strip = fig.add_subplot(grid[1])
    ordered_groups = [str(groups[i]) for i in leaf_order]
    strip = np.array([[list(colors).index(g)] for g in ordered_groups], dtype=float)
    ax_strip.imshow(strip, aspect="auto",
                    cmap=matplotlib.colors.ListedColormap([colors[k] for k in colors]),
                    vmin=0, vmax=max(len(colors) - 1, 1))
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    ax_strip.set_title("group", fontsize=8)

    ax_heat = fig.add_subplot(grid[2])
    limit = float(np.abs(sub).max())
    im = ax_heat.imshow(sub, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax_heat.set_yticks(range(len(leaf_order)))
    ax_heat.set_yticklabels([sample_ids[i] for i in leaf_order], fontsize=6)
    ax_heat.yaxis.tick_right()
    ax_heat.set_xticks(range(len(top)))
    ax_heat.set_xticklabels([feature_ids[i][:28] for i in top], rotation=90, fontsize=6)
    ax_heat.set_title(f"{n_features} most variable taxa (CLR)", fontsize=10)
    fig.colorbar(im, ax=ax_heat, fraction=0.015, pad=0.10, label="CLR (0 = sample average)")

    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=c, label=g)
               for g, c in colors.items()]
    ax_dendro.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath


def effect_volcano(records: list[dict], outpath: Path, alpha: float = 0.05,
                   group_a: str = "A", group_b: str = "B") -> Path:
    """Effect size against significance, using the Monte-Carlo adjusted p-value."""
    diff = np.array([r["mean_clr_diff"] for r in records])
    padj = np.array([r["p_adj_monte_carlo"] for r in records])
    significant = padj < alpha
    # Floor at the smallest representable positive value so -log10 is finite.
    y = -np.log10(np.clip(padj, 1e-300, None))

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    ax.scatter(diff[~significant], y[~significant], s=18, c="#bbbbbb",
               label=f"not significant (n={int((~significant).sum())})")
    ax.scatter(diff[significant], y[significant], s=30, c="#D55E00", alpha=0.9,
               label=f"FDR < {alpha} (n={int(significant.sum())})")
    ax.axhline(-np.log10(alpha), color="#555555", ls="--", lw=0.9)
    ax.axvline(0, color="#cccccc", lw=0.8)
    ax.set_xlabel(f"Mean CLR difference  ({group_a} − {group_b})")
    ax.set_ylabel("−log10 adjusted p (Monte-Carlo)")
    ax.set_title("Differential abundance")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath


def top_hits_boxplot(clr_values: np.ndarray, feature_ids: list[str], groups: np.ndarray,
                     records: list[dict], outpath: Path, n: int = 6) -> Path:
    """Box plots of the top differentially abundant taxa, with the points shown.

    Points over boxes on purpose: with tens of samples per group, a box plot
    alone hides whether a difference is a consistent shift or two outliers
    dragging a median.
    """
    index = {fid: i for i, fid in enumerate(feature_ids)}
    chosen = [r for r in records if r["feature_id"] in index][:n]
    if not chosen:
        raise ValueError("No plottable features among the top records")

    labels = sorted(set(map(str, groups)))
    colors = _color_map(groups)
    group_strings = np.asarray(list(map(str, groups)))

    ncols = min(3, len(chosen))
    nrows = int(np.ceil(len(chosen) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.4 * nrows), squeeze=False)

    for ax, record in zip(axes.ravel(), chosen):
        col = index[record["feature_id"]]
        data = [clr_values[group_strings == label, col] for label in labels]
        # Tick labels are set afterwards rather than passed in: the keyword
        # for this was renamed `labels` -> `tick_labels` in matplotlib 3.9,
        # and neither spelling works across both sides of that release.
        bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                        showfliers=False, medianprops=dict(color="black"))
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        for patch, label in zip(bp["boxes"], labels):
            patch.set_facecolor(colors[label])
            patch.set_alpha(0.35)
        rng = np.random.default_rng(0)
        for i, (values, label) in enumerate(zip(data, labels), start=1):
            ax.scatter(rng.normal(i, 0.06, values.size), values, s=14,
                       c=colors[label], alpha=0.85, edgecolor="white", linewidth=0.4)
        ax.set_ylabel("CLR")
        ax.set_title(f"{record['feature_id'][:34]}\nFDR={record['p_adj_monte_carlo']:.1e}",
                     fontsize=9)

    for ax in axes.ravel()[len(chosen):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath


def depth_bar(depths: np.ndarray, sample_ids: list[str], groups: np.ndarray,
              outpath: Path) -> Path:
    """Reads retained per sample after denoising — the first QC figure to read.

    Included because sequencing depth is the confounder that compositional
    methods are designed to neutralise, and the reader deserves to see how
    uneven it was before being told that it did not matter.
    """
    order = np.argsort(-depths)
    colors = _color_map(groups)
    fig, ax = plt.subplots(figsize=(max(7.0, 0.16 * len(depths)), 4.0))
    ax.bar(range(len(depths)), depths[order],
           color=[colors[str(groups[i])] for i in order])
    ax.set_ylabel("Reads after denoising")
    ax.set_xlabel("Sample (sorted by depth)")
    ax.set_xticks(range(len(depths)))
    ax.set_xticklabels([sample_ids[i] for i in order], rotation=90, fontsize=5)
    ax.set_title(f"Per-sample depth: median {np.median(depths):,.0f}, "
                 f"range {depths.min():,.0f}–{depths.max():,.0f} "
                 f"({depths.max() / max(depths.min(), 1):.0f}× spread)")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=c, label=g)
               for g, c in colors.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", outpath)
    return outpath
