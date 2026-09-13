#!/usr/bin/env python3
"""Compositional analysis of a QIIME2 feature table: CLR -> ordination -> tests -> report.

Runs as a Nextflow process or straight from a shell. `--project-root` is
passed explicitly rather than derived from __file__, for the same reason as
../scrnaseq-nextflow/bin/run_perturbation_pipeline.py: the caller already
knows where the pipeline root is, and reverse-engineering it from this
script's own path breaks once Nextflow has staged it through PATH.

Outputs written to --outdir:
    differential_abundance.csv    per-taxon results, both ASV and collapsed
    permanova.json                community-level test
    summary.json                  every number quoted in the report
    *.png                         figures
    report.md                     the written interpretation
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True, help="feature-table.tsv from `biom convert`")
    p.add_argument("--taxonomy", required=True, help="taxonomy.tsv from the classifier export")
    p.add_argument("--metadata", required=True, help="QIIME2 sample metadata TSV")
    p.add_argument("--outdir", required=True)
    p.add_argument("--project-root", required=True,
                   help="Pipeline root containing the compositional/ package")
    p.add_argument("--group-column", default="vegetation",
                   help="Metadata column defining the two groups to compare")
    p.add_argument("--group-a", default="yes")
    p.add_argument("--group-b", default="no")
    p.add_argument("--rank", default="genus", help="Taxonomic rank to collapse to")
    p.add_argument("--min-prevalence", type=float, default=0.15)
    p.add_argument("--min-total-count", type=int, default=25)
    p.add_argument("--min-depth", type=int, default=1000,
                   help="Drop samples with fewer reads than this after denoising. "
                        "Shallow samples are mostly zeros, so their CLR values end up "
                        "dominated by zero imputation rather than measurement.")
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--monte-carlo", type=int, default=128)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def analyse(counts_df, groups, label, args, modules, outdir):
    """Run the full compositional workflow on one feature table.

    Called twice — once on ASVs, once on the rank-collapsed table — because
    the two answer different questions and disagreeing is informative, not a
    bug. See README.md, "Why both ASV and genus level".
    """
    transforms, stats, ordination, plots = modules
    counts = counts_df.to_numpy(dtype=float)
    feature_ids = list(counts_df.columns)
    sample_ids = list(counts_df.index)

    clr_values = transforms.clr(counts)
    distances = transforms.aitchison_distance(counts)

    permanova_result = stats.permanova(
        distances, groups, permutations=args.permutations, random_state=args.seed
    )

    pca_result = ordination.pca(clr_values, feature_ids, n_components=5)
    link, leaf_order = ordination.hierarchical_clustering(distances)
    cluster_labels = ordination.cut_clusters(link, n_clusters=len(set(groups)))
    ari = ordination.cluster_group_agreement(cluster_labels, groups)

    records = stats.differential_abundance(
        counts, groups, args.group_a, args.group_b, feature_ids,
        n_monte_carlo=args.monte_carlo, random_state=args.seed,
    )
    significant = [r for r in records if r["p_adj_monte_carlo"] < args.alpha]

    prefix = f"{label}_"
    plots.pca_scatter(pca_result, groups, outdir / f"{prefix}pca.png",
                      title=f"Aitchison PCA — {label}-level ({counts.shape[1]} features)")
    plots.scree_bar(pca_result, outdir / f"{prefix}scree.png")
    plots.clustered_heatmap(clr_values, sample_ids, feature_ids, link, leaf_order,
                            groups, outdir / f"{prefix}heatmap_dendrogram.png")
    plots.effect_volcano(records, outdir / f"{prefix}volcano.png", alpha=args.alpha,
                         group_a=args.group_a, group_b=args.group_b)
    if significant:
        plots.top_hits_boxplot(clr_values, feature_ids, groups, significant,
                               outdir / f"{prefix}top_hits_boxplot.png")

    pd.DataFrame(records).to_csv(outdir / f"{prefix}differential_abundance.csv", index=False)

    return {
        "level": label,
        "n_samples": int(counts.shape[0]),
        "n_features": int(counts.shape[1]),
        "total_reads": int(counts.sum()),
        "permanova": {
            "pseudo_f": permanova_result.pseudo_f,
            "p_value": permanova_result.p_value,
            "r_squared": permanova_result.r_squared,
            "permutations": permanova_result.permutations,
            "group_sizes": permanova_result.group_sizes,
        },
        "pca_explained_variance": [float(v) for v in pca_result.explained_variance_ratio],
        "pc1_top_loadings": pca_result.top_loadings(0, 8),
        "cluster_vs_group_adjusted_rand": ari,
        "n_significant": len(significant),
        "n_significant_naive": int(sum(r["p_adj"] < args.alpha for r in records)),
        "top_hits": [
            {
                "feature_id": r["feature_id"],
                "mean_clr_diff": r["mean_clr_diff"],
                "effect_size": r["effect_size"],
                "p_adj": r["p_adj"],
                "p_adj_monte_carlo": r["p_adj_monte_carlo"],
                "enriched_in": r["enriched_in"],
            }
            for r in significant[:15]
        ],
    }


def write_report(summary: dict, args, path: Path) -> None:
    """Manuscript-style written interpretation, not just a dump of numbers."""
    lines = [
        "# Compositional analysis of the Atacama soil microbiome",
        "",
        f"Comparison: `{args.group_column}` = **{args.group_a}** vs **{args.group_b}**.",
        "",
        "## Method",
        "",
        "Counts were treated as compositions throughout. Zeros were imputed by",
        "multiplicative replacement, abundances centered-log-ratio transformed, and all",
        "distances computed as Euclidean distance in CLR space (Aitchison distance).",
        f"Community-level differences were tested by PERMANOVA with {args.permutations}",
        "permutations; per-taxon differences by Welch's t-test on CLR values with",
        "Benjamini-Hochberg FDR control, repeated across",
        f"{args.monte_carlo} Dirichlet draws from each sample's count posterior so that",
        "sequencing-depth uncertainty is carried into the reported p-values.",
        "",
    ]
    for level in summary["levels"]:
        pm = level["permanova"]
        verdict = "does" if pm["p_value"] < args.alpha else "does not"
        lines += [
            f"## {level['level'].capitalize()} level",
            "",
            f"- {level['n_samples']} samples, {level['n_features']} features, "
            f"{level['total_reads']:,} reads after filtering.",
            f"- PERMANOVA: pseudo-F = {pm['pseudo_f']:.2f}, R² = {pm['r_squared']:.3f}, "
            f"p = {pm['p_value']:.3f} — `{args.group_column}` {verdict} explain a "
            f"significant share of community variation "
            f"({100 * pm['r_squared']:.1f}% of the total).",
            f"- PC1 explains {100 * level['pca_explained_variance'][0]:.1f}% of variance; "
            f"PC2 {100 * level['pca_explained_variance'][1]:.1f}%.",
            f"- Unsupervised clustering vs. the known grouping: "
            f"adjusted Rand index {level['cluster_vs_group_adjusted_rand']:.2f}.",
            f"- Differentially abundant taxa at FDR < {args.alpha}: "
            f"**{level['n_significant']}** "
            f"(naive CLR test alone would report {level['n_significant_naive']}; "
            f"the difference is the cost of taking sequencing-depth uncertainty seriously).",
            "",
        ]
        if level["top_hits"]:
            lines += ["Top hits:", "",
                      "| Taxon | CLR difference | Effect size | FDR (Monte-Carlo) | Enriched in |",
                      "|---|---|---|---|---|"]
            for hit in level["top_hits"][:10]:
                lines.append(
                    f"| `{hit['feature_id'][:48]}` | {hit['mean_clr_diff']:+.2f} | "
                    f"{hit['effect_size']:+.2f} | {hit['p_adj_monte_carlo']:.2e} | "
                    f"{hit['enriched_in']} |"
                )
            lines.append("")
        else:
            lines += ["No taxon passed the FDR threshold at this level.", ""]

    lines += [
        "## Reading the caveats",
        "",
        "- PERMANOVA responds to differences in dispersion as well as in location. A",
        "  significant result means the groups differ, not necessarily that they are",
        "  centered in different places.",
        "- Every result is a statement about the *filtered subcomposition*. Features below",
        f"  the prevalence ({100 * args.min_prevalence:.0f}% of samples) and count "
f"({args.min_total_count} reads) thresholds were removed before transformation, "
f"and samples below {args.min_depth:,} reads were dropped entirely.",
        "- CLR values are relative to each sample's own geometric mean. A taxon 'enriched",
        "  in " + args.group_a + "' is enriched relative to the rest of that community; this",
        "  design cannot speak to absolute cell counts per gram of soil.",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    sys.path.insert(0, args.project_root)
    from compositional import io, ordination, plots, stats, transforms, utils

    logger = utils.get_logger("run_compositional_analysis")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = io.read_feature_table(args.table)
    taxonomy = io.read_taxonomy(args.taxonomy)
    metadata = io.read_metadata(args.metadata)
    table, metadata = io.align(table, metadata)

    if args.group_column not in metadata.columns:
        raise SystemExit(
            f"Metadata has no column {args.group_column!r}. Available: {list(metadata.columns)}"
        )
    group_series = metadata[args.group_column].astype(str)
    keep = group_series.isin([args.group_a, args.group_b])
    if keep.sum() < 4:
        raise SystemExit(
            f"Only {int(keep.sum())} samples in groups {args.group_a}/{args.group_b}; "
            f"values present: {sorted(group_series.unique())}"
        )
    logger.info("Restricting to %d of %d samples in the two compared groups",
                int(keep.sum()), len(keep))
    table = table.loc[keep.values]
    groups = group_series[keep.values].to_numpy()

    depths = table.sum(axis=1).to_numpy(dtype=float)
    plots.depth_bar(depths, list(table.index), groups, outdir / "sample_depth.png")

    asv_table = io.drop_empty_samples(table, min_depth=args.min_depth)
    groups = group_series.loc[asv_table.index].to_numpy()
    genus_table = io.collapse_to_rank(asv_table, taxonomy, rank=args.rank)

    modules = (transforms, stats, ordination, plots)
    levels = []
    for label, counts_df in [("asv", asv_table), (args.rank, genus_table)]:
        filtered = io.filter_features(
            counts_df, min_prevalence=args.min_prevalence, min_total_count=args.min_total_count
        )
        # Re-check depth after filtering: a sample can pass the raw threshold
        # and still be left near-empty once rare features are dropped.
        filtered = io.drop_empty_samples(filtered, min_depth=max(args.min_depth // 10, 1))
        level_groups = group_series.loc[filtered.index].to_numpy()
        levels.append(analyse(filtered, level_groups, label, args, modules, outdir))

    summary = {
        "group_column": args.group_column,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "alpha": args.alpha,
        "filters": {
            "min_prevalence": args.min_prevalence,
            "min_total_count": args.min_total_count,
            "min_depth": args.min_depth,
        },
        "permutations": args.permutations,
        "monte_carlo_draws": args.monte_carlo,
        "seed": args.seed,
        "levels": levels,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    (outdir / "permanova.json").write_text(json.dumps(
        {lv["level"]: lv["permanova"] for lv in levels}, indent=2))
    write_report(summary, args, outdir / "report.md")

    logger.info("Done. %s",
                ", ".join(f"{lv['level']}: {lv['n_significant']} significant taxa" for lv in levels))


if __name__ == "__main__":
    main()
