#!/usr/bin/env python3
"""QC -> normalize -> cluster -> differential expression, on a perturbation screen.

Called as a Nextflow process / Snakemake rule / directly from the shell —
`--project-root` is passed explicitly (rather than resolved from __file__)
so the import of the local `perturbation` package works the same way
whether this runs on the host, under conda, or under Docker: the caller
already has to know where the pipeline root is (Nextflow's `${projectDir}`,
Snakemake's `workflow.basedir`), so it's more robust to just pass it than
to reverse-engineer it from this script's own location once it's been
staged/symlinked/PATH-resolved by three different tools.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Raw perturbation screen .h5ad")
    p.add_argument("--outdir", required=True)
    p.add_argument("--project-root", required=True, help="Pipeline root containing the perturbation/ package")
    p.add_argument("--perturbation-col", default="perturbation")
    p.add_argument("--dose-col", default="dose_value")
    p.add_argument("--control-label", default="control")
    p.add_argument("--mt-prefix", default="MT-", help="'MT-' for human, 'mt-' for mouse")
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--max-pct-mt", type=float, default=20.0)
    p.add_argument("--resolution", type=float, default=0.5)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    sys.path.insert(0, args.project_root)
    from perturbation import analysis, metadata, preprocessing, qc, utils

    logger = utils.get_logger("run_perturbation_pipeline")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = utils.load_adata(args.input)
    logger.info("Loaded %s: %d cells x %d genes", args.input, adata.n_obs, adata.n_vars)

    adata = metadata.standardize_condition_columns(
        adata, args.perturbation_col, args.dose_col, args.control_label
    )
    validation = metadata.validate_conditions(adata)
    logger.info("Conditions found: %s", validation["conditions"])
    if not validation["has_control"]:
        raise RuntimeError("No control cells found — cannot run condition-vs-control DE")

    adata = qc.compute_qc_metrics(adata, mt_prefix=args.mt_prefix)
    adata = qc.filter_cells_and_genes(adata, min_genes=args.min_genes, max_pct_mt=args.max_pct_mt)

    adata = preprocessing.normalize_and_reduce(adata)
    adata = preprocessing.cluster(adata, resolution=args.resolution)

    de_results = analysis.differential_expression(adata, reference="control", groupby="condition")
    de_results.to_csv(outdir / "differential_expression.csv", index=False)

    summary = {
        "n_cells_final": adata.n_obs,
        "n_genes_final": adata.n_vars,
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "conditions": validation["conditions"],
        "n_significant_per_condition": {
            cond: int(((de_results["condition"] == cond) & (de_results["pvals_adj"] < 0.05)).sum())
            for cond in de_results["condition"].unique()
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", summary)

    utils.save_adata(adata, outdir / "processed.h5ad")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sc.pl.umap(adata, color="condition", ax=axes[0], show=False, title="Colored by drug")
    sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, title="Colored by Leiden cluster")
    plt.tight_layout()
    fig.savefig(outdir / "umap_condition_cluster.png", dpi=150)


if __name__ == "__main__":
    main()
