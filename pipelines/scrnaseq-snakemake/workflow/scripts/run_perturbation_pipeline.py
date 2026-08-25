"""Runs inside Snakemake's `script:` directive — `snakemake` object is injected.

Same underlying `perturbation` package and same logic as
../../scrnaseq-nextflow/bin/run_perturbation_pipeline.py — this is the
Snakemake-calling-convention adaptation of that script (snakemake@ input/
output/params instead of argparse), same relationship deseq2.R has to its
Nextflow counterpart in the bulk pipeline. Kept as a real second
implementation, not a symlink, so a match between the two isn't just
"one script running twice."
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perturbation import analysis, metadata, preprocessing, qc, utils

logger = utils.get_logger("run_perturbation_pipeline")

outdir = Path(snakemake.output["processed_h5ad"]).parent
outdir.mkdir(parents=True, exist_ok=True)

params = snakemake.params

adata = utils.load_adata(snakemake.input[0])
logger.info("Loaded %s: %d cells x %d genes", snakemake.input[0], adata.n_obs, adata.n_vars)

adata = metadata.standardize_condition_columns(
    adata, params["perturbation_col"], params["dose_col"], params["control_label"]
)
validation = metadata.validate_conditions(adata)
logger.info("Conditions found: %s", validation["conditions"])
if not validation["has_control"]:
    raise RuntimeError("No control cells found — cannot run condition-vs-control DE")

adata = qc.compute_qc_metrics(adata, mt_prefix=params["mt_prefix"])
adata = qc.filter_cells_and_genes(adata, min_genes=params["min_genes"], max_pct_mt=params["max_pct_mt"])

adata = preprocessing.normalize_and_reduce(adata)
adata = preprocessing.cluster(adata, resolution=params["leiden_resolution"])

de_results = analysis.differential_expression(adata, reference="control", groupby="condition")
de_results.to_csv(snakemake.output["de_csv"], index=False)

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
Path(snakemake.output["summary"]).write_text(json.dumps(summary, indent=2))
logger.info("Summary: %s", summary)

utils.save_adata(adata, snakemake.output["processed_h5ad"])

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sc.pl.umap(adata, color="condition", ax=axes[0], show=False, title="Colored by drug")
sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, title="Colored by Leiden cluster")
plt.tight_layout()
fig.savefig(snakemake.output["umap_plot"], dpi=150)
