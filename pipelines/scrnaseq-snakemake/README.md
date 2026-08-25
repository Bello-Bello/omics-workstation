# scrnaseq-snakemake

[![scrnaseq-snakemake CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml)

Snakemake port of [scrnaseq-nextflow's `perturbation.nf`](../scrnaseq-nextflow/perturbation.nf) —
same public drug-screen data, same QC → normalize → cluster → differential-expression logic,
different workflow manager. Same relationship as `rnaseq-snakemake`/`rnaseq-nextflow`, applied to
the perturbation-screen analysis specifically rather than the STARsolo quantification stage (see
below for why).

## Why not port the STARsolo quantification stage too?

`scrnaseq-nextflow/main.nf` quantifies test-scale chr19 data and its whole reason for existing is
validating that machinery, including a real macOS-specific STAR bug found and worked around via
Docker (documented in [that pipeline's README](../scrnaseq-nextflow/README.md)). Reimplementing
that in Snakemake would mostly duplicate the same STAR invocation and the same Docker requirement
without adding new signal. The perturbation analysis is where the actual interesting content
is — real biology, real differential expression, real dose-response — so that's what's built
twice for cross-verification, not the quantification mechanics.

## Structure

```
config/config.yaml               # dataset name + perturbation analysis params
workflow/
  Snakefile                       # includes rules/, defines `rule all`
  rules/perturbation.smk          # fetch_perturbation_data, perturbation_analysis rules
  scripts/
    fetch_perturbation_data.py    # snakemake@ version of the Nextflow fetch script
    run_perturbation_pipeline.py  # snakemake@ version of the Nextflow analysis script
    perturbation/                 # same package as scrnaseq-nextflow/perturbation/, copied
                                   # (not imported cross-directory) — same reasoning as
                                   # deseq2.R existing as two independent copies in the bulk
                                   # pipelines: each pipeline stays self-contained, the
                                   # calling convention differs (CLI args vs snakemake@
                                   # object), the underlying logic doesn't
envs/scanpy.yaml                  # identical env spec to scrnaseq-nextflow's, for a fair comparison
```

## Setup

```bash
conda activate snakemake   # or: conda create -n snakemake -c bioconda -c conda-forge snakemake -y
cd pipelines/scrnaseq-snakemake
```

## Run

```bash
snakemake --use-conda --cores 4 --snakefile workflow/Snakefile
```

Same dataset, same environment spec as the Nextflow version — first run downloads the ~470 MB
sci-Plex dataset and builds the conda env; subsequent runs reuse both.

## Output

`results/perturbation/`: `processed.h5ad`, `differential_expression.csv`, `summary.json`,
`umap_condition_cluster.png` — identical shape to `scrnaseq-nextflow`'s output.

## Cross-verification

Both pipelines pin the same package versions (`envs/scanpy.yaml` is byte-identical between the
two) and Scanpy's clustering/PCA/UMAP functions default to fixed random seeds, so a real
comparison is meaningful rather than coincidental. Confirmed by diffing both runs'
`summary.json`: **exact match** — 23,966 cells x 35,916 genes post-QC, 7 Leiden clusters, and
identical significant-DE-gene counts per drug (Dex: 5,415, Nutlin: 6,689, BMS: 6,291,
SAHA: 8,037). Same relationship as the bulk pipelines' "same real yeast RNA-seq experiment...
producing identical differential expression calls" result — cross-verified, not just
individually run.

## Status

- [x] Runs end-to-end
- [x] Verified against the Nextflow version — exact match on cell/gene/cluster counts and DE hit counts per drug
- [ ] CI
