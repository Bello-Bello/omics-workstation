# scrnaseq-snakemake

[![scrnaseq-snakemake CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml)

Snakemake port of [scrnaseq-nextflow's `perturbation.nf`](../scrnaseq-nextflow/perturbation.nf).
Same public drug-screen data, same QC → normalize → cluster →
differential-expression logic, different workflow manager. Same relationship
as `rnaseq-snakemake` and `rnaseq-nextflow`, except applied to the
perturbation analysis rather than to the STARsolo quantification stage.

## Why not port the STARsolo stage too

`scrnaseq-nextflow/main.nf` quantifies test-scale chr19 data, and its whole
reason for existing is to validate that machinery, including a
macOS-specific STAR bug I hit and worked around via Docker (documented in
[that pipeline's README](../scrnaseq-nextflow/README.md)). Rebuilding it in
Snakemake would mostly duplicate the same STAR invocation and the same Docker
requirement without adding signal. The perturbation analysis is where the
interesting content is: real biology, differential expression, dose response.
So that is the half built twice.

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
                                   # rather than imported cross-directory. Same reasoning as
                                   # deseq2.R existing as two independent copies in the bulk
                                   # pipelines: each pipeline stays self-contained, the
                                   # calling convention differs (CLI args vs snakemake@
                                   # object), the underlying logic does not.
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

Same dataset and same environment spec as the Nextflow version. The first run
downloads the 470 MB sci-Plex dataset and builds the conda env; later runs
reuse both.

## Output

`results/perturbation/`: `processed.h5ad`, `differential_expression.csv`,
`summary.json`, `umap_condition_cluster.png`. Identical in shape to
`scrnaseq-nextflow`'s output.

## Cross-verification

Both pipelines pin the same package versions (`envs/scanpy.yaml` is
byte-identical between them) and Scanpy's clustering, PCA and UMAP functions
default to fixed random seeds, so comparing the two runs is meaningful rather
than coincidental.

Diffing both runs' `summary.json` gives an **exact match**: 23,966 cells by
35,916 genes post-QC, 7 Leiden clusters, and identical significant-DE-gene
counts per drug (Dex 5,415, Nutlin 6,689, BMS 6,291, SAHA 8,037). Same kind of
result as the bulk pipelines' matching differential expression calls on the
yeast experiment.

## Status

- [x] Runs end-to-end
- [x] Verified against the Nextflow version, exact match on cell/gene/cluster counts and DE hit counts per drug
- [x] CI
