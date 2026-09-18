# Omics Workstation

[![rnaseq-snakemake CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-snakemake-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-snakemake-ci.yml)
[![rnaseq-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml)
[![scrnaseq-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-nextflow-ci.yml)
[![scrnaseq-snakemake CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-snakemake-ci.yml)
[![amplicon-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/amplicon-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/amplicon-nextflow-ci.yml)

Five pipelines covering bulk RNA-seq, single-cell RNA-seq and 16S amplicon
data. Each one runs on a published dataset, and each has CI that re-downloads
the data and re-executes the whole workflow on every push instead of linting
it.

The result worth starting with: the same yeast RNA-seq experiment (GSE110004,
wild-type against Rap1 transcription-factor depletion) analysed by two
independent workflow managers, returning the same differential expression
calls.

![Sample PCA of WT vs RAP1_IAA](assets/pca_wt_vs_rap1iaa.png)

## What's here

**[pipelines/rnaseq-snakemake/](pipelines/rnaseq-snakemake/)**
Bulk RNA-seq in Snakemake: FastQC, fastp, Salmon, MultiQC, DESeq2. Runs under
conda, with a documented (though locally untested) Apptainer path.

**[pipelines/rnaseq-nextflow/](pipelines/rnaseq-nextflow/)**
The same pipeline rewritten in Nextflow DSL2, verified under both conda and
Docker. Same data, matching results.

**[stats/rnaseq-stats-notebook/](stats/rnaseq-stats-notebook/)**
Statistical modelling on that pipeline's DESeq2 output: hypothesis testing,
the negative binomial GLM that DESeq2 actually fits, multiple testing
correction, PCA and clustering. Every method is implemented from scratch and
checked against DESeq2's own numbers.

**[pipelines/scrnaseq-nextflow/](pipelines/scrnaseq-nextflow/)**
Two entry points. STARsolo quantification from raw 10x reads, which is where
a macOS-specific STAR bug turned up and had to be diagnosed and worked
around. And a perturbation-screen case study on a public drug-dose-response
dataset (Srivatsan et al. 2020, 24k cells, 4 drugs), taken from QC through
clustering, differential expression and dose response. That analysis surfaced
a cytotoxicity confound, which the writeup explains instead of glossing.

**[pipelines/scrnaseq-snakemake/](pipelines/scrnaseq-snakemake/)**
The perturbation analysis rebuilt in Snakemake. Cell, gene and cluster counts
and per-drug DE hit counts match the Nextflow version exactly.

**[pipelines/amplicon-nextflow/](pipelines/amplicon-nextflow/)**
16S amplicon on Atacama Desert soil data. QIIME2 in a pinned container handles
demultiplexing, DADA2 denoising and taxonomic classification. The
compositional statistics downstream (CLR, Aitchison distance, PERMANOVA,
depth-aware differential abundance) are written from first principles and
cross-checked against scikit-bio and statsmodels. The biology comes out where
you would expect: *Bradyrhizobium* and *Candidatus* Nitrososphaera enriched in
vegetated soil, desiccation-resistant *Rubrobacter* in barren. Vegetation
accounts for 15.2% of community variation (PERMANOVA p = 0.001), and two
independent executions produced byte-identical results.

## Quick start

Each component is self-contained and carries its own setup and run
instructions. Start wherever is relevant:
[rnaseq-snakemake](pipelines/rnaseq-snakemake/README.md) ·
[rnaseq-nextflow](pipelines/rnaseq-nextflow/README.md) ·
[stats notebook](stats/rnaseq-stats-notebook/README.md) ·
[scrnaseq-nextflow](pipelines/scrnaseq-nextflow/README.md) ·
[scrnaseq-snakemake](pipelines/scrnaseq-snakemake/README.md) ·
[amplicon-nextflow](pipelines/amplicon-nextflow/README.md)

## Why two workflow managers for the same pipeline

Building an analysis twice, in Snakemake and in Nextflow, and getting the same
answer is a stronger claim than running either one carefully. It shows the
result belongs to the analysis and not to one tool's behaviour. I have done it
twice: once on the bulk RNA-seq pipeline, once on the perturbation screen.

The bulk pair turned up something I did not expect. Conda and Docker builds of
the same tool version can differ at the platform level, because conda resolves
an `osx-64` build of Salmon on macOS while any container is `linux-64`. Two
borderline genes (padj around 0.05 to 0.06) crossed the significance line
between the two runs, while everything clearly significant matched. Full
writeup in the
[Nextflow pipeline's README](pipelines/rnaseq-nextflow/README.md#conda-and-docker-can-disagree-at-the-margin).

The perturbation pair ran clean. Both are conda-only, no Docker involved, and
they landed on an exact match. That matters for reading the first finding: the
caveat is specifically about conda against Docker, not a general weakness in
reproducing these pipelines, and the second pair is the control case that
confirms it.

Where a result surprised me or a tool misbehaved, it is written up in the
relevant pipeline's README.
