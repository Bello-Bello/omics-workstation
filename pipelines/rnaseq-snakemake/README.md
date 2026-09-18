# rnaseq-snakemake

[![rnaseq-snakemake CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-snakemake-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-snakemake-ci.yml)

A bulk RNA-seq pipeline: **FastQC → fastp → Salmon → MultiQC → DESeq2**.

This is the Snakemake implementation. It is also built in Nextflow at
[pipelines/rnaseq-nextflow/](../rnaseq-nextflow/), with results cross-verified
between the two. The statistical modelling built on this pipeline's output
lives in [stats/rnaseq-stats-notebook/](../../stats/rnaseq-stats-notebook/).

Runs on yeast RNA-seq data (GSE110004): 3 wild-type (`WT`) against 3
Rap1-transcription-factor-depleted (`RAP1_IAA`) biological replicates.

## Why this dataset

Short runtime and well-characterised biology. That combination makes it a
tractable basis for cross-verification between two workflow managers, because
any divergence between the implementations can be pinned on the tooling
rather than on ambiguity in the underlying result.

## Structure

```
config/
  config.yaml     # paths, tool params, thread counts
  samples.tsv     # sample sheet: sample, condition, fq1, fq2
workflow/
  Snakefile       # entry point, wires the rules together
  rules/          # one .smk file per pipeline stage
  scripts/        # deseq2.R — the actual differential expression analysis
  envs/           # per-rule conda environments (keeps tool versions pinned)
```

## Setup

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html), then
   create the Snakemake environment:
   ```bash
   conda create -n snakemake -c bioconda -c conda-forge snakemake -y
   conda activate snakemake
   ```
2. Get the test data (paired-end yeast reads plus reference transcriptome,
   about 26 MB total, from the public `nf-core/test-datasets` GitHub repo)
   into `resources/reads/` and `resources/reference/`. File names and paths
   have to match `config/samples.tsv` and `config/config.yaml`.

## Run

```bash
cd pipelines/rnaseq-snakemake
snakemake --use-conda --cores 4 --snakefile workflow/Snakefile
```

Add `-n` first to dry-run and check the DAG before executing:
```bash
snakemake --use-conda --cores 4 -n --snakefile workflow/Snakefile
```

Visualize the DAG:
```bash
snakemake --dag --snakefile workflow/Snakefile | dot -Tpng > dag.png
```

### Running via containers instead of conda

Every rule also carries a `container:` directive (pinned BioContainers images,
same tags as the Nextflow version, which lists them). That path needs
[Apptainer](https://apptainer.org/), formerly Singularity, which is Linux and
HPC tooling with real friction on native macOS. It is **not tested locally in
this repo**; conda is the verified path here. On Linux or HPC:

```bash
snakemake --software-deployment-method apptainer --cores 4 --snakefile workflow/Snakefile
```

The `deseq2` rule's container needs the custom image built first. See
[rnaseq-nextflow/docker/deseq2.Dockerfile](../rnaseq-nextflow/docker/deseq2.Dockerfile),
which is the same image shared between both pipelines. It is referenced as
`docker-daemon://omics-deseq2:1.0` so that Apptainer uses the local Docker
cache instead of pulling from a registry.

**Reproducibility note.** The Nextflow version of this pipeline was run both
via conda and via Docker on the same machine. Results matched exactly for
every clearly significant gene, but 2 borderline genes (padj around 0.05 to
0.06) crossed the `padj < 0.05` line between the two runs. The cause is that
conda resolves an `osx-64` build of `salmon` on macOS, while Docker and
Apptainer containers are always `linux-64`. Those are different compiled
binaries for the same nominal version, and the tiny floating-point
differences are enough to move a marginal gene. Matching a version *number*
does not guarantee bit-identical results; that needs the same platform. Full
writeup in
[the Nextflow pipeline's README](../rnaseq-nextflow/README.md#conda-and-docker-can-disagree-at-the-margin).

## Output

- `results/fastqc/` — raw read QC reports
- `results/trimmed/` — trimmed reads + fastp reports
- `results/salmon/{sample}/quant.sf` — per-sample transcript quantification
- `results/multiqc/multiqc_report.html` — aggregated QC across all steps
- `results/deseq2/results.csv` — differential expression results
- `results/deseq2/pca.png` — sample PCA plot

## Status

- [x] Runs end-to-end on yeast RNA-seq data (6 samples, `WT` vs `RAP1_IAA`)
- [x] Ported to Nextflow (`pipelines/rnaseq-nextflow/`), cross-verified, matching DE results
- [x] Statistical modelling module built on this pipeline's output (`stats/rnaseq-stats-notebook/`)
- [x] `container:` directives added for all rules (Apptainer/HPC path, documented, not locally tested)
- [ ] Real screenshot of the PCA/MultiQC output in this README
