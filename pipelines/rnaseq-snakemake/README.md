# rnaseq-snakemake

A bulk RNA-seq pipeline: **FastQC → fastp → Salmon → MultiQC → DESeq2**.

This is flagship pipeline #1 (Snakemake implementation) in the
[skills roadmap](../../ROADMAP.md), reimplemented in Nextflow in
[pipelines/rnaseq-nextflow/](../rnaseq-nextflow/) — same biology, same real
data, different workflow manager, results cross-verified between the two.
The statistical modeling built on this pipeline's real output lives in
[stats/rnaseq-stats-notebook/](../../stats/rnaseq-stats-notebook/).

Runs on real yeast RNA-seq data (GSE110004): 3 wild-type (`WT`) vs 3
Rap1-transcription-factor-depleted (`RAP1_IAA`) biological replicates.

## Why this pipeline

It's the standard entry point for demonstrating reproducible-workflow skills:
short runtime, well-understood biology, and every step maps to a concept
worth understanding (QC, trimming, pseudo-alignment/quantification,
differential expression via a negative binomial GLM).

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
2. Get the real test data (paired-end yeast reads + reference transcriptome,
   ~26 MB total, sourced from the public `nf-core/test-datasets` GitHub repo)
   into `resources/reads/` and `resources/reference/` — file names and paths
   must match `config/samples.tsv` and `config/config.yaml`.

## Run

```bash
cd pipelines/rnaseq-snakemake
snakemake --use-conda --cores 4 --snakefile workflow/Snakefile
```

Add `-n` first to dry-run and check the DAG before actually executing:
```bash
snakemake --use-conda --cores 4 -n --snakefile workflow/Snakefile
```

Visualize the DAG (useful for your portfolio README/screenshots):
```bash
snakemake --dag --snakefile workflow/Snakefile | dot -Tpng > dag.png
```

### Running via containers instead of conda

Every rule also has a `container:` directive (pinned BioContainers images,
same tags as the Nextflow version — see that pipeline's README for the exact
list) as an alternative to conda. This path needs
[Apptainer](https://apptainer.org/) (formerly Singularity), which is Linux/HPC
tooling with real friction on native macOS — **not tested locally in this
repo**; conda is the verified, working path here. On Linux/HPC:

```bash
snakemake --software-deployment-method apptainer --cores 4 --snakefile workflow/Snakefile
```

The `deseq2` rule's container needs the custom image built first (see
[rnaseq-nextflow/docker/deseq2.Dockerfile](../rnaseq-nextflow/docker/deseq2.Dockerfile) —
same image, shared between both pipelines), referenced as
`docker-daemon://omics-deseq2:1.0` so Apptainer uses the local Docker cache
rather than pulling from a registry.

**Reproducibility note:** the Nextflow version of this pipeline was run both
via conda and via Docker on the same machine — results matched exactly for
every clearly-significant gene, but 2 borderline genes (padj ≈ 0.05–0.06)
flipped across the `padj < 0.05` line between the two runs. Cause: conda
resolves an `osx-64` build of `salmon` on macOS, while Docker/Apptainer
containers are always `linux-64` — genuinely different compiled binaries for
the same nominal version, producing tiny floating-point differences. Not a
bug — a real illustration that matching a version *number* doesn't guarantee
bit-identical results; true reproducibility needs the same platform, not just
the same tool version. Full writeup in the Nextflow pipeline's README.

## Output

- `results/fastqc/` — raw read QC reports
- `results/trimmed/` — trimmed reads + fastp reports
- `results/salmon/{sample}/quant.sf` — per-sample transcript quantification
- `results/multiqc/multiqc_report.html` — aggregated QC across all steps
- `results/deseq2/results.csv` — differential expression results
- `results/deseq2/pca.png` — sample PCA plot

## Status

- [x] Runs end-to-end on real yeast RNA-seq data (6 samples, `WT` vs `RAP1_IAA`)
- [x] Ported to Nextflow (`pipelines/rnaseq-nextflow/`) — cross-verified, matching DE results
- [x] Statistical modeling module built on this pipeline's real output (`stats/rnaseq-stats-notebook/`)
- [x] `container:` directives added for all rules (Apptainer/HPC path, documented, not locally tested)
- [ ] Real screenshot of the PCA/MultiQC output in this README
