# rnaseq-snakemake

A bulk RNA-seq pipeline: **FastQC → fastp → Salmon → MultiQC → DESeq2**.

This is flagship pipeline #1 (Snakemake implementation) in the
[skills roadmap](../../ROADMAP.md). It will be reimplemented in Nextflow
in `pipelines/rnaseq-nextflow/` once this version runs cleanly end to end.

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

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) and
   `mamba` (faster solver), then create the Snakemake environment:
   ```bash
   mamba create -n snakemake -c bioconda -c conda-forge snakemake
   conda activate snakemake
   ```
2. Get test data — use the small subset from
   [nf-core/test-datasets (rnaseq branch)](https://github.com/nf-core/test-datasets/tree/rnaseq)
   (chr21-only genome/transcriptome + a handful of paired-end FASTQs, runs in
   minutes on a laptop). Place reads under `resources/reads/` and the
   reference under `resources/reference/`, matching the paths in
   `config/config.yaml` and `config/samples.tsv`.
3. Edit `config/samples.tsv` to match your actual sample names/paths.

## Run

```bash
cd pipelines/rnaseq-snakemake
snakemake --use-conda --cores 4 --directory . --snakefile workflow/Snakefile
```

Add `-n` first to dry-run and check the DAG before actually executing:
```bash
snakemake --use-conda --cores 4 -n --snakefile workflow/Snakefile
```

Visualize the DAG (useful for your portfolio README/screenshots):
```bash
snakemake --dag --snakefile workflow/Snakefile | dot -Tpng > dag.png
```

## Output

- `results/fastqc/` — raw read QC reports
- `results/trimmed/` — trimmed reads + fastp reports
- `results/salmon/{sample}/quant.sf` — per-sample transcript quantification
- `results/multiqc/multiqc_report.html` — aggregated QC across all steps
- `results/deseq2/results.csv` — differential expression results
- `results/deseq2/pca.png` — sample PCA plot

## Status

- [ ] Runs end-to-end on nf-core test data
- [ ] Sample sheet + config generalized beyond the 4-sample example
- [ ] README includes a real screenshot of the PCA/MultiQC output
- [ ] Ported to Nextflow (`pipelines/rnaseq-nextflow/`)
