# rnaseq-stats-notebook

Statistical modelling built on **real DESeq2 output** from
[rnaseq-snakemake](../../pipelines/rnaseq-snakemake/) (yeast RNA-seq,
GSE110004, `WT` vs `RAP1_IAA`). No synthetic data anywhere in it.

Each method is reimplemented from scratch in Python and validated against
DESeq2's own numbers. The point is to make the modelling assumptions explicit
and auditable instead of delegating them to a package's defaults.

## Contents (`hypothesis_testing.ipynb`)

1. **Hypothesis testing.** t-test and Wilcoxon rank-sum against DESeq2's own
   per-gene p-values, and why RNA-seq needs a purpose-built model rather than
   off-the-shelf tests on small-n discrete count data.
2. **Negative binomial GLM.** The model DESeq2 actually fits, implemented from
   scratch via `statsmodels`, isolating what empirical Bayes shrinkage does
   and does not contribute. This section includes a correction I made
   mid-analysis, when the initial results contradicted my first hypothesis;
   the notebook's own markdown has the full account.
3. **Multiple testing correction.** Bonferroni against Benjamini-Hochberg/FDR,
   with a from-scratch BH implementation that reproduces DESeq2's `padj`
   column exactly.
4. **PCA and clustering.** Sample-level QC, PCA built from scratch via SVD,
   plus a correlation clustermap. The two methods disagreed here, and the
   notebook says so and works through why.
5. **Standard DE visualizations.** MA plot, volcano plot, and a heatmap of the
   significant genes.

## Setup

```bash
conda create -n stats -c conda-forge jupyter pandas scipy statsmodels matplotlib seaborn -y
conda activate stats
```

## Run

```bash
cd stats/rnaseq-stats-notebook
jupyter lab
```

Open `hypothesis_testing.ipynb` and run all cells. It reads directly from
`../../pipelines/rnaseq-snakemake/results/deseq2/`, so run that pipeline first
if those files are not there yet.

## Status

- [x] All 4 planned sections complete, executed outputs committed
- [x] Every from-scratch method validated against DESeq2's real output
