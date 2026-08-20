# rnaseq-stats-notebook

Statistical modeling module (Week 5-6 in the [skills roadmap](../../ROADMAP.md)) —
built entirely on **real DESeq2 output** from [rnaseq-snakemake](../../pipelines/rnaseq-snakemake/)
(real yeast RNA-seq, GSE110004, `WT` vs `RAP1_IAA`), not synthetic data.

Goal: understand what DESeq2 is actually computing well enough to explain it,
not just read its output — every method here is implemented from scratch in
Python and checked against DESeq2's own numbers.

## Contents (`hypothesis_testing.ipynb`)

1. **Hypothesis testing refresher** — t-test and Wilcoxon rank-sum vs DESeq2's
   own p-values per gene; why RNA-seq needs a purpose-built model rather than
   off-the-shelf tests on small-n discrete count data.
2. **Negative binomial GLM** — the actual model DESeq2 fits, implemented from
   scratch via `statsmodels`, isolating what empirical Bayes shrinkage does
   and doesn't contribute (a real correction made mid-analysis when initial
   results contradicted the first hypothesis — see the notebook's own
   markdown for the full account).
3. **Multiple testing correction** — Bonferroni vs Benjamini-Hochberg/FDR,
   with a from-scratch BH implementation validated to exactly reproduce
   DESeq2's `padj` column.
4. **PCA & clustering** — sample-level QC, PCA built from scratch via SVD,
   plus a correlation clustermap; includes an honest finding where the two
   methods disagreed and why.
5. **Standard DE visualizations** — MA plot, volcano plot, and a heatmap of
   the significant genes.

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

Open `hypothesis_testing.ipynb` and run all cells — it reads directly from
`../../pipelines/rnaseq-snakemake/results/deseq2/` (run that pipeline first
if those files aren't there yet).

## Status

- [x] All 4 planned sections complete, executed outputs committed
- [x] Every from-scratch method validated against DESeq2's real output
