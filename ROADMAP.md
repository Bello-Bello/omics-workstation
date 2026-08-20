# Skills Roadmap — Bioinformatics Scientist/Analyst (2-3 month target)

Goal: land a bioinformatics scientist/analyst role. Background: wet-lab, new to coding.
Strategy: depth over breadth. One flagship pipeline built twice (Snakemake, then Nextflow),
a stats module derived from real output, then a second pipeline if time allows.
Proteomics is intentionally deferred — see note at the bottom.

## Week 1-3 — Foundations (parallel track, don't gate on finishing this)
- [ ] Unix/shell: navigating, pipes, grep/awk/sed basics, `screen`/`tmux`, job basics
- [ ] Git: init/add/commit/branch/push, writing a decent README
- [ ] Python: variables, functions, loops, `pandas` (read_csv, groupby, merge, filtering)
- [ ] R: vectors, data frames, `tidyverse` (dplyr, ggplot2) — needed for DESeq2 later
- Resources: see [resources.md](resources.md)

## Week 2-5 — Flagship Pipeline #1: Bulk RNA-seq
Build the SAME pipeline twice — this is the point, not a detour:
1. **Snakemake** version first (pure Python, gentler debugging) — scaffold in
   [pipelines/rnaseq-snakemake/](pipelines/rnaseq-snakemake/)
2. **Nextflow** version second (industry standard, nf-core familiarity matters a lot
   in job postings) — scaffold in [pipelines/rnaseq-nextflow/](pipelines/rnaseq-nextflow/)

Pipeline steps (same in both): FastQC → fastp (trim) → Salmon (quant) → MultiQC → DESeq2 (R)

- [x] Get test data (see rnaseq-snakemake README)
- [x] Snakemake version runs end-to-end on test data
- [x] Nextflow version runs end-to-end on test data
- [x] Both containerized (conda env) and documented — Docker still open, see Week 9-12

## Week 5-6 — Statistical Modeling Module
Use the DESeq2 output from the pipeline above as real data — not a toy dataset.
- [x] Hypothesis testing refresher (t-test, Wilcoxon, when each applies)
- [x] Negative binomial GLM — the actual model behind DESeq2, understand what it's doing
- [x] Multiple testing correction (Bonferroni vs BH/FDR) and why RNA-seq needs it
- [x] PCA / clustering for QC (sample outlier detection)
- [ ] Write this up as a notebook: [stats/rnaseq-stats-notebook/](stats/)

## Week 6-9 — Flagship Pipeline #2 and #3
Original plan said pick one; user opted to do both, scRNA-seq first. Noted as
a scope expansion beyond the original "depth over breadth" default — revisit
pace if it starts crowding out the polish/application-sending phase.

- [ ] **scRNA-seq** (up first): quantification → Scanpy (Python), QC/filtering,
      clustering, marker genes. Build in Nextflow (you already know it from
      Pipeline #1).
- [ ] **Variant calling (DNA-seq)** (after scRNA-seq): GATK best-practices style —
      align (BWA) → mark duplicates → call variants (GATK HaplotypeCaller) →
      filter/annotate. Look at nf-core/sarek as a reference implementation,
      don't just copy it.

## Week 9-12 — Polish for job applications
- [x] Push everything to GitHub — https://github.com/Bello-Bello/omics-workstation
- [x] Docker containers for both flagship pipelines — Nextflow fully runs via
      Docker (tested, `-profile docker`); Snakemake has `container:`
      directives added but runs via Apptainer/HPC tooling, not locally tested
      on macOS — conda remains its verified path (documented honestly in
      both pipeline READMEs, including a real conda-vs-Docker platform
      reproducibility caveat found while testing)
- [x] GitHub Actions CI: real end-to-end runs (not just lint) on every push —
      both pipelines download real test data fresh and run for real on
      GitHub's runners; caught and fixed one genuine CI-only bug (a strict
      conda channel-priority setting that broke the deseq2 env solve)
- [x] Write real READMEs: top-level + all 3 component READMEs rewritten,
      CI badges added, a real result image embedded in the top-level README
- [ ] Pin the repo, link from CV/LinkedIn
- [ ] Optional: write a short blog post per pipeline (signals communication skill)

## Proteomics — deferred
Different toolchain (DIA-NN/MaxQuant/Skyline) and would dilute the 3 domains above
in a 2-3 month window. If a specific job posting requires it, revisit then rather
than building it speculatively now.

## Sequencing rule of thumb
Don't start Pipeline #2 until Pipeline #1 runs cleanly end-to-end in BOTH workflow
managers. A finished, well-documented single pipeline beats two half-finished ones
in an application review.
