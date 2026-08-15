# Omics Workstation

Working directory for building a bioinformatics scientist/analyst portfolio:
automation & reproducible workflows (Nextflow/Snakemake), statistical modeling,
and data manipulation.

Start here: **[ROADMAP.md](ROADMAP.md)** — the staged 12-week plan.
Learning materials: **[resources.md](resources.md)**.

## Layout

```
ROADMAP.md            # staged plan, check items off as you go
resources.md          # curated courses/docs per topic
pipelines/
  rnaseq-snakemake/    # flagship pipeline #1, Snakemake implementation (in progress)
  rnaseq-nextflow/      # flagship pipeline #1, Nextflow port (not started)
stats/                 # statistical modeling notebook, built on real pipeline output (not started)
```

## Current status

Working on: `pipelines/rnaseq-snakemake/` — bulk RNA-seq pipeline scaffold is in
place (FastQC → fastp → Salmon → MultiQC → DESeq2). Next step: pull test data
and get it running end to end. See that pipeline's own README for setup.
