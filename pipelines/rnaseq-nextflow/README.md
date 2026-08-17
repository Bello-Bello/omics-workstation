# rnaseq-nextflow

Nextflow (DSL2) port of [rnaseq-snakemake](../rnaseq-snakemake/) — same biology,
same tools, same real yeast dataset (GSE110004, WT vs RAP1_IAA), different
workflow manager. See the [skills roadmap](../../ROADMAP.md).

## Structure

```
main.nf              # workflow: wires processes together via channels
nextflow.config       # params + engine settings (conda, cpu defaults)
modules/               # one process per pipeline step
bin/deseq2.R           # auto-added to PATH by Nextflow's bin/ convention
envs/                   # same conda environments as the Snakemake version
config/samples.tsv      # same sample sheet format as the Snakemake version
resources/               # real yeast reads + reference (gitignored, self-contained copy)
```

## Setup

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow -y
conda activate nextflow
```

## Run

```bash
cd pipelines/rnaseq-nextflow
nextflow run main.nf
```

Nextflow builds each process's conda environment on first run (same one-time
cost as Snakemake's `--use-conda`), cached under `work/` afterward.

## Output

Same shape as the Snakemake version, under `results/`:
`fastqc/`, `trimmed/`, `salmon/`, `multiqc/`, `deseq2/results.csv` + `pca.png`.

## Status

- [ ] Runs end-to-end
- [ ] Compare against Snakemake version: same DE results?
