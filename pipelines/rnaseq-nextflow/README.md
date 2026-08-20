# rnaseq-nextflow

[![rnaseq-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml)

Nextflow (DSL2) port of [rnaseq-snakemake](../rnaseq-snakemake/) — same biology,
same tools, same real yeast dataset (GSE110004, WT vs RAP1_IAA), different
workflow manager. See the [skills roadmap](../../ROADMAP.md).

## Structure

```
main.nf              # workflow: wires processes together via channels
nextflow.config       # params + engine settings (conda/docker profiles, cpu defaults)
modules/               # one process per pipeline step, each with conda + container directives
bin/deseq2.R           # auto-added to PATH by Nextflow's bin/ convention
envs/                   # conda environments (used by default / -profile conda)
docker/deseq2.Dockerfile  # custom image for the 4-package R/Bioconductor step (build locally, see below)
config/samples.tsv      # same sample sheet format as the Snakemake version
resources/               # real yeast reads + reference (gitignored, self-contained copy)
```

## Setup

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow -y
conda activate nextflow
```

Docker is optional (see **Run** below) but if you want it, install
[Docker Desktop](https://www.docker.com/products/docker-desktop/) first, then
build the one custom image (everything else pulls pre-built BioContainers
images automatically):

```bash
docker build -t omics-deseq2:1.0 -f docker/deseq2.Dockerfile .
```

## Run

**Conda (default):**
```bash
cd pipelines/rnaseq-nextflow
nextflow run main.nf
```

**Docker:**
```bash
cd pipelines/rnaseq-nextflow
nextflow run main.nf -profile docker
```

Both paths run the identical pipeline logic against identical containers/envs
pinned to explicit versions — the only process not sourced from a public,
pre-built image is `DESEQ2`, which uses the locally-built `omics-deseq2:1.0`
(see **Setup**) since its 4 R/Bioconductor packages aren't bundled together
in any single BioContainers image.

Nextflow builds each process's conda environment (or pulls its container
image) on first run — a one-time cost, cached under `work/` afterward.

## A reproducibility caveat worth knowing

Running this pipeline via conda vs Docker on the same machine can produce
**very slightly different results for genes sitting right at a significance
threshold** — in one comparison, 2 borderline genes (padj ≈ 0.05–0.06) flipped
across `padj < 0.05` between the two runs, while every clearly-significant
gene matched exactly. The cause: conda resolves an **osx-64** build of
`salmon` on macOS, while every Docker container is necessarily **linux-64**
(Docker Desktop runs a Linux VM even on Mac) — two genuinely different
compiled binaries for the same nominal software version, which can produce
tiny floating-point differences in salmon's EM-based multi-mapping-read
resolution. This isn't a bug or a misconfiguration: it's a real illustration
that **matching a software version number doesn't guarantee bit-identical
results — true reproducibility needs the same platform/architecture, not
just the same tool version.** It's also a genuine argument in Docker's
favor over conda: a container's results are reproducible *across any
machine*, since the platform itself is fixed by the image, whereas a
conda-resolved environment inherits whatever platform conda is solving for
locally.

Each container tag here is explicit and pinned (not `latest`), so within
either path — conda or Docker — results are fully deterministic and
reproducible on their own; the two paths just aren't guaranteed to match
each other exactly at the margin.

## Output

Same shape as the Snakemake version, under `results/`:
`fastqc/`, `trimmed/`, `salmon/`, `multiqc/`, `deseq2/results.csv` + `pca.png`.

## Status

- [x] Runs end-to-end (conda)
- [x] Runs end-to-end (Docker)
- [x] Compared against Snakemake version — same DE results (conda-to-conda); see the reproducibility caveat above for the conda-vs-Docker comparison
