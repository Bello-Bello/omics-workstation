# rnaseq-nextflow

[![rnaseq-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/rnaseq-nextflow-ci.yml)

Nextflow (DSL2) port of [rnaseq-snakemake](../rnaseq-snakemake/). Same biology,
same tools, same yeast dataset (GSE110004, WT vs RAP1_IAA), different workflow
manager. Results are cross-verified against the Snakemake implementation.

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

Docker is optional (see **Run** below). If you want it, install
[Docker Desktop](https://www.docker.com/products/docker-desktop/) first, then
build the one custom image. Everything else pulls pre-built BioContainers
images automatically.

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

Both paths run identical pipeline logic against containers and environments
pinned to explicit versions. The one process not sourced from a public
pre-built image is `DESEQ2`, which uses the locally-built `omics-deseq2:1.0`
(see **Setup**), because its four R/Bioconductor packages are not bundled
together in any single BioContainers image.

Nextflow builds each process's conda environment, or pulls its container
image, on first run. That is a one-time cost, cached under `work/` afterward.

## Conda and Docker can disagree at the margin

Running this pipeline via conda and via Docker on the same machine can give
**slightly different results for genes sitting right at a significance
threshold**. In one comparison, 2 borderline genes (padj around 0.05 to 0.06)
crossed the `padj < 0.05` line between the two runs. Every clearly significant
gene matched exactly.

The cause is the build, not the configuration. Conda resolves an **osx-64**
build of `salmon` on macOS, while every Docker container is **linux-64**
(Docker Desktop runs a Linux VM even on a Mac). Those are different compiled
binaries for the same nominal software version, and they can produce tiny
floating-point differences in salmon's EM-based resolution of multi-mapping
reads.

So matching a software version number does not guarantee bit-identical
results. Reproducibility needs the same platform and architecture, not just
the same tool version. This is also a point in Docker's favour over conda for
this kind of work: a container's results reproduce on any machine because the
image fixes the platform, whereas a conda-resolved environment inherits
whatever platform conda is solving for locally.

Each container tag here is pinned rather than `latest`, so within either path
results are deterministic and reproducible on their own. The two paths just
are not guaranteed to match each other exactly at the margin.

## Output

Same shape as the Snakemake version, under `results/`:
`fastqc/`, `trimmed/`, `salmon/`, `multiqc/`, `deseq2/results.csv` + `pca.png`.

## Status

- [x] Runs end-to-end (conda)
- [x] Runs end-to-end (Docker)
- [x] Compared against Snakemake version, same DE results (conda-to-conda). See the section above for the conda-vs-Docker comparison.
