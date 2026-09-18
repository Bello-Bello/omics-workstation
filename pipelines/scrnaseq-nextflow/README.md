# scrnaseq-nextflow

[![scrnaseq-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/scrnaseq-nextflow-ci.yml)

Two workflows sharing one pipeline directory:

- **`main.nf`** — STARsolo quantification from raw 10x reads (test-scale chr19 data, see [status](#status))
- **`perturbation.nf`** — QC → normalize → cluster → differential expression on a public drug-screen count matrix (Srivatsan et al. 2020 sci-Plex, 24k A549 cells, 4 drugs)

## Why two entry points instead of one workflow

Real drug-screen sequencing produces hundreds of GB of raw FASTQ per screen.
Reprocessing one from scratch is impractical at that scale, and in practice a
lot of perturbation analysis starts from a delivered count matrix anyway
rather than re-running alignment per project.

`perturbation.nf` reflects that. It starts from a public, already-quantified
dataset and does the biological analysis. `main.nf` exists separately to prove
the FASTQ-to-matrix machinery works. The two answer different questions, and
keeping them apart makes it obvious which one runs on toy data and which runs
on real data.

## Structure

```
main.nf                          # STARsolo quantification (test data)
perturbation.nf                  # Perturbation analysis (real drug-screen data)
nextflow.config                  # Params for both workflows
modules/
  star_index.nf, starsolo.nf         # Quantification
  fetch_perturbation_data.nf         # Downloads the public dataset via pertpy
  perturbation_analysis.nf           # QC -> normalize -> cluster -> DE
bin/
  fetch_perturbation_data.py         # Called by fetch_perturbation_data.nf
  run_perturbation_pipeline.py       # Called by perturbation_analysis.nf
perturbation/                        # Reusable Python package (not dataset-specific)
  metadata.py                        # Standardizes condition/dose columns
  qc.py                              # Per-cell QC metrics + filtering
  preprocessing.py                   # Normalization, HVG, PCA, Leiden, UMAP
  analysis.py                        # Condition-vs-control differential expression
  utils.py                           # Logging, AnnData I/O
envs/
  star.yaml                          # STAR only (main.nf)
  scanpy.yaml                        # scanpy + pertpy (perturbation.nf)
config/samples.tsv                   # Sample sheet, main.nf only
analysis/
  scanpy_analysis.ipynb              # QC/clustering/marker-gene mechanics on main.nf's output
  perturbation_case_study.ipynb      # Full analysis of the drug-screen data
resources/perturbation/              # Downloaded dataset cache (gitignored, ~470 MB)
```

## Perturbation workflow: setup and run

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow -y
conda activate nextflow
cd pipelines/scrnaseq-nextflow
nextflow run perturbation.nf
```

There is no Docker profile for this workflow. It is pure Python and Scanpy
with no platform-specific binary behaviour, unlike STARsolo below, so plain
conda is reproducible enough on its own.

The first run downloads the dataset (about 470 MB, cached under
`resources/perturbation/` afterward) and then runs the analysis. On this
machine, QC through DE took roughly 2 minutes once the data was local.

## Quantification workflow (`main.nf`): setup and run

**On macOS, use `-profile docker` rather than conda. This is required, not
optional.** STAR 2.7.11b's conda and native build has a documented macOS
Sonoma+ bug ([alexdobin/STAR#2142](https://github.com/alexdobin/STAR/issues/2142))
in its shared-memory handling, which breaks `--soloFeatures Gene`
quantification outright (`Transcriptome.cpp: could not open input file
/geneInfo.tab`, regardless of how the genome index is set up). I tested for a
flag-based workaround and there isn't one; the only fix is patching the source
and recompiling. The Linux container sidesteps it, because the bug is in
macOS-specific shared-memory syscall behaviour rather than in STAR itself.

```bash
nextflow run main.nf -profile docker
```

On Linux, including CI, plain conda works:
```bash
nextflow run main.nf
```

## Output

**`perturbation.nf`** (`results/perturbation/`):
- `processed.h5ad` — filtered, normalized, clustered AnnData (23,966 cells x 35,916 genes on the real run)
- `differential_expression.csv` — one row per gene per drug, Wilcoxon test vs vehicle control (logfoldchange, p-value, adjusted p-value)
- `summary.json` — cell/gene/cluster counts and significant-DE-gene counts per drug
- `umap_condition_cluster.png` — UMAP coloured by drug and by Leiden cluster

**`main.nf`** (`results/starsolo/`):
- `{sample}.Solo.out/Gene/raw/` — sparse count matrix triplet per sample
- `{sample}.Log.final.out` — STAR mapping statistics

## The `perturbation/` package

Shared by `perturbation.nf`, the Snakemake port
([../scrnaseq-snakemake/](../scrnaseq-snakemake/)) and the notebook, so the
same QC, normalization, clustering and DE logic runs from all three instead of
being copy-pasted into each.

Moving to a new dataset only means changing column names: which column holds
the perturbation label, which holds the dose, what the control is called.
Those are passed as CLI args, Nextflow params or Snakemake config, not
hardcoded in the package.

## Dataset

Srivatsan et al. 2020, *Science*, sci-Plex. A549 human lung adenocarcinoma
cells exposed to one of four compounds (dexamethasone, nutlin-3a, BMS-345541,
vorinostat/SAHA) across seven doses in triplicate, plus vehicle controls.
Fetched via [pertpy](https://pertpy.readthedocs.io/)'s curated copy of the
[scPerturb](http://projects.sanderlab.org/scperturb/) release.

## Status

- [x] STARsolo quantification runs end-to-end (Docker; macOS-specific STAR bug found, diagnosed and worked around, documented above)
- [x] Downstream Scanpy analysis on `main.nf` output (QC, filtering, clustering, marker genes), [analysis/scanpy_analysis.ipynb](analysis/scanpy_analysis.ipynb)
- [x] Perturbation screen analysis on public drug-screen data (QC, clustering, differential expression), [analysis/perturbation_case_study.ipynb](analysis/perturbation_case_study.ipynb)
- [x] Cross-verified against an independent Snakemake implementation, [../scrnaseq-snakemake/](../scrnaseq-snakemake/). Both produce **exactly matching** results on the same data: 23,966 cells x 35,916 genes post-QC, 7 Leiden clusters, and identical significant-DE-gene counts per drug (Dex 5,415, Nutlin 6,689, BMS 6,291, SAHA 8,037)
- [x] CI
