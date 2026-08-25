# scrnaseq-nextflow

Two workflows sharing one pipeline directory:

- **`main.nf`** — STARsolo quantification from raw 10x reads (test-scale chr19 data, see [status](#status))
- **`perturbation.nf`** — QC → normalize → cluster → differential expression on a real public drug-screen count matrix (Srivatsan et al. 2020 sci-Plex, 24k A549 cells, 4 drugs)

## Why two entrypoints instead of one workflow?

Real drug-screen sequencing (sci-Plex, sci-RNA-seq) produces hundreds of GB of raw FASTQ per screen — reprocessing one from scratch isn't practical for a portfolio pipeline, and in practice a lot of real perturbation-analysis work starts from a delivered count matrix rather than re-running alignment per project anyway. `perturbation.nf` reflects that: it starts from a public, already-quantified dataset and does the actual biological analysis. `main.nf` is kept separately to prove the FASTQ→matrix quantification machinery itself works (see its own section below) — the two entrypoints answer different questions and are honest about which one is running on toy data vs. real data.

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
config/samples.tsv                   # Sample sheet — main.nf only
analysis/
  scanpy_analysis.ipynb              # QC/clustering/marker-gene mechanics on main.nf's output
  perturbation_case_study.ipynb      # Full walkthrough on the real drug-screen data
resources/perturbation/              # Downloaded dataset cache (gitignored, ~470 MB)
```

## Perturbation workflow — setup & run

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow -y
conda activate nextflow
cd pipelines/scrnaseq-nextflow
nextflow run perturbation.nf
```

No Docker profile for this workflow — it's pure Python/Scanpy with no platform-specific binary behavior (unlike STARsolo below), so plain conda is reproducible enough on its own. First run downloads the dataset (~470 MB, cached under `resources/perturbation/` afterward) then runs the analysis; on this machine, QC through DE end-to-end took about 2 minutes once the data was local.

## Quantification workflow (`main.nf`) — setup & run

**On macOS: use `-profile docker`, not conda — required, not optional.**
STAR 2.7.11b's conda/native build has a real, documented macOS Sonoma+ bug
([alexdobin/STAR#2142](https://github.com/alexdobin/STAR/issues/2142)) in
its shared-memory handling that breaks `--soloFeatures Gene` quantification
entirely (`Transcriptome.cpp: could not open input file /geneInfo.tab`,
regardless of genome index setup — confirmed by testing with no working
flag-based fix, only a source patch + recompile). The Linux container
sidesteps it completely since the bug is in macOS-specific shared-memory
syscall behavior, not the STAR version itself:

```bash
nextflow run main.nf -profile docker
```

On Linux (including CI), plain conda works fine and is untested-but-expected
to work given the bug is macOS-specific:
```bash
nextflow run main.nf
```

## Output

**`perturbation.nf`** (`results/perturbation/`):
- `processed.h5ad` — filtered, normalized, clustered AnnData (23,966 cells x 35,916 genes on the real run)
- `differential_expression.csv` — one row per gene per drug, Wilcoxon test vs. vehicle control (logfoldchange, p-value, adjusted p-value)
- `summary.json` — cell/gene/cluster counts and significant-DE-gene counts per drug
- `umap_condition_cluster.png` — UMAP colored by drug and by Leiden cluster

**`main.nf`** (`results/starsolo/`):
- `{sample}.Solo.out/Gene/raw/` — sparse count matrix triplet per sample
- `{sample}.Log.final.out` — real STAR mapping statistics

## The `perturbation/` package

Shared by `perturbation.nf`, the Snakemake port ([../scrnaseq-snakemake/](../scrnaseq-snakemake/)), and the notebook — the same QC/normalization/clustering/DE logic runs from all three rather than being copy-pasted into each. Column names for a new dataset (which column holds the perturbation label, the dose, what the control is called) are the only thing that changes between datasets — passed as CLI args / Nextflow params / Snakemake config, not hardcoded in the package itself.

## Dataset

Srivatsan et al. 2020, *Science* — sci-Plex: A549 (human lung adenocarcinoma) cells exposed to one of four compounds (dexamethasone, nutlin-3a, BMS-345541, vorinostat/SAHA) across seven doses in triplicate, plus vehicle controls. Fetched via [pertpy](https://pertpy.readthedocs.io/)'s curated copy of the [scPerturb](http://projects.sanderlab.org/scperturb/) release of this dataset.

## Status

- [x] STARsolo quantification runs end-to-end (Docker; real macOS-specific STAR bug found, diagnosed, and worked around — documented above)
- [x] Downstream Scanpy analysis on `main.nf` output (QC, filtering, clustering, marker genes) — [analysis/scanpy_analysis.ipynb](analysis/scanpy_analysis.ipynb)
- [x] Perturbation screen analysis on real public drug-screen data (QC, clustering, differential expression) — [analysis/perturbation_case_study.ipynb](analysis/perturbation_case_study.ipynb)
- [x] Cross-verified against an independent Snakemake implementation — [../scrnaseq-snakemake/](../scrnaseq-snakemake/): both produce **exactly matching** results on the same real data — 23,966 cells x 35,916 genes post-QC, 7 Leiden clusters, and identical significant-DE-gene counts per drug (Dex: 5,415, Nutlin: 6,689, BMS: 6,291, SAHA: 8,037)
- [ ] CI
