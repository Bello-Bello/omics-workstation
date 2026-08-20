# scrnaseq-nextflow

Flagship pipeline #2 in the [skills roadmap](../../ROADMAP.md) — single-cell
RNA-seq quantification (STARsolo) → downstream analysis in Scanpy (QC,
clustering, marker genes). Real 10x Genomics test data (mouse, chromosome 19
only, from nf-core/scrnaseq's own test fixtures) — not a designed
condition-vs-condition experiment like the bulk pipeline, this data validates
pipeline mechanics rather than testing a specific biological hypothesis.

## Structure

```
main.nf                      # workflow: groups multi-lane samples, runs STARsolo
nextflow.config                # params, conda default, docker profile
modules/
  star_index.nf                 # builds the STAR genome index
  starsolo.nf                    # per-sample barcode/UMI-aware quantification
envs/star.yaml                   # conda env (STAR only)
config/samples.tsv                # sample sheet — one row per (sample, lane)
resources/                         # real 10x reads + chr19 reference + barcode whitelist (gitignored)
```

## Setup

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow -y
conda activate nextflow
```

## Run

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
cd pipelines/scrnaseq-nextflow
nextflow run main.nf -profile docker
```

On Linux (including CI), plain conda works fine and is untested-but-expected
to work given the bug is macOS-specific:
```bash
nextflow run main.nf
```

## Output

`results/starsolo/{sample}.Solo.out/Gene/raw/` — the standard sparse count
matrix triplet (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) per sample,
plus `{sample}.Log.final.out` with real STAR mapping statistics.

## Status

- [x] STARsolo quantification runs end-to-end (Docker; real macOS-specific
      STAR bug found, diagnosed, and worked around — documented above)
- [ ] Downstream Scanpy analysis (QC, filtering, clustering, marker genes)
- [ ] CI
