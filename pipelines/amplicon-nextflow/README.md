# amplicon-nextflow

[![amplicon-nextflow CI](https://github.com/Bello-Bello/omics-workstation/actions/workflows/amplicon-nextflow-ci.yml/badge.svg)](https://github.com/Bello-Bello/omics-workstation/actions/workflows/amplicon-nextflow-ci.yml)

16S amplicon pipeline on **real public soil data**: raw multiplexed reads →
ASVs → taxonomy → compositional statistics.

The dataset is the Atacama Desert soil survey (Neilson et al. 2017, *mSystems*),
hosted publicly by QIIME2 — 75 samples of some of the most oligotrophic soil on
Earth, split between vegetated and barren sites along two transects. The
biological question the pipeline answers is the obvious one: **does the
presence of vegetation change the soil bacterial community, and which taxa
carry the difference?**

## The two halves

```
FETCH ─┬─> IMPORT_EMP ─> DEMUX ─> DENOISE_DADA2 ─┬─> CLASSIFY ─┐
       │                                         │             ├─> EXPORT ─> COMPOSITIONAL_ANALYSIS
       └─> FETCH_CLASSIFIER ─────────────────────┘─────────────┘
       └──────────── QIIME2 container ──────────────┘   └── conda env ──┘
```

**Left of `EXPORT_TABLE`** is QIIME2, in a pinned container, doing what QIIME2
is genuinely best at: EMP-protocol demultiplexing, DADA2 denoising, and
naive-Bayes taxonomic classification. This half is orchestration and
troubleshooting — install the framework, get the flags right, notice when a
step silently produces nothing.

**Right of it** is ordinary Python in a small conda env, with the compositional
methods written out in [`compositional/`](compositional/) rather than called
from a black box.

## Why the statistics are not `qiime composition ancombc`

QIIME2 ships a perfectly good differential abundance plugin. Calling it would
be one line and would teach nothing.

The reason this repository implements CLR, PERMANOVA and Benjamini-Hochberg
directly is the same reason [`stats/rnaseq-stats-notebook/`](../../stats/rnaseq-stats-notebook/)
reimplements DESeq2's negative binomial GLM: the interesting skill is not
knowing which command to run, it is knowing what the command assumes and when
those assumptions stop holding. A CLR that returns numbers is easy. Knowing
that the zero-replacement step can quietly go negative on a shallow sample —
and that the resulting error surfaces two functions later, pointing nowhere
near the cause — requires having written the zero-replacement step.

[`bin/validate_compositional.py`](bin/validate_compositional.py) then checks
this implementation against scikit-bio and statsmodels, so "written from
scratch" does not become "written wrong". It agrees with scikit-bio's
PERMANOVA to machine precision and with statsmodels' BH to 1e-16.

## Why both ASV and genus level

The pipeline runs the whole analysis twice, on the ASV table and on the
rank-collapsed table, and reports both.

ASVs are exact sequences, so two ASVs one base apart are separate features
even when they are the same organism for every practical purpose — which
splits real signal across columns and costs statistical power. Collapsing to
genus pools them, at the cost of resolution. Neither is the right answer in
general, and a result that appears at one level but not the other is
information about how robust it is, not a contradiction to be hidden.

The collapse keeps unclassified features in an explicit `Unassigned_genus`
bin rather than dropping them. In soil this bin is large — much of the Atacama
community has no cultured representative in Greengenes — and discarding it
would silently change every ratio in the table.

## What compositional means here, in one paragraph

A 16S table does not measure abundance. Sequencing returns a fixed number of
reads per sample, so a column of counts carries information about *ratios
between taxa within that sample*, nothing more. Two consequences bite
immediately: the constant-sum constraint manufactures negative correlations
between taxa that have no biological relationship, and Euclidean distance on
proportions treats "1% → 2%" (a doubling) as the same change as "50% → 51%"
(a 2% shift). The centered log-ratio transform divides each part by its own
sample's geometric mean before taking logs, which cancels the arbitrary total
and lands the data back in real space where PCA, clustering and t-tests are
valid again. Everything in [`compositional/transforms.py`](compositional/transforms.py)
follows from that.

## Structure

```
main.nf                          # the workflow
nextflow.config                  # params, container/conda wiring, test profile
modules/
  fetch_data.nf                    # data + pretrained classifier, storeDir-cached
  import_and_demux.nf              # EMPPairedEndSequences -> per-sample reads
  denoise.nf                       # DADA2 -> ASV table + rep seqs + stats
  classify.nf                      # naive-Bayes taxonomy
  export_table.nf                  # leave the .qza world: BIOM/TSV
  compositional_analysis.nf        # the statistics
compositional/                   # reusable package, not dataset-specific
  transforms.py                    # closure, zero replacement, CLR, Aitchison distance
  stats.py                         # PERMANOVA, differential abundance, BH
  ordination.py                    # PCA on CLR, hierarchical clustering, ARI
  plots.py                         # scatter / bar / box / heatmap / dendrogram
  io.py                            # QIIME2 exports -> aligned frames
bin/
  run_compositional_analysis.py    # the driver Nextflow calls
  validate_compositional.py        # cross-checks vs scikit-bio & statsmodels
envs/compositional.yaml          # the small analysis env (not QIIME2)
```

## Quick start

New to this? Start with **[TUTORIAL.md](TUTORIAL.md)** — a hands-on walk
through setup, running it, and reading the output, with what you should see
at each step. [WALKTHROUGH.md](WALKTHROUGH.md) explains the concepts behind
each stage once you have it running.

Needs Docker (for the QIIME2 container), conda, and Nextflow.

```bash
nextflow run main.nf                  # full run, ~1.35M reads
nextflow run main.nf -profile test    # same data, fewer permutations
python bin/validate_compositional.py  # statistics self-check, no data needed
```

Results land in `results/`:

| Path | What it is |
|---|---|
| `qiime/` | the `.qza`/`.qzv` artifacts, viewable at [view.qiime2.org](https://view.qiime2.org) |
| `qiime/denoising-stats.tsv` | read attrition through filter → merge → chimera removal |
| `exported/` | the plain-TSV feature table and taxonomy |
| `compositional/report.md` | the written interpretation |
| `compositional/summary.json` | every number quoted in the report |
| `compositional/*_differential_abundance.csv` | per-taxon results, ASV and genus |
| `compositional/*.png` | ordination, dendrogram+heatmap, volcano, box plots, depth |

## Parameters worth knowing about

| Param | Default | Why it matters |
|---|---|---|
| `trunc_len_f` / `trunc_len_r` | 150 / 150 | Truncate where quality falls off — but the two must still leave enough overlap to merge pairs. Too aggressive and the merge step silently returns almost nothing. |
| `trim_left_f` / `trim_left_r` | 13 / 13 | Removes the primer. |
| `min_depth` | 1000 | Samples below this are dropped. Shallow samples are mostly zeros, so their CLR values are dominated by imputation rather than measurement. |
| `min_prevalence` / `min_total_count` | 0.15 / 25 | Rare-feature filter. Keeping thousands of near-empty columns lets the zero replacement, not the data, drive the geometric mean. |
| `group_column` / `group_a` / `group_b` | `vegetation` / `yes` / `no` | Any two-level metadata column works — `transect-name` is the other obvious one in this dataset. |
| `monte_carlo` | 128 | Dirichlet draws per feature for the depth-aware p-value. Set to 0 cost by using `-profile test`. |

## Things this pipeline gets wrong on purpose, and says so

- **PERMANOVA responds to dispersion, not only location.** A significant
  result means the groups differ; it does not establish that they differ *in
  the way you assumed*. The generated report states this rather than leaving
  it to the reader.
- **Every result describes the filtered subcomposition.** Aitchison geometry
  is subcompositionally coherent, so this is permitted — but it is not free,
  and the thresholds are recorded in `summary.json` alongside the results.
- **CLR is relative, always.** "Enriched in vegetated soil" means enriched
  relative to the rest of that community. No 16S design of this kind can speak
  to absolute cells per gram; that needs spike-ins or qPCR.
- **A single classifier, one region.** Greengenes 13_8 trained on 515F/806R
  matches this dataset's primers. Swapping the reference without matching the
  region produces confident nonsense rather than an error.

## A real bug this found

The conventional multiplicative zero replacement assigns each zero
`0.65 / depth` of the composition. That is derived for samples with a sensible
number of reads, and it breaks silently when they do not have one: the imputed
mass is `n_zeros × delta`, which on the 1%-subsampled table reached **1.625**
for a composition that must sum to 1. Re-closure then hands the observed parts
a negative proportion, and the failure surfaces in `clr()` — two functions
downstream of the actual cause.

The fix caps total imputed mass per sample and warns when a row hits the cap,
because a row that hits it is telling you something true: more than half of
that sample is now imputation, and the right response is to raise `min_depth`,
not to trust the numbers. The regression test is in
`bin/validate_compositional.py` ("Shallow sample does not produce negative
parts").

That is also why the `test` profile reduces permutations rather than switching
to the 1% subsample: cutting the statistical work is honest, cutting the data
until it stops being data is not.

## Data

- **Reads and metadata** — Atacama soil microbiome, via QIIME2's public tutorial
  mirror. Not committed; fetched at run time and cached in `.data_cache/`.
  Source: Neilson, J.W. *et al.* (2017) "Significant impacts of increasing
  aridity on the arid soil microbiome", *mSystems* 2:e00195-16.
- **Classifier** — Greengenes 13_8, 99% OTUs, 515F/806R region, pretrained for
  scikit-learn 1.4.2 to match the pinned container.
