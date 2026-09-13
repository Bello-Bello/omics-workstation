# Walkthrough: 16S amplicon data, end to end

A step-by-step read-along for this pipeline. Each step says what runs, what it
produces, what can go wrong, and — where it matters — why the obvious
alternative is worse. You can follow it with the pipeline running or just read
it as notes.

Run the whole thing first if you want outputs to look at:

```bash
nextflow run main.nf
```

---

## Step 0 — What the data actually is

75 soil samples from the Atacama Desert, split between sites with vegetation
and sites without. Each sample was amplified at the **V4 region of the 16S
rRNA gene** using the 515F/806R primer pair, then sequenced.

Three things follow from "amplicon" that shape everything downstream:

1. **You are sequencing one gene, not genomes.** 16S is present in all
   bacteria and archaea and has conserved regions (for primers to bind) next
   to variable regions (to tell organisms apart). You get community
   composition, not gene content — you learn *who is there*, not *what they
   can do*. Shotgun metagenomics answers the second question and costs
   roughly an order of magnitude more.

2. **The reads are pooled.** All 75 samples were sequenced in one run, each
   tagged with a short barcode. Demultiplexing is the step that separates
   them again.

3. **The counts are not abundances.** This is the one that matters most, and
   Step 7 is entirely about it.

The three input files:

| File | What it holds |
|---|---|
| `forward.fastq.gz` | read 1 of every pair, from every sample, pooled |
| `reverse.fastq.gz` | read 2 of every pair, same order |
| `barcodes.fastq.gz` | the barcode for each pair, same order |

That third parallel file is the **Earth Microbiome Project (EMP) layout**.
Other protocols put the barcode inline in the read header, which needs a
different import path in QIIME2 — this is the first place a pipeline can go
wrong by assuming the wrong format.

---

## Step 1 — Fetch

`modules/fetch_data.nf` → `FETCH_DATA`, `FETCH_CLASSIFIER`

Downloads the three FASTQs, the metadata sheet, and the pretrained classifier.

**Worth noticing:** the process uses Nextflow's `storeDir`, not `publishDir`.
`storeDir` makes the download a *persistent cache* — rerun the pipeline and
Nextflow skips the process entirely because the output already exists. That
matters here because you will rerun this pipeline many times while tuning
DADA2, and re-downloading 300 MB each time is the kind of friction that makes
people stop iterating.

**The metadata sheet** has a quirk worth knowing. Row 2 is not data:

```
sample-id   barcode-sequence   elevation   ...
#q2:types   categorical        numeric     ...
BAQ1370.1.2 GCCCAAGTTCAC       1370        ...
```

That `#q2:types` row declares column types to QIIME2. Read the file with
pandas naively and every numeric column becomes text, which surfaces much
later as a baffling comparison error. `compositional/io.py::read_metadata`
drops it at the door.

---

## Step 2 — Import

`modules/import_and_demux.nf` → `IMPORT_EMP`

```bash
qiime tools import \
    --type EMPPairedEndSequences \
    --input-path emp-paired-end-sequences \
    --output-path emp-paired-end-sequences.qza
```

This does almost nothing computationally — it wraps the FASTQs in a `.qza`.

**What a `.qza` is:** a zip file containing your data plus a UUID, the type,
and the full provenance of every command that produced it. Rename one to
`.zip` and look inside. The provenance is the genuinely useful part: months
later you can ask an artifact what parameters made it, and it answers.

**What it costs:** you are now inside QIIME2's type system, and everything has
to enter and leave through `qiime tools import` / `export`. That is why this
pipeline has an explicit `EXPORT_TABLE` step — the boundary is real, so it is
drawn where you can see it in the DAG.

---

## Step 3 — Demultiplex

`modules/import_and_demux.nf` → `DEMUX`

```bash
qiime demux emp-paired \
    --m-barcodes-file sample_metadata.tsv \
    --m-barcodes-column barcode-sequence \
    --p-rev-comp-mapping-barcodes \
    --i-seqs emp-paired-end-sequences.qza \
    --o-per-sample-sequences demux.qza \
    --o-error-correction-details demux-details.qza
```

Reads get matched to samples by barcode, with error correction (Golay codes
tolerate a small number of miscalls).

**`--p-rev-comp-mapping-barcodes` is the trap.** For this dataset the barcodes
as sequenced are the reverse complement of the ones written in the metadata.
Omit the flag and QIIME2 does **not** error — it demultiplexes almost nothing
and hands you a nearly empty table. That looks like a data quality disaster
and sends you investigating the wrong thing entirely.

The general lesson, which recurs: **the dangerous failures are the ones that
return valid-looking output.** When a step produces far less than expected,
suspect an orientation or format flag before suspecting the data.

Check `results/qiime/demux.qzv` at [view.qiime2.org](https://view.qiime2.org) —
per-sample read counts and the quality profiles you need for the next step.

---

## Step 4 — Denoise with DADA2

`modules/denoise.nf` → `DENOISE_DADA2`

The heart of the pipeline, and the slowest step.

```bash
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --p-trim-left-f 13 --p-trim-left-r 13 \
    --p-trunc-len-f 150 --p-trunc-len-r 150 \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza
```

**What DADA2 does:** learns the run's own error model — how often this
sequencer turns an A into a G at this quality score — then uses it to decide
whether a rare sequence is a real biological variant or a misread of an
abundant neighbour. Then it merges read pairs and removes chimeras (artifacts
where two different templates fused during PCR).

**ASVs vs OTUs, and why it matters.** The older approach clustered reads at
97% identity into Operational Taxonomic Units, on the reasoning that you
cannot distinguish a real variant from an error, so bin them together. DADA2's
error model *can* distinguish them, so it returns Amplicon Sequence Variants:
exact sequences. Two practical consequences — you can resolve organisms
differing by a single base, and ASVs are comparable across studies without
re-clustering, which OTUs are not.

**The four parameters:**

- `trim_left` cuts the primer off the 5′ end.
- `trunc_len` cuts the 3′ end where quality degrades.

`trunc_len` is where runs get ruined. Truncate too little and low-quality
tails poison the error model; truncate too much and **the forward and reverse
reads no longer overlap enough to merge** — and DADA2 will not tell you that
is the problem. It just returns very few features.

**So read `denoising-stats.tsv` every time.** It is a per-sample table of
attrition:

```
input → filtered → denoised → merged → non-chimeric
```

A big drop at `merged` means your truncation lengths left insufficient
overlap. A big drop at `non-chimeric` suggests over-amplification. This one
file answers "why is my table so empty?" faster than anything else, which is
why the pipeline exports it as a TSV unconditionally instead of leaving it
inside a `.qza`.

---

## Step 5 — Classify

`modules/classify.nf` → `CLASSIFY_TAXONOMY`

```bash
qiime feature-classifier classify-sklearn \
    --i-classifier classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza
```

A naive Bayes classifier over k-mer profiles assigns each ASV a taxonomic
string like `k__Bacteria; p__Actinobacteria; c__Thermoleophilia; ...`.

**Two things must match, and neither failure is loud:**

1. **The amplified region.** The classifier is trained on 515F/806R sequences
   because that is what this dataset amplified. A classifier trained on the
   full-length gene, or a different variable region, returns confident
   nonsense rather than an error.

2. **The scikit-learn version.** The classifier is a pickled sklearn model.
   The container here has sklearn 1.4.2, so the pipeline pins the
   `sklearn-1.4.2` classifier URL. Mismatch and you get an opaque unpickling
   error — annoying, but at least it fails loudly, unlike (1).

**Expect a lot of `Unassigned` in soil.** Much of the Atacama community has no
cultured representative in the reference database. That is a real property of
the ecosystem, not a pipeline defect — which is why Step 7's genus collapse
keeps unassigned features in an explicit bin instead of dropping them.

---

## Step 6 — Export

`modules/export_table.nf` → `EXPORT_TABLE`

```bash
qiime tools export --input-path table.qza --output-path table_export
biom convert -i table_export/feature-table.biom -o feature-table.tsv --to-tsv
```

Two unwrappings: `.qza` → BIOM, then BIOM → TSV. BIOM is a sparse-matrix
format for count tables; `biom convert` gives you something pandas can read.

**Orientation matters here.** QIIME2 writes **features as rows, samples as
columns**. Everything downstream in this pipeline wants the transpose. Get it
backwards and the CLR still computes — it just silently means nothing, because
you have normalised across the wrong axis. `io.read_feature_table` transposes
once, at the boundary, so no other function has to think about it.

---

## Step 7 — The compositional step

`compositional/transforms.py`

This is the conceptual core. Everything before it was mechanics.

### The problem

Your sample has 12,000 reads. Not because the soil had 12,000 cells — because
that is how the flow cell divided its capacity that day. Sequence the same
soil again and you might get 30,000.

So a count column tells you about **ratios between taxa within that sample**
and nothing else. The data live on a *simplex* (everything sums to a constant),
not in ordinary Euclidean space. Two things break as a result:

**Spurious negative correlation.** Parts must sum to a constant, so one taxon
rising forces every other proportion down. Correlate raw proportions and you
manufacture negative associations with no biological basis. Every taxon looks
like it competes with every other taxon.

**Differences are not distances.** "1% → 2%" and "50% → 51%" are both +1
percentage point. The first is a doubling; the second is a 2% change.
Euclidean distance on proportions treats them identically. It should not.

### The fix

Work with log-ratios, which do not care about the total.

**Closure** — rescale each sample to sum to 1. This does not destroy
information; it makes explicit what sequencing already did implicitly.

**Zero replacement** — `log(0)` is undefined, so every zero needs a decision.
In amplicon data a zero is almost always a *rounded* zero: below detection at
this depth, not truly absent. Multiplicative replacement substitutes a small
δ and shrinks the observed parts proportionally so the row still sums to 1.

> The "multiplicative" part is not pedantry. Adding a flat pseudocount
> (`counts + 1`) distorts the ratios between the non-zero parts, and distorts
> them hardest in the shallowest samples — exactly backwards, since those are
> the samples you already trust least.

**CLR** — divide each part by the **geometric mean of its own sample**, then
take the log:

```
clr(x)_i = log( x_i / g(x) )
```

Scale every count in a sample by *k*: the geometric mean scales by *k* too,
the ratio is unchanged, the CLR is unchanged. **That invariance is the entire
point** — sequencing depth has been cancelled out by construction rather than
corrected for by a normalisation factor you have to justify.

**How to read a CLR value:** "how much more of this taxon than the typical
taxon in this sample", on a log scale. Positive = enriched relative to that
sample's own average. Zero = average. **Not** "absent" — that is why the
heatmap uses a diverging colour scale centered at zero.

A structural consequence: CLR rows sum to zero by construction, so the data
are rank-deficient by one. PCA is fine with that. Anything that inverts a
covariance matrix is not.

### The bug this pipeline actually hit

δ defaults to `0.65 / depth`. Total imputed mass is `n_zeros × δ`. On the 1%
subsample, one sample had **4 reads across 12 features** — 10 zeros × 0.1625 =
**1.625** of mass, for a composition that must sum to 1. Re-closure then
assigned the observed parts a *negative* proportion, and the error surfaced
two functions later inside `clr()`, pointing nowhere near the cause.

The fix caps total imputed mass and warns when a row hits the cap — because a
row that hits it is telling you something true: more than half that sample is
imputation now, and the answer is to raise `min_depth`, not to trust the
numbers. See `README.md`, "A real bug this found".

---

## Step 8 — Ordination and clustering

`compositional/ordination.py`

Once CLR has moved the data into real space, ordinary Euclidean methods become
*correct*, not merely convenient.

**Aitchison distance** = Euclidean distance between CLR vectors. Its advantage
over Bray-Curtis is **subcompositional coherence**: drop half the taxa and the
distances among the samples you kept are computed from the remaining ratios
only. They do not lurch around because the total changed. Bray-Curtis on
proportions offers no such guarantee — filtering rare taxa can move your
ordination, which is unsettling once you notice it.

**PCA on CLR values** (an *Aitchison biplot*) is equivalent to PCoA on the
Aitchison distance matrix, but keeps the loadings — so you can ask which taxa
drive each axis, not just where the samples land. `pca_result.top_loadings(0)`
answers that.

> Implementation note: computed by SVD, not by eigendecomposition of the
> covariance matrix. CLR rows sum to zero, so that matrix is rank-deficient by
> one and its smallest eigenvalue should be 0 — eigensolvers sometimes return
> a small *negative* number there, and the square root becomes NaN. SVD on the
> centered matrix never forms the covariance matrix, so the problem cannot
> arise.

**Hierarchical clustering** builds a dendrogram from the same distances. The
figure orders heatmap rows by the dendrogram, not by group, with the group
shown only as a colour strip — so you see what the clustering *found* rather
than what the design *asserts*. The adjusted Rand index in `summary.json` puts
one number on the agreement: 0 is chance, 1 is perfect.

---

## Step 9 — PERMANOVA

`compositional/stats.py::permanova`

*"Does vegetation explain community structure at all?"* — one test on the whole
community, before any per-taxon testing.

```
SS_total  = (1/N) Σ_{i<j} d_ij²
SS_within = Σ_g (1/n_g) Σ_{i<j ∈ g} d_ij²
pseudo-F  = (SS_between / (a-1)) / (SS_within / (N-a))
```

It partitions variance into between- and within-group parts using **only the
pairwise distances** — samples never need coordinates, which is what lets it
work on any metric.

**Why permutation.** The statistic has the shape of an F, but not its null
distribution: distances among the same N points are not independent. So you
shuffle the group labels several hundred times and ask how often chance
produces a pseudo-F this large.

**Why the p-value denominator is `permutations + 1`.** The observed value
counts as one of the draws. That keeps the test valid, and it means the
smallest achievable p-value with 999 permutations is 0.001, never 0. If you
ever see `p = 0` reported from a permutation test, the implementation is
wrong.

**The caveat to carry into any report:** PERMANOVA responds to differences in
**dispersion** as well as in location. A significant result means the groups
differ — not necessarily that they are *centered* in different places. One
group being more variable is enough. The generated `report.md` states this
rather than leaving the reader to know it.

---

## Step 10 — Differential abundance

`compositional/stats.py::differential_abundance`

*"Which taxa differ?"* — Welch's t-test per feature on CLR values, then
Benjamini-Hochberg across all of them.

**Why not a negative binomial GLM,** as in the RNA-seq work in
`stats/rnaseq-stats-notebook/`? Because that framing models counts as draws
whose mean scales with library size, and the simplex constraint here couples
taxa within a sample in a way that violates it. Once CLR has done its job, a
location test on the transformed values is both simpler and better matched to
what the data are.

**Benjamini-Hochberg**, briefly: testing 400 taxa at α = 0.05 buys you ~20
false positives by construction. BH sorts the p-values, scales each by
`n / rank`, then takes a cumulative minimum from the largest downwards. That
last pass is the step-up part and the easy thing to omit — without it the
adjusted values are not monotone, and you can reject a larger p-value while
accepting a smaller one.

**The Monte-Carlo pass** is the part worth understanding.

A point-estimate CLR pretends the observed counts *are* the composition. For a
taxon seen 3 times in 8,000 reads, that is a confident claim from very little
evidence. So the pipeline draws the composition from its posterior instead —
`Dirichlet(counts + 0.5)`, the Jeffreys prior — reruns the entire test on each
draw, and averages the adjusted p-values. This is the idea behind ALDEx2.

**Expect the Monte-Carlo column to be less significant than the naive one, and
expect the gap to be widest for the rarest taxa.** That gap is the honest
answer to "how much of this hit is signal, and how much is a small number that
landed well". Both columns are reported — reporting only the better-looking one
would be the easiest possible way to mislead with this pipeline.

**One more thing to expect:** taxa you did not plant can still come up
significant. That is not necessarily a false positive in the usual sense —
under the simplex constraint, genuinely enriching some taxa depresses the CLR
of others. The constraint is real, so the shift is real. It just is not
independent evidence.

---

## Step 11 — Read the output

```
results/compositional/
  report.md                          written interpretation, with caveats
  summary.json                       every number in the report
  asv_differential_abundance.csv     per-ASV results
  genus_differential_abundance.csv   per-genus results
  sample_depth.png                   QC: how uneven was depth?
  *_pca.png                          ordination
  *_scree.png                        variance per component
  *_heatmap_dendrogram.png           clustering + top variable taxa
  *_volcano.png                      effect size vs significance
  *_top_hits_boxplot.png             the individual hits, with points shown
```

**Read them in this order:**

1. `denoising-stats.tsv` — did enough reads survive? If not, nothing below matters.
2. `sample_depth.png` — how uneven was depth? This is the confounder CLR neutralises; look at it before being told it did not matter.
3. `*_pca.png` — does anything separate at all?
4. `permanova.json` — is that separation more than chance?
5. `*_volcano.png` and the CSVs — which taxa, and how confidently?
6. `report.md` — the written version, including what it cannot claim.

---

## Verify the statistics

```bash
python bin/validate_compositional.py
```

Three tiers: mathematical invariants, agreement with scikit-bio and
statsmodels, and recovery of planted signal against a null dataset.

The third tier is the one that matters. Agreeing with scikit-bio proves both
implementations do the same thing; recovering planted truth and staying quiet
on null data proves the thing is the right thing.

---

## What this pipeline cannot tell you

Worth being able to say out loud in an interview:

- **Nothing about absolute abundance.** "Enriched in vegetated soil" means
  relative to the rest of that community. Cells per gram needs spike-ins or
  qPCR.
- **Nothing about function.** 16S is taxonomy. What the community *does*
  needs shotgun metagenomics — assembly, gene calling, pangenomes, BLAST/HMM
  annotation.
- **Nothing causal.** Vegetated and barren sites differ in pH, moisture and
  temperature too. Any of those could drive the community difference; the
  metadata has those columns, and testing them is the obvious next analysis.
- **Nothing about strains.** ASVs resolve exact V4 sequences, not genomes.
  Two organisms with identical V4 are one feature here.
