# Tutorial: run this pipeline yourself

A hands-on guide. Every step has a command to run and what you should see
afterwards, so you can tell at each point whether it worked.

This is the *doing* companion to [WALKTHROUGH.md](WALKTHROUGH.md), which
explains the concepts. Do this one first; read that one when you want to know
why a step exists.

**Time:** about 30 minutes, most of it waiting for DADA2.
**You need:** a terminal. No prior QIIME2 experience.

---

## Part 0 — Set up (once)

### 0.1 Check what you already have

```bash
for t in docker conda nextflow java; do
  printf "%-10s %s\n" "$t:" "$(command -v $t || echo MISSING)"
done
```

You need all four. Install whatever is missing:

| Tool | Install | What it's for |
|---|---|---|
| Docker | [Docker Desktop](https://www.docker.com/products/docker-desktop/) | runs QIIME2 without installing it |
| conda | [Miniforge](https://github.com/conda-forge/miniforge) | the Python analysis environment |
| Java | `brew install openjdk` | Nextflow runs on the JVM |
| Nextflow | `curl -s https://get.nextflow.io \| bash` | the workflow engine |

If `java` is missing after `brew install openjdk`, Homebrew installed it
"keg-only", meaning it deliberately did not put it on your PATH. Add it:

```bash
echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
java -version
```

### 0.2 Start Docker

Open Docker Desktop and wait for the whale icon in your menu bar to stop
animating. Then confirm the daemon is actually up:

```bash
docker info > /dev/null && echo "Docker is running"
```

> **The single most common failure in this tutorial** is Docker not being
> started. Nextflow's error for it is long and mentions permissions, not
> Docker — so check here first, always.

### 0.3 Get the code

```bash
git clone https://github.com/Bello-Bello/omics-workstation.git
cd omics-workstation/pipelines/amplicon-nextflow
```

---

## Part 1 — Check the statistics before running anything

The analysis code can be tested on its own, with no sequencing data at all.
Start here: it takes a minute and tells you the maths is sound before you
spend half an hour on the pipeline.

```bash
conda env create -f envs/compositional.yaml -p ./env
./env/bin/python bin/validate_compositional.py
```

**You should see** 19 checks, all `[PASS]`, ending with `All checks passed.`
Among them:

```
  [PASS] CLR is invariant to sequencing depth — a sample scaled 1e6x gives identical CLR values
  [PASS] Benjamini-Hochberg matches statsmodels — max |diff| = 2.22e-16
  [PASS] PERMANOVA pseudo-F matches scikit-bio (null) — 1.105920313 vs 1.105920313
  [PASS] All planted taxa recovered — 12/12 recovered, 14 total hits
```

**Read those four lines carefully — they are the whole argument of this
pipeline in miniature:**

1. The CLR cancels sequencing depth. Multiply a sample by a million, get the
   same numbers back.
2. The hand-written statistics agree with the standard libraries.
3. On data with a *known* answer, the pipeline finds it — and finds nothing
   on data with no answer to find.

---

## Part 2 — Run the pipeline

```bash
nextflow run main.nf
```

Then go make coffee. It takes 10–15 minutes, and DADA2 is most of it.

**What's happening, in order:**

| Stage | Time | What it does |
|---|---|---|
| `FETCH_DATA` | 1–3 min | downloads ~300 MB of soil reads |
| `FETCH_CLASSIFIER` | ~30 s | downloads the taxonomy classifier |
| `IMPORT_EMP` | seconds | wraps the reads as a QIIME2 artifact |
| `DEMUX` | 1–2 min | splits pooled reads back into 75 samples |
| `DENOISE_DADA2` | **8–10 min** | error-corrects reads into ASVs |
| `CLASSIFY_TAXONOMY` | 1–2 min | names each ASV |
| `EXPORT_TABLE` | seconds | converts out of QIIME2 format |
| `COMPOSITIONAL_ANALYSIS` | ~1 min | the statistics |

**Success looks like:**

```
[SUCCESS] completed=8 failed=0 cached=0
```

> **If it fails**, don't rerun from scratch. Fix the problem, then add
> `-resume` — Nextflow reuses every step that already succeeded:
> ```bash
> nextflow run main.nf -resume
> ```
> See Part 6 for the common failures.

---

## Part 3 — Read the results, in the right order

The order matters. Each step below is only worth looking at if the one before
it looked sane.

### 3.1 Did enough reads survive?

```bash
head -5 results/qiime/denoising-stats.tsv | column -t
```

Columns go left to right: `input → filtered → denoised → merged → non-chimeric`.
Reads drop at each stage. That is normal — what matters is *where*.

```bash
./env/bin/python -c "
import pandas as pd
d = pd.read_csv('results/qiime/denoising-stats.tsv', sep='\t', skiprows=[1])
print('samples:', len(d))
print('median reads in: ', int(d['input'].median()))
print('median reads out:', int(d['non-chimeric'].median()))
print('median retained: {:.0f}%'.format(100*(d['non-chimeric']/d['input']).median()))
"
```

**Diagnosis:** a big drop at **merged** means your forward and reverse reads
no longer overlap — `trunc_len` was too aggressive. A big drop at
**non-chimeric** suggests over-amplification during PCR. This one file
answers "why is my table empty?" faster than anything else.

### 3.2 How uneven was sequencing depth?

```bash
open results/compositional/sample_depth.png
```

You will see a steep curve — some samples have hundreds of times more reads
than others. **This is the problem the whole compositional approach exists to
solve.** Look at it before you are told it doesn't matter.

### 3.3 Does anything separate?

```bash
open results/compositional/genus_pca.png
```

Each dot is a soil sample. Orange = vegetated, blue = barren. You should see
them separating left-to-right along PC1.

That axis is not "abundance" — it's built from *log-ratios*, so a sample's
position is unaffected by how many reads it happened to get. That is the
payoff of the CLR.

### 3.4 Is the separation more than chance?

```bash
cat results/compositional/permanova.json
```

Expect roughly:

```json
"genus": { "pseudo_f": 8.93, "p_value": 0.001, "r_squared": 0.152 }
```

Reading it: vegetation explains **15%** of the variation between these
communities, and a p of 0.001 means that in 999 random reshuffles of the
group labels, *none* produced a split this strong.

> **A detail worth internalising:** 0.001 is the smallest p this test can
> produce with 999 permutations — it's `1/(999+1)`. A permutation test can
> never return p = 0. If you ever see one reported, the implementation is
> broken.

### 3.5 Which organisms?

```bash
open results/compositional/genus_volcano.png
column -s, -t < results/compositional/genus_differential_abundance.csv | head -8
```

Then the biology:

| Taxon | Where | Why it makes sense |
|---|---|---|
| *Bradyrhizobium* | vegetated | nitrogen-fixing symbiont that lives with plant roots |
| *Candidatus* Nitrososphaera | vegetated | ammonia-oxidising archaeon — nitrogen cycling |
| DA101 | vegetated | classic rhizosphere organism |
| *Rubrobacter* | barren | famously resistant to drying and radiation |

**This is the moment the pipeline stops being plumbing.** Nothing told it
about plants. It read sequences and recovered the nitrogen-cycling community
that lives around roots, plus the extremophile that survives where nothing
grows.

### 3.6 The written version

```bash
cat results/compositional/report.md
```

Generated by the pipeline, caveats included. Note that it reports two numbers
for significance — 10 taxa by the depth-aware test, 23 by the naive one. The
honest number is the smaller one.

---

## Part 4 — Poke at it

Now the parts you can only understand by touching.

### 4.1 See the CLR cancel sequencing depth

```bash
./env/bin/python -c "
import numpy as np, sys; sys.path.insert(0,'.')
from compositional.transforms import clr

sample = np.array([[100., 50., 30., 20.]])
deep   = sample * 1000

print('shallow counts:', sample[0])
print('deep counts:   ', deep[0])
print()
print('shallow CLR:', np.round(clr(sample)[0], 6))
print('deep CLR:   ', np.round(clr(deep)[0], 6))
print()
print('identical:', np.allclose(clr(sample), clr(deep)))
"
```

One sample sequenced a thousand times deeper. **Identical CLR values.** No
normalisation factor, no correction step — the transform is built so depth
cancels out.

### 4.2 See why raw proportions lie

```bash
./env/bin/python -c "
import numpy as np
before = np.array([100., 100., 100., 100.])
after  = np.array([400., 100., 100., 100.])   # only taxon 1 actually grew

print('taxon 2 raw counts: ', before[1], '->', after[1], ' (unchanged)')
print('taxon 2 proportion: {:.0%} -> {:.0%}'.format(
      before[1]/before.sum(), after[1]/after.sum()))
"
```

Taxon 2 did not change. Its *proportion* fell from 25% to 14%, because taxon 1
grew. Analyse proportions directly and you will report a decline that never
happened. **This is the trap the CLR avoids**, and it is the single most
important idea in microbiome statistics.

### 4.3 Look inside a `.qza` file

```bash
cp results/qiime/table.qza /tmp/table.zip
unzip -o -q /tmp/table.zip -d /tmp/qza_peek
find /tmp/qza_peek -name "*.yaml" | head -3
cat $(find /tmp/qza_peek -name "metadata.yaml" | head -1)
```

A `.qza` is a zip file. Inside is your data plus **complete provenance** —
every command and parameter that produced it. Months later you can ask a file
how it was made, and it answers.

### 4.4 Run one step by hand

```bash
docker run --rm -v "$PWD/results:/data" quay.io/qiime2/amplicon:2024.10 \
  qiime feature-table summarize \
    --i-table /data/qiime/table.qza \
    --o-visualization /data/qiime/table-summary.qzv
```

Then upload `results/qiime/table-summary.qzv` to
[view.qiime2.org](https://view.qiime2.org) — an interactive summary, in your
browser, nothing installed.

---

## Part 5 — Change something

This is where it becomes yours.

### 5.1 Ask a different question

The metadata has more than vegetation in it. Compare the two transects:

```bash
nextflow run main.nf -resume \
  --group_column transect-name --group_a Baquedano --group_b Yungay
```

Only the analysis step reruns — `-resume` reuses the denoising. About a
minute. Then compare `results/compositional/report.md` with what you saw
before. Is geography a stronger signal than vegetation?

### 5.2 Change the filtering and watch it matter

```bash
nextflow run main.nf -resume --min_prevalence 0.5 --min_total_count 100
```

Far stricter: a taxon must appear in half the samples. Fewer features
survive, so the geometric mean is computed over a different set — **every CLR
value changes.** Results describe whatever subcomposition you kept, which is
why the pipeline records the thresholds in `summary.json` next to the results.

### 5.3 Break it deliberately

```bash
nextflow run main.nf --min_depth 50000
```

No sample has 50,000 reads, so everything is filtered out. Read the error.
Good pipelines fail with a message that tells you what to change — decide
whether this one does.

---

## Part 6 — When it goes wrong

| Symptom | Cause | Fix |
|---|---|---|
| Long error mentioning permissions or sockets | Docker isn't running | Start Docker Desktop, wait for the whale to settle |
| `Cannot find Java` / `NXF_JAVA_HOME` | Java not on PATH | `brew install openjdk`, then the PATH line in 0.1 |
| Hangs on `FETCH_DATA` | slow network, 300 MB | wait; it caches, so it only happens once |
| DADA2 gives almost no features | `trunc_len` too short to merge pairs | check `denoising-stats.tsv` for a drop at `merged` |
| `No module named compositional` | ran the script directly | run through Nextflow, or pass `--project-root .` |
| Conda solve fails | channel priority | `conda config --set channel_priority flexible` |
| Out of memory in DADA2 | Docker Desktop's RAM cap | Docker Desktop → Settings → Resources → raise to 8 GB |

**The general rule from this pipeline:** the dangerous failures are the ones
that return valid-looking output. An empty table is more often a wrong flag
than bad data. Check `denoising-stats.tsv` before you doubt the sequencing.

---

## Part 7 — Make it yours

Ordered by effort:

1. **Rerun with `--rank family`** instead of genus. Does the biology hold at a
   coarser level? Fewer, better-populated features — does that help or hurt?
2. **Test a numeric variable.** The metadata has pH, elevation and soil
   temperature. PERMANOVA as written takes categories, so you would need to
   bin them, or extend it to handle continuous predictors.
3. **Add a rarefaction curve** — plot features observed against reads
   sampled, to show whether sequencing went deep enough.
4. **Bring your own data.** Any 16S dataset with a metadata sheet works.
   Change the URLs in `nextflow.config` and the `group_column`.
5. **Reimplement it in Snakemake** and check the two agree, as the RNA-seq
   pipelines in this repo do. That cross-verification is a stronger
   reproducibility claim than either alone.

---

## What to take away

Three things are worth being able to say out loud:

1. **Sequencing counts are ratios, not abundances.** Everything else follows.
2. **The CLR makes depth irrelevant by construction** — not corrected for,
   cancelled.
3. **A test that cannot fail is not evidence.** That's why the validation
   script checks against a null dataset, and why the pipeline reports the
   more conservative p-value rather than the prettier one.

Then read [WALKTHROUGH.md](WALKTHROUGH.md) for the mechanism behind each
stage, and [README.md](README.md) for the design decisions and the real bug
this pipeline found.
