# Curated Resources

## Foundations
- Unix/shell: [The Unix Workbench](https://seankross.com/the-unix-workbench/) (free) or Software Carpentry's [Shell novice](https://swcarpentry.github.io/shell-novice/)
- Git: [Software Carpentry Git novice](https://swcarpentry.github.io/git-novice/)
- Python for biologists: [Python for Biologists](https://pythonforbiologists.com/) or Software Carpentry [Python novice-gapminder](https://swcarpentry.github.io/python-novice-gapminder/)
- R/tidyverse: [R for Data Science](https://r4ds.hadley.nz/) (free online)

## Reproducible workflows
- Snakemake: [official tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) — do this end to end before touching Nextflow
- Nextflow: [nf-core training](https://training.nextflow.io/) — the single best resource, also teaches nf-core conventions employers recognize
- Browse [nf-core pipelines](https://nf-co.re/pipelines) as reference implementations (rnaseq, sarek, scrnaseq) — read the code, don't just run it

## Bulk RNA-seq specific
- [RNA-seq analysis in R (Combine Australia workshop)](https://combine-australia.github.io/RNAseq-R/) — covers DESeq2 well
- Salmon docs: https://salmon.readthedocs.io/
- Test data: [nf-core/test-datasets (rnaseq branch)](https://github.com/nf-core/test-datasets/tree/rnaseq) — small, fast to download and run

## Statistical modeling
- [Modern Statistics for Modern Biology](https://www.huber.embl.de/msmb/) (free, EMBL) — directly relevant to omics stats
- DESeq2 vignette (explains the negative binomial GLM in plain terms): https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## scRNA-seq (if chosen for Pipeline #2)
- [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/index.html)
- [Single-cell best practices book](https://www.sc-best-practices.org/) (free, excellent)

## Variant calling (if chosen for Pipeline #2)
- [GATK best practices workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- nf-core/sarek source as reference

## Polish & portfolio presentation (Docker, CI, GitHub)
- Docker: [Docker's official "Get Started" guide](https://docs.docker.com/get-started/) — enough to write and run a Dockerfile for a pipeline environment
- Docker for research/bioinformatics specifically: [Carpentries Incubator: Docker Introduction](https://carpentries-incubator.github.io/docker-introduction/) — free, hands-on, aimed at exactly this use case
- GitHub Actions: [official quickstart](https://docs.github.com/en/actions/quickstart) — enough to write a simple "run the pipeline on test data on every push" workflow
- See how a real project does CI: browse [nf-core/rnaseq's GitHub Actions workflows](https://github.com/nf-core/rnaseq/tree/master/.github/workflows) — don't copy wholesale, but see the shape of a real CI setup for a Nextflow pipeline
- README writing: [makeareadme.com](https://www.makeareadme.com/) — a short, practical guide to what a README should actually contain
- What "good" looks like in practice: browse a few pinned repos on bioinformatics-focused GitHub profiles (search "bioinformatics" on GitHub, sort by stars) for README structure ideas — not to copy, but to calibrate what a reviewer expects to see in the first 30 seconds

## Job search
- Search "bioinformatics scientist" / "bioinformatics analyst" on LinkedIn/Indeed and note recurring tool requirements — feed back into what you prioritize
- GitHub is your portfolio: pin the 2 flagship pipeline repos, make sure the READMEs are the first thing a reviewer sees
