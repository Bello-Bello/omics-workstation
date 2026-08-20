# BioContainers' single-package images (like bioconductor-deseq2) are built
# via a stripped multi-stage process with no conda/mamba left inside them —
# fine for one tool, but can't be extended with `mamba install` afterward.
# Building our own environment from a full conda base is more reliable when
# multiple R/Bioconductor packages need to be installed together.
FROM condaforge/miniforge3:24.9.2-0

RUN mamba create -y -n deseq2 -c bioconda -c conda-forge \
    bioconductor-deseq2=1.42.* \
    bioconductor-tximport=1.30.* \
    bioconductor-rtracklayer=1.62.* \
    r-ggplot2=3.5.* \
    && mamba clean -afy

# Put the environment's binaries first on PATH so `Rscript` etc. resolve to
# this env without needing `conda activate` in every RUN/CMD.
ENV PATH=/opt/conda/envs/deseq2/bin:$PATH
