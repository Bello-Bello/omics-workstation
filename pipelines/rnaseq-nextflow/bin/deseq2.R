#!/usr/bin/env Rscript
# Called as: deseq2.R <samples.tsv> <annotation.gtf>
# Nextflow stages each sample's salmon output as a directory named after the
# sample (e.g. WT_REP1/quant.sf) in this process's work directory — no
# `snakemake@` object here, just plain command-line args and relative paths.
args <- commandArgs(trailingOnly = TRUE)
samples_file <- args[1]
gtf_file <- args[2]

library(tximport)
library(DESeq2)
library(ggplot2)

samples <- read.delim(samples_file, stringsAsFactors = FALSE)
quant_files <- file.path(samples$sample, "quant.sf")
names(quant_files) <- samples$sample

tx2gene <- rtracklayer::import(gtf_file)
tx2gene <- as.data.frame(tx2gene)
tx2gene <- unique(tx2gene[tx2gene$type == "transcript", c("transcript_id", "gene_id")])

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)
dds <- DESeq(dds)

res <- results(dds)
write.csv(as.data.frame(res), "results.csv")

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("Sample PCA")
ggsave("pca.png", p, width = 6, height = 5)
