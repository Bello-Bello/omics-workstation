# Runs inside Snakemake's `script:` directive — `snakemake` object is injected.
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tximport)
library(DESeq2)
library(ggplot2)

samples <- read.delim(snakemake@input[["samples"]], stringsAsFactors = FALSE)
quant_files <- snakemake@input[["quants"]]
names(quant_files) <- samples$sample

# tx2gene mapping from the GTF — swap for a real annotation package mapping
# if your reference uses a different transcript ID convention.
tx2gene <- rtracklayer::import(snakemake@input[["gtf"]])
tx2gene <- as.data.frame(tx2gene)
tx2gene <- unique(tx2gene[tx2gene$type == "transcript", c("transcript_id", "gene_id")])

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)
dds <- DESeq(dds)

res <- results(dds)
write.csv(as.data.frame(res), snakemake@output[["results"]])

# Per-sample normalized counts (not just the per-gene summary stats above) —
# needed downstream for the stats module, which runs its own hypothesis
# tests per gene and compares them against DESeq2's.
write.csv(counts(dds, normalized = TRUE), snakemake@output[["counts"]])

# vst() needs ~1000+ genes to subsample from; with 124 real genes here,
# varianceStabilizingTransformation() (the exact, non-subsampled version)
# is the appropriate choice instead.
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("Sample PCA")
ggsave(snakemake@output[["pca"]], p, width = 6, height = 5)

sink()
sink(type = "message")
