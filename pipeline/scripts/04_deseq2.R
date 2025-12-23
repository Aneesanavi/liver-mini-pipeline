#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
counts_file <- ifelse(length(args) >= 1, args[1], "data/raw_counts.tsv")
meta_file   <- ifelse(length(args) >= 2, args[2], "data/metadata.tsv")
outdir      <- ifelse(length(args) >= 3, args[3], "results")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(counts_file), file.exists(meta_file))

counts <- read.delim(counts_file, check.names = FALSE)
meta <- read.delim(meta_file, check.names = FALSE, colClasses = c("character","character"))

# Require columns
stopifnot(all(c("sample","group") %in% colnames(meta)))

# Keep sample IDs consistent (avoid losing leading zeros)
meta$sample <- sprintf("%04d", as.integer(meta$sample))
colnames(counts)[-1] <- sprintf("%04d", as.integer(colnames(counts)[-1]))

gene_col <- colnames(counts)[1]
genes <- counts[[gene_col]]

mat <- as.matrix(counts[, -1, drop = FALSE])
rownames(mat) <- genes

# Align metadata to count columns
stopifnot(all(colnames(mat) %in% meta$sample))
meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
stopifnot(all(meta$sample == colnames(mat)))

# Control vs disease only
keep <- meta$group %in% c("control","disease")
meta2 <- meta[keep, , drop = FALSE]
mat2  <- mat[, keep, drop = FALSE]

meta2$group <- factor(meta2$group, levels = c("control","disease"))
stopifnot(all(colnames(mat2) == meta2$sample))

# Packages
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(countData = round(mat2), colData = meta2, design = ~ group)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("group","disease","control"))

# Save results
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id","log2FoldChange","lfcSE","stat","pvalue","padj")]
res_df <- res_df[order(res_df$padj, res_df$pvalue), ]

write.table(res_df, file.path(outdir,"deseq2_results.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
write.table(sig, file.path(outdir,"deseq2_significant_padj_lt_0.05.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# MA plot
png(file.path(outdir,"deseq2_MAplot.png"), width=900, height=600)
plotMA(res, main="DESeq2 MA Plot (disease vs control)", ylim=c(-6,6))
dev.off()

# Volcano
vol <- res_df
vol$neglog10padj <- -log10(pmax(vol$padj, 1e-300))
vol$is_sig <- !is.na(vol$padj) & vol$padj < 0.05 & abs(vol$log2FoldChange) >= 1

png(file.path(outdir,"deseq2_volcano.png"), width=900, height=600)
plot(vol$log2FoldChange, vol$neglog10padj,
     xlab="log2 Fold Change (disease vs control)",
     ylab="-log10(adjusted p-value)",
     main="Volcano Plot (DESeq2)",
     pch=19)
abline(v=c(-1,1), lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()

# Normalized counts
norm_counts <- counts(dds, normalized=TRUE)
norm_df <- data.frame(gene_id=rownames(norm_counts), norm_counts, check.names=FALSE)
write.table(norm_df, file.path(outdir,"deseq2_normalized_counts.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Save metadata used (critical for later steps)
write.table(meta2, file.path(outdir,"metadata_used.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Summary
summary_lines <- c(
  paste0("Samples used: ", sum(meta2$group=="control"), " control; ", sum(meta2$group=="disease"), " disease"),
  paste0("Genes tested (after filter): ", nrow(dds)),
  paste0("Significant (padj < 0.05): ", sum(!is.na(res$padj) & res$padj < 0.05)),
  paste0("Significant (padj < 0.05 & |log2FC| >= 1): ", sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1))
)
writeLines(summary_lines, file.path(outdir,"deseq2_summary.txt"))

cat("OK: DESeq2 outputs written to", outdir, "\n")
