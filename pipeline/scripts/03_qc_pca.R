#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
counts_file <- ifelse(length(args) >= 1, args[1], "data/raw_counts.tsv")
meta_file   <- ifelse(length(args) >= 2, args[2], "data/metadata.tsv")
outdir      <- ifelse(length(args) >= 3, args[3], "results")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(counts_file), file.exists(meta_file))

counts <- read.delim(counts_file, check.names = FALSE)
meta <- read.delim(meta_file, check.names = FALSE, colClasses = c("character","character"))

gene_col <- colnames(counts)[1]
genes <- counts[[gene_col]]
mat <- as.matrix(counts[, -1])
rownames(mat) <- genes

stopifnot(all(colnames(mat) %in% meta$sample))
meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
stopifnot(all(meta$sample == colnames(mat)))

libsize <- colSums(mat)
png(file.path(outdir, "qc_library_size.png"), width=900, height=500)
plot(libsize, ylab="Library size (total counts)", xlab="Sample index"); abline(h=median(libsize), lty=2)
dev.off()

detected_genes <- colSums(mat > 0)
png(file.path(outdir, "qc_detected_genes.png"), width=900, height=500)
plot(detected_genes, ylab="Genes with count > 0", xlab="Sample index"); abline(h=median(detected_genes), lty=2)
dev.off()

cpm <- sweep(mat, 2, libsize / 1e6, "/")
log_cpm <- log1p(cpm)

gene_var <- apply(log_cpm, 1, var)
log_cpm2 <- log_cpm[gene_var > 0, , drop = FALSE]

pca <- prcomp(t(log_cpm2), scale. = TRUE)
pc_df <- data.frame(sample=colnames(log_cpm2), PC1=pca$x[,1], PC2=pca$x[,2], group=meta$group, stringsAsFactors=FALSE)

png(file.path(outdir, "pca_logcpm.png"), width=900, height=600)
plot(pc_df$PC1, pc_df$PC2, xlab="PC1", ylab="PC2", main="PCA (log1p-CPM; colored by group)", pch=19, col=as.factor(pc_df$group))
legend("topright", legend=levels(as.factor(pc_df$group)), pch=19, col=1:length(levels(as.factor(pc_df$group))))
dev.off()

keep <- pc_df$group %in% c("control","disease")
pc_df2 <- pc_df[keep, , drop=FALSE]
png(file.path(outdir, "pca_logcpm_control_vs_disease.png"), width=900, height=600)
plot(pc_df2$PC1, pc_df2$PC2, xlab="PC1", ylab="PC2", main="PCA (control vs disease only)", pch=19, col=as.factor(pc_df2$group))
legend("topright", legend=levels(as.factor(pc_df2$group)), pch=19, col=1:length(levels(as.factor(pc_df2$group))))
dev.off()

cor_mat <- cor(log_cpm2, method="pearson")
png(file.path(outdir, "sample_correlation.png"), width=900, height=800)
image(cor_mat, main="Sample correlation (log1p-CPM)", axes=FALSE)
dev.off()

qc_summary <- data.frame(sample=colnames(mat), group=meta$group, library_size=libsize, detected_genes=detected_genes, stringsAsFactors=FALSE)
write.table(qc_summary, file.path(outdir, "qc_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("OK: QC + PCA outputs written to", outdir, "\n")
