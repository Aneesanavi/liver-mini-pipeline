#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
res_file  <- ifelse(length(args) >= 1, args[1], "results/deseq2_results.tsv")
norm_file <- ifelse(length(args) >= 2, args[2], "results/deseq2_normalized_counts.tsv")
meta_file <- ifelse(length(args) >= 3, args[3], "results/metadata_used.tsv")
outdir    <- ifelse(length(args) >= 4, args[4], "results")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(res_file), file.exists(norm_file), file.exists(meta_file))

# ---- 1) Load DESeq2 results + metadata used ----
res_df <- read.delim(res_file, check.names = FALSE)
meta <- read.delim(meta_file, check.names = FALSE, colClasses = c("character", "character"))

stopifnot(all(c("sample","group") %in% colnames(meta)))
meta <- meta[meta$group %in% c("control","disease"), , drop = FALSE]

# Keep sample IDs consistent (4 digits)
meta$sample <- sprintf("%04d", as.integer(meta$sample))

# ---- 2) Load normalized counts ----
norm_df <- read.delim(norm_file, check.names = FALSE)
stopifnot("gene_id" %in% colnames(norm_df))

genes <- norm_df$gene_id
mat <- as.matrix(norm_df[, setdiff(colnames(norm_df), "gene_id"), drop = FALSE])
rownames(mat) <- genes

# Ensure sample colnames match metadata style
colnames(mat) <- sprintf("%04d", as.integer(colnames(mat)))

# Align samples to metadata
stopifnot(all(meta$sample %in% colnames(mat)))
mat <- mat[, meta$sample, drop = FALSE]
stopifnot(all(colnames(mat) == meta$sample))

# ---- 3) Choose top 50 genes (padj + effect size) ----
res_df <- res_df[!is.na(res_df$padj), ]
res_df$absLFC <- abs(res_df$log2FoldChange)

res_top <- subset(res_df, padj < 0.05 & absLFC >= 1)
if (nrow(res_top) < 50) res_top <- subset(res_df, padj < 0.05)

res_top <- res_top[order(res_top$padj, -res_top$absLFC), ]
top_genes <- head(res_top$gene_id, 50)

# ---- 4) Build scaled heatmap matrix ----
hm <- mat[top_genes, , drop = FALSE]
hm_log <- log2(hm + 1)
hm_scaled <- t(scale(t(hm_log)))
hm_scaled[is.na(hm_scaled)] <- 0

# Save exact matrix used (genes x samples)
heatmap_export <- data.frame(gene_id = rownames(hm_scaled), hm_scaled, check.names = FALSE)
write.table(heatmap_export, file.path(outdir, "top50_heatmap_matrix_scaled.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(data.frame(gene_id = top_genes),
            file.path(outdir, "top50_genes_for_heatmap.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 5) Heatmap PNG (base R) ----
ord <- order(meta$group)
hm_scaled_plot <- hm_scaled[, ord, drop = FALSE]
meta_ord <- meta[ord, , drop = FALSE]

png(file.path(outdir, "top50_heatmap.png"), width = 1200, height = 900)
par(mar = c(7, 7, 4, 2))
image(
  x = 1:ncol(hm_scaled_plot),
  y = 1:nrow(hm_scaled_plot),
  z = t(hm_scaled_plot[nrow(hm_scaled_plot):1, ]),
  axes = FALSE,
  xlab = "", ylab = "",
  main = "Top 50 DE genes (scaled log2 normalized counts)"
)
axis(1, at = 1:ncol(hm_scaled_plot), labels = meta_ord$group, las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(hm_scaled_plot), labels = rev(rownames(hm_scaled_plot)), las = 2, cex.axis = 0.5)
box()

grp_changes <- which(meta_ord$group[-1] != meta_ord$group[-nrow(meta_ord)])
if (length(grp_changes) > 0) abline(v = grp_changes + 0.5, lty = 2)
dev.off()

# ---- 6) Add SYMBOL annotation for top 50 genes ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)

suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))

top_ensembl <- sub("\\..*$", "", top_genes)

ann_top <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(top_ensembl),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)
ann_top <- ann_top[!duplicated(ann_top$ENSEMBL), ]

symbol_map <- setNames(ann_top$SYMBOL, ann_top$ENSEMBL)
top_symbols <- symbol_map[top_ensembl]
top_symbols[is.na(top_symbols)] <- top_ensembl[is.na(top_symbols)]

write.table(
  data.frame(gene_id = top_genes, ensembl = top_ensembl, symbol = top_symbols),
  file.path(outdir, "top50_gene_symbol_map.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("OK: Heatmap + exports saved in:", outdir, "\n")
cat("Counts:", ncol(hm_scaled), "samples,", nrow(hm_scaled), "genes\n")
print(table(meta$group))
