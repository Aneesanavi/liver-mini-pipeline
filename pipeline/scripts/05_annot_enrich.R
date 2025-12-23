#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
res_file <- ifelse(length(args) >= 1, args[1], "results/deseq2_results.tsv")
outdir   <- ifelse(length(args) >= 2, args[2], "results")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(res_file))

res_df <- read.delim(res_file, check.names = FALSE)
stopifnot(all(c("gene_id","padj","log2FoldChange") %in% colnames(res_df)))

# ---- MINIMAL annotation only ----
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("org.Hs.eg.db","AnnotationDbi")) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

res_df$ensembl <- sub("\\..*$", "", res_df$gene_id)

ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_df$ensembl),
  keytype = "ENSEMBL",
  columns = c("ENTREZID","SYMBOL")
)
ann <- ann[!duplicated(ann$ENSEMBL), ]
res_annot <- merge(res_df, ann, by.x="ensembl", by.y="ENSEMBL", all.x=TRUE, sort=FALSE)

write.table(res_annot, file.path(outdir, "deseq2_results_annotated.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Top lists
res_rank <- res_df[!is.na(res_df$padj), ]
res_rank$absLFC <- abs(res_rank$log2FoldChange)
res_rank <- res_rank[order(res_rank$padj, -res_rank$absLFC), ]

write.table(head(res_rank[, c("gene_id","padj","log2FoldChange")], 10),
            file.path(outdir, "top10_genes_by_padj.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

write.table(head(res_rank[, c("gene_id","padj","log2FoldChange")], 50),
            file.path(outdir, "top50_genes_by_padj.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create empty enrichment placeholders so Nextflow outputs still match
write.table(data.frame(), file.path(outdir, "enrichGO_BP.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(), file.path(outdir, "enrichKEGG.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("OK: annotation done; enrichment skipped (placeholders written) ->", outdir, "\n")
