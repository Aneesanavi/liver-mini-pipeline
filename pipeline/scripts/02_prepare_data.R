#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
counts_path <- ifelse(length(args) >= 1, args[1], "data/GSE126848/GSE126848_Gene_counts_raw.txt")
series_path <- ifelse(length(args) >= 2, args[2], "data/GSE126848_series_matrix.txt")
out_counts  <- ifelse(length(args) >= 3, args[3], "data/raw_counts.tsv")
out_meta    <- ifelse(length(args) >= 4, args[4], "data/metadata.tsv")
out_map     <- ifelse(length(args) >= 5, args[5], "data/GSE126848/sample_map.tsv")

stopifnot(file.exists(counts_path))
stopifnot(file.exists(series_path))

dir.create("data/GSE126848", showWarnings = FALSE, recursive = TRUE)

counts <- read.delim(counts_path, check.names = FALSE)
write.table(counts, out_counts, sep="\t", row.names=FALSE, quote=FALSE)

sample_names <- colnames(counts)[-1]
sample_names <- sprintf("%04d", as.integer(sample_names))

metadata <- data.frame(sample = sample_names, group = "unknown", stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

eset <- getGEO(filename = series_path, GSEMatrix = TRUE)
stopifnot(inherits(eset, "ExpressionSet"))
pd <- pData(eset)

stopifnot(nrow(pd) == length(sample_names))

disease_text <- pd[["characteristics_ch1.2"]]
group <- ifelse(
  grepl("control|healthy|normal", disease_text, ignore.case = TRUE), "control",
  ifelse(grepl("NAFL|NAFLD|NASH|MASH", disease_text, ignore.case = TRUE), "disease", "unknown")
)

map <- data.frame(
  sample = sample_names,
  geo_accession = pd$geo_accession,
  title = pd$title,
  disease_text = disease_text,
  group = group,
  stringsAsFactors = FALSE
)

write.table(map, out_map, sep="\t", row.names=FALSE, quote=FALSE)

metadata$group <- map$group
write.table(metadata, out_meta, sep="\t", row.names=FALSE, quote=FALSE)

meta_check <- read.delim(out_meta, check.names = FALSE, colClasses = c("character","character"))
mat_sample_names <- sprintf("%04d", as.integer(colnames(counts)[-1]))
stopifnot(all(mat_sample_names == meta_check$sample))

cat("OK: wrote\n")
cat("-", out_counts, "\n")
cat("-", out_meta, "\n")
cat("-", out_map, "\n")
cat("Group counts:\n")
print(table(meta_check$group))
