#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
meta_file <- ifelse(length(args) >= 1, args[1], "data/metadata.tsv")
hm_file   <- ifelse(length(args) >= 2, args[2], "results/top50_heatmap_matrix_scaled.tsv")
outdir    <- ifelse(length(args) >= 3, args[3], "results")
USE_PCA   <- TRUE
N_PC      <- 10
K_FOLD    <- 5
THR       <- 0.5

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(meta_file), file.exists(hm_file))

safe_id <- function(x) sprintf("%04d", as.integer(as.character(x)))

metrics_binary <- function(y_true, prob, thr = 0.5) {
  stopifnot(all(y_true %in% c(0,1)))
  pred <- ifelse(prob >= thr, 1, 0)
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  
  acc <- (tp + tn) / max(1, (tp + tn + fp + fn))
  sens <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
  spec <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  
  auc <- NA
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(response = y_true, predictor = prob, quiet = TRUE)
    auc <- as.numeric(pROC::auc(roc_obj))
  }
  
  list(acc=acc, bal_acc=bal_acc, auc=auc, tp=tp, tn=tn, fp=fp, fn=fn)
}

meta <- read.delim(meta_file, check.names = FALSE, colClasses = c("character","character"))
stopifnot(all(c("sample","group") %in% colnames(meta)))
meta <- meta[meta$group %in% c("control","disease"), , drop = FALSE]
meta$sample <- safe_id(meta$sample)

hm_df <- read.delim(hm_file, check.names = FALSE)
stopifnot("gene_id" %in% colnames(hm_df))

genes <- hm_df$gene_id
Xg <- as.matrix(hm_df[, setdiff(colnames(hm_df), "gene_id"), drop = FALSE])
rownames(Xg) <- genes

X <- t(Xg)
rownames(X) <- safe_id(rownames(X))
storage.mode(X) <- "double"

common <- intersect(meta$sample, rownames(X))
stopifnot(length(common) >= 10)

meta <- meta[match(common, meta$sample), , drop=FALSE]
X <- X[common, , drop=FALSE]
stopifnot(all(meta$sample == rownames(X)))

y <- ifelse(meta$group == "disease", 1, 0)

if (USE_PCA) {
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  max_pc <- min(ncol(pca$x), nrow(pca$x) - 1)
  k <- min(N_PC, max_pc)
  X_feat <- pca$x[, 1:k, drop=FALSE]
  colnames(X_feat) <- paste0("PC", seq_len(k))
  feat_name <- paste0("PCA(", k, " PCs)")
} else {
  X_feat <- X
  feat_name <- "Top50 genes (scaled)"
}

feat_out <- data.frame(sample=rownames(X_feat), group=meta$group, X_feat, check.names=FALSE)
write.table(feat_out, file.path(outdir,"ml_features_used.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- Holdout split ----
set.seed(42)
idx <- sample(seq_len(nrow(X_feat)))
train_n <- floor(0.8 * nrow(X_feat))
train_idx <- idx[1:train_n]
test_idx  <- idx[(train_n + 1):nrow(X_feat)]

X_train <- X_feat[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X_feat[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

df_train <- data.frame(y=y_train, X_train)
df_test  <- data.frame(y=y_test,  X_test)

fit <- glm(y ~ ., data=df_train, family=binomial())
prob_test <- predict(fit, newdata=df_test, type="response")
m_hold <- metrics_binary(y_test, prob_test, thr=THR)

writeLines(c(
  paste0("Features: ", feat_name),
  "Model: logistic regression",
  paste0("Train samples: ", nrow(X_train)),
  paste0("Test samples: ", nrow(X_test)),
  paste0("Holdout Accuracy: ", round(m_hold$acc, 3)),
  paste0("Holdout Balanced Accuracy: ", round(m_hold$bal_acc, 3)),
  paste0("Holdout AUC: ", ifelse(is.na(m_hold$auc), "NA (install pROC)", round(m_hold$auc, 3))),
  paste0("Confusion tp/tn/fp/fn: ", m_hold$tp, "/", m_hold$tn, "/", m_hold$fp, "/", m_hold$fn)
), file.path(outdir,"ml_holdout_logreg_summary.txt"))

# ---- K-fold CV ----
set.seed(42)
fold_id <- sample(rep(1:K_FOLD, length.out = nrow(X_feat)))

accs <- numeric(K_FOLD)
balaccs <- numeric(K_FOLD)
aucs <- rep(NA_real_, K_FOLD)

for (k in 1:K_FOLD) {
  te <- which(fold_id == k)
  tr <- setdiff(seq_len(nrow(X_feat)), te)
  
  df_tr <- data.frame(y=y[tr], X_feat[tr,,drop=FALSE])
  df_te <- data.frame(y=y[te], X_feat[te,,drop=FALSE])
  
  fit_k <- glm(y ~ ., data=df_tr, family=binomial())
  p_k <- predict(fit_k, newdata=df_te, type="response")
  
  mk <- metrics_binary(df_te$y, p_k, thr=THR)
  accs[k] <- mk$acc
  balaccs[k] <- mk$bal_acc
  aucs[k] <- mk$auc
}

cv_lines <- c(
  paste0("Features: ", feat_name),
  "Model: logistic regression",
  paste0("K-fold: ", K_FOLD),
  paste0("CV Accuracy mean: ", round(mean(accs), 3), " | SD: ", round(sd(accs), 3)),
  paste0("CV Balanced Accuracy mean: ", round(mean(balaccs), 3), " | SD: ", round(sd(balaccs), 3)),
  paste0("CV AUC mean: ", ifelse(all(is.na(aucs)), "NA (install pROC)", paste0(round(mean(aucs, na.rm=TRUE), 3), " | SD: ", round(sd(aucs, na.rm=TRUE), 3)))),
  paste0("Fold Accuracy: ", paste(round(accs, 3), collapse=", ")),
  paste0("Fold Balanced Acc: ", paste(round(balaccs, 3), collapse=", ")),
  paste0("Fold AUC: ", paste(ifelse(is.na(aucs), "NA", round(aucs, 3)), collapse=", "))
)
writeLines(cv_lines, file.path(outdir,"ml_cv_logreg_summary.txt"))

# ---- SHAP (fastshap) on training set ----
if (!requireNamespace("fastshap", quietly=TRUE)) {
  install.packages("fastshap")
}
if (requireNamespace("fastshap", quietly=TRUE)) {
  pred_wrapper <- function(object, newdata) predict(object, newdata=newdata, type="response")
  
  set.seed(42)
  shap_values <- fastshap::explain(
    object = fit,
    X = as.data.frame(X_train),
    pred_wrapper = pred_wrapper,
    nsim = 200,
    adjust = TRUE
  )
  
  shap_df <- cbind(sample=rownames(X_train), shap_values)
  write.table(shap_df, file.path(outdir,"shap_values_train.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  
  shap_imp <- data.frame(
    feature=colnames(shap_values),
    mean_abs_shap=apply(abs(shap_values), 2, mean),
    stringsAsFactors=FALSE
  )
  shap_imp <- shap_imp[order(-shap_imp$mean_abs_shap), , drop=FALSE]
  write.table(shap_imp, file.path(outdir,"shap_importance.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  
  png(file.path(outdir,"shap_importance_bar.png"), width=1000, height=650)
  par(mar=c(8,6,4,2))
  barplot(shap_imp$mean_abs_shap, names.arg=shap_imp$feature, las=2,
          main="SHAP importance (mean |SHAP|) — Logistic Regression", ylab="mean(|SHAP|)")
  dev.off()
}

coef_df <- data.frame(term=names(coef(fit)), coef=as.numeric(coef(fit)), stringsAsFactors=FALSE)
write.table(coef_df, file.path(outdir,"logreg_coefficients.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("OK: ML outputs written to", outdir, "\n")
