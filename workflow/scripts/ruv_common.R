## Shared helpers for RUVg-based repeat DE analyses.
## Sourced from:
##   paper/scripts/ruv_gse230647_bulk_report.Rmd
##   paper/scripts/ruv_gse126543_bulk_report.Rmd
##   workflow/scripts/gse230647_sc_report.Rmd
##
## Expects edgeR, RUVSeq, and matrixStats to be attached by the caller.

load_counts_tsv <- function(path, samples, require_all = TRUE) {
  stopifnot(file.exists(path))
  m <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(m)
  if (require_all) {
    stopifnot(setequal(colnames(m), samples))
  } else {
    stopifnot(all(samples %in% colnames(m)))
  }
  m <- m[, samples, drop = FALSE]
  storage.mode(m) <- "integer"
  stopifnot(all(is.finite(m)), all(m >= 0))
  m
}

rle_score <- function(mat) {
  med <- matrixStats::rowMedians(log2(mat + 1))
  rle <- sweep(log2(mat + 1), 1, med, "-")
  mean(apply(rle, 2, function(x) mean(abs(x - median(x)))))
}

## Fit RUVg across a range of k values. Returns the per-k RUVg fits and an
## RLE-style mean-abs-deviation diagnostic table (k = 0 is the unnormalized
## baseline).
fit_ruv_range <- function(genes_expr, controls, pheno_df, ks = 1:4) {
  stopifnot(is.matrix(genes_expr), ncol(genes_expr) == nrow(pheno_df),
            length(controls) > 0, all(controls %in% rownames(genes_expr)))
  set1 <- RUVSeq::newSeqExpressionSet(counts = genes_expr, phenoData = pheno_df)
  fits <- list()
  for (k in ks) {
    fits[[as.character(k)]] <- RUVSeq::RUVg(set1, controls, k = k)
  }
  rle_before <- rle_score(genes_expr)
  rle_after <- vapply(fits, function(s) rle_score(RUVSeq::normCounts(s)),
                      numeric(1))
  diag_rle <- data.frame(k = c(0, ks),
                         rle_mean_abs_dev = c(rle_before, rle_after))
  list(fits = fits, diag_rle = diag_rle)
}

pick_W_matrix <- function(ruv_fits, k_selected, samples) {
  key <- as.character(k_selected)
  stopifnot(key %in% names(ruv_fits))
  pdata <- Biobase::pData(ruv_fits[[key]])
  W <- as.matrix(pdata[, grep("^W_", colnames(pdata)), drop = FALSE])
  stopifnot(nrow(W) == length(samples), identical(rownames(W), samples),
            ncol(W) == k_selected)
  W
}

## Rank empirical control genes by naive-DE p-value: genes with the highest
## p-values are the least associated with the contrast. `test` is passed to
## glmQLFTest as either `coef = ...` or `contrast = ...`.
rank_empirical_controls <- function(counts, design, test, n_controls) {
  stopifnot(qr(design)$rank == ncol(design), n_controls <= nrow(counts))
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- do.call(edgeR::glmQLFTest, c(list(fit), test))
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  controls <- rownames(tt)[order(tt$PValue, decreasing = TRUE)][seq_len(n_controls)]
  stopifnot(length(controls) == n_controls, !anyDuplicated(controls))
  list(naive_de = tt, controls = controls)
}

## Fit edgeR with RUV covariates on a single repeat count matrix and return a
## ranked topTags data frame with feature_id column.
run_ruv_repeat_de <- function(counts, design, test) {
  stopifnot(qr(design)$rank == ncol(design))
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- do.call(edgeR::glmQLFTest, c(list(fit), test))
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tt$feature_id <- rownames(tt)
  tt[order(tt$PValue), c("feature_id","logFC","logCPM","F","PValue","FDR")]
}
