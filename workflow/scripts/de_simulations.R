#!/usr/bin/env Rscript
## de_simulations.R
##
## Count-level power benchmark for repeat differential expression as a function
## of gene library size and repeat library size. Reads two pre-computed count
## matrices (rows = features, cols = samples), fits negative-binomial mean and
## dispersion from each, then sweeps a 2D grid of multiplicative library-size
## scalers. At each grid cell it simulates new counts with a known fold change
## planted at a fixed set of repeat features, runs several normalization
## strategies, and reports power and false positive rate over n_iter replicates.
##
## See docs/de_simulations.md for the plain-English description.

suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(Biobase)
  library(matrixStats)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(viridis)
})


load_count_matrix <- function(path, samples = NULL) {
  stopifnot(file.exists(path))
  m <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE)
  m <- as.matrix(m)
  if (!is.null(samples)) {
    stopifnot(all(samples %in% colnames(m)))
    m <- m[, samples, drop = FALSE]
  }
  storage.mode(m) <- "integer"
  stopifnot(all(is.finite(m)), all(m >= 0))
  m
}


fit_nb_params <- function(counts, design) {
  stopifnot(is.matrix(counts), ncol(counts) == nrow(design))
  min_present <- max(2L, ncol(counts) %/% 4L)
  keep <- rowSums(counts >= 1) >= min_present
  counts <- counts[keep, , drop = FALSE]
  stopifnot(nrow(counts) >= 2)
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  mu <- pmax(rowMeans(counts), 1e-6)
  disp <- dge$tagwise.dispersion
  if (is.null(disp) || length(disp) != length(mu)) {
    disp <- rep_len(dge$common.dispersion, length(mu))
  }
  list(mu = mu, dispersion = disp, feature_id = rownames(counts))
}


simulate_count_matrix <- function(mu, dispersion, n_samples, lib_scale = 1,
                                  w_loadings = NULL, w_factor = NULL,
                                  de_logfc = NULL, condition = NULL,
                                  seed = 1L) {
  set.seed(seed)
  n_feat <- length(mu)
  stopifnot(length(dispersion) == n_feat, n_samples >= 2)
  log_mu <- matrix(log(mu * lib_scale), nrow = n_feat, ncol = n_samples)
  if (!is.null(w_loadings) && !is.null(w_factor)) {
    stopifnot(length(w_loadings) == n_feat, length(w_factor) == n_samples)
    log_mu <- log_mu + outer(w_loadings, w_factor)
  }
  if (!is.null(de_logfc) && !is.null(condition)) {
    stopifnot(is.factor(condition), length(condition) == n_samples,
              length(de_logfc) == n_feat, length(levels(condition)) == 2)
    treated <- as.numeric(condition == levels(condition)[2])
    log_mu <- log_mu + outer(de_logfc, treated)
  }
  size <- 1 / pmax(dispersion, 1e-4)
  mu_mat <- exp(log_mu)
  counts <- matrix(rnbinom(n_feat * n_samples,
                           size = rep(size, n_samples),
                           mu = as.vector(mu_mat)),
                   nrow = n_feat)
  storage.mode(counts) <- "integer"
  counts
}


pick_empirical_controls <- function(counts, design, n_controls) {
  stopifnot(qr(design)$rank == ncol(design),
            n_controls >= 1, n_controls <= nrow(counts))
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  rownames(tt)[order(tt$PValue, decreasing = TRUE)][seq_len(n_controls)]
}


run_de <- function(counts, design) {
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  topTags(qlf, n = Inf, sort.by = "none")$table
}


run_de_no_norm <- function(counts, design) {
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- rep(1, ncol(counts))
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  topTags(qlf, n = Inf, sort.by = "none")$table
}


## Power and FPR scoring. planted_features is the FULL planted-truth set,
## including features that might have been filtered out before testing.
## Filtered planted features count as FN, which is the correct behavior when
## comparing across grid cells with very different repeat library sizes.
score_calls <- function(tt_df, planted_features, fdr = 0.05) {
  called <- rownames(tt_df)[!is.na(tt_df$FDR) & tt_df$FDR <= fdr]
  truth <- planted_features
  null <- setdiff(rownames(tt_df), truth)
  tp <- sum(called %in% truth)
  fn <- length(truth) - tp
  fp <- sum(called %in% null)
  tn <- length(null) - fp
  list(power = tp / max(1L, tp + fn),
       fpr = fp / max(1L, fp + tn),
       tp = tp, fn = fn, fp = fp, tn = tn)
}


run_one_grid_cell <- function(gene_params, repeat_params, samples, condition,
                              gene_lib_scale, repeat_lib_scale,
                              fc, n_de_repeats, n_de_genes,
                              sigma_w, n_controls, k_max, fdr,
                              n_iter, seed_base) {
  n_samples <- length(samples)
  stopifnot(length(condition) == n_samples)
  cond_factor <- factor(condition, levels = sort(unique(condition)))
  stopifnot(length(levels(cond_factor)) == 2)
  design <- model.matrix(~ cond_factor)

  set.seed(seed_base)
  w_factor <- rnorm(n_samples)
  gene_w_loadings <- rnorm(length(gene_params$mu), 0, sigma_w)
  repeat_w_loadings <- rnorm(length(repeat_params$mu), 0, sigma_w)

  set.seed(seed_base + 1L)
  n_de_repeats <- min(n_de_repeats, length(repeat_params$mu))
  de_repeat_idx <- sample.int(length(repeat_params$mu), n_de_repeats)
  de_repeat_features <- repeat_params$feature_id[de_repeat_idx]
  repeat_logfc <- numeric(length(repeat_params$mu))
  repeat_logfc[de_repeat_idx] <- log(fc)

  n_de_genes <- min(n_de_genes, length(gene_params$mu))
  de_gene_idx <- sample.int(length(gene_params$mu), n_de_genes)
  gene_logfc <- numeric(length(gene_params$mu))
  gene_logfc[de_gene_idx] <- log(2) * sample(c(-1, 1), n_de_genes, replace = TRUE)

  iter_results <- list()
  for (i in seq_len(n_iter)) {
    seed_i <- seed_base + 100L + i
    gene_counts <- simulate_count_matrix(
      mu = gene_params$mu, dispersion = gene_params$dispersion,
      n_samples = n_samples, lib_scale = gene_lib_scale,
      w_loadings = gene_w_loadings, w_factor = w_factor,
      de_logfc = gene_logfc, condition = cond_factor, seed = seed_i)
    rownames(gene_counts) <- gene_params$feature_id
    colnames(gene_counts) <- samples

    repeat_counts <- simulate_count_matrix(
      mu = repeat_params$mu, dispersion = repeat_params$dispersion,
      n_samples = n_samples, lib_scale = repeat_lib_scale,
      w_loadings = repeat_w_loadings, w_factor = w_factor,
      de_logfc = repeat_logfc, condition = cond_factor, seed = seed_i + 1L)
    rownames(repeat_counts) <- repeat_params$feature_id
    colnames(repeat_counts) <- samples

    keep_g <- rowSums(gene_counts >= 1) >= 2L
    gene_counts_f <- gene_counts[keep_g, , drop = FALSE]
    keep_r <- rowSums(repeat_counts >= 1) >= 2L
    repeat_counts_f <- repeat_counts[keep_r, , drop = FALSE]

    method_scores <- list()

    if (nrow(repeat_counts_f) >= 2 && qr(design)$rank == ncol(design)) {
      method_scores[["none"]] <- score_calls(
        run_de_no_norm(repeat_counts_f, design), de_repeat_features, fdr)
      method_scores[["TMM"]] <- score_calls(
        run_de(repeat_counts_f, design), de_repeat_features, fdr)
    } else {
      method_scores[["none"]] <- list(power = NA, fpr = NA,
                                      tp = NA, fn = NA, fp = NA, tn = NA)
      method_scores[["TMM"]] <- method_scores[["none"]]
    }

    n_ctrl <- min(n_controls, max(0L, nrow(gene_counts_f) - 1L))
    if (n_ctrl >= 1 && nrow(gene_counts_f) >= 2) {
      controls <- pick_empirical_controls(gene_counts_f, design, n_ctrl)
      pheno <- AnnotatedDataFrame(data.frame(
        condition = cond_factor, row.names = samples))
      set1 <- newSeqExpressionSet(counts = gene_counts_f, phenoData = pheno)
      for (k in seq_len(k_max)) {
        score <- list(power = NA, fpr = NA, tp = NA, fn = NA, fp = NA, tn = NA)
        ruv_fit <- tryCatch(RUVg(set1, controls, k = k), error = function(e) NULL)
        if (!is.null(ruv_fit)) {
          W <- as.matrix(pData(ruv_fit)[, grep("^W_", colnames(pData(ruv_fit))),
                                         drop = FALSE])
          design_ruv <- model.matrix(~ W + cond_factor)
          if (qr(design_ruv)$rank == ncol(design_ruv) &&
              nrow(repeat_counts_f) >= 2) {
            tt_ruv <- run_de(repeat_counts_f, design_ruv)
            score <- score_calls(tt_ruv, de_repeat_features, fdr)
          }
        }
        method_scores[[paste0("RUVg_k", k)]] <- score
      }
    } else {
      for (k in seq_len(k_max)) {
        method_scores[[paste0("RUVg_k", k)]] <- list(
          power = NA, fpr = NA, tp = NA, fn = NA, fp = NA, tn = NA)
      }
    }

    iter_df <- do.call(rbind, lapply(names(method_scores), function(m) {
      r <- method_scores[[m]]
      data.frame(method = m, iter = i,
                 power = r$power, fpr = r$fpr,
                 tp = r$tp, fn = r$fn, fp = r$fp, tn = r$tn,
                 stringsAsFactors = FALSE)
    }))
    iter_results[[i]] <- iter_df
  }
  do.call(rbind, iter_results)
}


run_grid <- function(gene_counts_path, repeat_counts_path, metadata,
                     gene_lib_grid, repeat_lib_grid, fc,
                     n_de_repeats, n_de_genes, sigma_w,
                     n_controls, k_max, fdr, n_iter, seed) {
  samples <- metadata$sample
  condition <- metadata$condition
  gene_mat <- load_count_matrix(gene_counts_path, samples)
  repeat_mat <- load_count_matrix(repeat_counts_path, samples)
  cond_factor <- factor(condition, levels = sort(unique(condition)))
  design <- model.matrix(~ cond_factor)
  gene_params <- fit_nb_params(gene_mat, design)
  repeat_params <- fit_nb_params(repeat_mat, design)

  cells <- expand.grid(gene_lib_scale = gene_lib_grid,
                       repeat_lib_scale = repeat_lib_grid,
                       KEEP.OUT.ATTRS = FALSE)
  results <- vector("list", nrow(cells))
  for (ci in seq_len(nrow(cells))) {
    g <- cells$gene_lib_scale[ci]
    r <- cells$repeat_lib_scale[ci]
    cell_seed <- as.integer(seed) * 1000L + ci * 10L
    cat(sprintf("[grid %d/%d] gene_lib_scale=%g repeat_lib_scale=%g\n",
                ci, nrow(cells), g, r))
    df <- run_one_grid_cell(
      gene_params, repeat_params, samples, condition,
      gene_lib_scale = g, repeat_lib_scale = r,
      fc = fc, n_de_repeats = n_de_repeats, n_de_genes = n_de_genes,
      sigma_w = sigma_w, n_controls = n_controls, k_max = k_max,
      fdr = fdr, n_iter = n_iter, seed_base = cell_seed)
    df$gene_lib_scale <- g
    df$repeat_lib_scale <- r
    results[[ci]] <- df
  }
  list(results = do.call(rbind, results),
       gene_params = gene_params,
       repeat_params = repeat_params)
}


summarize_grid <- function(results_df) {
  results_df %>%
    group_by(method, gene_lib_scale, repeat_lib_scale) %>%
    summarise(power_mean = mean(power, na.rm = TRUE),
              power_sd = sd(power, na.rm = TRUE),
              fpr_mean = mean(fpr, na.rm = TRUE),
              fpr_sd = sd(fpr, na.rm = TRUE),
              n_iter = n(),
              .groups = "drop")
}


plot_heatmap <- function(summary_df, metric = c("power_mean", "fpr_mean"),
                         title = NULL) {
  metric <- match.arg(metric)
  ggplot(summary_df, aes(x = factor(gene_lib_scale),
                         y = factor(repeat_lib_scale),
                         fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", .data[[metric]])), size = 2.5) +
    facet_wrap(~ method) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = "gene library size scaler",
         y = "repeat library size scaler",
         fill = metric, title = title) +
    theme_minimal()
}


# Snakemake glue. Only runs when invoked via the snakemake script: directive.
if (exists("snakemake")) {
  log_path <- snakemake@log[[1]]
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_path, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  params <- snakemake@params
  outdir <- params$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  gene_counts_path <- snakemake@input$gene_counts
  repeat_counts_path <- snakemake@input$repeat_counts

  if (nzchar(params$metadata)) {
    md_raw <- read.table(params$metadata, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
    metadata <- data.frame(
      sample = md_raw[[params$sample_column]],
      condition = md_raw[[params$condition_column]],
      stringsAsFactors = FALSE)
  } else {
    header <- read.table(gene_counts_path, header = TRUE, sep = "\t",
                         nrows = 1, check.names = FALSE)
    samples <- colnames(header)[-1]
    half <- length(samples) %/% 2
    metadata <- data.frame(
      sample = samples,
      condition = c(rep("A", half), rep("B", length(samples) - half)),
      stringsAsFactors = FALSE)
  }

  grid_out <- run_grid(
    gene_counts_path = gene_counts_path,
    repeat_counts_path = repeat_counts_path,
    metadata = metadata,
    gene_lib_grid = as.numeric(unlist(params$gene_lib_grid)),
    repeat_lib_grid = as.numeric(unlist(params$repeat_lib_grid)),
    fc = as.numeric(params$fc),
    n_de_repeats = as.integer(params$n_de_repeats),
    n_de_genes = as.integer(params$n_de_genes),
    sigma_w = as.numeric(params$sigma_w),
    n_controls = as.integer(params$n_controls),
    k_max = as.integer(params$k_max),
    fdr = as.numeric(params$fdr),
    n_iter = as.integer(params$n_iter),
    seed = as.integer(params$seed))

  results <- grid_out$results
  summary_df <- summarize_grid(results)

  write.table(results, snakemake@output$results_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_df, snakemake@output$summary_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(list(results = results, summary = summary_df,
               metadata = metadata,
               gene_params = grid_out$gene_params,
               repeat_params = grid_out$repeat_params,
               params = as.list(params)),
          snakemake@output$results_rds)

  pdf(snakemake@output$heatmap_pdf, width = 11, height = 7)
  print(plot_heatmap(summary_df, "power_mean",
                     sprintf("Power at FC=%g, FDR=%g",
                             as.numeric(params$fc), as.numeric(params$fdr))))
  print(plot_heatmap(summary_df, "fpr_mean",
                     sprintf("False positive rate at FDR=%g",
                             as.numeric(params$fdr))))
  dev.off()

  sink(type = "message")
  sink(type = "output")
  close(log_con)
}
