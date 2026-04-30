#!/usr/bin/env Rscript
## Unit tests for workflow/scripts/de_simulations.R
##
## Sources the helpers from the production script and uses them on small
## synthetic matrices. Exits non-zero on any failure. Designed to require only
## edgeR, RUVSeq, EDASeq, Biobase, and matrixStats from the edger conda env.

suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
  library(EDASeq)
  library(Biobase)
  library(matrixStats)
})

script_path <- Sys.getenv("DE_SIM_SCRIPT", unset = "")
if (!nzchar(script_path)) {
  script_path <- file.path(getwd(), "workflow", "scripts", "de_simulations.R")
}
stopifnot(file.exists(script_path))

env <- new.env(parent = globalenv())
source(script_path, local = env)


report <- function(name, ok) {
  cat(sprintf("%s %s\n", if (ok) "PASS" else "FAIL", name))
  if (!ok) quit(status = 1)
}


## simulate_count_matrix: shape and non-negativity
m <- env$simulate_count_matrix(mu = c(10, 50, 200, 1000),
                               dispersion = rep(0.1, 4),
                               n_samples = 6, lib_scale = 1, seed = 1)
report("simulate_count_matrix shape and >= 0",
       is.matrix(m) && nrow(m) == 4 && ncol(m) == 6 && all(m >= 0))

## lib_scale increases mean counts roughly proportionally
m1 <- env$simulate_count_matrix(rep(100, 200), rep(0.1, 200), 10, 1, seed = 2)
m5 <- env$simulate_count_matrix(rep(100, 200), rep(0.1, 200), 10, 5, seed = 2)
report("lib_scale of 5 raises mean counts above 3x of lib_scale 1",
       mean(m5) > 3 * mean(m1))

## planted DE recovered by score_calls when FDR is small only on planted rows
features <- paste0("f", seq_len(100))
fake_tt <- data.frame(FDR = c(rep(0.01, 10), rep(0.5, 90)), row.names = features)
sc <- env$score_calls(fake_tt, planted_features = features[1:10], fdr = 0.05)
report("score_calls perfect recovery",
       sc$power == 1 && sc$fpr == 0 && sc$tp == 10 && sc$fn == 0)

## complete miss: significant rows are entirely outside the planted set
sc2 <- env$score_calls(fake_tt, planted_features = features[91:100], fdr = 0.05)
report("score_calls complete miss",
       sc2$power == 0 && sc2$fp == 10 && sc2$tn == 80)

## planted features that are not in tt_df count as FN, not silently dropped
truth_with_missing <- c(features[1:10], "missing_feature")
sc3 <- env$score_calls(fake_tt, planted_features = truth_with_missing, fdr = 0.05)
report("score_calls penalizes filtered-out planted features",
       sc3$power == 10 / 11 && sc3$fn == 1)

## pick_empirical_controls returns the requested number of unique features
set.seed(7)
counts <- matrix(rpois(50 * 8, lambda = 100), nrow = 50)
rownames(counts) <- paste0("g", seq_len(50))
colnames(counts) <- paste0("s", seq_len(8))
design <- model.matrix(~ factor(rep(c("A", "B"), each = 4)))
ctrl <- env$pick_empirical_controls(counts, design, n_controls = 10)
report("pick_empirical_controls returns 10 unique gene IDs",
       length(ctrl) == 10 && !anyDuplicated(ctrl) && all(ctrl %in% rownames(counts)))

## end-to-end run on a tiny synthetic problem
gene_params <- list(
  mu = rep(c(50, 200), length.out = 60),
  dispersion = rep(0.1, 60),
  feature_id = paste0("g", seq_len(60)))
repeat_params <- list(
  mu = rep(c(20, 80), length.out = 30),
  dispersion = rep(0.2, 30),
  feature_id = paste0("r", seq_len(30)))
samples <- paste0("s", seq_len(8))
condition <- rep(c("A", "B"), each = 4)
res <- env$run_one_grid_cell(
  gene_params, repeat_params, samples, condition,
  gene_lib_scale = 1, repeat_lib_scale = 1,
  fc = 5, n_de_repeats = 10, n_de_genes = 15,
  sigma_w = 0.3, n_controls = 20, k_max = 1, fdr = 0.05,
  n_iter = 2, seed_base = 42)
report("run_one_grid_cell returns data frame with expected columns",
       is.data.frame(res)
       && all(c("power", "fpr", "method", "iter") %in% colnames(res))
       && all(c("none", "TMM", "RUVg_k1") %in% unique(res$method))
       && nrow(res) >= 6)

## summarize_grid aggregates iterations correctly
res$gene_lib_scale <- 1
res$repeat_lib_scale <- 1
sm <- env$summarize_grid(res)
report("summarize_grid produces one row per (method, gene_lib, repeat_lib)",
       nrow(sm) == length(unique(res$method))
       && all(c("power_mean", "power_sd", "fpr_mean", "fpr_sd") %in% colnames(sm)))

cat("OK: all de_simulations unit tests passed\n")
