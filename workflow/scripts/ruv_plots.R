## Shared plotting and table helpers for RUVg edgeR top-tags tables.
## Callers must have attached ggplot2, ggrepel, dplyr, and knitr, and must have
## sourced ggtheme.R so that theme_ng() is on the search path.

ruv_volcano <- function(tt, title_tag, k,
                        label_top = 20, label_fdr_cutoff = 0.01,
                        sig_cutoff = 0.05) {
  stopifnot(all(c("logFC", "FDR", "feature_id") %in% colnames(tt)))
  tt$sig <- tt$FDR < sig_cutoff
  labels <- tt %>%
    dplyr::filter(FDR < label_fdr_cutoff) %>%
    dplyr::arrange(FDR) %>%
    utils::head(label_top)
  ggplot2::ggplot(tt, ggplot2::aes(x = logFC, y = -log10(FDR), colour = sig)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggrepel::geom_text_repel(data = labels,
                             ggplot2::aes(label = feature_id),
                             size = 2.5, max.overlaps = 15) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey60",
                                            "TRUE" = "firebrick")) +
    ggplot2::labs(title = paste(title_tag, "(RUVg k =", k, ")"),
                  x = "log2 FC (RUVg-normalized)",
                  y = "-log10 FDR",
                  colour = paste0("FDR < ", sig_cutoff)) +
    theme_ng()
}

## Print a kable of the top `n` features below `fdr_cutoff`, or a placeholder
## line if none pass. Intended for results='asis' chunks.
ruv_top_table <- function(tt, n = 30, fdr_cutoff = 0.05) {
  stopifnot(all(c("FDR", "feature_id") %in% colnames(tt)))
  top <- tt %>%
    dplyr::filter(FDR < fdr_cutoff) %>%
    dplyr::arrange(FDR) %>%
    utils::head(n)
  if (nrow(top) > 0) {
    print(knitr::kable(top, digits = 4, format = "html"))
  } else {
    cat("no features at FDR <", fdr_cutoff, "\n\n")
  }
  invisible(top)
}

## Scatter of a library-level QC metric (e.g. mismatch rate) per perturbation.
## One point per (perturbation, library) pair, with a median crossbar per
## perturbation. The metric is library-level, so differences across
## perturbations reflect which libraries each one drew cells from rather than
## per-perturbation alignment quality. Use as a batch-confound check.
##
## Required columns in `df`: pert_id, value. Optional: n_cells (sets point
## size). Override the y-axis title via `value_label`.
perturbation_qc_scatter <- function(df,
                                    value_label = "library mismatch rate, %",
                                    title_tag = NULL,
                                    point_alpha = 0.6,
                                    median_colour = "firebrick") {
  stopifnot(all(c("pert_id", "value") %in% colnames(df)))
  size_col <- if ("n_cells" %in% colnames(df)) "n_cells" else NULL
  medians <- df %>%
    dplyr::group_by(pert_id) %>%
    dplyr::summarise(value = stats::median(value, na.rm = TRUE),
                     .groups = "drop")
  pert_order <- medians %>%
    dplyr::arrange(value) %>%
    dplyr::pull(pert_id)
  df$pert_id <- factor(df$pert_id, levels = pert_order)
  medians$pert_id <- factor(medians$pert_id, levels = pert_order)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pert_id, y = value))
  if (!is.null(size_col)) {
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(size = .data[[size_col]]), width = 0.25, height = 0,
      alpha = point_alpha, colour = "grey40")
  } else {
    p <- p + ggplot2::geom_jitter(width = 0.25, height = 0,
                                  alpha = point_alpha, colour = "grey40")
  }
  p <- p +
    ggplot2::geom_crossbar(data = medians,
                           ggplot2::aes(ymin = value, ymax = value),
                           width = 0.6, colour = median_colour, fatten = 2) +
    ggplot2::labs(x = NULL, y = value_label, title = title_tag,
                  size = if (is.null(size_col)) NULL else "cells per (pert, lib)") +
    theme_ng()
  p
}

## Per-perturbation logFC distribution as a horizontal boxplot per perturbation.
## `df` needs pert_id and logFC. `facet_col` (optional) facets by
## polyA-competence or any other stratifier. Use above the per-perturbation
## volcano stack.
class_logfc_box <- function(df,
                            facet_col = NULL,
                            title_tag = NULL,
                            highlight_fdr = 0.05) {
  stopifnot(all(c("pert_id", "logFC") %in% colnames(df)))
  df$pert_id <- factor(df$pert_id,
                       levels = df %>%
                         dplyr::group_by(pert_id) %>%
                         dplyr::summarise(m = stats::median(logFC, na.rm = TRUE),
                                          .groups = "drop") %>%
                         dplyr::arrange(m) %>%
                         dplyr::pull(pert_id))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pert_id, y = logFC)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_boxplot(outlier.size = 0.6, fill = "grey90") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "log2 FC", title = title_tag) +
    theme_ng()
  if (!is.null(facet_col) && facet_col %in% colnames(df)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)))
  }
  if ("FDR" %in% colnames(df)) {
    sig <- df %>% dplyr::filter(.data$FDR < highlight_fdr)
    if (nrow(sig) > 0) {
      p <- p + ggplot2::geom_jitter(data = sig, width = 0.2, height = 0,
                                    size = 0.8, alpha = 0.7,
                                    colour = "firebrick")
    }
  }
  p
}

## Heatmap of per-perturbation x per-class (or per-family) logFC. Rows are
## perturbations, columns are class_id (or family_id). Shows per-class
## response patterns that volcanoes hide at coarser granularities.
class_logfc_heatmap <- function(df, group_col = "class_id",
                                title_tag = NULL,
                                fc_clip = 3) {
  stopifnot(all(c("pert_id", "logFC", group_col) %in% colnames(df)))
  agg <- df %>%
    dplyr::group_by(.data$pert_id, .data[[group_col]]) %>%
    dplyr::summarise(median_lfc = stats::median(.data$logFC, na.rm = TRUE),
                     .groups = "drop")
  agg$median_lfc <- pmax(pmin(agg$median_lfc, fc_clip), -fc_clip)
  ggplot2::ggplot(agg,
                  ggplot2::aes(x = .data[[group_col]], y = .data$pert_id,
                               fill = .data$median_lfc)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white",
                                  high = "firebrick", midpoint = 0,
                                  limits = c(-fc_clip, fc_clip)) +
    ggplot2::labs(x = group_col, y = NULL,
                  title = title_tag,
                  fill = "median log2 FC") +
    theme_ng()
}
