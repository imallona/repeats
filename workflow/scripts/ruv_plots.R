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
