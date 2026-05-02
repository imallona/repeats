## Roll up STARsolo gene_id-level count matrices to family_id or class_id.
##
## STARsolo Solo.out/Gene/raw/matrix.mtx (and UniqueAndMult-EM.mtx) emits
## counts at the locus/gene_id level. RepeatMasker family and class roll-ups
## come from the per-feature_set locus_map TSV at
##   <indices_base>/indices/<genome_tag>/<feature_set>_locus_map.tsv
## with columns (no header): transcript_id, gene_id, family_id, class_id.
##
## All operations keep the input sparse (dgCMatrix in, dgCMatrix out) so they
## are cheap on cell-level matrices. The caller densifies (as.matrix) only
## when needed, e.g. before DGEList construction.
##
## Sourceable by any sc Rmd that has Matrix and dplyr attached.

load_locus_map <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  cols <- c("transcript_id", "gene_id", "family_id", "class_id")
  df <- readr::read_tsv(path, col_names = cols, show_col_types = FALSE,
                        progress = FALSE)
  df <- dplyr::distinct(df, gene_id, family_id, class_id)
  df
}

## Aggregate matrix rows by a vector of group labels aligned to rownames(m).
## Returns a Matrix sparse object (dgCMatrix) with one row per unique group.
## Rows whose group label is NA are dropped before aggregation.
aggregate_rows_by_group <- function(m, group) {
  stopifnot(length(group) == nrow(m))
  ok <- !is.na(group)
  if (!all(ok)) {
    m <- m[ok, , drop = FALSE]
    group <- group[ok]
  }
  group_f <- factor(group)
  indicator <- Matrix::sparseMatrix(
    i = as.integer(group_f),
    j = seq_along(group_f),
    x = 1,
    dims = c(nlevels(group_f), length(group_f)),
    dimnames = list(levels(group_f), NULL))
  indicator %*% m
}

## Roll a gene_id-level matrix to a coarser granularity using the locus map.
## Input and output are sparse (dgCMatrix); the gene_id pass-through returns
## the input unchanged. NULL is returned when the locus map is unavailable
## and the granularity is not gene_id.
rollup_matrix <- function(m, locus_map, granularity) {
  if (granularity == "gene_id") return(m)
  if (is.null(locus_map)) return(NULL)
  stopifnot(granularity %in% c("family_id", "class_id"))
  group_for <- setNames(locus_map[[granularity]], locus_map$gene_id)
  group <- unname(group_for[rownames(m)])
  aggregate_rows_by_group(m, group)
}

## Convenience: roll a list of pseudobulk objects (each with $counts and
## $meta) to a target granularity. $counts is rolled while preserving its
## sparsity class; $meta is unchanged. NULL entries pass through unmodified.
## Names of `pb_list` are expected to follow the "<mode>_<feature_set>"
## convention used by the repogle and gse230647 reports so the right locus
## map can be picked from `locus_maps[[feature_set]]`.
rollup_pseudobulks <- function(pb_list, locus_maps, granularity) {
  out <- vector("list", length(pb_list))
  names(out) <- names(pb_list)
  for (key in names(pb_list)) {
    pb <- pb_list[[key]]
    if (is.null(pb)) next
    fs <- sub("^[^_]+_", "", key)
    lmap <- locus_maps[[fs]]
    rolled <- rollup_matrix(pb$counts, lmap, granularity)
    if (is.null(rolled)) next
    out[[key]] <- list(counts = rolled, meta = pb$meta)
  }
  out
}
