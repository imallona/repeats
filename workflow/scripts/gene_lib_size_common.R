## Gene-library-size helpers for repeat DE analyses.
## Sourced alongside ruv_common.R from the paper-side Rmds.
##
## Motivation: TMM on the repeat count matrix assumes that most repeat
## subfamilies are not differentially expressed, which is not true under
## TDP-43 perturbation (coordinated repeat derepression). TMM norm factors
## derived from the gene matrix sidestep that assumption, since genes are the
## bulk of the library and the "most features not DE" assumption holds for
## them. This file provides the gene-derived library-size path; the RUV path
## lives in ruv_common.R.
##
## Expects edgeR to be attached by the caller.

## Fit edgeR on a repeat count matrix using library-size and norm factors
## derived from a gene count matrix (no RUV covariates). Sample ordering
## between the two matrices must match; the stopifnot guards against silent
## misalignment.
run_gene_lib_repeat_de <- function(repeat_counts, gene_counts, design, test) {
  stopifnot(qr(design)$rank == ncol(design),
            identical(colnames(repeat_counts), colnames(gene_counts)))
  gene_dge <- edgeR::DGEList(counts = gene_counts)
  gene_dge <- edgeR::calcNormFactors(gene_dge, method = "TMM")
  dge <- edgeR::DGEList(counts = repeat_counts,
                        lib.size = gene_dge$samples$lib.size,
                        norm.factors = gene_dge$samples$norm.factors)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- do.call(edgeR::glmQLFTest, c(list(fit), test))
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tt$feature_id <- rownames(tt)
  tt[order(tt$PValue), c("feature_id","logFC","logCPM","F","PValue","FDR")]
}
