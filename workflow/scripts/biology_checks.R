# Shared biology helpers for repeat DGE/detection reports.
#
# Sourced by gse126543_bulk_report.Rmd, gse230647_bulk_report.Rmd, and
# gse230647_sc_report.Rmd to expose a common set of label/join helpers used
# to facet plots and stratify tables by:
#
#   1. PolyA competence of the repeat class (polya_label).
#   2. Family age bracket from RepeatMasker subfamily name or milliDiv
#      (family_age_label).
#   3. Strand concordance with overlapping host genes
#      (strand_vs_gene_label).
#   4. Internal A-tract proxy for polyA-priming risk (a_run_label).
#
# The age, strand, and A-tract columns are sourced from a per-gene_id
# annotation TSV produced by the reference-building workflow
# (repeat_biology_annotation.tsv). load_biology_annotation returns NULL
# when the file is absent so reports degrade gracefully to polyA-only
# faceting.
#
# Pure classification logic (regex on subfamily names, milliDiv buckets,
# longest-A-run scan) lives on the Python side in
# workflow/scripts/biology_annotations.py and is exercised by the pytest
# suite at test/unit/test_biology_annotations.py.

POLYA_CLASSES <- c("LINE", "SINE", "LTR", "Retroposon")

polya_label <- function(class_ids) {
  out <- ifelse(is.na(class_ids), "unknown",
                ifelse(class_ids %in% POLYA_CLASSES, "polyA-competent",
                       "not polyA-competent"))
  factor(out, levels = c("polyA-competent", "not polyA-competent", "unknown"))
}

family_age_label <- function(age_strings) {
  out <- ifelse(is.na(age_strings) | age_strings == "", "unknown", age_strings)
  factor(out, levels = c("young", "medium", "old", "unknown"))
}

strand_vs_gene_label <- function(strand_strings) {
  out <- ifelse(is.na(strand_strings) | strand_strings == "",
                "unknown", strand_strings)
  factor(out, levels = c("sense", "antisense", "both", "intergenic", "unknown"))
}

# Bins a per-subfamily fraction of loci with a downstream A-run >= the
# upstream threshold (default 12 bp) into a low/medium/high bucket.
#
# class_ids is an optional vector aligned to `fractions`. Subfamilies in
# classes listed in `exclude_classes` (default: Simple_repeat,
# Low_complexity, Satellite) are returned as NA rather than binned. The
# downstream A-run flag is intended to proxy the risk that oligo-dT
# internal-primes on genomic DNA instead of on a real RNA polyA tail, which
# only matters for classes whose RNA is supposed to be polyadenylated and
# recovered by polyA-selected library prep. Simple_repeat rows have repeat
# units that are themselves A-containing, so the flag fires tautologically
# and carries no information about internal-priming of real transcripts;
# Low_complexity and Satellite classes do not produce polyadenylated
# transcripts either and fall into the same "flag is uninformative" bin.
# Excluding them here keeps the a_run facets focused on subfamilies where
# the priming question is biologically meaningful (LINEs, SINEs, LTRs,
# Retroposons, and non-polyA-competent transcribed units like sn/scRNAs).
a_run_label <- function(fractions,
                        class_ids = NULL,
                        low_max = 0.1,
                        high_min = 0.5,
                        exclude_classes = c("Simple_repeat",
                                            "Low_complexity",
                                            "Satellite")) {
  out <- ifelse(is.na(fractions), "unknown",
                ifelse(fractions <= low_max, "low",
                       ifelse(fractions >= high_min, "high", "medium")))
  if (!is.null(class_ids)) {
    drop <- !is.na(class_ids) & class_ids %in% exclude_classes
    out[drop] <- NA_character_
  }
  factor(out, levels = c("low", "medium", "high", "unknown"))
}

load_biology_annotation <- function(path) {
  # Per-gene_id annotation TSV with columns:
  #   gene_id, family_id, class_id, n_loci,
  #   family_age, median_millidiv,
  #   frac_sense_gene, frac_antisense_gene, frac_intergenic, frac_both,
  #   frac_downstream_a_run
  # Returns NULL when the file is missing so reports can no-op.
  if (is.null(path) || !file.exists(path)) return(NULL)
  readr::read_tsv(path, show_col_types = FALSE)
}

# Joins a DE top-table with the biology annotation by repeat identifier
# (rownames of tt against the chosen key column). Adds the columns used by
# the facet helpers above. Missing rows or a NULL annotation produce NA
# columns so callers can branch on availability.
annotate_biology <- function(tt, biology_annotation, key = "gene_id") {
  if (is.null(biology_annotation)) {
    tt$class_id <- NA_character_
    tt$family_age <- NA_character_
    tt$median_millidiv <- NA_real_
    tt$strand_vs_gene <- NA_character_
    tt$frac_downstream_a_run <- NA_real_
    return(tt)
  }
  ids <- rownames(tt)
  idx <- match(ids, biology_annotation[[key]])
  tt$class_id              <- biology_annotation$class_id[idx]
  tt$family_age            <- biology_annotation$family_age[idx]
  tt$median_millidiv       <- biology_annotation$median_millidiv[idx]
  tt$strand_vs_gene        <- biology_annotation$strand_vs_gene[idx]
  tt$frac_downstream_a_run <- biology_annotation$frac_downstream_a_run[idx]
  tt
}
