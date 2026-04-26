# Methods

## Overview

A single Snakefile drives five pipeline modes selected by `pipeline_type` in the config.

The simulation pipeline generates reads from repeat loci (SmartSeq2 or 10x Chromium), aligns them with four tools, and evaluates quantification accuracy against ground truth.

The bulk pipeline downloads FASTQs from SRA via fasterq-dump. Reads are aligned with STAR, kallisto, and salmon. edgeR produces differential expression reports at gene_id and family_id granularity. Read layout is controlled by `real_data.library_layout`: `paired` (R1 + R2) or `single` (R1 only). Paired-end runs use `bulk_paired.snmk`. Single-end runs use `bulk_single.snmk`, which passes `--single --fragment-length --sd` to kallisto and `--unmatedReads` to salmon.

The sc pipeline downloads 10x Chromium FASTQs from SRA with `fasterq-dump --include-technical`. File numbering follows fasterq-dump conventions. File `_1` is the index read and is discarded. The CB+UMI and cDNA reads are selected via `sc_cbumi_index` and `sc_cdna_index` (defaults 2 and 3, matching 10x Chromium v2/v3). STARsolo (CB_UMI_Simple) and kallisto|bustools quantify genic and intergenic repeat annotations.

The de_simulation pipeline is a count-level power benchmark. It loads pre-computed gene and repeat count TSVs from a real run, fits NB parameters, and resamples count matrices over a 2D grid of gene and repeat library size scalers, planting a known fold change at a fixed feature set. It scores no-norm, TMM, and RUVg(k=1..k_max) for power and FPR. See [docs/de_simulations.md](de_simulations.md) for the full description.

## Reference preparation

Genome sequence (GRCh38), gene annotation (Ensembl GTF), and repeat annotation (RepeatMasker, UCSC RMSK format converted to GTF via makeTEgtf) are downloaded or provided by the user. All shared reference files and aligner indices sit under a single `indices_base` directory. Multiple runs with different simulation parameters (e.g. noise sweeps) reuse the same indices without rebuilding.

Internal files use bare chromosome names without the `chr` prefix (Ensembl convention).  UCSC-sourced files with `chr`-prefixed chromosome names are stripped at decompression time.

When `filter_genic` is enabled, three annotation subsets are produced using bedtools intersect:

- `genic_repeats`: repeat loci overlapping at least one gene body
- `intergenic_repeats`: repeat loci with no overlap to any gene body
- `repeat_harboring_genes`: gene records overlapping at least one repeat locus

These subsets are saved as GTF files for downstream filtering and stratified evaluation.

## Simulation

Reads are generated directly from repeat element loci defined in a RepeatMasker GTF annotation. The genome FASTA is streamed one chromosome at a time. Only chromosomes with at least one sampled locus are loaded into memory, which avoids holding the full genome in RAM.

For each simulated cell a random subset of repeat loci is drawn without replacement. The number of expressed loci per cell varies around a user-specified mean (`n_expressed_per_cell`, default 1000). The per-cell count is sampled uniformly from [0.7x, 1.3x] of that mean. The count assigned to each expressed locus comes from a geometric distribution with mean 5, truncated at 50. The geometric distribution gives the heavy-tailed, low-count profile typical of scRNA-seq without external dependencies.

Each read is a subsequence sampled uniformly at random from the repeat locus sequence. Substitution errors are introduced independently at each base with probability `mutation_rate` (default 0.001). A noise sweep runs four mutation rates (0%, 1%, 5%, 10%) to test how aligners hold up under sequencing errors.

In SmartSeq2 mode each cell is written to a separate gzipped FASTQ file alongside a STARsolo-compatible manifest.

In Chromium mode cell barcodes are drawn at random (16 bp, matching 10x v3). A single paired-end FASTQ pair is generated: read 1 carries the cell barcode and UMI (12 bp), read 2 carries the cDNA.

Ground truth is a long-format TSV with one row per expressed (cell, repeat) pair. `family_id` and `class_id` columns come from the GTF. Only loci whose genome sequence is usable (< 10% N bases, >= 50 bp) appear in the ground truth.

## Quantification granularity

All repeat-aware aligners produce count matrices at configurable aggregation levels:

| Level     | RepeatMasker field | Example value |
|-----------|--------------------|---------------|
| locus     | transcript_id      | AluSz6_dup1   |
| gene_id   | gene_id            | AluSz6        |
| family_id | family_id          | Alu           |
| class_id  | class_id           | SINE          |

The `granularities` config key specifies which levels to compute.  Ground truth is aggregated to the same level before evaluation so metrics are always computed at matching resolution.

## Aligners and quantifiers

### STARsolo

STAR 2.7+ runs in SmartSeq (`soloType SmartSeq`) or Chromium (`soloType CB_UMI_Simple`) mode. A custom `sjdbGTFfile` selects the repeat or gene annotation. Multimapper handling is controlled by `--soloMultiMappers`: `Unique` (only uniquely mapping reads contribute) or `EM` (multimapping reads distributed by the EM algorithm). Both modes are evaluated for all technologies.

### Kallisto

For SmartSeq2 (simulation) and bulk single-end, `kallisto quant` runs per sample against a repeat-sequence pseudo-transcriptome. Single-end bulk passes `--single --fragment-length --sd`. Fragment length and standard deviation are set via `real_data.fragment_length` and `real_data.fragment_sd` (defaults 200 and 20). Per-sample `abundance.tsv` files are merged and aggregated at the requested granularity.

For bulk paired-end, `kallisto quant` is run without `--single`, passing both R1 and R2.

For 10x Chromium, `kallisto bus` runs on the combined R1/R2 files. The resulting BUS file is sorted and corrected with bustools. Per-cell count matrices are produced with `bustools count`.

### Alevin (Salmon)

For SmartSeq2 and bulk paired-end, `salmon quant` runs per sample in quasi-mapping mode. `-1` and `-2` pass the two mates. `--minAssignedFrags 1` stops Salmon from exiting when few reads map.

For bulk single-end, `salmon quant` is run with `--unmatedReads` instead of `-1`/`-2`.

For 10x Chromium, `salmon alevin` is run in Chromium v3 mode using the repeat pseudo-transcriptome index.  Output is converted by `normalize_alevin_chromium.py`.

### Bowtie2: pseudo-genome approach

Bowtie2 is indexed against a pseudo-genome FASTA where each repeat locus is its own reference sequence, named `transcript_id::chrom:start-end(strand)`. Reads align with `bowtie2 -a` (all alignments). Per-locus counts come from `samtools idxstats`, which counts alignments per reference without coordinate lookup.

For SmartSeq2, one BAM is produced per cell; featureCounts is then run on configurable cell chunks to limit peak memory.

For 10x Chromium, read 2 (cDNA) aligns to the repeat index. Cell barcode and UMI tags from read 1 are attached via `umi_tools extract`. UMI deduplication runs via `umi_tools dedup` before counting.

## Normalization

Each aligner's native output is converted to a common feature x cell TSV (rows = repeat features, columns = cell identifiers) by a per-aligner `normalize_*.py` script.  Only non-zero rows are written.  This format is the sole input to the evaluation step.

For 10x Chromium, kallisto and alevin normalization uses sparse accumulators. Only non-zero `{cell_index: count}` pairs are stored per feature group. Memory stays proportional to expressed pairs rather than O(features x cells).

The bowtie2 Chromium counting script streams the CB-tagged deduplicated BAM once via `samtools view`, accumulating `(barcode, feature)` counts in a sparse dict without per-cell BAM splitting.

## Evaluation

`evaluate.py` compares each aligner's count matrix against the simulation ground truth at three levels:

Global (one row per aligner x granularity): Pearson r, Spearman r, log1p RMSE, precision, recall, F1, Jaccard index, specificity.

Per cell: Pearson r and Spearman r computed separately for each cell across all repeat loci.

Per repeat class: All global metrics stratified by `class_id` (LINE, SINE, LTR, DNA transposon, satellite, etc.), revealing which repeat types are accurately quantified.

Compute resource metrics (wall time, CPU time, peak RSS, I/O) are read from Snakemake benchmark files and included in the summary table.

The `aggregate_global_metrics` rule concatenates all per-aligner global metric TSVs into a single `summary_global_metrics.tsv`.

An HTML report (`evaluation_report.html`) is rendered by `evaluation_report.Rmd` using `rmarkdown::render` with ggplot2 + patchwork. The report also writes paper-ready CSV tables to `{eval_dir}/paper_csv/`: `global_metrics_long.csv` (tidy long format), `global_metrics_wide.csv`, `aligner_ranking_per_metric.csv` (rank per feature_set x granularity x metric, direction-aware), `best_aligner_per_metric.csv` (best, second, delta), `per_cell_metrics_summary.csv` (n, mean, median, sd, q25, q75 per aligner x granularity x feature_set), `per_class_metrics.csv`, and `resources_wide.csv`.

A separate noise sweep report (`noise_sweep_report.Rmd`) loads `summary_global_metrics.tsv` from each noise-level run and plots metric degradation as a function of mutation rate. It writes paper-ready CSV tables to `{csv_outdir}/paper_csv/` (by default the directory of the HTML output): `noise_metrics_long.csv`, `noise_metrics_wide.csv` (columns per mutation rate), `noise_degradation_slopes.csv` (linear fit `value ~ mutation_rate` with slope, intercept, R^2, n_points per aligner x feature_set x granularity x metric), `noise_robustness_ranking.csv` (ranks aligners by shallowest slope for accuracy metrics), `per_cell_noise_summary.csv`, and `per_class_noise_metrics.csv`.

## Conda environments

| Environment | File                 | Key packages                              |
|-------------|----------------------|-------------------------------------------|
| star        | envs/star.yaml       | STAR, pigz, bedtools, gffread, samtools   |
| kallisto    | envs/kallisto.yaml   | kallisto, bustools                        |
| alevin      | envs/alevin.yaml     | salmon                                    |
| bowtie2     | envs/bowtie2.yaml    | bowtie2, samtools                         |
| umi_tools   | envs/umi_tools.yaml  | umi_tools, samtools (>= 1.12), pysam       |
| evaluation  | envs/evaluation.yaml | python, scipy                             |
| rmarkdown   | envs/rmarkdown.yaml  | R, rmarkdown, ggplot2, patchwork          |
| edger       | envs/edger.yaml      | R, edgeR, rmarkdown, ggplot2              |
| ruvseq      | paper/envs/ruvseq.yaml | R, edgeR, RUVSeq, EDASeq, rmarkdown     |

## Paper analyses

The `paper/` directory holds a separate Snakefile for dataset-specific reanalyses used in the manuscript. It is kept out of `workflow/Snakefile` for two reasons. First, paper analyses are cross-dataset, with a single config referencing outputs of multiple method runs. Second, paper rules should not retrigger method jobs via Snakemake provenance tracking.

Entry point: `paper/Snakefile`, driven by `paper/configs/tdp43.yaml`. Counts are consumed from the method-pipeline output directory declared in each dataset entry (`counts_dir` for bulk, `starsolo_base` for sc), so the relevant method pipeline must be run first for each dataset. Dataset keys use last author plus study label plus assay kind; each entry also records the source GEO accession (`gse`) and the method-pipeline config (`source_config`) that produced its inputs.

Paper datasets currently wired in `paper/configs/tdp43.yaml`:

| Paper dataset key | GSE | Assay | Source method config |
|---|---|---|---|
| `lee_nuclear_tdp_retention_bulk` | GSE126543 | bulk | `workflow/configs/gse126543_bulk.yaml` |
| `polymenidou_tdp_kd_oe_cultures_bulk` | GSE230647 | bulk | `workflow/configs/gse230647_bulk.yaml` |
| `polymenidou_tdp_oe_cultures_scrnaseq` | GSE230647 | scRNA-seq | `workflow/configs/gse230647_sc.yaml` |

### TDP-43 bulk differential expression

The `render_ruv_bulk` rule runs a per-dataset RUVg-normalized edgeR analysis on STAR unique-count matrices. A dataset wildcard lets `polymenidou_tdp_kd_oe_cultures_bulk` and `lee_nuclear_tdp_retention_bulk` share a single rule definition. Inputs are all repeats, genic repeats, and intergenic repeats at gene_id and family_id granularity, plus the gene matrix. STAR unique counts outperform kallisto and salmon on repeats in the simulation benchmark. The kallisto and salmon matrices are not consumed by the paper rule. The baseline bulk edgeR reports are kept; the paper rule adds a second, RUVg-normalized analysis rather than replacing them.

Per-dataset procedure, encoded in `paper/scripts/ruv_{dataset}_report.Rmd`:

1. Load the gene matrix and the six repeat matrices. Filter the gene matrix to rows with count >= 10 in at least half the samples; repeat matrices are not filtered at this stage.
2. Fit a naive edgeR quasi-likelihood model on the filtered gene matrix with TMM normalization. The design is two-group for GSE126543 and four-group for GSE230647, with F-test across the three non-reference coefficients. Rank genes by p-value. Take the top `n_controls` with the largest p-values (default 5000) as empirical controls. These are the genes least associated with the contrasts of interest.
3. Run RUVg with k = 1..4 on the gene matrix. Score each k by the mean-absolute-deviation of the relative log expression distribution. The full k = 0..4 table is written to `ruv_rle_diagnostic.tsv` for audit. The operating k is set via `k_selected` (default 2), not picked automatically from the RLE. RUVg produces per-sample W factors that do not depend on feature granularity or featureset, so the same W is reused across all repeat contrasts.
4. For each repeat featureset and contrast: build a DGEList on the raw repeat counts, compute TMM factors, and fit `~ W + group` with `glmQLFit` and `glmQLFTest`. The ranked TSV contains feature_id, logFC, logCPM, F, PValue, and FDR. GSE126543 drops the two unsorted samples and runs a single tdp43pos vs tdp43neg contrast. GSE230647 runs two contrasts: TDP-43 knockdown vs non-targeting and TDP-43 overexpression vs DOX-off control.

Outputs per dataset are written under `results/paper/tdp43/{dataset}/`. The report lives alongside one ranked DE TSV per featureset x contrast named `ruv_k{k}_{featureset}_{label}.tsv`. RDS snapshots (`ruv_state.rds`, `ruv_naive_gene_de.rds`, `ruv_repeat_de.rds`) are saved next to `ruv_empirical_controls.txt`, `ruv_rle_diagnostic.tsv`, and `ruv_summary.tsv`.

Shared helper functions used by both bulk Rmds live in `workflow/scripts/ruv_common.R`: `load_counts_tsv`, `rle_score`, `rank_empirical_controls`, `fit_ruv_range`, `pick_W_matrix`, and `run_ruv_repeat_de`. The same helpers are reused by the single-cell pseudo-bulk RUVg analysis so that control selection, RUVg fitting, and edgeR contrast code are identical across bulk and sc.

A parallel gene-library-size repeat DE is run alongside RUV in the same Rmd. Library-size factors come from TMM on the filtered gene matrix and are transferred to each repeat DGEList; no RUV covariates enter the design. TMM computed directly on the repeat matrix is avoided because TDP-43 perturbation moves many repeat subfamilies in concert, which violates TMM's "most features not DE" assumption and biases repeat logFC toward zero. The gene-library-size path serves as a transparent baseline against which the RUV-normalized analysis is compared in the integrative Rmd. Outputs are `gene_lib_repeat_de.rds` (list of topTags by featureset x contrast) and per-(featureset, contrast) TSVs named `gene_lib_{featureset}_{label}.tsv`. The helper is `run_gene_lib_repeat_de` in `workflow/scripts/gene_lib_size_common.R`, sourced alongside `ruv_common.R`.

### Single-cell pseudo-bulk RUVg (polymenidou_tdp_oe_cultures_scrnaseq)

The `render_ruv_sc` rule dispatches on the `kind: sc` field in the dataset config and renders `paper/scripts/ruv_polymenidou_tdp_oe_cultures_scrnaseq_report.Rmd`. Inputs are the STARsolo raw matrices produced by the sc method pipeline (`starsolo_base`), the GEO cell metadata file, and a sample metadata TSV. No counts are re-derived from raw reads. Pseudo-bulking is per (sample, cluster) pair: each GEO-annotated cluster in each of the three samples contributes its own pseudo-bulk column, dropping any pair with fewer than `min_cells_per_col` cells (default 10). With three samples and 15+ annotated clusters the gene matrix ends up with tens of columns rather than 6, which gives the design comfortable residual degrees of freedom. The binary test contrast stays cluster 12 versus any non-12 cluster, encoded as `group = factor(is_cluster12)` with `TRUE` for cluster-12 columns. Sample of origin enters the design through `sample_block = factor(sample_id)`. The pseudo-bulk gene matrix is filtered to genes with count >= 10 in at least half the columns. A naive `~ sample_block + group` model ranks genes by p-value; the top `n_controls` (default 3000) with the largest p-values become empirical controls. RUVg is fit for k = 1..4. The operating k is fixed via `k_selected` (default 1). Repeat DE uses `~ sample_block + W + group` with the `groupTRUE` coefficient as the test. The same `workflow/scripts/ruv_common.R` helpers used by the bulk reports are reused here. Outputs mirror the bulk layout: ranked TSVs `ruv_k{k}_{mode}_{featureset}_cluster12_vs_other.tsv`, RDS snapshots (`ruv_state.rds`, `ruv_naive_gene_de.rds`, `ruv_repeat_de.rds`, `pseudobulk_counts.rds`, `geo_meta_annotated.rds`), `ruv_empirical_controls_{mode}.txt`, `ruv_rle_diagnostic_{mode}.tsv`, `ruv_summary.tsv`, and CSV exports under `<report_outdir>/csv/` (`cluster_sizes_per_sample.csv`, `group_sizes_per_sample.csv`, `pseudobulk_col_meta.csv`, `pseudobulk_gene_lib_sizes.csv`, `ruv_W_factors.csv`, `ruv_rle_diagnostic.csv`, `ruv_de_long.csv`, `ruv_summary.csv`). `ruv_W_factors.csv` and `pseudobulk_gene_lib_sizes.csv` carry `sample_id`, `cluster_num`, `cluster_name`, `n_cells`, and `is_cluster12` columns so downstream code can map pseudo-bulk columns back to their (sample, cluster) origin.

RUVg is paper-side only, matching the bulk convention: `workflow/scripts/gse230647_sc_report.Rmd` runs gene-library-normalized edgeR plus QC only, and all RUVg-normalized DE for this dataset lives in the paper-side Rmd. Control selection, RUVg fitting, and edgeR contrast code are executed under the paper Snakefile and consume only the method pipeline's STARsolo outputs.

The sc Rmd adds two analyses on top of the RUV cluster-12-vs-other contrast. First, a parallel gene-library-size repeat DE is run for the same cluster-12-vs-other contrast, with library-size factors transferred from the pseudo-bulk gene matrix via `run_gene_lib_repeat_de`. Output is `gene_lib_repeat_de.rds` and per-(mode, featureset) TSVs named `gene_lib_{mode}_{featureset}_cluster12_vs_other.tsv`. Second, a per-cluster-vs-rest sweep runs the same gene-library-size analysis for every cluster with pseudo-bulk support in at least two samples. The sweep exists to empirically test the cluster-12-as-OE-counterpart hypothesis; the integrative Rmd correlates each cluster's repeat logFC profile against Polymenidou bulk OE and ranks clusters by that correlation. Output is `per_cluster_gene_lib_repeat_de.rds`. RUV is not run per cluster because RUV-g empirical controls are selected as genes least associated with the contrast under test, which means a per-cluster RUV analysis would need fresh control ranking per contrast and would mix noise-suppression with the question under study; the gene-library-size path is preferred here.

Run:

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd paper
snakemake --use-conda --cores 4 --configfile configs/tdp43.yaml
```

### Cross-dataset TDP-43 synthesis

The paper aims at a single TDP-43 view combining three inputs, all rendered by `paper/Snakefile`:

- `polymenidou_tdp_kd_oe_cultures_bulk` (GSE230647 bulk): knockdown (`tdp43_kd vs nt_kd_ctrl`) and overexpression (`tdp43_oe vs oe_ctrl`).
- `lee_nuclear_tdp_retention_bulk` (GSE126543 bulk): sorted nuclei, `tdp43pos vs tdp43neg`.
- `polymenidou_tdp_oe_cultures_scrnaseq` (GSE230647 scRNA-seq): cluster-12 vs other cells.

Cluster 12 in the sc dataset matches the TDP-43 overexpression arm of the bulk datasets. Cells in cluster 12 express the HA-TDP43 transgene under DOX-on conditions, while the rest of the DOX-on sc population shows little to no transgene. So the cluster-12 contrast is a neuron-resolved analogue of the bulk `tdp43_oe vs oe_ctrl` comparison in `polymenidou_tdp_kd_oe_cultures_bulk`. The other DOX-on and DOX-off sc samples, and the KD/NT bulk arms, are treated separately.

The synthesis rule `render_integrative_tdp43` in `paper/Snakefile` consumes the three per-dataset `ruv_rds/` directories and produces `results/paper/tdp43/integrative_tdp43/integrative_tdp43_report.html` plus ranked CSVs and RDS snapshots under `integrative_rds/`. Two concordance matrices are reported: one on RUV-normalized repeat DE tables (`integrative_repeat_concordance_summary.csv`), one on gene-library-size repeat DE tables (`integrative_gene_lib_repeat_concordance_summary.csv`). The two are also written side by side in `integrative_ruv_vs_gene_lib_concordance_side_by_side.csv`. Comparing the two separates cross-study signal robust to normalization choice from signal that depends on RUV factor estimation.

Pairs reported across both matrices: `lee_tdppos vs poly_kd`, `lee_tdppos vs poly_oe`, `lee_tdppos vs sc_cluster12`, `poly_kd vs poly_oe`, `poly_oe vs sc_cluster12`, `poly_kd vs sc_cluster12`. Each pair is reported separately for genic and intergenic repeats. Columns: n (features inner-joined on feature_id), pearson, spearman, n_sig_both (significant in both contrasts at `fdr_cutoff`), n_concordant_sig (same sign), n_discordant_sig (opposite sign). Biology-stratified scatter panels and correlation tables follow: polyA competence, downstream A-run (polyA priming proxy), family age. The biology annotation is loaded from any per-dataset `biology_annotation.rds` and joined on `feature_id = gene_id`.

`polya_class` and `a_run_class` are orthogonal axes. `polya_class` is class-level RNA biology (LINE, SINE, LTR, Retroposon generate polyA-tailed RNA during their life cycle; other classes do not). `a_run_class` is locus-level genomic DNA context (fraction of loci with a downstream A-run >= 12 bp) and proxies the risk of internal oligo-dT priming. Signal concentrated in polyA-competent + low-a_run subfamilies is the cleanest; polyA-competent + medium/high-a_run subfamilies (L1HS, SVA_B-F, young LTRs) are biologically interesting but cannot be disambiguated from priming artifact using RNA-seq alone. `a_run_label` in `workflow/scripts/biology_checks.R` excludes `Simple_repeat`, `Low_complexity`, and `Satellite` classes (returns NA) because their repeat units are themselves A-containing so the flag fires tautologically and carries no information about internal-priming of real polyadenylated transcripts; the exclusion keeps the a_run facets focused on classes where the priming question is meaningful.

A per-cluster ranking panel tests the cluster-12 choice empirically. Per-cluster gene-library-size repeat DE tables from the sc Rmd are inner-joined against the Polymenidou bulk OE gene-library-size repeat DE table, and Pearson and Spearman correlations of logFC are computed per cluster. Output is `integrative_per_cluster_vs_bulk_oe_ranking.csv` plus a bar chart of Pearson r per cluster, sorted; cluster 12 is highlighted. If cluster 12 does not top the ranking, the choice of sc counterpart of bulk OE deserves re-examination before further inference.
- Bulk per-contrast ranked TSVs at `results/paper/tdp43/{lee_nuclear_tdp_retention_bulk,polymenidou_tdp_kd_oe_cultures_bulk}/ruv_k{k}_{featureset}_{contrast}.tsv`.
- sc per-(mode, featureset) ranked TSVs at `results/paper/tdp43/polymenidou_tdp_oe_cultures_scrnaseq/ruv_k{k}_{mode}_{featureset}_cluster12_vs_other.tsv` and the long-format CSV `csv/ruv_de_long.csv`.

