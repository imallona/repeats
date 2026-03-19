# Methods

## Overview

The pipeline quantifies repeat element expression from single-cell RNA-seq data
(SmartSeq2 or 10x Chromium) using four aligners/quantifiers and evaluates them
against simulation ground truth.

## Reference preparation

Genome sequence (GRCh38), gene annotation (Ensembl GTF), and repeat annotation
(RepeatMasker, UCSC RMSK format converted to GTF via `gtfToGenePred`/`makeTEgtf`)
are downloaded or provided by the user.

All internal files use bare chromosome names without the `chr` prefix (Ensembl convention).
UCSC-sourced files with `chr`-prefixed chromosome names are stripped at decompression time
(`decompress_genome`, `decompress_repeats_gtf`, `decompress_genes_gtf` rules in
`reference.snmk`), ensuring consistency across all downstream steps.

Optional genic/intergenic repeat subsets are produced with `bedtools intersect` when
`filter_genic: true` in the config.

## Simulation

SmartSeq2 reads are simulated with `simulate_reads.py` (in-house script).  For each
simulated cell, a set of expressed repeat loci is drawn from the RepeatMasker annotation;
reads are sampled from the corresponding genomic sequence using ART or a simple
uniform-coverage model.  The exact per-cell per-locus read counts are recorded in
`ground_truth.tsv` and used as the reference for all evaluation metrics.

## Aligners and quantifiers

### STARsolo

STAR 2.7+ is run in `SmartSeq` (SS2) or `CB_UMI_Simple` (Chromium) solo mode with a
custom `sjdbGTFfile` that selects either the repeat or gene annotation.  Multimapper
handling is controlled by `--soloMultiMappers` (`Unique` or `EM`); SmartSeq mode only
supports `Unique`.  Output is a 10x-style sparse matrix in `Solo.out/Gene/raw/`.

### Kallisto

For SmartSeq2, `kallisto quant` is run independently per cell in bulk mode against a
repeat-sequence pseudo-transcriptome (`indices/{genome_tag}/repeats.fa`).  Per-cell
`abundance.tsv` files are merged and optionally aggregated by transcript-to-gene mapping
at the requested granularity level (see **Quantification granularity** below).

### Alevin (Salmon)

For SmartSeq2, `salmon quant` is run per cell in quasi-mapping mode.  Outputs are
merged in the normalization step.  The `--minAssignedFrags 1` flag is used to prevent
Salmon from exiting with an error when too few reads map (common with small simulations).

### Bowtie2 — pseudo-genome approach

Standard alignment to the full genome followed by featureCounts is problematic for
repeat elements because:
- Repeat locus names in the BAM (derived from the pseudo-FASTA header) do not match
  chromosome names in the genomic GTF.
- featureCounts assigns reads based on genomic coordinates, conflating loci that share
  the same repeat gene_id.

Instead, Bowtie2 is indexed against a **pseudo-genome FASTA** extracted from the genome
using the repeat GTF (`extract_repeat_fasta` rule).  Each repeat locus becomes its own
reference sequence, named `transcript_id::chrom:start-end(strand)`.

Reads are aligned with `bowtie2 -a` (report all alignments) and per-locus read counts
are obtained from `samtools idxstats`, which counts alignments per reference sequence
without coordinate lookup.  The `count_pseudo_genome.py` script:

1. Parses `samtools idxstats` output per cell BAM.
2. Looks up each transcript_id in the locus map (`repeats_locus_map.tsv`) to resolve
   gene_id, family_id, and class_id.
3. Aggregates counts at the requested granularity and writes a feature × cell TSV.

This approach ensures exact per-locus counting and seamless multi-granularity output.

## Quantification granularity

All repeat-aware aligners (Bowtie2, Kallisto, Alevin for SS2) produce count matrices at
configurable aggregation levels:

| Level | RepeatMasker field | Example value |
|---|---|---|
| `locus` | `transcript_id` | `AluSz6::chr10:1234-1534(+)` |
| `gene_id` | `gene_id` | `AluSz6` |
| `family_id` | `family_id` | `Alu` |
| `class_id` | `class_id` | `SINE` |

The `granularities` config key specifies which levels to compute (default: `[family_id]`).

Ground truth is aggregated to the same level before evaluation, so metrics are always
computed at the same resolution as the counts.

## Normalization

Each aligner's native output is converted to a common feature × cell TSV by a
per-aligner `normalize_*.py` script called from `normalize.snmk`.  All downstream
evaluation rules consume this format.

## Evaluation

`evaluate.py` compares the normalized count matrix against the simulation ground truth:

- **Global metrics** (one row per aligner × granularity): Pearson r, Spearman r,
  log1p RMSE, precision, recall, F1, Jaccard index, specificity.
- **Per-cell metrics**: Pearson r and Spearman r for each simulated cell.
- **Per-repeat-class metrics**: all global metrics stratified by RepeatMasker `class_id`.

The `aggregate_global_metrics` rule concatenates all per-aligner global metric TSVs
into a single `summary_global_metrics.tsv`.

An HTML report (`evaluation_report.html`) is rendered by `evaluation_report.Rmd` using
`rmarkdown::render` and presents all metrics as interactive plots (ggplot2 + patchwork).

## Conda environments

| Environment | File | Key packages |
|---|---|---|
| `star` | `envs/star.yaml` | STAR, pigz, bedtools, gffread, samtools |
| `kallisto` | `envs/kallisto.yaml` | kallisto, bustools |
| `alevin` | `envs/alevin.yaml` | salmon |
| `bowtie2` | `envs/bowtie2.yaml` | bowtie2, samtools |
| `evaluation` | `envs/evaluation.yaml` | python, scipy |
| `rmarkdown` | `envs/rmarkdown.yaml` | R, rmarkdown, ggplot2, patchwork |
