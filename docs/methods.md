# Methods

## Overview

The pipeline quantifies repeat element expression from single-cell RNA-seq data
(SmartSeq2 or 10x Chromium) using four aligners/quantifiers and evaluates them
against simulation ground truth.

## Reference preparation

Genome sequence (GRCh38), gene annotation (Ensembl GTF), and repeat annotation
(RepeatMasker, UCSC RMSK format converted to GTF via makeTEgtf) are downloaded
or provided by the user.  All shared reference files and aligner indices are
placed under a single `indices_base` directory so that multiple runs with
different simulation parameters (e.g. noise sweeps) reuse the same indices
without rebuilding them.

Internal files use bare chromosome names without the `chr` prefix (Ensembl
convention).  UCSC-sourced files with `chr`-prefixed chromosome names are
stripped at decompression time.

When `filter_genic` is enabled, three annotation subsets are produced using
bedtools intersect:

- `genic_repeats`: repeat loci overlapping at least one gene body
- `intergenic_repeats`: repeat loci with no overlap to any gene body
- `repeat_harboring_genes`: gene records overlapping at least one repeat locus

These subsets are saved as GTF files for downstream filtering and stratified
evaluation.

## Simulation

Reads are generated directly from repeat element loci defined in a RepeatMasker
GTF annotation.  The genome FASTA is streamed one chromosome at a time; only
chromosomes containing at least one sampled locus are loaded into memory, which
avoids holding the full genome sequence in RAM simultaneously.

For each simulated cell a random subset of repeat loci is drawn without
replacement.  The number of expressed loci per cell varies around a
user-specified mean (`n_expressed_per_cell`, default 1000) by sampling uniformly
in the range [0.7×, 1.3×] of that mean.  The count assigned to each expressed
locus is drawn from a geometric distribution with mean 5, truncated at 50.  The
geometric distribution produces the heavy-tailed, low-count profile typical of
scRNA-seq data without requiring external dependencies.

Each read is a subsequence sampled uniformly at random from within the repeat
locus sequence.  Substitution errors are introduced independently at each base
with probability `mutation_rate` (configurable; default 0.001).  A noise sweep
is provided across four mutation rates — 0%, 1%, 5%, 10% — to assess aligner
robustness to sequencing errors.

In **SmartSeq2** mode each cell is written to a separate gzipped FASTQ file
alongside a STARsolo-compatible manifest.

In **Chromium** mode cell barcodes are drawn at random (16 bp, matching 10x v3)
and a single paired-end FASTQ pair is generated.  Read 1 carries the synthetic
cell barcode and UMI (12 bp); read 2 carries the cDNA sequence.

The ground truth is stored as a long-format TSV with one row per expressed
(cell, repeat) pair, including `family_id` and `class_id` columns from the GTF.
Only loci for which the genome contained a usable sequence (< 10% N bases,
≥ 50 bp) appear in the ground truth.

## Quantification granularity

All repeat-aware aligners produce count matrices at configurable aggregation
levels:

| Level     | RepeatMasker field | Example value |
|-----------|--------------------|---------------|
| locus     | transcript_id      | AluSz6_dup1   |
| gene_id   | gene_id            | AluSz6        |
| family_id | family_id          | Alu           |
| class_id  | class_id           | SINE          |

The `granularities` config key specifies which levels to compute.  Ground truth
is aggregated to the same level before evaluation so metrics are always computed
at matching resolution.

## Aligners and quantifiers

### STARsolo

STAR 2.7+ is run in SmartSeq (`soloType SmartSeq`) or Chromium (`soloType
CB_UMI_Simple`) mode with a custom `sjdbGTFfile` selecting the repeat or gene
annotation.  Multimapper handling is controlled by `--soloMultiMappers`:
`Unique` (only uniquely mapping reads contribute) or `EM` (multimapping reads
are distributed by the EM algorithm).  Both modes are evaluated for all
technologies.

### Kallisto

For **SmartSeq2**, `kallisto quant` is run per cell in bulk single-end mode
against a repeat-sequence pseudo-transcriptome.  Per-cell `abundance.tsv` files
are merged and aggregated at the requested granularity.

For **10x Chromium**, `kallisto bus` is run on the combined R1/R2 files; the
resulting BUS file is sorted and corrected with bustools, and per-cell count
matrices are produced with `bustools count`.

### Alevin (Salmon)

For **SmartSeq2**, `salmon quant` is run per cell in quasi-mapping mode.
`--minAssignedFrags 1` prevents Salmon from exiting when few reads map.

For **10x Chromium**, `salmon alevin` is run in Chromium v3 mode using the
repeat pseudo-transcriptome index.  Output is converted by
`normalize_alevin_chromium.py`.

### Bowtie2 — pseudo-genome approach

Bowtie2 is indexed against a pseudo-genome FASTA where each repeat locus is its
own reference sequence, named `transcript_id::chrom:start-end(strand)`.  Reads
are aligned with `bowtie2 -a` (all alignments) and per-locus counts are obtained
from `samtools idxstats`, which counts alignments per reference sequence without
coordinate lookup.

For **SmartSeq2**, one BAM is produced per cell; featureCounts is then run on
configurable cell chunks to limit peak memory.

For **10x Chromium**, read 2 (cDNA) is aligned to the repeat index and cell
barcode/UMI tags from read 1 are attached via `umi_tools extract`.  UMI
deduplication is performed with `umi_tools dedup` before counting.

## Normalization

Each aligner's native output is converted to a common feature × cell TSV (rows =
repeat features, columns = cell identifiers) by a per-aligner `normalize_*.py`
script.  Only non-zero rows are written.  This format is the sole input to the
evaluation step.

For 10x Chromium, kallisto and alevin normalization uses sparse accumulators:
only non-zero `{cell_index: count}` pairs are stored per feature group, keeping
memory proportional to expressed pairs rather than O(features × cells).

The bowtie2 Chromium counting script streams the CB-tagged deduplicated BAM once
via `samtools view`, accumulating `(barcode, feature)` counts in a sparse dict
without per-cell BAM splitting.

## Evaluation

`evaluate.py` compares each aligner's count matrix against the simulation ground
truth at three levels:

**Global** (one row per aligner × granularity): Pearson r, Spearman r, log1p
RMSE, precision, recall, F1, Jaccard index, specificity.

**Per cell**: Pearson r and Spearman r computed separately for each cell across
all repeat loci.

**Per repeat class**: All global metrics stratified by `class_id` (LINE, SINE,
LTR, DNA transposon, satellite, etc.), revealing which repeat types are
accurately quantified.

Compute resource metrics (wall time, CPU time, peak RSS, I/O) are read from
Snakemake benchmark files and included in the summary table.

The `aggregate_global_metrics` rule concatenates all per-aligner global metric
TSVs into a single `summary_global_metrics.tsv`.

An HTML report (`evaluation_report.html`) is rendered by
`evaluation_report.Rmd` using `rmarkdown::render` with ggplot2 + patchwork.

A separate **noise sweep report** (`noise_sweep_report.Rmd`) loads
`summary_global_metrics.tsv` from each noise-level run and plots metric
degradation as a function of mutation rate.

## Conda environments

| Environment | File                 | Key packages                              |
|-------------|----------------------|-------------------------------------------|
| star        | envs/star.yaml       | STAR, pigz, bedtools, gffread, samtools   |
| kallisto    | envs/kallisto.yaml   | kallisto, bustools                        |
| alevin      | envs/alevin.yaml     | salmon                                    |
| bowtie2     | envs/bowtie2.yaml    | bowtie2, samtools                         |
| umi_tools   | envs/umi_tools.yaml  | umi_tools, samtools (≥ 1.12), pysam       |
| evaluation  | envs/evaluation.yaml | python, scipy                             |
| rmarkdown   | envs/rmarkdown.yaml  | R, rmarkdown, ggplot2, patchwork          |
