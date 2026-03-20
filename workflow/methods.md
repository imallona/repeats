# Methods

## overview

The pipeline quantifies repeat element expression from single-cell RNA-seq data
(SmartSeq2 or 10x Chromium) using four aligners/quantifiers and evaluates them
against simulation ground truth.

## Reference preparation

Genome sequence (GRCh38), gene annotation (Ensembl GTF), and repeat annotation
(RepeatMasker, UCSC RMSK format converted to GTF via makeTEgtf) are downloaded
or provided by the user.

All internal files use bare chromosome names without the chr prefix (Ensembl convention).
UCSC-sourced files with chr-prefixed chromosome names are stripped at decompression time.

Repeats can be classified in genic/intergenic on user request.

## Simulation

Reads are simulated with simulate_reads.py (in-house script) for both SmartSeq2 and
10x Chromium technologies.

For SmartSeq2, each simulated cell is an independent FASTQ file.  A set of expressed
repeat loci is drawn from the RepeatMasker annotation and reads are sampled from the
corresponding genomic sequence.  Cell barcodes are not used. This can be hacked for bulk RNA-seq analysis.

For 10x Chromium, a single paired-end FASTQ pair is generated.  Read 1 carries a
synthetic cell barcode (16 nt) and UMI (12 nt); read 2 carries the cDNA sequence.
Cell barcodes are drawn from the 10x v3 whitelist.  UMIs are unique per (cell, locus)
pair to allow UMI deduplication after alignment.

In both cases the exact per-cell per-locus read counts are recorded in
ground_truth.tsv and used as the reference for all evaluation metrics.

## Quantification granularity

All repeat-aware aligners produce count matrices at configurable aggregation levels:

| Level     | RepeatMasker field | Example value              |
|-----------|--------------------|-----------------------------|
| gene_id   | gene_id            | AluSz6                      |
| family_id | family_id          | Alu                         |
| class_id  | class_id           | SINE                        |

The granularities config key specifies which levels to compute.

Ground truth is aggregated to the same level before evaluation, so metrics are always
computed at the same resolution as the counts.


## Aligners and quantifiers

### STARsolo

STAR 2.7+ is run in SmartSeq (SS2) or CB_UMI_Simple (Chromium) solo mode with a
custom sjdbGTFfile that selects either the repeat or gene annotation.  Multimapper
handling is controlled by --soloMultiMappers (Unique or EM); SmartSeq mode only
supports Unique.  Output is a 10x-style sparse matrix in `Solo.out/Gene/raw/`.

### Kallisto

For SmartSeq2, kallisto quant is run independently per cell in bulk mode against a
repeat-sequence pseudo-transcriptome (indices/{genome_tag}/repeats.fa).  Per-cell
abundance.tsv files are merged and optionally aggregated by transcript-to-gene mapping
at the requested granularity level (see Quantification granularity below).

For 10x Chromium, kallisto bus is run in single-cell mode to produce a BUS file,
which is sorted and corrected with bustools.  Per-cell count matrices are produced with
bustools count and then converted to the common feature x cell TSV format by
normalize_kallisto_chromium.py.  The locus_map (col 0 = locus_id) is used to map
bustools locus identifiers to gene_id / family_id / class_id.

### Alevin (Salmon)

For SmartSeq2, salmon quant is run per cell in quasi-mapping mode.  Outputs are
merged in the normalization step.  The --minAssignedFrags 1 flag is used to prevent
Salmon from exiting with an error when too few reads map (common with small simulations).

For 10x Chromium, salmon alevin is run in Chromium v3 mode (--chromium) using the
repeat pseudo-transcriptome as the index.  The alevin output (quants_mat.gz,
quants_mat_cols.txt, quants_mat_rows.txt) is converted by normalize_alevin_chromium.py.
The locus_map (col 0 = locus_id) is used to map alevin feature identifiers to the
requested granularity.

### Bowtie2 - pseudo-genome approach

Standard alignment to the full genome followed by featureCounts is problematic for
repeat elements because repeat locus names in the BAM do not match chromosome names
in the genomic GTF, and featureCounts assigns reads by genomic coordinates, conflating
loci that share the same repeat gene_id.

Instead, Bowtie2 is indexed against a pseudo-genome FASTA extracted from the genome
using the repeat GTF (extract_repeat_fasta rule).  Each repeat locus becomes its own
reference sequence, named transcript_id::chrom:start-end(strand).

Reads are aligned with bowtie2 -a (report all alignments) and per-locus read counts
are obtained from samtools idxstats, which counts alignments per reference sequence
without coordinate lookup.  The count_pseudo_genome.py script:

1. Parses samtools idxstats output per cell BAM.
2. Looks up each transcript_id in the locus map (repeats_locus_map.tsv) to resolve
   gene_id, family_id, and class_id.
3. Aggregates counts at the requested granularity and writes a feature x cell TSV.

For 10x Chromium, UMI deduplication is performed with umi_tools dedup (--per-cell)
before counting.  This step uses a dedicated conda environment (umi_tools.yaml) that
bundles a recent samtools (>= 1.12) to support the --no-PG flag.

This approach keeps the granularity.

## Normalization (= harmonization)

Each aligner's native output is converted to a common feature x cell TSV by a
per-aligner normalize_*.py script called from normalize.snmk.  All downstream
evaluation rules consume this format.

## evaluation

`evaluate.py` compares the normalized count matrix against the simulation ground truth:

- Global metrics (one row per aligner x granularity): Pearson r, Spearman r,
  log1p RMSE, precision, recall, F1, Jaccard index, specificity.
- Per-cell metrics: Pearson r and Spearman r for each simulated cell.
- Per-repeat-class metrics: all global metrics stratified by RepeatMasker class_id.

Compute resource metrics (wall time, CPU time, peak RSS, I/O) are read from Snakemake
benchmark files and included in the summary table.

The aggregate_global_metrics rule concatenates all per-aligner global metric TSVs
into a single summary_global_metrics.tsv.

An HTML report (evaluation_report.html) is rendered by evaluation_report.Rmd using
`rmarkdown::render` and presents all metrics as plots (ggplot2 + patchwork).

## Conda environments

| Environment  | File                   | Key packages                          |
|--------------|------------------------|---------------------------------------|
| star         | envs/star.yaml         | STAR, pigz, bedtools, gffread, samtools |
| kallisto     | envs/kallisto.yaml     | kallisto, bustools                    |
| alevin       | envs/alevin.yaml       | salmon                                |
| bowtie2      | envs/bowtie2.yaml      | bowtie2, samtools                     |
| umi_tools    | envs/umi_tools.yaml    | umi_tools, samtools (>= 1.12), pysam  |
| evaluation   | envs/evaluation.yaml   | python, scipy                         |
| rmarkdown    | envs/rmarkdown.yaml    | R, rmarkdown, ggplot2, patchwork      |
