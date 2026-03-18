# Methods

## Simulation

Reads are generated directly from repeat element loci defined in a RepeatMasker GTF annotation (TEtranscripts format). The genome FASTA is streamed one chromosome at a time; only chromosomes containing at least one sampled locus are loaded into memory. This avoids holding the full genome sequence in memory simultaneously.

For each simulated cell a random subset of repeat loci is drawn without replacement. The number of expressed loci per cell varies around a user-specified mean by sampling uniformly in the range [0.7x, 1.3x] of that mean. The count assigned to each expressed locus is drawn from a geometric distribution with mean 5, truncated at 50. The geometric distribution was chosen because it produces the heavy-tailed, low-count profile typical of scRNA-seq data without requiring external dependencies. Each read is a subsequence sampled uniformly at random from within the repeat locus sequence. Substitution errors are introduced independently at each base with probability 0.001 to approximate typical Illumina error rates.

In SmartSeq2 mode each cell is written to a separate gzipped FASTQ file. A STARsolo-compatible manifest is written alongside.

In Chromium mode cell barcodes are drawn at random (16 bp by default, matching 10x v3). For each read a unique UMI is sampled and placed in read 1 together with the cell barcode. Read 2 carries the cDNA sequence. The barcodes file can be passed as a whitelist to any aligner that accepts one.

The ground truth is stored as a long-format TSV with one row per expressed (cell, repeat) pair, including the repeat family_id and class_id columns from the GTF. Only loci for which the genome contained a usable sequence (less than 10% N bases, at least 50 bp) appear in the ground truth.

## Reference building

Repeat sequences for pseudo-alignment tools (kallisto, salmon) are extracted from the genome using the repeat GTF coordinates and bedtools getfasta. Each repeat locus is treated as its own transcript and gene, so the transcript-to-gene map is the identity. This is the simplest assumption and avoids collapsing loci from the same family, which would obscure locus-level quantification differences between aligners.

When `filter_genic` is enabled, three annotation subsets are produced using bedtools intersect:
- `genic_repeats`: repeat loci overlapping at least one gene body
- `intergenic_repeats`: repeat loci with no overlap to any gene body
- `repeat_harboring_genes`: gene records overlapping at least one repeat locus

These subsets are saved as GTF files for downstream filtering and stratified evaluation.

## Aligners

**STARsolo** aligns reads to the full genome and counts using the repeat GTF. For SmartSeq2 data it uses soloType SmartSeq and the per-cell manifest. For Chromium data it uses soloType CB_UMI_Simple with barcode and UMI length parameters matching the simulation or chemistry. Two multimapping strategies are available: `unique` (only reads mapping to a single locus contribute to counts) and `multimapping` (reads mapping to multiple loci are distributed using the EM algorithm). SmartSeq2 mode only supports `unique` because STARsolo does not implement EM for that soloType.

**kallisto** aligns reads to a pseudo-genome of repeat sequences using pseudoalignment. For SmartSeq2 data, `kallisto quant` is run in single-end mode for each cell individually. For Chromium data, `kallisto bus` is run on the combined R1/R2 files; the resulting BUS file is sorted and counted with bustools, collapsing transcript-level counts to gene-level counts using the identity transcript-to-gene map.

**salmon alevin** uses selective alignment against the repeat sequence index. For SmartSeq2 data, `salmon quant` is run per cell in single-end mode. For Chromium data, `salmon alevin` is run with the chemistry flag corresponding to the barcode and UMI lengths.

**bowtie2** aligns reads to the repeat sequence index. For SmartSeq2 data, one BAM is produced per cell. featureCounts is then run on groups of cells (chunks) to limit peak memory use. The chunk size is configurable. The per-chunk count tables are merged row-wise into a single feature x cell matrix. For Chromium data, read 2 (the cDNA read) is aligned to the repeat index; cell barcode and UMI tags from read 1 are attached to the BAM using umi_tools extract. UMI deduplication is performed with umi_tools dedup before featureCounts.

## Normalization

Each aligner's native output is converted to a common feature x cell TSV where rows are repeat gene_ids and columns are cell identifiers. Only rows with at least one non-zero count are written. This format is aligner-agnostic and is the only input consumed by the evaluation step.

## Evaluation

The evaluation compares each aligner's count matrix against the simulation ground truth. Metrics are computed at three levels.

**Global**: Pearson and Spearman correlation between ground truth and observed counts across all (cell, repeat) pairs. Log1p-RMSE (root mean squared error on log-transformed counts with pseudocount 1). Detection metrics using the set of (cell, repeat) pairs with count greater than zero: precision, recall, F1, Jaccard index, and specificity.

**Per cell**: Pearson and Spearman correlation computed separately for each cell across all repeat loci.

**Per repeat class**: All global metrics stratified by the class_id field from the GTF (typically LINE, SINE, LTR, DNA transposon, and satellite). Because different repeat classes differ substantially in copy number and sequence conservation, per-class metrics reveal which repeat types are accurately quantified by each aligner.

**Computational metrics**: Wall time, CPU time, peak resident set size, and I/O volume are recorded for each rule using Snakemake's built-in benchmark directive and collated into a summary table alongside the accuracy metrics.

A final summary TSV concatenates the global metrics row from every aligner for direct comparison.
