# workflow/modules/

Each `.snmk` file is a self-contained Snakemake module included by the main `Snakefile`.
Modules are conditionally included based on `config['mode']` and `config['aligners']`.

## Module index

| `download_references.snmk` | always | Download genome FASTA, genes GTF, and RepeatMasker annotation from Ensembl/UCSC |
| `reference.snmk` | always | Build indices (STAR, Kallisto, Salmon, Bowtie2), extract repeat/gene FASTAs, genic/intergenic filtering |
| `data_acquisition.snmk` | always | SRA download rules for real-data mode |
| `simulations.snmk` | `mode: simulation` | Simulate SmartSeq2 per-cell FASTQs or Chromium R1/R2 from repeat loci |
| `starsolo.snmk` | `starsolo` in aligners | STARsolo alignment (SmartSeq2 manifest or Chromium CB_UMI_Simple) |
| `kallisto.snmk` | `kallisto` in aligners | Kallisto pseudoalignment (bulk per-cell for SS2, kallisto bus + bustools for Chromium) |
| `alevin.snmk` | `alevin` in aligners | Salmon quant (bulk per-cell for SS2) or salmon alevin (Chromium, sketch mode for simulation) |
| `bowtie2.snmk` | `bowtie2` in aligners | Bowtie2 alignment + featureCounts with optional cell/feature chunking |
| `normalize.snmk` | always | Convert each aligner's native output to a common feature × cell TSV |
| `evaluation.snmk` | `mode: simulation` | Compare normalized counts against simulation ground truth; produce metric tables |

## Chromosome subsetting and `genome_tag`

All reference indices (STAR, Kallisto, Salmon, Bowtie2) are built under
`{base}/references/{genome_tag}/` where `genome_tag` is derived from
`config.reference.chromosomes`:

| Config value | `genome_tag` | path (index) |
|---|---|---|
| *(empty)* | `all` | `references/all/star/` |
| `["chr10"]` | `chr10` | `references/chr10/star/` |
| `["chr10", "chrX"]` | `chr10_chrX` | `references/chr10_chrX/star/` |

## Adding a new aligner

1. Create `modules/<aligner>.snmk` with rules gated on `sim_technology`.
2. Add a normalization block in `normalize.snmk` producing `counts/<aligner>_{feature_set}.tsv`.
3. The `evaluation.snmk` `evaluate_single_mode_aligner` rule will pick it up automatically
   once the aligner name is added to `config['aligners']`.
