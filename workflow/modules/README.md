# workflow/modules/

Each `.snmk` file is a self-contained Snakemake module included by the main `Snakefile`. Modules are conditionally included based on `config['pipeline_type']` and `config['aligners']`.

## Module index

| module | active when | purpose |
|---|---|---|
| `download_references.snmk` | always | Download genome FASTA, genes GTF, and RepeatMasker annotation from Ensembl/UCSC |
| `reference.snmk` | always | Build indices (STAR, Kallisto, Salmon, Bowtie2), extract repeat/gene FASTAs, genic/intergenic filtering |
| `data_acquisition.snmk` | `pipeline_type: simulation` | Legacy SRA download rule reused by the simulation pipeline (colon_cancer_cell_lines layout) |
| `download_sra.snmk` | `pipeline_type: bulk` or `sc` | SRA download for real-data bulk and sc pipelines; honours `real_data.sc_cbumi_index` / `sc_cdna_index` for 10x Chromium splitting |
| `simulations.snmk` | `pipeline_type: simulation` | Simulate SmartSeq2 per-cell FASTQs or Chromium R1/R2 from repeat loci |
| `starsolo.snmk` | `starsolo` in aligners, simulation pipeline | STARsolo alignment (SmartSeq2 manifest or Chromium CB_UMI_Simple) |
| `starsolo_sc.snmk` | `starsolo` in aligners, sc pipeline | STARsolo CB_UMI_Simple for real 10x Chromium data |
| `kallisto.snmk` | `kallisto` in aligners, simulation pipeline | Kallisto pseudoalignment (per-cell for SS2, `kallisto bus` + bustools for Chromium) |
| `kallisto_sc.snmk` | `kallisto` in aligners, sc pipeline | `kallisto bus` + bustools for real 10x Chromium data |
| `alevin.snmk` | `alevin` in aligners, simulation pipeline | Salmon quant (per-cell for SS2) or salmon alevin (Chromium) |
| `bowtie2.snmk` | `bowtie2` in aligners, simulation pipeline | Bowtie2 alignment to repeat pseudo-genome + `samtools idxstats` counting |
| `bulk_paired.snmk` | `pipeline_type: bulk`, `library_layout: paired` | STAR, kallisto, salmon for paired-end bulk RNA-seq |
| `bulk_single.snmk` | `pipeline_type: bulk`, `library_layout: single` | STAR, kallisto (--single), salmon (--unmatedReads) for single-end bulk RNA-seq |
| `normalize.snmk` | always | Convert each aligner's native output to a common feature x cell TSV |
| `evaluation.snmk` | `pipeline_type: simulation` | Compare normalized counts against simulation ground truth; produce metric tables and HTML report |
| `noise_report.snmk` | `pipeline_type: noise_report` | Render the cross-mutation-rate noise sweep HTML report from existing evaluation dirs |
| `de_simulations.snmk` | `pipeline_type: de_simulation` | Count-level power benchmark for repeat DE: resamples gene/repeat counts over a 2D library-size grid, plants a known fold change, scores no-norm/TMM/RUVg, and renders the HTML report |

## Chromosome subsetting and `genome_tag`

All indices are built under `{base}/indices/{genome_tag}/` where `genome_tag` is derived from `config.reference.chromosomes` (the user-supplied value is used as-is for labelling):

| config value (user-provided) | `genome_tag` | example path |
|---|---|---|
| *(empty)* | `all` | `indices/all/star/` |
| `["chr10"]` | `chr10` | `indices/chr10/star/` |
| `["chr10", "chrX"]` | `chr10_chrX` | `indices/chr10_chrX/star/` |

Internally, all GTF and FASTA files use bare chromosome names without the `chr` prefix (Ensembl convention), regardless of the source annotation.  UCSC `chr`-prefixed files are stripped at the `decompress_*` rules.

## Quantification granularity

Bowtie2, Kallisto (SmartSeq2), and Alevin (SmartSeq2) all support multi-level aggregation of repeat counts via the `granularities` config key:

```yaml
granularities:
  - locus      # individual repeat locus (transcript_id)
  - gene_id    # repeat gene / family member
  - family_id  # RepeatMasker family (e.g. L1HS, AluSz)
  - class_id   # RepeatMasker class (e.g. LINE, SINE, LTR)
```

Default when not set: `[family_id]`.

The aggregation is performed at the normalization step using a locus map (`indices/{genome_tag}/repeats_locus_map.tsv`) built from the repeats GTF by `reference.snmk::build_repeat_locus_map`.

## Bowtie2 pseudo-genome approach

Rather than aligning to the full genome and using featureCounts (which requires chromosome names in the BAM to match those in the annotation GTF), Bowtie2 is indexed against a **pseudo-genome FASTA** of extracted repeat sequences (`indices/{genome_tag}/repeats.fa`).  Each repeat locus becomes its own reference "chromosome" named `transcript_id::chrom:start-end(strand)`.

Read counts per locus are obtained with `samtools idxstats` and aggregated by the `count_pseudo_genome.py` script.  This approach has two key advantages:

1. Counting is exact per-locus without GTF/chromosome name reconciliation.
2. Multi-granularity outputs (locus -> family -> class) are produced in a single pass.

## Technology and library layout

Single-cell runs set `simulation.technology` or `real_data.technology`. Bulk runs set `real_data.library_layout`. The values are:

Single-cell (`technology`):

| value | reads | quantification |
|---|---|---|
| `smartseq2` | one FASTQ per cell | kallisto `quant`, salmon `quant`, STARsolo manifest |
| `chromium` | R1 (CB+UMI) + R2 (cDNA) | STARsolo CB_UMI_Simple, kallisto `bus`, salmon `alevin` |

Bulk (`library_layout` under `real_data`):

| value | reads | quantification |
|---|---|---|
| `single` | R1 only per sample | kallisto `quant`, salmon `quant`, STAR |
| `paired` | R1 + R2 per sample | kallisto `quant`, salmon `quant`, STAR |

SmartSeq2 and bulk single both produce one FASTQ per sample or cell and share the same per-file normalization scripts. The difference is that SmartSeq2 runs the simulation pipeline with per-cell ground truth, while bulk single processes real SRA downloads without cell barcodes.

## Adding a new aligner

1. Create `modules/<aligner>.snmk` with rules gated on `sim_technology`.
2. Add a normalization block in `normalize.snmk` producing `counts/<aligner>_{feature_set}[_{granularity}].tsv`.
3. Add the aligner name to `granular_aligners` in `evaluation.snmk` if it supports granularity, or add a dedicated `evaluate_<aligner>` rule otherwise.
