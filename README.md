# Repeat element quantification in bulk and single cells

## Status

Active development (2024-2026). See [TODO](TODO.md) for open items.
Previous exploratory version: [v0.1](https://github.com/imallona/repeats/releases/tag/v0.1).

## Acknowledgements

We are extremely grateful to SNF for their funding in ca. 2020.

## Overview

The pipeline:

1. Prepares references (genome, Ensembl gene annotation, RepeatMasker repeat annotation).
2. Simulates single-cell RNA-seq reads from repeat loci (SmartSeq2 or 10x Chromium).
3. Quantifies repeat expression with STARsolo, Kallisto (bustools), Alevin (Salmon),
   and Bowtie2 (pseudo-genome approach).
4. Evaluates each aligner against simulation ground truth at three aggregation levels:
   locus (gene_id), repeat family (family_id), and repeat class (class_id).
5. Produces an HTML evaluation report with accuracy, precision/recall, and compute
   resource plots.

For method details see [workflow/methods.md](workflow/methods.md).
Workflow diagrams (mermaid, renders on GitHub) are in [docs/diagrams.md](docs/diagrams.md).

## Requirements

- Snakemake 8+
- conda / mamba (for --use-conda)

All other dependencies are installed automatically via per-rule conda environments.

Exact package pins (platform-locked explicit exports) are in `workflow/envs/explicit/`:

| File | Environment | Key packages |
|---|---|---|
| `repeats_star.txt` | repeats_star | STAR, samtools, bedtools, pigz, gffread |
| `repeats_kallisto.txt` | repeats_kallisto | kallisto, bustools |
| `repeats_alevin.txt` | repeats_alevin | salmon |
| `repeats_bowtie2.txt` | repeats_bowtie2 | bowtie2, samtools, subread (featureCounts) |
| `repeats_umi_tools.txt` | repeats_umi_tools | umi_tools, samtools >= 1.12, pysam |
| `repeats_evaluation.txt` | repeats_evaluation | python, scipy, pysam |
| `rmarkdown.txt` | rmarkdown | R, rmarkdown, ggplot2, patchwork |

## Running

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd workflow
snakemake --configfile configs/simulation_smartseq2.yaml --use-conda --cores 10
snakemake --configfile configs/simulation_chromium.yaml  --use-conda --cores 10
```

The report is written to `{base}/evaluation/evaluation_report.html` as defined in
the config file.

## Configuration

Two reference configs are provided under `workflow/configs/`:

| File | Technology | Cells | Chr subset |
|---|---|---|---|
| `simulation_smartseq2.yaml` | SmartSeq2 | 20 | chr10 |
| `simulation_chromium.yaml`  | 10x Chromium | 50 | chr10 |

Key parameters:

- `base`: output directory for all results.
- `indices_base`: directory for shared aligner indices (can be shared across runs).
- `feature_sets`: which repeat subsets to quantify (`repeats`, `genic_repeats`, `intergenic_repeats`).
- `granularities`: aggregation levels (`gene_id`, `family_id`, `class_id`).
- `aligner_params.{aligner}.multimapper_modes`: `unique` (best hit only) or `multi` (EM/all-alignments).

## Implementation notes

Chromium normalization (kallisto, alevin) uses sparse accumulators to avoid
allocating a dense cells x features matrix.  Only non-zero (cell_index, count)
pairs are stored per feature group, keeping memory proportional to expressed
pairs rather than O(features x cells).

Bowtie2 Chromium counting streams the CB-tagged BAM once and accumulates
per-(barcode, locus) counts directly, without splitting into per-cell BAM files.

See workflow/methods.md for full details.

## Contact

izaskun dot mallona dot work at gmail.com

## License

GNU General Public License (GPL)
