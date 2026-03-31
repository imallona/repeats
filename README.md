# Repetitive element quantification in bulk and single-cell RNA-seq

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

For method details see [docs/methods.md](docs/methods.md).
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

### Via Makefile (recommended)

A `Makefile` at the project root orchestrates all configs and renders reports.

```
# Run everything (all simulations + noise sweep + reports)
make

# Individual targets
make simulation_smartseq2      # SmartSeq2 base simulation
make simulation_chromium       # Chromium base simulation
make noise_smartseq2           # SmartSeq2 noise sweep (0%, 1%, 5%, 10% mutation rate)
make noise_chromium            # Chromium noise sweep
make report_noise_smartseq2    # Render SmartSeq2 noise sweep HTML report
make report_noise_chromium     # Render Chromium noise sweep HTML report
make help                      # List all targets
```

Tune parallelism with `make CORES=20`. The Makefile activates the `snakemake` conda
environment automatically.

### Manually

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd workflow
snakemake --configfile configs/simulation_smartseq2.yaml --use-conda --cores 10 --rerun-triggers mtime
snakemake --configfile configs/simulation_chromium.yaml  --use-conda --cores 10 --rerun-triggers mtime
```

The per-run evaluation report is written to `{base}/evaluation/evaluation_report.html`.

## Configuration

Simulation configs are under `workflow/configs/`:

| File | Technology | Cells | Expressed loci/cell | Mutation rate |
|---|---|---|---|---|
| `simulation_smartseq2.yaml` | SmartSeq2 | 500 | 1000 | 0.1% |
| `simulation_chromium.yaml`  | 10x Chromium | 500 | 1000 | 0.1% |
| `simulation_smartseq2_noise_0pct.yaml` | SmartSeq2 | 500 | 1000 | 0% |
| `simulation_smartseq2_noise_1pct.yaml` | SmartSeq2 | 500 | 1000 | 1% |
| `simulation_smartseq2_noise_5pct.yaml` | SmartSeq2 | 500 | 1000 | 5% |
| `simulation_smartseq2_noise_10pct.yaml` | SmartSeq2 | 500 | 1000 | 10% |
| `simulation_chromium_noise_0pct.yaml` | 10x Chromium | 500 | 1000 | 0% |
| `simulation_chromium_noise_1pct.yaml` | 10x Chromium | 500 | 1000 | 1% |
| `simulation_chromium_noise_5pct.yaml` | 10x Chromium | 500 | 1000 | 5% |
| `simulation_chromium_noise_10pct.yaml` | 10x Chromium | 500 | 1000 | 10% |

Unused real-data configs have been moved to `workflow/configs/old/`.

### Bulk RNA-seq reanalysis

Two publicly available bulk RNA-seq datasets are reanalysed for repeat element quantification:

| GSE | Description | Design | Snakefile | Config |
|---|---|---|---|---|
| GSE126543 | NeuN+ neuronal nuclei, human frontal cortex (FTD/ALS), TDP-43 nuclear positive vs negative, 7 donors | paired by donor | Snakefile_gse126543 | configs/gse126543_bulk.yaml |
| GSE230647 | iPSC-derived neural cultures, TDP-43 OE (HA-TDP43 DOX on vs off, n=4) and KD (TDP-43 shRNA vs NT shRNA, n=4) | two independent experiments | Snakefile_gse230647 | configs/gse230647_bulk.yaml |

Both pipelines download paired-end FASTQs from SRA, align with STAR, kallisto, and salmon,
quantify repeats at gene_id and family_id granularities, and render an edgeR differential
expression report.

Key parameters:

- `base`: run-specific output directory.
- `indices_base`: shared directory for aligner indices and decompressed references.
  All noise-sweep configs share the same `indices_base` so indices are built once.
- `simulation.mutation_rate`: per-base substitution rate applied to simulated reads.
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

See [docs/methods.md](docs/methods.md) for full details.

## Testing

Unit tests, integration tests, and Snakemake dry-run tests are in the
`test/` directory. See [test/README.md](test/README.md) for details on
design, coverage, and how to run the tests.

## Contact

izaskun dot mallona dot work at gmail.com

## License

GNU General Public License (GPL)
