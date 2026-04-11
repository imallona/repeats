# Repetitive element quantification in bulk and single-cell RNA-seq

## Status

Active development (2024-2026). See [TODO](TODO.md) for open items.
Previous exploratory version: [v0.1](https://github.com/imallona/repeats/releases/tag/v0.1).

## Acknowledgements

We are extremely grateful to SNF for their funding in ca. 2020.

## Overview

A single `Snakefile` drives four pipeline modes selected by `pipeline_type` in the config:

- `simulation` - simulate reads from repeat loci (SmartSeq2 or 10x Chromium), align,
  and evaluate against ground truth
- `bulk` - download FASTQs from SRA (paired or single, set via `real_data.library_layout`),
  align with STAR/kallisto/salmon, quantify repeats, render a differential expression report
- `sc` - download 10x Chromium FASTQs from SRA, run STARsolo and kallisto|bustools,
  render a single-cell report
- `noise_report` - render a noise sweep HTML report across multiple simulation runs

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

### Memory-aware scheduling

Alignment rules and the report rule declare `resources: mem_mb` from
`resources.max_mem_mb` (override the report budget with `resources.report_mem_mb`).
Pass `--resources mem_mb=N` to cap the concurrent sum, for example
`--resources mem_mb=16000` on a 16 GB host.

### Manually

```
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
cd workflow

# simulation
snakemake --configfile configs/simulation_smartseq2.yaml --use-conda --cores 10
snakemake --configfile configs/simulation_chromium.yaml  --use-conda --cores 10

# bulk reanalysis
snakemake --configfile configs/gse230647_bulk.yaml --use-conda --cores 10
snakemake --configfile configs/gse126543_bulk.yaml --use-conda --cores 10

# single-cell reanalysis
snakemake --configfile configs/gse230647_sc.yaml --use-conda --cores 10
```

The simulation evaluation report is written to `{base}/evaluation/evaluation_report.html`.
Bulk and SC reports are written to `{base}/{report.output}`.

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

| GSE | Description | Design | Config |
|---|---|---|---|
| GSE126543 | NeuN+ neuronal nuclei, human frontal cortex (FTD/ALS), TDP-43 nuclear positive vs negative, 7 donors | paired by donor | `configs/gse126543_bulk.yaml` |
| GSE230647 | iPSC-derived neural cultures, TDP-43 OE (HA-TDP43 DOX on vs off, n=4) and KD (TDP-43 shRNA vs NT shRNA, n=4) | two independent experiments | `configs/gse230647_bulk.yaml` |

Both pipelines download paired-end FASTQs from SRA (`library_layout: paired`), align with STAR, kallisto, and salmon,
quantify repeats at gene_id and family_id granularities, and render an edgeR differential
expression report.

### Single-cell reanalysis

A single-cell analysis covers the TDP-43-HA overexpression scRNA-seq experiment from
GSE230647. Three DOX-on 10x Chromium 3' v3 samples are included (two sequencing lanes
each): GSM7230453 (2-week culture), GSM7230454 and GSM7230455 (4-week culture,
replicates A and B). The DOX-off control and unlabeled samples are excluded.

| Item | Value |
|---|---|
| Config | `configs/gse230647_sc.yaml` |
| Sample metadata | `configs/gse230647_sc_sample_metadata.tsv` |
| Report script | `scripts/gse230647_sc_report.Rmd` |
| Output | `../results/gse230647_sc/gse230647_sc_report.html` |

The pipeline downloads each SRR with `fasterq-dump --include-technical`. For 10x
Chromium, file `_2` carries the CB+UMI read and file `_3` carries the cDNA read
(configured via `sc_cbumi_index` and `sc_cdna_index` in the config). STARsolo
(CB_UMI_Simple, 10x v3 geometry) and kallisto|bustools are run against genic and
intergenic repeat annotations. The GEO supplementary file
`GSE230647_single_cell_metadata_tdp43ha.txt.gz` maps each cell barcode to its
cluster assignment (17 clusters; cluster 12 marks TDP-43-HA-expressing cells).

The report builds pseudo-bulk count matrices by aggregating STARsolo raw counts
for cluster-12 cells vs all other annotated cells, then runs edgeR (quasi-likelihood,
sample as blocking factor) to identify differentially expressed repeat elements.

Indices are shared with the bulk GSE230647 pipeline via `indices_base: ../results/shared`,
so building them once is sufficient for both analyses.

Key config parameters (shared across pipeline types):

- `pipeline_type`: `simulation`, `bulk`, `sc`, or `noise_report`.
- `base`: run-specific output directory.
- `indices_base`: shared directory for aligner indices. Multiple runs can share one index.
- `feature_sets`: which repeat subsets to quantify (`repeats`, `genic_repeats`, `intergenic_repeats`).
- `granularities`: aggregation levels (`gene_id`, `family_id`, `class_id`).

Simulation-specific:
- `simulation.mutation_rate`: per-base substitution rate applied to simulated reads.
- `aligner_params.{aligner}.multimapper_modes`: `unique` or `multi` (EM).

Bulk-specific (`real_data`):
- `library_layout`: `paired` (default) or `single`.

SC-specific (`real_data`):
- `sc_cbumi_index`: fasterq-dump file number containing CB+UMI (default `2`).
- `sc_cdna_index`: fasterq-dump file number containing cDNA (default `3`).
- `cb_length`, `umi_length`: barcode and UMI lengths. Single source of truth
  for both STARsolo (`--soloCBlen`/`--soloUMIlen`) and kallisto bus
  (`-x 0,0,cb:0,cb,cb+umi:1,0,0` custom geometry), so both aligners stay in
  sync if you change the chemistry.

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
