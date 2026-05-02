# paper/

TDP-43 paper analyses. Consumes count matrices from the bulk method pipeline (`workflow/Snakefile`) and STARsolo raw matrices from the sc method pipeline. Every dataset gets a standalone RUVg-normalized edgeR Rmd rendered by `paper/Snakefile`.

Dataset keys are set in `configs/tdp43.yaml`. Each entry records the GEO accession and the source method-pipeline config so results are traceable back to the pipeline that produced them. See `docs/methods.md` for the full procedure.

## Layout

```
paper/
  Snakefile              # rule render_ruv_bulk and render_ruv_sc dispatched by kind
  configs/tdp43.yaml     # dataset registry
  envs/ruvseq.yaml       # R + edgeR + RUVSeq + EDASeq + ggrepel
  scripts/               # one Rmd per dataset
```

Outputs land under `../results/paper/tdp43/{dataset}/` with an HTML report, RDS snapshots under `ruv_rds/`, and paper-ready CSVs under `ruv_rds/csv/`.

## Running

Always activate the snakemake env first:

```
source ~/miniconda3/bin/activate
conda activate snakemake
cd paper
```

### All datasets

```
snakemake --use-conda --cores 4 --configfile configs/tdp43.yaml
```

This only works on a machine that holds both the bulk count TSVs and the sc STARsolo outputs. In practice bulk inputs live on `barbara` and sc inputs live on `vm`, so target individual datasets instead.

### Polymenidou TDP-43 OE cultures scRNA-seq (`polymenidou_tdp_oe_cultures_scrnaseq`)

Runs on `vm` under `/mnt/src/repeats`. Requires the sc method pipeline to have finished so that `../results/gse230647_sc/starsolo/` and the GEO cell metadata are in place.

```
snakemake --use-conda --cores 4 \
  ../results/paper/tdp43/polymenidou_tdp_oe_cultures_scrnaseq/ruv_polymenidou_tdp_oe_cultures_scrnaseq_report.html
```

Each (sample, GEO cluster) pair becomes its own pseudo-bulk column (pairs with fewer than 10 cells are dropped). The contrast is cluster 12 versus any non-12 cluster under `~ sample_block + W + group`. Output HTML includes RLE diagnostics, volcano plots, and top-table listings; the numeric results land as TSVs plus tidy CSVs for paper figures.

### Polymenidou TDP-43 KD/OE cultures bulk (`polymenidou_tdp_kd_oe_cultures_bulk`)

Runs on `barbara`. Requires the bulk method pipeline outputs under `../results/gse230647_bulk/counts/`.

```
snakemake --use-conda --cores 4 \
  ../results/paper/tdp43/polymenidou_tdp_kd_oe_cultures_bulk/ruv_polymenidou_tdp_kd_oe_cultures_bulk_report.html
```

### Lee nuclear TDP-43 retention bulk (`lee_nuclear_tdp_retention_bulk`)

Runs on `barbara`. Requires the bulk method pipeline outputs under `../results/gse126543_bulk/counts/`.

```
snakemake --use-conda --cores 4 \
  ../results/paper/tdp43/lee_nuclear_tdp_retention_bulk/ruv_lee_nuclear_tdp_retention_bulk_report.html
```

## Useful logs and re-runs

Per-dataset logs: `../results/paper/tdp43/logs/render_ruv_{dataset}.log`.

Force a single report to re-render (e.g., after editing its Rmd) without re-running siblings: delete the target HTML and rerun the same `snakemake ...` command, or add `--forcerun render_ruv_sc` / `--forcerun render_ruv_bulk`.
## repogle_2022 use case

A second use case sits alongside the TDP-43 analyses: reanalysis of the K562 essential Perturb-seq screen from Replogle et al. 2022. Documented end-to-end in [../docs/repogle_2022.md](../docs/repogle_2022.md). It uses its own Snakefile, config, and selection script so the TDP-43 entry point stays unchanged.

| Item | Value |
|---|---|
| Entry point | `Snakefile_repogle_2022` |
| Config | `configs/repogle_2022.yaml` |
| Selection script | `scripts/select_repogle_2022_cells.py` |
| Conda env | `envs/repogle_metadata.yaml` (Python, h5py, numpy) |
| Output base | `../results/paper/repogle_2022/` |

Phase 1 (implemented) downloads three figshare files (KD6 manifest, K562 essential raw bulk h5ad, K562 essential raw singlecell h5ad), verifies size and md5, and emits three slim TSVs under `data/`:

- `selected_perturbations.tsv` (top 30 by Replogle's `TE_ratio` plus their published value for validation)
- `cells_to_perturbation.tsv` (cell barcode, gemgroup, gene_transcript, sgID_AB, UMI count for ~6,200 cells)
- `selected_libraries.tsv` (per-gemgroup fastq and BAM filenames parsed from the KD6 manifest)

Run on a machine with the figshare and SRA bandwidth (the singlecell h5ad is 10 GB, fetched once):

```
source ~/miniconda3/bin/activate
conda activate snakemake
cd paper
snakemake -s Snakefile_repogle_2022 --use-conda --cores 2 \
  --configfile configs/repogle_2022.yaml
```

Phase 2 prep (SRA accession resolution) and Phase 3 (per-perturbation pseudobulk RUVg DE) are wired in `Snakefile_repogle_2022`. The remaining manual step is the workflow-side SRA download + STARsolo alignment, run separately on the compute machine via `workflow/configs/repogle_kd6_sc.yaml`.

