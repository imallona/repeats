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
