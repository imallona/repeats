# de_simulations

Count-level power benchmark for repeat DE.

## Question

If a repeat is truly N times more abundant in one condition than the other, what is the probability the pipeline calls it DE, and how does that probability change with gene and repeat library size?

## How

1. Load gene and repeat count TSVs from a real run.
2. Fit per-feature NB mean and dispersion with edgeR.
3. For each (gene_lib_scale, repeat_lib_scale) on the grid, simulate `n_iter` count matrices: scale the mean by the cell's scaler, add a shared latent W loading per feature, plant a fixed logFC at a fixed set of features.
4. Run no-norm, TMM, and RUVg(k=1..k_max) on each replicate. RUVg picks empirical controls from the gene matrix.
5. Score each call set at FDR <= `fdr` against the planted truth.

Planted features that fail the count filter still count as false negatives.

## Output

- `de_simulations_results.tsv`: one row per (cell, replicate, method) with power, FPR, TP/FN/FP/TN.
- `de_simulations_summary.tsv`: per-cell mean and SD across replicates.
- `de_simulations_results.rds`: full results plus fitted NB params and the planted indices.
- `de_simulations_heatmap.pdf`: power and FPR heatmaps faceted by method.
- `de_simulations_report.html`: rendered report with the heatmaps embedded plus precision, calls volume, F1, confusion outcome bars, per-replicate spread, and a sanity NB-params check. Each plot has a "How to read this" paragraph next to it. The Rmd source is `workflow/scripts/de_simulations_report.Rmd`.

## How to read the heatmap

Rows are repeat library scalers, columns are gene library scalers, color is the metric (0 to 1). Bright cells mean the method recovers the planted FC at that depth combination. Dark bottom rows mean repeat depth is the bottleneck. Dark left columns mean gene depth starves RUVg's W estimation.

## Inputs

- `gene_counts`: TSV, feature_id rownames, sample columns.
- `repeat_counts`: TSV, same sample columns. Granularity (locus, family, class) is whatever the input rows are.
- `metadata` (optional): TSV with `sample` and `condition`. Empty splits samples in halves (A then B).

The two count TSVs come from a real bulk run. This pipeline does not align reads.

## Config keys

Under `config['de_simulation']`:

- `gene_counts`, `repeat_counts`: required input paths.
- `metadata`, `sample_column`, `condition_column`: optional sample table.
- `fc`: planted fold change (default 3).
- `n_iter`: replicates per cell (default 20).
- `n_de_repeats`, `n_de_genes`: planted DE counts.
- `gene_lib_grid`, `repeat_lib_grid`: scaler vectors.
- `sigma_w`: SD of per-feature W loading.
- `k_max`: max RUVg factors.
- `n_controls`: empirical controls per replicate.
- `fdr`: FDR threshold (default 0.05).
- `seed`: master seed.

## Run

Via the Makefile from the project root:

```
make de_polymenidou_bulk CORES=N
```

This invokes the main `workflow/Snakefile` with `pipeline_type: de_simulation` (declared inside the configfile), runs the `de_simulations` rule, and then renders the HTML report in the same DAG.

Or directly:

```
snakemake --use-conda --cores N --configfile configs/de_simulations_polymenidou_bulk.yaml
```

Outputs go to `{base}/de_simulations/`. The report is rendered by `rule render_de_simulations_report` inside `workflow/modules/de_simulations.snmk` using the `rmarkdown` conda env, and depends on the four artifacts above so it is rebuilt whenever the simulation outputs change.

## Scope

This is count-level only. It does not test alignment or multi-mapper assignment. For those, see `workflow/modules/simulations.snmk` and the noise sweep configs.
