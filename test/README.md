# Tests

This directory contains the test suite for the repeats pipeline. Tests are
run with pytest and cover the Python scripts in workflow/scripts/ as well as
the correctness of the evaluation logic end-to-end.

See also: the main [README.md](../README.md) for pipeline usage.

## Running the tests

All commands below assume the repository root as the working directory unless
stated otherwise.

### Unit and integration tests

Install dependencies:

```
pip install pytest pytest-cov scipy
```

Run unit and integration tests:

```
cd /path/to/repeats
pytest test/unit/ test/integration/ -v
```

Run with coverage:

```
cd /path/to/repeats
pytest test/unit/ test/integration/ --cov=workflow/scripts --cov-report=term-missing
```

### Snakemake dry-run tests

These require snakemake in PATH and are tagged with the `workflow` marker.
Run via pytest (from the repo root):

```
cd /path/to/repeats
pytest test/workflow/test_snakemake_dryrun.py -v -m workflow
```

Or run the dry-runs directly with snakemake. The working directory must be
`workflow/` for the production Snakefile, and paths to the test Snakefile and
configs are given relative to that directory:

```
cd /path/to/repeats/workflow

snakemake --configfile configs/simulation_smartseq2.yaml --dry-run --quiet

snakemake --configfile configs/simulation_chromium.yaml --dry-run --quiet

snakemake -s ../test/workflow/Snakefile_test \
  --configfile ../test/workflow/configs/test_negative_control.yaml \
  --dry-run --quiet
```

### Negative control workflow (full run)

Requires reference data accessible from the runner and snakemake with conda
support. Run from the `workflow/` directory:

```
cd /path/to/repeats/workflow

snakemake -s ../test/workflow/Snakefile_test \
  --configfile ../test/workflow/configs/test_negative_control.yaml \
  --use-conda \
  --cores 4 \
  --conda-frontend mamba
```


## Directory layout

```
test/
  unit/                    unit tests for individual script functions
  integration/             end-to-end tests using synthetic input files
  workflow/                snakemake dry-run tests and the negative control workflow
    Snakefile_test         test workflow that reuses the production snmk modules
    configs/               config for the negative control run
    envs/                  conda environment for test workflow rules
```


## Unit tests (test/unit/)

Each file targets one script in workflow/scripts/.

test_evaluate.py covers:
- pearson_r, spearman_r: known values, edge cases (too few points, constant vector)
- log1p_rmse: perfect recovery gives 0.0, symmetry, empty input
- detection_metrics: all-correct, no-overlap, partial overlap, specificity
- build_aligned_vectors: zero-filling for missing features, cell ordering
- load_ground_truth: all four granularities (locus, gene_id, family_id, class_id)
  and the valid_locus_ids filter that was added to fix the genic/intergenic
  partitioning bug
- load_count_matrix: feature tracking including zero-count rows
- compute_metrics_for_subset: perfect observer and null observer

test_simulate_reads.py covers:
- reverse_complement: known sequences, N passthrough, involution property
- parse_gtf_attribute: quoted and unquoted values, missing key
- extract_repeat_sequence: N content filter, strand handling, length filter
- sample_subseq: output length, N padding for short sequences
- sample_count_geometric: always >= 1, respects max_count, approximate mean
- build_cell_plan and build_locus_to_cells: output structure
- write_ground_truth: TSV format and row count
- parse_gtf_repeats_by_chrom: minimal GTF parsing, length filter, chrom filter

test_build_rmsk_gtf.py covers:
- make_gtf_attributes: format and dup index incrementing
- convert_rmsk_to_gtf: basic output, dup index, length filter, chrom filter,
  1-based GTF coordinates, gzip output

test_parse_gtf_t2g.py covers:
- parse_attrs: quoted values, missing keys, empty string
- main() with sys.argv patching: 2-col and 4-col output, correct mapping,
  deduplication, feature-type filtering


## Integration tests (test/integration/)

test_feature_set_consistency.py documents and tests a bug that caused lower
evaluation metrics for genic_repeats and intergenic_repeats at gene_id
granularity compared to the full repeats set, even for a hypothetically
perfect observer.

Root cause: evaluate.py was aggregating ground truth counts across all loci
for a given gene_id (including both genic and intergenic copies), then
filtering at the gene_id level. For a gene_id like AluSz6 that has copies
in both contexts, the genic truth count included intergenic loci, inflating
truth relative to observed and degrading Pearson correlation and recall.

The fix adds locus-level filtering before aggregation (see valid_locus_ids
in load_ground_truth). The tests in this file verify:

1. gene_ids overlap between genic and intergenic locus maps (structural
   precondition for the bug, expected in real data)
2. full repeats at gene_id with a perfect observer gives metrics of 1.0
   (baseline sanity check)
3. genic_repeats at locus granularity with a perfect observer gives 1.0
   (locus level is unaffected by the partition issue)
4. genic_repeats at gene_id with a perfect genic observer now gives 1.0
   (regression test for the fix; would fail if the fix were reverted)
5. intergenic_repeats at gene_id with a perfect observer gives 1.0
6. subset metrics are not lower than full-set metrics for a perfect observer
   (monotonicity property)

test_metamorphic.py tests properties that must hold across input transformations
without needing exact expected values:

- scale invariance: multiplying observed counts by a constant does not change
  Pearson or Spearman
- scaling does change log1p_rmse
- monotone recall: adding true positives never decreases recall
- noise degrades metrics: progressively noisier observers give lower Pearson
- perfect beats null: a perfect observer strictly outperforms a zero observer
  on all metrics
- granularity aggregation consistency: a perfect locus-level count matrix
  summed to family_id also gives metrics of 1.0
- random counts give low correlation (below 0.7 for 50-feature random data)


## Workflow tests (test/workflow/)

test_snakemake_dryrun.py verifies that snakemake --dry-run succeeds for both
production configs (simulation_smartseq2.yaml and simulation_chromium.yaml)
and for the test negative control config. These tests are skipped if snakemake
is not in the PATH and are marked with the workflow marker.

Snakefile_test is a separate Snakemake workflow that includes all production
snmk modules and adds the negative control rules:

simulate_from_genes runs simulate_reads.py with the Ensembl gene GTF instead
of the repeat GTF. Reads come from gene body regions rather than repeat
elements. This gives a ground truth over gene features.

check_negative_control_recall runs evaluate.py comparing that gene-body
ground truth against the repeat quantification output and asserts that recall
is below the threshold in testing.negative_control_max_recall (default 0.10).
A low recall means the pipeline correctly does not attribute gene-body reads
to repeat elements.

The conda environment for these rules is defined in
test/workflow/envs/test_evaluation.yaml (python, scipy, pandas, numpy).


## GitHub Actions

.github/workflows/tests.yml defines three jobs:

unit-and-integration runs on every push and pull request. It uses plain
Python and does not need any bioinformatics tools installed.

snakemake-dryrun runs on every push and pull request. It installs snakemake
via micromamba and runs dry-run checks on both production configs and the test
Snakefile. Snakemake installs conda dependencies itself when --use-conda is
passed, but dry-runs do not trigger conda installs.

negative-control-run is only triggered manually (workflow_dispatch) or when
the repository variable ENABLE_FULL_INTEGRATION is set to true. It runs the
full negative control workflow with --use-conda and requires reference data
to be accessible from the runner.
