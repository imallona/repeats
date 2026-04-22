"""
Integration tests for evaluation consistency across feature sets.

Background
----------
The pipeline evaluates each (feature_set x granularity x aligner) combination.
An initially puzzling result was that genic_repeats and intergenic_repeats at
gene_id granularity showed lower metrics than the full repeats set, even for a
hypothetically perfect observer.

Root cause (now fixed in evaluate.py)
--------------------------------------
At gene_id granularity, load_ground_truth() previously aggregated ALL loci for
a given gene_id (summing genic + intergenic copies), then filtered the already-
aggregated dict by checking whether the gene_id appeared in the genic locus_map.

For a gene_id like AluSz6 that has copies in BOTH genic and intergenic regions:

  - genic locus_map contains AluSz6 (it has genic copies)
  - OLD truth['c1']['AluSz6'] = 5 (genic) + 3 (intergenic) = 8   [inflated]
  - observed from genic normalization = 5 (only genic copy)
  - -> truth > observed -> degraded Pearson / RMSE, false-negative detections

The fix moves locus_map loading to BEFORE load_ground_truth() and passes a
valid_locus_ids set so filtering happens AT LOCUS LEVEL before any aggregation.

Scenario used by these tests
-----------------------------
We construct a minimal ground truth with:
  - genic_alu_1      -> AluSz6, Alu, SINE, count=5   (GENIC)
  - intergenic_alu_1 -> AluSz6, Alu, SINE, count=3   (INTERGENIC - same gene_id!)
  - genic_alu_2      -> AluSx,  Alu, SINE, count=7   (GENIC)
  - intergenic_alu_2 -> AluSx,  Alu, SINE, count=6   (INTERGENIC - same gene_id!)
  - genic_L1_1       -> L1PA2,  L1,  LINE, count=2   (GENIC)
  - genic_mir_1      -> MIR,    MIR, SINE, count=4   (GENIC)
  - genic_L2_1       -> L2,     L2,  LINE, count=3   (GENIC)
  - intergenic_dna_1 -> TcMar,  TcMar, DNA, count=1  (INTERGENIC)
  - [TcMar has a genic reference locus genic_tcmar_ref not expressed in
     simulation, but present in the genic locus_map - this triggers the
     false-negative detection failure in the old code.]

Genic locus_map contains:
  genic_alu_1, genic_alu_2, genic_L1_1, genic_mir_1, genic_L2_1, genic_tcmar_ref

Intergenic locus_map contains:
  intergenic_alu_1, intergenic_alu_2, intergenic_dna_1

AluSz6 and AluSx appear in BOTH locus_maps (gene_id level).

With the fix:
  genic truth at gene_id: AluSz6=5, AluSx=7, L1PA2=2, MIR=4, L2=3, TcMar=0
  genic observed at gene_id (from genic normalization): AluSz6=5, AluSx=7, L1PA2=2, MIR=4, L2=3
  -> perfect match -> pearson=1.0, recall=1.0, precision=1.0
"""
import csv
import math
import os
import subprocess
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'workflow', 'scripts')
)
EVALUATE = os.path.join(SCRIPTS_DIR, 'evaluate.py')

# ---------------------------------------------------------------------------
# Shared fixture: build the scenario files
# ---------------------------------------------------------------------------

GROUND_TRUTH_ROWS = [
    # cell_id, locus_id, repeat_id (=gene_id), family_id, class_id, true_count
    ('cell_001', 'genic_alu_1',       'AluSz6', 'Alu',   'SINE', '5'),
    ('cell_001', 'intergenic_alu_1',  'AluSz6', 'Alu',   'SINE', '3'),
    ('cell_001', 'genic_alu_2',       'AluSx',  'Alu',   'SINE', '7'),
    ('cell_001', 'intergenic_alu_2',  'AluSx',  'Alu',   'SINE', '6'),
    ('cell_001', 'genic_L1_1',        'L1PA2',  'L1',    'LINE', '2'),
    ('cell_001', 'genic_mir_1',       'MIR',    'MIR',   'SINE', '4'),
    ('cell_001', 'genic_L2_1',        'L2',     'L2',    'LINE', '3'),
    ('cell_001', 'intergenic_dna_1',  'TcMar',  'TcMar', 'DNA',  '1'),
]

# Genic locus_map includes a reference locus for TcMar that was NOT expressed
# in the simulation.  In the old code, TcMar's intergenic count (1) would
# bleed into the genic truth because TcMar is in the genic locus_map.
GENIC_LOCUS_MAP = [
    ('genic_alu_1',      'AluSz6', 'Alu',   'SINE'),
    ('genic_alu_2',      'AluSx',  'Alu',   'SINE'),
    ('genic_L1_1',       'L1PA2',  'L1',    'LINE'),
    ('genic_mir_1',      'MIR',    'MIR',   'SINE'),
    ('genic_L2_1',       'L2',     'L2',    'LINE'),
    # reference locus for TcMar - exists in genome but NOT expressed
    ('genic_tcmar_ref',  'TcMar',  'TcMar', 'DNA'),
]

INTERGENIC_LOCUS_MAP = [
    ('intergenic_alu_1',  'AluSz6', 'Alu',   'SINE'),
    ('intergenic_alu_2',  'AluSx',  'Alu',   'SINE'),
    ('intergenic_dna_1',  'TcMar',  'TcMar', 'DNA'),
]

FULL_LOCUS_MAP = GENIC_LOCUS_MAP + INTERGENIC_LOCUS_MAP

# "Perfect genic observer" at gene_id level - exactly what genic normalization produces:
# only counts from genic loci, using GENIC_LOCUS_MAP.
GENIC_OBSERVED_GENE_ID = {
    'AluSz6': 5,
    'AluSx':  7,
    'L1PA2':  2,
    'MIR':    4,
    'L2':     3,
    # TcMar: 0 (its only genic reference locus was not expressed)
}

# "Perfect intergenic observer" at gene_id level
INTERGENIC_OBSERVED_GENE_ID = {
    'AluSz6': 3,
    'AluSx':  6,
    'TcMar':  1,
}

# "Perfect full observer" at gene_id level
FULL_OBSERVED_GENE_ID = {
    'AluSz6': 8,
    'AluSx':  13,
    'L1PA2':  2,
    'MIR':    4,
    'L2':     3,
    'TcMar':  1,
}

# "Perfect genic observer" at locus level
GENIC_OBSERVED_LOCUS = {
    'genic_alu_1':  5,
    'genic_alu_2':  7,
    'genic_L1_1':   2,
    'genic_mir_1':  4,
    'genic_L2_1':   3,
}


# ---------------------------------------------------------------------------
# File-writing helpers
# ---------------------------------------------------------------------------

def write_ground_truth(tmp_path):
    p = tmp_path / 'ground_truth.tsv'
    with open(p, 'w') as fh:
        fh.write('cell_id\tlocus_id\trepeat_id\tfamily_id\tclass_id\ttrue_count\n')
        for row in GROUND_TRUTH_ROWS:
            fh.write('\t'.join(row) + '\n')
    return p


def write_locus_map(tmp_path, name, rows):
    p = tmp_path / name
    with open(p, 'w') as fh:
        for row in rows:
            fh.write('\t'.join(row) + '\n')
    return p


def write_count_matrix(tmp_path, name, feature_counts, cells=('cell_001',)):
    p = tmp_path / name
    features = sorted(feature_counts)
    with open(p, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(cells) + '\n')
        for feat in features:
            vals = [str(feature_counts[feat]) for _ in cells]
            fh.write(feat + '\t' + '\t'.join(vals) + '\n')
    return p


def run_evaluate(tmp_path, gt, counts, locus_map, granularity, feature_set):
    prefix = str(tmp_path / f'{feature_set}_{granularity}')
    cmd = [
        sys.executable, EVALUATE,
        '--ground-truth', str(gt),
        '--observed-counts', str(counts),
        '--aligner', 'test',
        '--multimapper-mode', 'unique',
        '--granularity', granularity,
        '--feature-set', feature_set,
        '--locus-map', str(locus_map),
        '--output-prefix', prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'evaluate.py failed:\n{result.stderr}')
    return read_global_metrics(prefix + '_global_metrics.tsv')


def read_global_metrics(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return next(reader)


def parse_float(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return float('nan')


# ---------------------------------------------------------------------------
# Test 1: Structural - gene_id overlap between genic and intergenic locus_maps
# ---------------------------------------------------------------------------

def test_gene_ids_overlap_between_genic_and_intergenic_locus_maps():
    """
    Documents the structural precondition for the bug:
    AluSz6 and AluSx appear as gene_ids in BOTH genic and intergenic locus_maps.
    This is expected in real data (repeat families span both contexts) and is
    the reason why locus-level filtering before aggregation is necessary.
    """
    genic_gene_ids = {row[1] for row in GENIC_LOCUS_MAP}
    intergenic_gene_ids = {row[1] for row in INTERGENIC_LOCUS_MAP}
    shared = genic_gene_ids & intergenic_gene_ids
    assert len(shared) > 0, (
        'Expected at least one gene_id in both genic and intergenic maps '
        f'but got genic={genic_gene_ids}, intergenic={intergenic_gene_ids}'
    )
    assert 'AluSz6' in shared
    assert 'AluSx' in shared


def test_full_locus_map_partitions_at_locus_level_not_gene_id():
    """
    At locus level, genic and intergenic locus_maps are disjoint.
    At gene_id level they overlap.  This is why the old post-aggregation
    filter could not correctly separate the two subsets.
    """
    genic_loci = {row[0] for row in GENIC_LOCUS_MAP}
    intergenic_loci = {row[0] for row in INTERGENIC_LOCUS_MAP}
    # Locus level: disjoint
    assert genic_loci & intergenic_loci == set()
    # Gene_id level: overlapping
    genic_gene_ids = {row[1] for row in GENIC_LOCUS_MAP}
    intergenic_gene_ids = {row[1] for row in INTERGENIC_LOCUS_MAP}
    assert genic_gene_ids & intergenic_gene_ids != set()


# ---------------------------------------------------------------------------
# Test 2: Full repeats, gene_id - perfect observer -> perfect metrics (baseline)
# ---------------------------------------------------------------------------

def test_full_repeats_gene_id_perfect_observer_gives_perfect_metrics(tmp_path):
    """
    Sanity baseline: if observed = truth for all repeats at gene_id level,
    all accuracy metrics should be 1.0 / 0.0.
    """
    gt = write_ground_truth(tmp_path)
    lm = write_locus_map(tmp_path, 'full_locus_map.tsv', FULL_LOCUS_MAP)
    counts = write_count_matrix(tmp_path, 'full_counts.tsv', FULL_OBSERVED_GENE_ID)
    m = run_evaluate(tmp_path, gt, counts, lm, 'gene_id', 'repeats')
    assert parse_float(m['pearson_r']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['recall']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['precision']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['log1p_rmse']) == pytest.approx(0.0, abs=1e-4)


# ---------------------------------------------------------------------------
# Test 3: Genic repeats, locus granularity - perfect observer -> perfect metrics
# ---------------------------------------------------------------------------

def test_genic_repeats_locus_granularity_perfect_observer(tmp_path):
    """
    At locus level, a perfect genic observer should give perfect metrics.
    No cross-partition contamination is possible at locus level because
    locus IDs are unique per repeat instance.
    """
    gt = write_ground_truth(tmp_path)
    lm = write_locus_map(tmp_path, 'genic_locus_map.tsv', GENIC_LOCUS_MAP)
    counts = write_count_matrix(tmp_path, 'genic_locus_counts.tsv', GENIC_OBSERVED_LOCUS)
    m = run_evaluate(tmp_path, gt, counts, lm, 'locus', 'genic_repeats')
    assert parse_float(m['pearson_r']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['recall']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['precision']) == pytest.approx(1.0, abs=1e-4)


# ---------------------------------------------------------------------------
# Test 4 (THE KEY TEST): Genic repeats, gene_id - perfect genic observer
#          should now give perfect metrics after the fix
# ---------------------------------------------------------------------------

def test_genic_repeats_gene_id_perfect_genic_observer_gives_perfect_metrics(tmp_path):
    """
    Core regression test for the locus-level filtering fix.

    With the OLD code (post-aggregation filter):
      - truth['cell_001']['AluSz6'] = 5 + 3 = 8  (sums ALL AluSz6 loci)
      - observed['cell_001']['AluSz6'] = 5  (only genic loci)
      - -> truth > observed -> degraded pearson, false-negative for TcMar

    With the FIXED code (locus-level filter BEFORE aggregation):
      - only genic loci contribute to truth
      - truth['cell_001']['AluSz6'] = 5, truth['cell_001']['AluSx'] = 7
      - TcMar: its intergenic count is excluded; genic_tcmar_ref is not expressed
        -> TcMar absent from both truth and observed -> no false negative
      - observed matches truth exactly -> all metrics = 1.0
    """
    gt = write_ground_truth(tmp_path)
    lm = write_locus_map(tmp_path, 'genic_locus_map.tsv', GENIC_LOCUS_MAP)
    counts = write_count_matrix(tmp_path, 'genic_gene_id_counts.tsv', GENIC_OBSERVED_GENE_ID)
    m = run_evaluate(tmp_path, gt, counts, lm, 'gene_id', 'genic_repeats')

    pearson = parse_float(m['pearson_r'])
    recall = parse_float(m['recall'])
    precision = parse_float(m['precision'])
    rmse = parse_float(m['log1p_rmse'])

    assert pearson == pytest.approx(1.0, abs=1e-4), (
        f'pearson={pearson} - expected 1.0; '
        'if this fails, the locus-level filtering fix may have been reverted'
    )
    assert recall == pytest.approx(1.0, abs=1e-4), f'recall={recall}'
    assert precision == pytest.approx(1.0, abs=1e-4), f'precision={precision}'
    assert rmse == pytest.approx(0.0, abs=1e-4), f'log1p_rmse={rmse}'


# ---------------------------------------------------------------------------
# Test 5: Intergenic repeats, gene_id - same fix applies
# ---------------------------------------------------------------------------

def test_intergenic_repeats_gene_id_perfect_intergenic_observer_gives_perfect_metrics(tmp_path):
    """
    Mirror of test 4 for the intergenic subset.

    After the fix, intergenic truth for AluSz6 = 3 (only intergenic_alu_1),
    not 8 (which would include genic_alu_1).
    """
    gt = write_ground_truth(tmp_path)
    lm = write_locus_map(tmp_path, 'intergenic_locus_map.tsv', INTERGENIC_LOCUS_MAP)
    counts = write_count_matrix(
        tmp_path, 'intergenic_counts.tsv', INTERGENIC_OBSERVED_GENE_ID)
    m = run_evaluate(tmp_path, gt, counts, lm, 'gene_id', 'intergenic_repeats')

    assert parse_float(m['pearson_r']) == pytest.approx(1.0, abs=1e-4), (
        f"pearson={m['pearson_r']} - expected 1.0 for perfect intergenic observer"
    )
    assert parse_float(m['recall']) == pytest.approx(1.0, abs=1e-4)
    assert parse_float(m['precision']) == pytest.approx(1.0, abs=1e-4)


# ---------------------------------------------------------------------------
# Test 6: Genic + intergenic metrics should not be lower than full-set metrics
#          for a perfect full observer (monotonicity property)
# ---------------------------------------------------------------------------

def test_subset_metrics_not_lower_than_full_for_perfect_observer(tmp_path):
    """
    Metamorphic property: if an aligner perfectly recovers ALL repeats,
    evaluating a clean subset should not give LOWER metrics than the full set.

    We use a perfect full observer and check that both genic and intergenic
    sub-evaluations also give perfect metrics (1.0).

    This test would have failed before the fix because the inflated truth
    counts caused pearson < 1.0 even for a perfect observer.
    """
    gt = write_ground_truth(tmp_path)

    # Full repeats evaluation
    lm_full = write_locus_map(tmp_path, 'full_lm.tsv', FULL_LOCUS_MAP)
    cnt_full = write_count_matrix(tmp_path, 'full_cnt.tsv', FULL_OBSERVED_GENE_ID)
    m_full = run_evaluate(tmp_path, gt, cnt_full, lm_full, 'gene_id', 'repeats')

    # Genic repeats evaluation (using genic counts only)
    lm_genic = write_locus_map(tmp_path, 'genic_lm.tsv', GENIC_LOCUS_MAP)
    cnt_genic = write_count_matrix(tmp_path, 'genic_cnt.tsv', GENIC_OBSERVED_GENE_ID)
    m_genic = run_evaluate(tmp_path, gt, cnt_genic, lm_genic, 'gene_id', 'genic_repeats')

    # Intergenic repeats evaluation
    lm_inter = write_locus_map(tmp_path, 'inter_lm.tsv', INTERGENIC_LOCUS_MAP)
    cnt_inter = write_count_matrix(
        tmp_path, 'inter_cnt.tsv', INTERGENIC_OBSERVED_GENE_ID)
    m_inter = run_evaluate(tmp_path, gt, cnt_inter, lm_inter, 'gene_id', 'intergenic_repeats')

    full_pearson = parse_float(m_full['pearson_r'])
    genic_pearson = parse_float(m_genic['pearson_r'])
    inter_pearson = parse_float(m_inter['pearson_r'])

    # All should be 1.0 (perfect observer) - none should be lower than full
    assert full_pearson == pytest.approx(1.0, abs=1e-4)
    assert genic_pearson >= full_pearson - 0.01, (
        f'genic pearson={genic_pearson} is lower than full={full_pearson}; '
        'subset metrics should not be worse than full-set metrics for a perfect observer'
    )
    assert inter_pearson >= full_pearson - 0.01, (
        f'intergenic pearson={inter_pearson} is lower than full={full_pearson}'
    )
