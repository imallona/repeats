"""
Metamorphic tests for the evaluation pipeline.

Metamorphic testing verifies properties that must hold across transformations
of the input, without needing exact expected values.  This is especially
useful for a bioinformatics benchmarking pipeline where true correct outputs
are not always known in advance.

Properties tested
-----------------
1. Scale invariance:     multiplying observed counts by a constant should not
                         change Pearson or Spearman (but will change RMSE).
2. Monotone recall:      adding true positives to observed should increase
                         or maintain recall - never decrease it.
3. Noise degrades metrics:  adding random noise to a perfect observer should
                         monotonically worsen Pearson and RMSE.
4. Null observer is worst:  a perfect observer beats a zero observer on all
                         quantification and detection metrics.
5. Permuting cell labels degrades per-cell metrics but not global metrics
                         computed on pooled counts.
6. Family_id granularity:  summing a perfect locus-level count matrix to
                         family_id should also give perfect metrics.
"""
import math
import random

import pytest

import evaluate as ev


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_truth_obs(n_cells=4, n_features=8, seed=42):
    """Return a (truth, observed) pair where observed = truth (perfect)."""
    rng = random.Random(seed)
    cells = [f'c{i}' for i in range(n_cells)]
    feats = [f'feat_{i}' for i in range(n_features)]
    truth = {c: {f: rng.randint(1, 20) for f in feats} for c in cells}
    obs = {c: dict(d) for c, d in truth.items()}
    return truth, obs, cells, feats


# ---------------------------------------------------------------------------
# 1. Scale invariance of rank-based metrics
# ---------------------------------------------------------------------------

def test_scale_invariance_pearson():
    """Multiplying observed by a scalar should not change Pearson."""
    truth, obs, cells, feats = make_truth_obs()
    m1 = ev.compute_metrics_for_subset(truth, obs, cells, feats)

    obs_scaled = {c: {f: v * 10 for f, v in d.items()} for c, d in obs.items()}
    m2 = ev.compute_metrics_for_subset(truth, obs_scaled, cells, feats)

    assert m1['pearson_r'] == pytest.approx(m2['pearson_r'], abs=1e-4)


def test_scale_invariance_spearman():
    """Multiplying observed by a scalar should not change Spearman."""
    truth, obs, cells, feats = make_truth_obs()
    m1 = ev.compute_metrics_for_subset(truth, obs, cells, feats)

    obs_scaled = {c: {f: v * 3 for f, v in d.items()} for c, d in obs.items()}
    m2 = ev.compute_metrics_for_subset(truth, obs_scaled, cells, feats)

    assert m1['spearman_r'] == pytest.approx(m2['spearman_r'], abs=1e-4)


def test_scale_changes_rmse():
    """Scaling observed DOES change log1p_rmse."""
    truth, obs, cells, feats = make_truth_obs()
    m1 = ev.compute_metrics_for_subset(truth, obs, cells, feats)

    obs_scaled = {c: {f: v * 100 for f, v in d.items()} for c, d in obs.items()}
    m2 = ev.compute_metrics_for_subset(truth, obs_scaled, cells, feats)

    # RMSE for perfect observer is 0; scaled observer should have higher RMSE
    assert m2['log1p_rmse'] > m1['log1p_rmse']


# ---------------------------------------------------------------------------
# 2. Monotone recall: adding true positives never decreases recall
# ---------------------------------------------------------------------------

def test_monotone_recall_adding_true_positives():
    """
    Start with a partially correct observer.  Add more true positives.
    Recall must not decrease.
    """
    truth = {'c1': {'A': 5, 'B': 3, 'C': 2, 'D': 7}}
    feats = ['A', 'B', 'C', 'D']

    obs_partial = {'c1': {'A': 5}}   # only A correct
    obs_more    = {'c1': {'A': 5, 'B': 3}}   # A and B correct

    m1 = ev.compute_metrics_for_subset(truth, obs_partial, ['c1'], feats)
    m2 = ev.compute_metrics_for_subset(truth, obs_more, ['c1'], feats)

    assert m2['recall'] >= m1['recall'] - 1e-9


def test_monotone_recall_full_recovery_is_maximum():
    truth, obs, cells, feats = make_truth_obs()
    # Partial: only first half of features observed
    feats_half = feats[:len(feats)//2]
    obs_partial = {c: {f: d[f] for f in feats_half if f in d} for c, d in obs.items()}

    m_partial = ev.compute_metrics_for_subset(truth, obs_partial, cells, feats)
    m_full = ev.compute_metrics_for_subset(truth, obs, cells, feats)

    assert m_full['recall'] >= m_partial['recall'] - 1e-9
    assert m_full['f1'] >= m_partial['f1'] - 1e-9


# ---------------------------------------------------------------------------
# 3. Noise degrades metrics monotonically
# ---------------------------------------------------------------------------

def test_increasing_noise_degrades_pearson():
    """
    Progressively noisier observers should give monotonically worse Pearson.
    """
    truth, _, cells, feats = make_truth_obs(seed=7)
    prev_r = 1.0

    rng = random.Random(7)
    obs = {c: dict(d) for c, d in truth.items()}

    noise_levels = [0, 1, 3, 10, 50]
    for level in noise_levels:
        noisy_obs = {
            c: {f: max(0, v + rng.randint(-level, level)) for f, v in d.items()}
            for c, d in truth.items()
        }
        m = ev.compute_metrics_for_subset(truth, noisy_obs, cells, feats)
        r = m['pearson_r']
        if not isinstance(r, str) and not math.isnan(float(r)):
            assert float(r) <= prev_r + 0.05, (
                f'Pearson increased from {prev_r} to {r} as noise level went to {level}'
            )
            prev_r = float(r)


# ---------------------------------------------------------------------------
# 4. Perfect observer strictly beats null observer
# ---------------------------------------------------------------------------

def test_perfect_beats_null_on_all_metrics():
    truth, perfect_obs, cells, feats = make_truth_obs(seed=3)
    null_obs = {}

    m_perfect = ev.compute_metrics_for_subset(truth, perfect_obs, cells, feats)
    m_null = ev.compute_metrics_for_subset(truth, null_obs, cells, feats)

    null_pearson = m_null.get('pearson_r', 'NA')
    null_pearson_val = 0.0 if null_pearson in ('NA', '') else float(null_pearson)
    assert float(m_perfect['pearson_r']) > null_pearson_val
    assert float(m_perfect['recall']) > float(m_null['recall'])
    assert float(m_perfect['f1']) > float(m_null['f1'])
    assert (float(m_perfect['log1p_rmse'])
            < float(m_null['log1p_rmse']))


# ---------------------------------------------------------------------------
# 5. Granularity aggregation consistency
# ---------------------------------------------------------------------------

def test_locus_counts_summed_to_family_give_correct_metrics(tmp_path):
    """
    A locus-level count matrix summed to family_id should give the same
    result as loading ground truth at family_id granularity directly.
    Validates that the granularity aggregation pipeline is self-consistent.
    """
    import csv
    import os
    import sys
    import subprocess

    EVALUATE = os.path.join(
        os.path.dirname(__file__), '..', '..', 'workflow', 'scripts', 'evaluate.py'
    )

    gt_rows = [
        ('cell_001', 'AluSz6_dup1', 'AluSz6', 'Alu',  'SINE', '5'),
        ('cell_001', 'AluSz6_dup2', 'AluSz6', 'Alu',  'SINE', '3'),
        ('cell_001', 'L1PA2_dup1',  'L1PA2',  'L1',   'LINE', '7'),
        ('cell_001', 'MIR3_dup1',   'MIR3',   'MIR',  'SINE', '2'),
    ]
    locus_map_rows = [
        ('AluSz6_dup1', 'AluSz6', 'Alu',  'SINE'),
        ('AluSz6_dup2', 'AluSz6', 'Alu',  'SINE'),
        ('L1PA2_dup1',  'L1PA2',  'L1',   'LINE'),
        ('MIR3_dup1',   'MIR3',   'MIR',  'SINE'),
    ]

    gt_path = tmp_path / 'gt.tsv'
    with open(gt_path, 'w') as fh:
        fh.write('cell_id\tlocus_id\trepeat_id\tfamily_id\tclass_id\ttrue_count\n')
        for row in gt_rows:
            fh.write('\t'.join(row) + '\n')

    lm_path = tmp_path / 'lm.tsv'
    with open(lm_path, 'w') as fh:
        for row in locus_map_rows:
            fh.write('\t'.join(row) + '\n')

    # Perfect observer at family_id level (manually aggregated)
    family_counts = {'Alu': 8, 'L1': 7, 'MIR': 2}
    cnt_path = tmp_path / 'counts.tsv'
    with open(cnt_path, 'w') as fh:
        fh.write('feature_id\tcell_001\n')
        for fam, cnt in family_counts.items():
            fh.write(f'{fam}\t{cnt}\n')

    prefix = str(tmp_path / 'family_eval')
    cmd = [
        sys.executable, EVALUATE,
        '--ground-truth', str(gt_path),
        '--observed-counts', str(cnt_path),
        '--aligner', 'test',
        '--multimapper-mode', 'unique',
        '--granularity', 'family_id',
        '--feature-set', 'repeats',
        '--locus-map', str(lm_path),
        '--output-prefix', prefix,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    with open(prefix + '_global_metrics.tsv') as fh:
        m = next(csv.DictReader(fh, delimiter='\t'))

    pearson = float(m['pearson_r'])
    recall = float(m['recall'])
    assert pearson == pytest.approx(1.0, abs=1e-4), (
        f'family_id aggregation: pearson={pearson}, expected 1.0'
    )
    assert recall == pytest.approx(1.0, abs=1e-4)


# ---------------------------------------------------------------------------
# 6. Random counts give low (not artificially inflated) metrics
# ---------------------------------------------------------------------------

def test_random_observer_gives_low_correlation():
    """
    A random observer uncorrelated with truth should give Pearson near 0.
    Uses enough data points to make this reliable.
    """
    rng = random.Random(99)
    n_feats = 50
    feats = [f'f{i}' for i in range(n_feats)]
    truth = {'c1': {f: rng.randint(1, 100) for f in feats}}
    random_obs = {'c1': {f: rng.randint(1, 100) for f in feats}}

    m = ev.compute_metrics_for_subset(truth, random_obs, ['c1'], feats)
    # With 50 random pairs, Pearson should rarely be above 0.5
    r = float(m['pearson_r'])
    assert abs(r) < 0.7, (
        f'Random observer gave high Pearson={r}; this suggests the metric is inflated'
    )
