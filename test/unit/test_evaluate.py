"""
Unit tests for workflow/scripts/evaluate.py

Covers ~80% of the pure-function logic: metric computations, data loading,
and the vector-alignment helpers. Does not require any external bioinformatics
tools - only scipy (already in the evaluation conda env).
"""
import csv
import io
import math
import os
import sys
import textwrap

import pytest

# conftest.py adds scripts/ to sys.path
import evaluate as ev


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_tsv(tmp_path, name, rows, fieldnames=None):
    p = tmp_path / name
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    return p


def make_ground_truth(tmp_path, rows):
    """Write a ground_truth.tsv and return its path."""
    fieldnames = ['cell_id', 'locus_id', 'repeat_id', 'family_id', 'class_id', 'true_count']
    return write_tsv(tmp_path, 'ground_truth.tsv', rows, fieldnames)


def make_count_matrix(tmp_path, feature_cell_dict, cell_ids):
    """
    Write a feature x cell TSV.
    feature_cell_dict: {feature_id: {cell_id: count}}
    """
    p = tmp_path / 'counts.tsv'
    features = sorted(feature_cell_dict)
    with open(p, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(cell_ids) + '\n')
        for feat in features:
            vals = [str(feature_cell_dict[feat].get(c, 0)) for c in cell_ids]
            fh.write(feat + '\t' + '\t'.join(vals) + '\n')
    return p


def make_locus_map(tmp_path, rows):
    """
    Write a locus_map TSV (no header):
    transcript_id  gene_id  family_id  class_id
    """
    p = tmp_path / 'locus_map.tsv'
    with open(p, 'w') as fh:
        for row in rows:
            fh.write('\t'.join(row) + '\n')
    return p


# ---------------------------------------------------------------------------
# pearson_r
# ---------------------------------------------------------------------------

def test_pearson_r_identical_vectors():
    r, p = ev.pearson_r([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])
    assert r == pytest.approx(1.0, abs=1e-9)


def test_pearson_r_perfectly_anti_correlated():
    r, _ = ev.pearson_r([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
    assert r == pytest.approx(-1.0, abs=1e-9)


def test_pearson_r_zero_correlation():
    # constant y -> variance = 0 -> scipy returns nan
    r, _ = ev.pearson_r([1, 2, 3], [5, 5, 5])
    assert math.isnan(r)


def test_pearson_r_too_few_points():
    r, p = ev.pearson_r([1, 2], [1, 2])
    assert math.isnan(r)


def test_pearson_r_known_value():
    # x = [1, 3], y = [2, 4] -> r = 1.0
    r, _ = ev.pearson_r([0, 1, 2, 3], [0, 2, 4, 6])
    assert r == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# spearman_r
# ---------------------------------------------------------------------------

def test_spearman_r_identical_vectors():
    r, _ = ev.spearman_r([10, 20, 30], [10, 20, 30])
    assert r == pytest.approx(1.0, abs=1e-9)


def test_spearman_r_rank_invariant_to_scale():
    # spearman should give 1.0 for any monotone transform
    r, _ = ev.spearman_r([1, 2, 3, 4], [1, 4, 9, 16])
    assert r == pytest.approx(1.0, abs=1e-9)


def test_spearman_r_anti_correlated():
    r, _ = ev.spearman_r([1, 2, 3], [3, 2, 1])
    assert r == pytest.approx(-1.0, abs=1e-9)


def test_spearman_r_too_few_points():
    r, _ = ev.spearman_r([1], [1])
    assert math.isnan(r)


# ---------------------------------------------------------------------------
# log1p_rmse
# ---------------------------------------------------------------------------

def test_log1p_rmse_perfect():
    """Identical vectors -> RMSE = 0."""
    assert ev.log1p_rmse([0, 5, 10, 100], [0, 5, 10, 100]) == pytest.approx(0.0)


def test_log1p_rmse_empty():
    assert math.isnan(ev.log1p_rmse([], []))


def test_log1p_rmse_known_single_pair():
    # log1p(3) - log1p(0) = log(4) ~ 1.386; RMSE = 1.386
    result = ev.log1p_rmse([3], [0])
    assert result == pytest.approx(math.log(4), rel=1e-6)


def test_log1p_rmse_symmetric():
    a, b = [1, 2, 3], [4, 5, 6]
    assert ev.log1p_rmse(a, b) == pytest.approx(ev.log1p_rmse(b, a), rel=1e-9)


# ---------------------------------------------------------------------------
# detection_metrics
# ---------------------------------------------------------------------------

def test_detection_metrics_perfect_recall_and_precision():
    expressed = {'A', 'B', 'C'}
    all_f = {'A', 'B', 'C', 'D'}
    p, r, f1, j, spec = ev.detection_metrics(expressed, expressed, all_f)
    assert p == pytest.approx(1.0)
    assert r == pytest.approx(1.0)
    assert f1 == pytest.approx(1.0)
    assert j == pytest.approx(1.0)


def test_detection_metrics_no_true_positives():
    truth = {'A', 'B'}
    obs = {'C', 'D'}
    all_f = {'A', 'B', 'C', 'D'}
    p, r, f1, j, spec = ev.detection_metrics(truth, obs, all_f)
    assert r == pytest.approx(0.0)   # no TP
    assert p == pytest.approx(0.0)   # no TP
    assert f1 == pytest.approx(0.0)
    assert j == pytest.approx(0.0)


def test_detection_metrics_all_features_observed_none_true():
    truth = set()
    obs = {'A', 'B', 'C'}
    all_f = {'A', 'B', 'C'}
    p, r, f1, j, spec = ev.detection_metrics(truth, obs, all_f)
    # TP=0, FP=3, FN=0, TN=0
    # truth empty -> TP=0, FN=0 -> TP/(TP+FN) = 0/0 -> code returns 0.0 (correct)
    assert r == pytest.approx(0.0)
    # precision = 0 (TP=0, FP=3)
    assert p == pytest.approx(0.0)


def test_detection_metrics_partial_overlap():
    truth = {'A', 'B', 'C', 'D'}
    obs = {'C', 'D', 'E', 'F'}
    all_f = {'A', 'B', 'C', 'D', 'E', 'F', 'G'}
    p, r, f1, j, spec = ev.detection_metrics(truth, obs, all_f)
    # TP=2, FP=2, FN=2, TN=1
    assert p == pytest.approx(2 / 4)
    assert r == pytest.approx(2 / 4)
    assert j == pytest.approx(2 / 6)


def test_detection_metrics_specificity():
    truth = {'A'}
    obs = {'A'}
    all_f = {'A', 'B', 'C', 'D'}
    # TN=3, FP=0 -> spec = 1.0
    p, r, f1, j, spec = ev.detection_metrics(truth, obs, all_f)
    assert spec == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# build_aligned_vectors
# ---------------------------------------------------------------------------

def test_build_aligned_vectors_zero_fill_missing():
    truth = {'c1': {'A': 5, 'B': 3}}
    obs = {'c1': {'A': 4}}
    t_vec, o_vec = ev.build_aligned_vectors(truth, obs, ['c1'], ['A', 'B'])
    assert t_vec == [5, 3]
    assert o_vec == [4, 0]   # B missing in obs -> 0


def test_build_aligned_vectors_ordering():
    truth = {'c1': {'A': 1}, 'c2': {'A': 2}}
    obs = {'c1': {'A': 1}, 'c2': {'A': 2}}
    t_vec, o_vec = ev.build_aligned_vectors(truth, obs, ['c1', 'c2'], ['A'])
    assert t_vec == [1, 2]
    assert o_vec == [1, 2]


def test_build_aligned_vectors_missing_cell():
    truth = {'c1': {'A': 5}}
    obs = {}
    t_vec, o_vec = ev.build_aligned_vectors(truth, obs, ['c1'], ['A'])
    assert t_vec == [5]
    assert o_vec == [0]


# ---------------------------------------------------------------------------
# load_ground_truth - granularity modes
# ---------------------------------------------------------------------------

GT_ROWS = [
    {'cell_id': 'c1', 'locus_id': 'AluSz6_dup1', 'repeat_id': 'AluSz6',
     'family_id': 'Alu', 'class_id': 'SINE', 'true_count': '5'},
    {'cell_id': 'c1', 'locus_id': 'AluSz6_dup2', 'repeat_id': 'AluSz6',
     'family_id': 'Alu', 'class_id': 'SINE', 'true_count': '3'},
    {'cell_id': 'c1', 'locus_id': 'L1PA2_dup1', 'repeat_id': 'L1PA2',
     'family_id': 'L1', 'class_id': 'LINE', 'true_count': '7'},
    {'cell_id': 'c2', 'locus_id': 'AluSz6_dup1', 'repeat_id': 'AluSz6',
     'family_id': 'Alu', 'class_id': 'SINE', 'true_count': '2'},
]


def test_load_ground_truth_locus_granularity(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, meta = ev.load_ground_truth(str(p), granularity='locus')
    # Each locus_id is a separate feature
    assert truth['c1']['AluSz6_dup1'] == 5
    assert truth['c1']['AluSz6_dup2'] == 3
    assert truth['c1']['L1PA2_dup1'] == 7
    assert truth['c2']['AluSz6_dup1'] == 2


def test_load_ground_truth_gene_id_granularity(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, meta = ev.load_ground_truth(str(p), granularity='gene_id')
    # AluSz6_dup1 and AluSz6_dup2 aggregate to AluSz6
    assert truth['c1']['AluSz6'] == 8
    assert truth['c1']['L1PA2'] == 7
    assert truth['c2']['AluSz6'] == 2


def test_load_ground_truth_family_id_granularity(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, meta = ev.load_ground_truth(str(p), granularity='family_id')
    # Both AluSz6 loci -> 'Alu' family
    assert truth['c1']['Alu'] == 8
    assert truth['c1']['L1'] == 7


def test_load_ground_truth_class_id_granularity(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, meta = ev.load_ground_truth(str(p), granularity='class_id')
    assert truth['c1']['SINE'] == 8   # Alu is SINE
    assert truth['c1']['LINE'] == 7


def test_load_ground_truth_valid_locus_ids_filters_before_aggregation(tmp_path):
    """
    With valid_locus_ids = {AluSz6_dup1 only}, AluSz6 gene_id should
    total 5, not 8 (excludes AluSz6_dup2 which is not in valid set).
    This exercises the core fix for the genic/intergenic partitioning bug.
    """
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, _ = ev.load_ground_truth(
        str(p), granularity='gene_id',
        valid_locus_ids={'AluSz6_dup1', 'L1PA2_dup1'}
    )
    assert truth['c1']['AluSz6'] == 5   # only dup1, not dup2
    assert truth['c1']['L1PA2'] == 7
    assert 'AluSz6_dup2' not in truth['c1']


def test_load_ground_truth_valid_locus_ids_at_locus_granularity(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, _ = ev.load_ground_truth(
        str(p), granularity='locus',
        valid_locus_ids={'AluSz6_dup1'}
    )
    assert 'AluSz6_dup1' in truth['c1']
    assert 'AluSz6_dup2' not in truth.get('c1', {})


# ---------------------------------------------------------------------------
# load_count_matrix
# ---------------------------------------------------------------------------

def test_load_count_matrix_basic(tmp_path):
    p = make_count_matrix(
        tmp_path,
        {'AluSz6': {'c1': 5, 'c2': 0}, 'L1PA2': {'c1': 3, 'c2': 7}},
        ['c1', 'c2']
    )
    obs, matrix_feats = ev.load_count_matrix(str(p))
    assert obs['c1']['AluSz6'] == 5
    assert obs['c2']['L1PA2'] == 7
    # zero counts are not stored
    assert 'AluSz6' not in obs.get('c2', {})
    # but the feature is in matrix_features (zero-count rows are tracked)
    assert 'AluSz6' in matrix_feats


def test_load_count_matrix_all_features_tracked(tmp_path):
    p = make_count_matrix(
        tmp_path,
        {'A': {'c1': 0}, 'B': {'c1': 1}},
        ['c1']
    )
    _, matrix_feats = ev.load_count_matrix(str(p))
    assert 'A' in matrix_feats
    assert 'B' in matrix_feats


# ---------------------------------------------------------------------------
# compute_metrics_for_subset - end-to-end
# ---------------------------------------------------------------------------

def test_compute_metrics_perfect_observer():
    truth = {'c1': {'A': 5, 'B': 3, 'C': 0}, 'c2': {'A': 1, 'B': 8, 'C': 2}}
    obs = truth
    metrics = ev.compute_metrics_for_subset(truth, obs, ['c1', 'c2'], ['A', 'B', 'C'])
    assert metrics['pearson_r'] == pytest.approx(1.0, abs=1e-6)
    assert metrics['recall'] == pytest.approx(1.0)
    assert metrics['precision'] == pytest.approx(1.0)
    assert metrics['log1p_rmse'] == pytest.approx(0.0, abs=1e-9)


def test_compute_metrics_null_observer():
    truth = {'c1': {'A': 5, 'B': 3}, 'c2': {'A': 1, 'B': 8}}
    obs = {}
    metrics = ev.compute_metrics_for_subset(truth, obs, ['c1', 'c2'], ['A', 'B'])
    assert metrics['recall'] == pytest.approx(0.0)
    assert metrics['precision'] == pytest.approx(0.0 if metrics['precision'] is not None else 0.0)


def test_compute_metrics_returns_all_expected_keys():
    truth = {'c1': {'A': 1, 'B': 2, 'C': 3}}
    obs = {'c1': {'A': 1, 'B': 2, 'C': 3}}
    m = ev.compute_metrics_for_subset(truth, obs, ['c1'], ['A', 'B', 'C'])
    for key in ('pearson_r', 'spearman_r', 'log1p_rmse',
                'precision', 'recall', 'f1', 'jaccard', 'specificity'):
        assert key in m, f"Missing key: {key}"


# ---------------------------------------------------------------------------
# compute_per_cell_metrics
# ---------------------------------------------------------------------------

def test_compute_per_cell_metrics_basic():
    truth = {'c1': {'A': 5, 'B': 3}, 'c2': {'A': 1}}
    obs   = {'c1': {'A': 4},          'c2': {'A': 1, 'B': 2}}
    rows  = ev.compute_per_cell_metrics(truth, obs, ['c1', 'c2'], ['A', 'B'])
    assert len(rows) == 2
    assert rows[0]['cell_id'] == 'c1'
    assert rows[0]['n_truth_expressed'] == 2
    assert rows[1]['n_observed_expressed'] == 2


def test_compute_per_cell_metrics_empty_cell():
    truth = {'c1': {}}
    obs   = {}
    rows  = ev.compute_per_cell_metrics(truth, obs, ['c1'], ['A'])
    assert rows[0]['n_truth_expressed'] == 0
    assert rows[0]['n_observed_expressed'] == 0


def test_compute_per_cell_metrics_returns_all_keys():
    truth = {'c1': {'A': 3, 'B': 1, 'C': 2}}
    obs   = {'c1': {'A': 3, 'B': 1, 'C': 2}}
    rows  = ev.compute_per_cell_metrics(truth, obs, ['c1'], ['A', 'B', 'C'])
    assert 'pearson_r' in rows[0]
    assert 'spearman_r' in rows[0]
    assert rows[0]['pearson_r'] == pytest.approx(1.0, abs=1e-4)


# ---------------------------------------------------------------------------
# load_benchmark
# ---------------------------------------------------------------------------

def test_load_benchmark_existing_file(tmp_path):
    bench = tmp_path / 'bench.txt'
    bench.write_text('s\tcpu_time\tmax_rss\tio_in\tio_out\n'
                     '12.3\t11.2\t500.0\t100.0\t200.0\n')
    result = ev.load_benchmark(str(bench))
    assert result['wall_time_s'] == '12.3'
    assert result['max_rss_mb'] == '500.0'
    assert result['io_in_mb'] == '100.0'


def test_load_benchmark_missing_path():
    result = ev.load_benchmark('/nonexistent/path.txt')
    assert result == {}


def test_load_benchmark_none():
    assert ev.load_benchmark(None) == {}


# ---------------------------------------------------------------------------
# write_tsv
# ---------------------------------------------------------------------------

def test_write_tsv_basic(tmp_path):
    rows = [{'a': 1, 'b': 2}, {'a': 3, 'b': 4}]
    out  = tmp_path / 'out.tsv'
    ev.write_tsv(rows, str(out))
    with open(out) as fh:
        lines = fh.readlines()
    assert lines[0].strip() == 'a\tb'
    assert lines[1].strip() == '1\t2'
    assert lines[2].strip() == '3\t4'


def test_write_tsv_empty_with_fallback_fields(tmp_path):
    out = tmp_path / 'out.tsv'
    ev.write_tsv([], str(out), fallback_fields=['x', 'y'])
    with open(out) as fh:
        content = fh.read()
    assert 'x\ty' in content


# ---------------------------------------------------------------------------
# load_ground_truth default granularity fallback
# ---------------------------------------------------------------------------

def test_load_ground_truth_unknown_granularity_falls_back(tmp_path):
    p = make_ground_truth(tmp_path, GT_ROWS)
    truth, _ = ev.load_ground_truth(str(p), granularity='unknown_gran')
    # falls back to repeat_id key (same as gene_id)
    assert 'AluSz6' in truth['c1']


# ---------------------------------------------------------------------------
# main() end-to-end integration tests
# ---------------------------------------------------------------------------

def _write_gt(path):
    path.write_text(
        'cell_id\tlocus_id\trepeat_id\tfamily_id\tclass_id\ttrue_count\n'
        'c1\tAluSz6_dup1\tAluSz6\tAlu\tSINE\t5\n'
        'c1\tL1PA2_dup1\tL1PA2\tL1\tLINE\t3\n'
        'c2\tAluSz6_dup1\tAluSz6\tAlu\tSINE\t2\n'
    )


def _write_counts(path):
    path.write_text(
        'feature_id\tc1\tc2\n'
        'AluSz6\t4\t2\n'
        'L1PA2\t3\t0\n'
    )


def test_evaluate_main_creates_output_files(tmp_path, monkeypatch):
    import sys
    gt     = tmp_path / 'gt.tsv'
    counts = tmp_path / 'counts.tsv'
    _write_gt(gt)
    _write_counts(counts)
    prefix = str(tmp_path / 'out')
    monkeypatch.setattr(sys, 'argv', [
        'evaluate.py',
        '--ground-truth', str(gt),
        '--observed-counts', str(counts),
        '--aligner', 'test_aligner',
        '--output-prefix', prefix,
    ])
    ev.main()
    assert os.path.exists(prefix + '_global_metrics.tsv')
    assert os.path.exists(prefix + '_per_cell_metrics.tsv')
    assert os.path.exists(prefix + '_per_family_metrics.tsv')


def test_evaluate_main_global_metrics_content(tmp_path, monkeypatch):
    import sys
    gt     = tmp_path / 'gt.tsv'
    counts = tmp_path / 'counts.tsv'
    _write_gt(gt)
    _write_counts(counts)
    prefix = str(tmp_path / 'out')
    monkeypatch.setattr(sys, 'argv', [
        'evaluate.py',
        '--ground-truth', str(gt),
        '--observed-counts', str(counts),
        '--aligner', 'starsolo',
        '--multimapper-mode', 'unique',
        '--granularity', 'gene_id',
        '--feature-set', 'repeats',
        '--output-prefix', prefix,
    ])
    ev.main()
    import csv as _csv
    with open(prefix + '_global_metrics.tsv') as fh:
        rows = list(_csv.DictReader(fh, delimiter='\t'))
    assert len(rows) == 1
    assert rows[0]['aligner'] == 'starsolo'
    assert rows[0]['granularity'] == 'gene_id'


def test_evaluate_main_with_locus_map(tmp_path, monkeypatch):
    import sys
    gt     = tmp_path / 'gt.tsv'
    counts = tmp_path / 'counts.tsv'
    lm     = tmp_path / 'locus_map.tsv'
    _write_gt(gt)
    _write_counts(counts)
    lm.write_text('AluSz6_dup1\tAluSz6\tAlu\tSINE\nL1PA2_dup1\tL1PA2\tL1\tLINE\n')
    prefix = str(tmp_path / 'out')
    monkeypatch.setattr(sys, 'argv', [
        'evaluate.py',
        '--ground-truth', str(gt),
        '--observed-counts', str(counts),
        '--aligner', 'alevin',
        '--locus-map', str(lm),
        '--output-prefix', prefix,
    ])
    ev.main()
    assert os.path.exists(prefix + '_global_metrics.tsv')


def test_evaluate_main_with_benchmark(tmp_path, monkeypatch):
    import sys
    gt     = tmp_path / 'gt.tsv'
    counts = tmp_path / 'counts.tsv'
    bench  = tmp_path / 'bench.txt'
    _write_gt(gt)
    _write_counts(counts)
    bench.write_text('s\tcpu_time\tmax_rss\tio_in\tio_out\n5.0\t4.5\t300.0\t50.0\t80.0\n')
    prefix = str(tmp_path / 'out')
    monkeypatch.setattr(sys, 'argv', [
        'evaluate.py',
        '--ground-truth', str(gt),
        '--observed-counts', str(counts),
        '--aligner', 'kallisto',
        '--benchmark', str(bench),
        '--output-prefix', prefix,
    ])
    ev.main()
    import csv as _csv
    with open(prefix + '_global_metrics.tsv') as fh:
        rows = list(_csv.DictReader(fh, delimiter='\t'))
    assert rows[0]['wall_time_s'] == '5.0'
