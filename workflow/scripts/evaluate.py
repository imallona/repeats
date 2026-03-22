#!/usr/bin/env python3
"""
Compare aligner output against simulation ground truth and produce metrics tables.

Outputs:
  {output_prefix}_global_metrics.tsv    - one row with global accuracy and compute metrics
  {output_prefix}_per_cell_metrics.tsv  - per-cell Pearson and Spearman correlations
  {output_prefix}_per_family_metrics.tsv - metrics broken down by repeat class_id
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from scipy import stats


def load_ground_truth(gt_path, granularity='gene_id', valid_locus_ids=None):
    """
    Returns:
      truth: {cell_id: {feature_id: count}}
      repeat_meta: {feature_id: (family_id, class_id)}

    granularity controls how ground truth counts are aggregated:
      locus     - use locus_id (finest level: individual repeat instance)
      gene_id   - group by repeat_id (= gene_id; sums across loci of same gene)
      family_id - sum true_count by cell_id + family_id
      class_id  - sum true_count by cell_id + class_id
    """
    from collections import defaultdict as _dd
    raw = []
    with open(gt_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        cell_col = reader.fieldnames[0]
        for row in reader:
            raw.append({
                'cell_id': row[cell_col],
                'locus_id': row.get('locus_id', row.get('repeat_id', '')),
                'repeat_id': row.get('repeat_id', row.get('locus_id', '')),
                'family_id': row['family_id'],
                'class_id': row['class_id'],
                'true_count': int(row['true_count'])
            })

    if granularity == 'locus':
        key_col = 'locus_id'
    elif granularity == 'gene_id':
        key_col = 'repeat_id'
    elif granularity == 'family_id':
        key_col = 'family_id'
    elif granularity == 'class_id':
        key_col = 'class_id'
    else:
        key_col = 'repeat_id'

    truth = _dd(lambda: _dd(int))
    repeat_meta = {}
    for r in raw:
        # Filter at locus level before aggregating so cross-partition
        # gene_ids (e.g. AluSz6 in both genic and intergenic) do not
        # inflate truth counts for a subset evaluation.
        if valid_locus_ids is not None and r['locus_id'] not in valid_locus_ids:
            continue
        feat = r[key_col]
        truth[r['cell_id']][feat] += r['true_count']
        repeat_meta[feat] = (r['family_id'], r['class_id'])

    return {c: dict(v) for c, v in truth.items()}, repeat_meta


def load_count_matrix(tsv_path):
    """
    Returns:
      observed: {cell_id: {repeat_id: count}}
      all_matrix_features: set of all feature IDs in the matrix (incl. zero-count rows)
    """
    observed = defaultdict(dict)
    all_matrix_features = set()
    with open(tsv_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        cell_ids = header[1:]
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            repeat_id = parts[0]
            all_matrix_features.add(repeat_id)
            for cell_idx, val in enumerate(parts[1:]):
                count = float(val)
                if count > 0:
                    observed[cell_ids[cell_idx]][repeat_id] = count
    return dict(observed), all_matrix_features


def pearson_r(x_vals, y_vals):
    if len(x_vals) < 3:
        return float('nan'), float('nan')
    r, p = stats.pearsonr(x_vals, y_vals)
    return r, p


def spearman_r(x_vals, y_vals):
    if len(x_vals) < 3:
        return float('nan'), float('nan')
    r, p = stats.spearmanr(x_vals, y_vals)
    return r, p


def log1p_rmse(x_vals, y_vals):
    n = len(x_vals)
    if n == 0:
        return float('nan')
    total = sum((math.log1p(x) - math.log1p(y)) ** 2 for x, y in zip(x_vals, y_vals))
    return math.sqrt(total / n)


def detection_metrics(truth_expressed_set, observed_expressed_set, all_features):
    tp = len(truth_expressed_set & observed_expressed_set)
    fp = len(observed_expressed_set - truth_expressed_set)
    fn = len(truth_expressed_set - observed_expressed_set)
    tn = len(all_features - truth_expressed_set - observed_expressed_set)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    jaccard = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    return precision, recall, f1, jaccard, specificity


def build_aligned_vectors(truth, observed, cells, features):
    truth_vec = []
    obs_vec = []
    for cell_id in cells:
        t_cell = truth.get(cell_id, {})
        o_cell = observed.get(cell_id, {})
        for feat in features:
            truth_vec.append(t_cell.get(feat, 0))
            obs_vec.append(o_cell.get(feat, 0))
    return truth_vec, obs_vec


def compute_metrics_for_subset(truth, observed, cells, features):
    truth_vec, obs_vec = build_aligned_vectors(truth, observed, cells, features)

    pearson, _ = pearson_r(truth_vec, obs_vec)
    spearman, _ = spearman_r(truth_vec, obs_vec)
    rmse_log = log1p_rmse(truth_vec, obs_vec)

    all_features_set = set(features)
    truth_expressed = {f for cell in cells for f in truth.get(cell, {}) if f in all_features_set}
    obs_expressed = {f for cell in cells for f in observed.get(cell, {}) if f in all_features_set}

    precision, recall, f1, jaccard, specificity = detection_metrics(
        truth_expressed, obs_expressed, all_features_set)

    return {
        'n_cells': len(cells),
        'n_features': len(features),
        'pearson_r': round(pearson, 4) if not math.isnan(pearson) else 'NA',
        'spearman_r': round(spearman, 4) if not math.isnan(spearman) else 'NA',
        'log1p_rmse': round(rmse_log, 4) if not math.isnan(rmse_log) else 'NA',
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1': round(f1, 4),
        'jaccard': round(jaccard, 4),
        'specificity': round(specificity, 4)
    }


def compute_per_cell_metrics(truth, observed, cells, all_features):
    rows = []
    for cell_id in cells:
        t_cell = truth.get(cell_id, {})
        o_cell = observed.get(cell_id, {})
        t_vals = [t_cell.get(f, 0) for f in all_features]
        o_vals = [o_cell.get(f, 0) for f in all_features]
        pearson, _ = pearson_r(t_vals, o_vals)
        spearman, _ = spearman_r(t_vals, o_vals)
        rows.append({
            'cell_id': cell_id,
            'pearson_r': round(pearson, 4) if not math.isnan(pearson) else 'NA',
            'spearman_r': round(spearman, 4) if not math.isnan(spearman) else 'NA',
            'n_truth_expressed': len(t_cell),
            'n_observed_expressed': len(o_cell)
        })
    return rows


def load_benchmark(benchmark_path):
    if not benchmark_path or not os.path.exists(benchmark_path):
        return {}
    with open(benchmark_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        row = next(reader, {})
    return {
        'wall_time_s': row.get('s', 'NA'),
        'cpu_time_s': row.get('cpu_time', 'NA'),
        'max_rss_mb': row.get('max_rss', 'NA'),
        'io_in_mb': row.get('io_in', 'NA'),
        'io_out_mb': row.get('io_out', 'NA')
    }


def write_tsv(rows, path, fallback_fields=None):
    fields = rows[0].keys() if rows else (fallback_fields or [])
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        if fields:
            writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ground-truth', required=True)
    ap.add_argument('--observed-counts', required=True,
                    help='Standard feature x cell TSV from a normalize_*.py script')
    ap.add_argument('--aligner', required=True)
    ap.add_argument('--multimapper-mode', default='unique')
    ap.add_argument('--granularity', default='gene_id',
                    choices=['locus', 'gene_id', 'family_id', 'class_id'])
    ap.add_argument('--feature-set', default='repeats',
                    help='Feature set being evaluated (repeats, genic_repeats, '
                         'intergenic_repeats). Used to filter ground truth to the '
                         'feature_set\'s loci when --locus-map is provided.')
    ap.add_argument('--benchmark', default=None)
    ap.add_argument('--locus-map', default=None,
                    help='4-col TSV (transcript_id, gene_id, family_id, class_id). '
                         'Expands the feature universe for specificity and restricts '
                         'ground truth to this feature_set\'s loci.')
    ap.add_argument('--output-prefix', required=True)
    args = ap.parse_args()

    # Load locus map FIRST so truth is filtered at locus level before aggregation.
    # Without this, at gene_id granularity a gene_id such as AluSz6 that has copies
    # in both genic and intergenic regions carries intergenic counts into a
    # genic_repeats evaluation (and vice-versa), inflating truth vs observed.
    locus_map_features = set()
    valid_locus_ids = None
    if args.locus_map and os.path.exists(args.locus_map):
        valid_locus_ids = set()
        col = {'locus': 0, 'gene_id': 1, 'family_id': 2, 'class_id': 3}.get(
            args.granularity, 1)
        with open(args.locus_map) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if not parts or not parts[0]:
                    continue
                valid_locus_ids.add(parts[0])  # column 0 is always transcript_id
                if len(parts) > col:
                    locus_map_features.add(parts[col])
        print(f'  Locus map loaded: {len(valid_locus_ids)} loci / '
              f'{len(locus_map_features)} features at {args.granularity} level',
              file=sys.stderr)

    print(f'Loading ground truth from {args.ground_truth} at granularity={args.granularity}',
          file=sys.stderr)
    truth, repeat_meta = load_ground_truth(
        args.ground_truth, granularity=args.granularity, valid_locus_ids=valid_locus_ids)

    print(f'Loading observed counts from {args.observed_counts}', file=sys.stderr)
    observed, matrix_features = load_count_matrix(args.observed_counts)

    feature_universe = (
        matrix_features |
        {f for cell in truth for f in truth[cell]} |
        {f for cell in observed for f in observed[cell]} |
        locus_map_features
    )
    all_features = sorted(feature_universe)

    common_cells = sorted(set(truth.keys()) & set(observed.keys()))
    print(f'  {len(common_cells)} common cells, {len(all_features)} features', file=sys.stderr)

    global_metrics = compute_metrics_for_subset(truth, observed, common_cells, all_features)
    global_metrics['aligner'] = args.aligner
    global_metrics['multimapper_mode'] = args.multimapper_mode
    global_metrics['feature_set'] = args.feature_set
    global_metrics['granularity'] = args.granularity

    bench = load_benchmark(args.benchmark)
    global_metrics.update(bench)

    write_tsv([global_metrics], args.output_prefix + '_global_metrics.tsv')

    per_cell_rows = compute_per_cell_metrics(truth, observed, common_cells, all_features)
    for row in per_cell_rows:
        row['aligner'] = args.aligner
        row['multimapper_mode'] = args.multimapper_mode
        row['feature_set'] = args.feature_set
        row['granularity'] = args.granularity
    write_tsv(per_cell_rows, args.output_prefix + '_per_cell_metrics.tsv',
              fallback_fields=['cell_id', 'pearson_r', 'spearman_r',
                               'n_truth_expressed', 'n_observed_expressed',
                               'aligner', 'multimapper_mode', 'feature_set', 'granularity'])

    class_to_features = defaultdict(list)
    for feat in all_features:
        if feat in repeat_meta:
            class_id = repeat_meta[feat][1]
        else:
            class_id = 'unknown'
        class_to_features[class_id].append(feat)

    per_family_rows = []
    for class_id, class_features in sorted(class_to_features.items()):
        row = compute_metrics_for_subset(truth, observed, common_cells, class_features)
        row['class_id'] = class_id
        row['aligner'] = args.aligner
        row['multimapper_mode'] = args.multimapper_mode
        row['feature_set'] = args.feature_set
        row['granularity'] = args.granularity
        per_family_rows.append(row)
    write_tsv(per_family_rows, args.output_prefix + '_per_family_metrics.tsv',
              fallback_fields=['class_id', 'aligner', 'multimapper_mode', 'feature_set',
                               'granularity', 'pearson_r', 'spearman_r', 'n_cells',
                               'n_features', 'sensitivity', 'specificity'])

    print(f'Metrics written to {args.output_prefix}_*.tsv', file=sys.stderr)


if __name__ == '__main__':
    main()
