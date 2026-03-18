#!/usr/bin/env python3
"""
Convert salmon alevin output to the standard feature x cell TSV.
Reads quants_mat.gz, quants_mat_rows.txt (cells), quants_mat_cols.txt (features).
Processes row by row to avoid loading the full dense matrix.
"""

import argparse
import gzip
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--alevin-dir', required=True,
                    help='Directory containing quants_mat.gz, quants_mat_rows.txt, quants_mat_cols.txt')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    alevin_dir = args.alevin_dir

    with open(os.path.join(alevin_dir, 'quants_mat_rows.txt')) as fh:
        cell_barcodes = [line.rstrip('\n') for line in fh]
    with open(os.path.join(alevin_dir, 'quants_mat_cols.txt')) as fh:
        feature_ids = [line.rstrip('\n') for line in fh]

    n_cells = len(cell_barcodes)
    n_features = len(feature_ids)
    print(f'{n_cells} cells, {n_features} features', file=sys.stderr)

    feature_cell_counts = {}

    with gzip.open(os.path.join(alevin_dir, 'quants_mat.gz'), 'rt') as fh:
        for cell_idx, line in enumerate(fh):
            values = line.rstrip('\n').split('\t')
            for feat_idx, val in enumerate(values):
                count = float(val)
                if count > 0:
                    feat_id = feature_ids[feat_idx]
                    if feat_id not in feature_cell_counts:
                        feature_cell_counts[feat_id] = {}
                    feature_cell_counts[feat_id][cell_barcodes[cell_idx]] = count

    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(cell_barcodes) + '\n')
        for feat_id in sorted(feature_cell_counts):
            row = [feature_cell_counts[feat_id].get(cb, 0) for cb in cell_barcodes]
            fh.write(feat_id + '\t' + '\t'.join(str(v) for v in row) + '\n')

    print(f'{len(feature_cell_counts)} expressed features written', file=sys.stderr)


if __name__ == '__main__':
    main()
