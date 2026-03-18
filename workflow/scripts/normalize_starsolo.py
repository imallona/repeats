#!/usr/bin/env python3
"""
Convert a STARsolo Solo.out/Gene/raw directory to the standard
feature x cell TSV used by evaluate.py.
Rows are features (repeat gene_ids), columns are cell barcodes.
Only rows with at least one non-zero count are written.
"""

import argparse
import gzip
import os
import sys
import scipy.io


def read_lines(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return [line.rstrip('\n') for line in fh]


def load_starsolo_raw(raw_dir):
    barcodes = read_lines(os.path.join(raw_dir, 'barcodes.tsv'))
    feature_lines = read_lines(os.path.join(raw_dir, 'features.tsv'))
    feature_ids = [line.split('\t')[0] for line in feature_lines]
    mat = scipy.io.mmread(os.path.join(raw_dir, 'matrix.mtx')).tocsc()
    return barcodes, feature_ids, mat


def write_count_matrix(barcodes, feature_ids, sparse_mat, output_path):
    with open(output_path, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(barcodes) + '\n')
        for feat_idx, feat_id in enumerate(feature_ids):
            row = sparse_mat.getrow(feat_idx)
            if row.nnz == 0:
                continue
            values = row.toarray().flatten()
            fh.write(feat_id + '\t' + '\t'.join(str(int(v)) for v in values) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-dir', required=True, help='Solo.out/Gene/raw directory')
    ap.add_argument('--output', required=True, help='Output TSV path')
    args = ap.parse_args()

    barcodes, feature_ids, mat = load_starsolo_raw(args.raw_dir)
    print(f'{len(feature_ids)} features, {len(barcodes)} barcodes', file=sys.stderr)
    write_count_matrix(barcodes, feature_ids, mat, args.output)
    print(f'Written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
