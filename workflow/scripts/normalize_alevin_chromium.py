#!/usr/bin/env python3
"""
Convert salmon alevin output (--dumpMtx) to the standard feature x cell TSV.
Reads quants_mat.mtx.gz (sparse MTX), quants_mat_rows.txt (cells),
quants_mat_cols.txt (features/loci).

The MTX from alevin --dumpMtx is barcodes x genes (cells x features); it is
transposed to features x cells before output.

Supports granularity aggregation for repeat feature sets via --locus-map
(4-col TSV: transcript_id, gene_id, family_id, class_id).  The feature IDs
from alevin (quants_mat_cols.txt) are locus IDs (= gene_id column of locus_map);
these are aggregated to the requested --granularity by summing counts.
Strategy mirrors normalize_starsolo.py.
"""

import argparse
import os
import sys
from collections import defaultdict

import scipy.io


def load_locus_to_group(locus_map_path, granularity):
    """Return dict: locus_id (gene_id col) -> group_id at the requested granularity.

    locus_map columns (no header): transcript_id, gene_id, family_id, class_id.
    """
    col_idx = {'gene_id': 1, 'family_id': 2, 'class_id': 3}
    if granularity not in col_idx:
        sys.exit(f'Unknown granularity: {granularity}')
    target_col = col_idx[granularity]
    mapping = {}
    with open(locus_map_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= target_col:
                continue
            locus = parts[0]  # locus_id (col 0) as used in alevin quants_mat_cols.txt
            target = parts[target_col]
            if locus and target:
                mapping[locus] = target
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--alevin-dir', required=True,
                    help='Directory containing quants_mat.mtx.gz, quants_mat_rows.txt, quants_mat_cols.txt')
    ap.add_argument('--output', required=True)
    ap.add_argument('--locus-map', default=None,
                    help='4-col locus map TSV for granularity aggregation (repeats only)')
    ap.add_argument('--granularity', default='gene_id',
                    choices=['gene_id', 'family_id', 'class_id'],
                    help='Aggregation granularity (default: gene_id)')
    args = ap.parse_args()

    alevin_dir = args.alevin_dir

    with open(os.path.join(alevin_dir, 'quants_mat_rows.txt')) as fh:
        cell_barcodes = [line.rstrip('\n') for line in fh]
    with open(os.path.join(alevin_dir, 'quants_mat_cols.txt')) as fh:
        feature_ids = [line.rstrip('\n') for line in fh]

    # MTX from alevin --dumpMtx: rows=cells, cols=genes → transpose to features x cells
    mtx_path = os.path.join(alevin_dir, 'quants_mat.mtx.gz')
    mat = scipy.io.mmread(mtx_path).T.tocsr()
    n_features, n_cells = mat.shape
    print(f'{n_cells} cells, {n_features} features', file=sys.stderr)

    # Build locus -> target mapping
    if args.locus_map:
        locus_to_group = load_locus_to_group(args.locus_map, args.granularity)
        print(f'Aggregating to {args.granularity} ({len(locus_to_group)} loci mapped)',
              file=sys.stderr)
    else:
        locus_to_group = None

    group_counts = defaultdict(lambda: [0.0] * n_cells)

    for feat_idx, locus in enumerate(feature_ids):
        row = mat.getrow(feat_idx)
        if row.nnz == 0:
            continue
        target = locus_to_group.get(locus) if locus_to_group is not None else locus
        if target is None:
            continue
        values = row.toarray().flatten()
        dest = group_counts[target]
        for ci, v in enumerate(values):
            if v:
                dest[ci] += v

    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(cell_barcodes) + '\n')
        for feat_id in sorted(group_counts):
            row = group_counts[feat_id]
            if all(v == 0 for v in row):
                continue
            fh.write(feat_id + '\t' + '\t'.join(str(v) for v in row) + '\n')

    print(f'{len(group_counts)} expressed features written', file=sys.stderr)


if __name__ == '__main__':
    main()
