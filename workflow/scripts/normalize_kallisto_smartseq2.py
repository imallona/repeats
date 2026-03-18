#!/usr/bin/env python3
"""
Merge per-cell kallisto abundance.tsv files into the standard feature x cell TSV.
Each abundance.tsv has columns: target_id, length, eff_length, est_counts, tpm.
Only est_counts is used.
"""

import argparse
import os
import sys
from collections import defaultdict


def read_abundance(tsv_path):
    counts = {}
    with open(tsv_path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            target_id = parts[0]
            est_counts = float(parts[3])
            if est_counts > 0:
                counts[target_id] = est_counts
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cell-dirs', nargs='+', required=True,
                    help='Directories containing abundance.tsv, one per cell')
    ap.add_argument('--cell-ids', nargs='+', required=True,
                    help='Cell IDs in the same order as --cell-dirs')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.cell_dirs) != len(args.cell_ids):
        sys.exit('Number of --cell-dirs and --cell-ids must match')

    cell_counts = {}
    all_features = set()
    for cell_id, cell_dir in zip(args.cell_ids, args.cell_dirs):
        abundance_path = os.path.join(cell_dir, 'abundance.tsv')
        counts = read_abundance(abundance_path)
        cell_counts[cell_id] = counts
        all_features.update(counts.keys())

    sorted_features = sorted(all_features)
    sorted_cells = args.cell_ids

    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(sorted_cells) + '\n')
        for feat in sorted_features:
            row = [cell_counts[c].get(feat, 0) for c in sorted_cells]
            if all(v == 0 for v in row):
                continue
            fh.write(feat + '\t' + '\t'.join(str(v) for v in row) + '\n')

    print(f'{len(sorted_features)} expressed features across {len(sorted_cells)} cells',
          file=sys.stderr)


if __name__ == '__main__':
    main()
