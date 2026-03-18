#!/usr/bin/env python3
"""
Merge per-cell salmon quant.sf files into the standard feature x cell TSV.
Uses NumReads column from quant.sf.
"""

import argparse
import os
import sys


def read_quant_sf(sf_path):
    counts = {}
    with open(sf_path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            target_id = parts[0]
            num_reads = float(parts[4])
            if num_reads > 0:
                counts[target_id] = num_reads
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cell-dirs', nargs='+', required=True)
    ap.add_argument('--cell-ids', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.cell_dirs) != len(args.cell_ids):
        sys.exit('Number of --cell-dirs and --cell-ids must match')

    cell_counts = {}
    all_features = set()
    for cell_id, cell_dir in zip(args.cell_ids, args.cell_dirs):
        sf_path = os.path.join(cell_dir, 'quant.sf')
        counts = read_quant_sf(sf_path)
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
