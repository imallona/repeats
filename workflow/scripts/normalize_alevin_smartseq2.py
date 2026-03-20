#!/usr/bin/env python3
"""
Merge per-cell salmon quant.sf files into the standard feature x cell TSV.
Uses NumReads column from quant.sf.
Supports optional --t2g to aggregate transcript-level counts to a higher granularity.
Target IDs may have '::chrom:start-end(strand)' suffixes which are stripped before t2g lookup.
"""

import argparse
import os
import sys
from collections import defaultdict


def load_t2g(t2g_path):
    """Returns dict: transcript_id -> group_id."""
    mapping = {}
    with open(t2g_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def read_quant_sf(sf_path, t2g):
    group_counts = defaultdict(float)
    with open(sf_path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            target_id = parts[0]
            num_reads = float(parts[4])
            if num_reads == 0:
                continue
            tid = target_id.split('::')[0]
            group = t2g.get(tid)
            if group is None:
                continue
            group_counts[group] += num_reads
    return dict(group_counts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cell-dirs', nargs='+', required=True)
    ap.add_argument('--cell-ids', nargs='+', required=True)
    ap.add_argument('--t2g', required=True,
                    help='Two-column TSV: transcript_id, group_id (no header)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.cell_dirs) != len(args.cell_ids):
        sys.exit('Number of --cell-dirs and --cell-ids must match')

    print('Loading t2g from {}'.format(args.t2g), file=sys.stderr)
    t2g = load_t2g(args.t2g)
    print('  {} transcript entries loaded'.format(len(t2g)), file=sys.stderr)

    cell_counts = {}
    all_features = set()
    for cell_id, cell_dir in zip(args.cell_ids, args.cell_dirs):
        sf_path = os.path.join(cell_dir, 'quant.sf')
        counts = read_quant_sf(sf_path, t2g)
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

    print('{} expressed features across {} cells'.format(len(sorted_features), len(sorted_cells)),
          file=sys.stderr)


if __name__ == '__main__':
    main()
