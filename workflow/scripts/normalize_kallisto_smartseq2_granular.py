#!/usr/bin/env python3
"""
Merge per-cell kallisto abundance.tsv files into a feature x cell TSV,
aggregating transcript-level est_counts by a t2g mapping (for repeat granularity).

The t2g TSV has two columns (no header): transcript_id, group_id.
For granularity=locus, this is transcript_id -> transcript_id (identity).
For granularity=gene_id/family_id/class_id, this maps to the higher-level group.

Each abundance.tsv has columns: target_id, length, eff_length, est_counts, tpm.
Only est_counts is used. Counts for all transcripts mapping to the same group are summed.
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


def read_abundance(tsv_path, t2g):
    """
    Read abundance.tsv, aggregate est_counts by group via t2g.
    Returns dict: group_id -> aggregated_count.
    """
    group_counts = defaultdict(float)
    with open(tsv_path) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            target_id = parts[0]
            est_counts = float(parts[3])
            if est_counts == 0:
                continue
            tid = target_id.split('::')[0]
            group = t2g.get(tid)
            if group is None:
                continue
            group_counts[group] += est_counts
    return dict(group_counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cell-dirs', nargs='+', required=True,
                    help='Directories containing abundance.tsv, one per cell')
    ap.add_argument('--cell-ids', nargs='+', required=True,
                    help='Cell IDs in the same order as --cell-dirs')
    ap.add_argument('--t2g', required=True,
                    help='Two-column TSV: transcript_id, group_id (no header)')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.cell_dirs) != len(args.cell_ids):
        sys.exit('Number of --cell-dirs and --cell-ids must match')

    print(f'Loading t2g from {args.t2g}', file=sys.stderr)
    t2g = load_t2g(args.t2g)
    print(f'  {len(t2g)} transcript entries loaded', file=sys.stderr)

    cell_counts = {}
    all_features = set()
    for cell_id, cell_dir in zip(args.cell_ids, args.cell_dirs):
        abundance_path = os.path.join(cell_dir, 'abundance.tsv')
        counts = read_abundance(abundance_path, t2g)
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
