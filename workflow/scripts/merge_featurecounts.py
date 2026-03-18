#!/usr/bin/env python3
"""
Merge featureCounts chunk TSVs into a single feature x cell matrix.

For SmartSeq2: chunk filenames contain cell{N}_feat{M}; rows are stacked
(feature axis from feat chunks) and columns are concatenated (cell axis from
cell chunks).  For Chromium (feat-only chunks) the cell axis is already
captured by featureCounts --byReadGroup, so we just stack row-wise.
"""
import argparse
import re
import sys
from collections import defaultdict


def read_tsv(path):
    rows = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if parts:
                rows[parts[0]] = parts[1:]
    return hdr[1:], rows


def parse_chunk_ids(path):
    """Return (cell_chunk_id, feat_chunk_id) from filename, or (0, feat_id)."""
    m = re.search(r'cell(\d+)_feat(\d+)', path)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.search(r'feat(\d+)', path)
    if m:
        return 0, int(m.group(1))
    return 0, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chunks', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    # group by cell_chunk; within each group merge feat chunks vertically
    cell_groups = defaultdict(dict)
    for path in args.chunks:
        cc, fc = parse_chunk_ids(path)
        cells, rows = read_tsv(path)
        cell_groups[cc][fc] = (cells, rows)

    merged_groups = {}
    for cc in sorted(cell_groups):
        all_cells = []
        all_rows = {}
        for fc in sorted(cell_groups[cc]):
            cells, rows = cell_groups[cc][fc]
            if not all_cells:
                all_cells = cells
            all_rows.update(rows)
        merged_groups[cc] = (all_cells, all_rows)

    all_features = list(dict.fromkeys(
        f for _, rows in merged_groups.values() for f in rows))
    all_cells = [c for cells, _ in merged_groups.values() for c in cells]

    with open(args.output, 'w') as out:
        out.write('feature_id\t' + '\t'.join(all_cells) + '\n')
        for feat in all_features:
            row_vals = []
            for cells, rows in merged_groups.values():
                row_vals.extend(rows.get(feat, ['0'] * len(cells)))
            out.write(feat + '\t' + '\t'.join(row_vals) + '\n')


if __name__ == '__main__':
    main()
