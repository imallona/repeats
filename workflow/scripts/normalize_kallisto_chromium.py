#!/usr/bin/env python3
"""
Convert bustools count output (MTX + barcodes + genes) to feature x cell TSV.

The MTX from bustools count is (cells x loci); it is transposed to (loci x cells)
before output.

For repeat feature sets, supports granularity aggregation via --locus-map
(4-col TSV: transcript_id, gene_id, family_id, class_id).  The genes file from
bustools count contains locus IDs (= gene_id column of locus_map); these are
aggregated to the requested --granularity by summing counts across loci that
share the same group.

For the genes feature set (no --locus-map), counts are written as-is at locus
level.  granularity=gene_id is always a no-op (identity mapping).
"""

import argparse
import sys
from collections import defaultdict

import scipy.io


def load_locus_to_group(locus_map_path, granularity):
    """Return dict: locus_id (gene_id col) -> group_id at the requested granularity.

    locus_map columns (no header): transcript_id, gene_id, family_id, class_id
    The bustools count genes file contains gene_id values (col 2, index 1).
    """
    col_idx = {'locus': 0, 'gene_id': 1, 'family_id': 2, 'class_id': 3}
    if granularity not in col_idx:
        sys.exit(f'Unknown granularity: {granularity}')
    target_col = col_idx[granularity]
    mapping = {}
    with open(locus_map_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= target_col:
                continue
            locus = parts[0]      # locus_id (col 0) as used in bustools counts.genes.txt
            target = parts[target_col]
            if locus and target:
                mapping[locus] = target
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mtx',      required=True, help='counts.mtx from bustools count')
    ap.add_argument('--barcodes', required=True, help='counts.barcodes.txt')
    ap.add_argument('--genes',    required=True, help='counts.genes.txt (locus IDs)')
    ap.add_argument('--output',   required=True)
    ap.add_argument('--locus-map', default=None,
                    help='4-col locus map TSV for granularity aggregation (repeats only)')
    ap.add_argument('--granularity', default='gene_id',
                    choices=['locus', 'gene_id', 'family_id', 'class_id'],
                    help='Aggregation granularity (default: gene_id)')
    args = ap.parse_args()

    with open(args.barcodes) as fh:
        barcodes = [line.rstrip('\n') for line in fh]
    with open(args.genes) as fh:
        loci = [line.rstrip('\n') for line in fh]

    # bustools count MTX: rows=cells, cols=loci, then transpose to loci x cells
    mat = scipy.io.mmread(args.mtx).T.tocsr()
    n_features, n_cells = mat.shape
    print(f'{n_features} loci, {n_cells} cells', file=sys.stderr)

    # Build locus -> target mapping
    if args.locus_map:
        locus_to_group = load_locus_to_group(args.locus_map, args.granularity)
        print(f'Aggregating to {args.granularity} ({len(locus_to_group)} loci mapped)',
              file=sys.stderr)
    else:
        locus_to_group = None  # identity: locus IS the target group

    # Accumulate per-group counts (sparse: only store non-zero cell indices)
    # sparse accumulator: {cell_index: count} per group avoids allocating
    # a dense n_cells list for every feature group (saves O(groups x cells) memory).
    group_counts = defaultdict(lambda: defaultdict(float))
    for feat_idx, locus_id in enumerate(loci):
        row = mat.getrow(feat_idx)
        if row.nnz == 0:
            continue
        target = locus_to_group.get(locus_id) if locus_to_group is not None else locus_id
        if target is None:
            continue
        dest = group_counts[target]
        for cell_idx, v in zip(row.indices, row.data):
            if v:
                dest[cell_idx] += v

    sorted_groups = sorted(group_counts)
    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(barcodes) + '\n')
        for group in sorted_groups:
            cell_vals = group_counts[group]
            if not cell_vals:
                continue
            vals = [str(int(cell_vals.get(ci, 0))) for ci in range(n_cells)]
            fh.write(group + '\t' + '\t'.join(vals) + '\n')

    print(f'{len(sorted_groups)} features written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
