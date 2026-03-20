#!/usr/bin/env python3
"""
Count reads per repeat feature from bowtie2 pseudo-genome BAMs.

For each BAM, runs samtools idxstats and parses the reference names.
Reference names have the form: transcript_id::chrom:start-end(strand)
Splits on '::' to extract transcript_id, then aggregates by the
requested granularity level using the locus map TSV.

Locus map TSV columns (tab-separated, no header):
    transcript_id  gene_id  family_id  class_id

Outputs a feature x cell count TSV (rows=features, cols=cells, first col=feature_id).
"""

import argparse
import subprocess
import sys
from collections import defaultdict


def load_locus_map(locus_map_path, granularity):
    """
    Returns dict: transcript_id -> group_id based on granularity.
    For 'locus': transcript_id -> transcript_id (identity).
    For 'gene_id': transcript_id -> gene_id
    For 'family_id': transcript_id -> family_id
    For 'class_id': transcript_id -> class_id
    """
    col_index = {'locus': 0, 'gene_id': 1, 'family_id': 2, 'class_id': 3}
    if granularity not in col_index:
        sys.exit(f'Unknown granularity: {granularity}. Choose from: locus, gene_id, family_id, class_id')

    tid_col = 0
    grp_col = col_index[granularity]

    mapping = {}
    with open(locus_map_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            tid = parts[tid_col]
            grp = parts[grp_col]
            mapping[tid] = grp
    return mapping


def count_bam(bam_path, tid_to_group):
    """
    Run samtools idxstats on bam_path and aggregate mapped reads by group.
    Returns dict: group_id -> count (int).
    """
    result = subprocess.run(
        ['samtools', 'idxstats', bam_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    group_counts = defaultdict(int)
    for line in result.stdout.decode().splitlines():
        parts = line.split('\t')
        if len(parts) < 4:
            continue
        ref_name = parts[0]
        if ref_name == '*':
            continue
        # mapped reads = column 3 (index 2)
        try:
            mapped = int(parts[2])
        except ValueError:
            continue
        if mapped == 0:
            continue
        # ref_name format: transcript_id::chrom:start-end(strand)
        tid = ref_name.split('::')[0]
        group = tid_to_group.get(tid)
        if group is None:
            continue
        group_counts[group] += mapped
    return dict(group_counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bams', nargs='+', required=True,
                    help='BAM files, one per cell')
    ap.add_argument('--cell-ids', nargs='+', required=True,
                    help='Cell IDs in the same order as --bams')
    ap.add_argument('--locus-map', required=True,
                    help='TSV: transcript_id, gene_id, family_id, class_id (no header)')
    ap.add_argument('--granularity', required=True,
                    choices=['locus', 'gene_id', 'family_id', 'class_id'])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    if len(args.bams) != len(args.cell_ids):
        sys.exit('Number of --bams and --cell-ids must match')

    print(f'Loading locus map from {args.locus_map}', file=sys.stderr)
    tid_to_group = load_locus_map(args.locus_map, args.granularity)
    print(f'  {len(tid_to_group)} transcript entries loaded', file=sys.stderr)

    cell_counts = {}
    all_features = set()
    for cell_id, bam_path in zip(args.cell_ids, args.bams):
        print(f'Counting {cell_id} from {bam_path}', file=sys.stderr)
        counts = count_bam(bam_path, tid_to_group)
        cell_counts[cell_id] = counts
        all_features.update(counts.keys())

    sorted_features = sorted(all_features)
    sorted_cells = args.cell_ids

    print(f'Writing {len(sorted_features)} features x {len(sorted_cells)} cells to {args.output}',
          file=sys.stderr)
    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(sorted_cells) + '\n')
        for feat in sorted_features:
            row = [cell_counts[c].get(feat, 0) for c in sorted_cells]
            fh.write(feat + '\t' + '\t'.join(str(v) for v in row) + '\n')

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()
