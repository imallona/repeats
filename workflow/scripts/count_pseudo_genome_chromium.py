#!/usr/bin/env python3
"""
Count reads per (cell barcode, repeat feature) from a CB-tagged bowtie2 BAM.

Reads must have a CB:Z tag set by tag_bam_chromium.py.  The script streams
the BAM once via samtools view, accumulates a sparse (cb, feature) -> count
dict, and writes a feature x cell TSV.

Reference names in the pseudo-genome BAM have the form:
    transcript_id::chrom:start-end(strand)

The transcript_id prefix is looked up in the locus map to resolve the
requested granularity (gene_id, family_id, class_id, or locus).

Memory: O(n_expressed_cb_feature_pairs); the full matrix is never allocated.
"""

import argparse
import subprocess
import sys
from collections import defaultdict


def load_locus_map(locus_map_path, granularity):
    col_index = {'locus': 0, 'gene_id': 1, 'family_id': 2, 'class_id': 3}
    if granularity not in col_index:
        sys.exit(f'unknown granularity: {granularity}')
    grp_col = col_index[granularity]
    mapping = {}
    with open(locus_map_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            tid = parts[0]
            grp = parts[grp_col]
            if tid and grp:
                mapping[tid] = grp
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam',       required=True, help='CB-tagged BAM (dedup.bam)')
    ap.add_argument('--locus-map', required=True, help='4-col locus map TSV (no header)')
    ap.add_argument('--granularity', required=True,
                    choices=['locus', 'gene_id', 'family_id', 'class_id'])
    ap.add_argument('--output',    required=True)
    args = ap.parse_args()

    print(f'loading locus map ({args.granularity})', file=sys.stderr)
    tid_to_group = load_locus_map(args.locus_map, args.granularity)
    print(f'  {len(tid_to_group)} loci loaded', file=sys.stderr)

    # sparse accumulator: cb -> feature -> count
    counts = defaultdict(lambda: defaultdict(int))
    all_features = set()
    all_cbs = set()
    n_reads = 0
    n_counted = 0

    proc = subprocess.Popen(
        ['samtools', 'view', args.bam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    for raw in proc.stdout:
        line = raw.decode()
        if line.startswith('@'):
            continue
        fields = line.split('\t')
        if len(fields) < 12:
            continue
        ref = fields[2]
        if ref == '*':
            continue
        n_reads += 1
        cb = None
        for tag in fields[11:]:
            if tag.startswith('CB:Z:'):
                cb = tag[5:].rstrip('\n')
                break
        if cb is None:
            continue
        tid = ref.split('::')[0]
        feat = tid_to_group.get(tid)
        if feat is None:
            continue
        counts[cb][feat] += 1
        all_features.add(feat)
        all_cbs.add(cb)
        n_counted += 1

    proc.wait()
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr.read().decode())
        sys.exit(f'samtools view failed (exit {proc.returncode})')

    print(f'  {n_reads} reads processed, {n_counted} counted', file=sys.stderr)
    print(f'  {len(all_cbs)} cells, {len(all_features)} features', file=sys.stderr)

    sorted_features = sorted(all_features)
    sorted_cbs = sorted(all_cbs)

    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(sorted_cbs) + '\n')
        for feat in sorted_features:
            row = [str(counts[cb].get(feat, 0)) for cb in sorted_cbs]
            fh.write(feat + '\t' + '\t'.join(row) + '\n')

    print(f'wrote {len(sorted_features)} features x {len(sorted_cbs)} cells to {args.output}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
