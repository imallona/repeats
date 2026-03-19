#!/usr/bin/env python3
"""
Convert bustools count output (MTX + barcodes + genes) to the standard feature x cell TSV.
"""

import argparse
import sys
import scipy.io


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtx', required=True, help='counts.mtx from bustools count')
    ap.add_argument('--barcodes', required=True, help='barcodes file from bustools count')
    ap.add_argument('--genes', required=True, help='genes/transcripts file from bustools count')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    with open(args.barcodes) as fh:
        barcodes = [line.rstrip('\n') for line in fh]
    with open(args.genes) as fh:
        genes = [line.rstrip('\n') for line in fh]

    mat = scipy.io.mmread(args.mtx).tocsr()
    print(f'{mat.shape[0]} features, {mat.shape[1]} cells', file=sys.stderr)

    with open(args.output, 'w') as fh:
        fh.write('feature_id\t' + '\t'.join(barcodes) + '\n')
        for feat_idx, gene_id in enumerate(genes):
            row = mat.getrow(feat_idx)
            if row.nnz == 0:
                continue
            values = row.toarray().flatten()
            fh.write(gene_id + '\t' + '\t'.join(str(int(v)) for v in values) + '\n')

    print(f'Written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
