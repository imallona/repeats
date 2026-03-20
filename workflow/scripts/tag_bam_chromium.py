#!/usr/bin/env python3
"""
Add CB (cell barcode) and UB (UMI) tags to a BAM by matching read names to R1 FASTQ.

R1 layout: <barcode><umi> (no adapter).  The barcode and UMI lengths are given
via --cb-length and --umi-length.
"""

import argparse
import gzip
import sys

import pysam


def parse_r1_barcodes(r1_path, cb_length, umi_length):
    """Return dict: read_name -> (barcode, umi)."""
    mapping = {}
    opener = gzip.open if r1_path.endswith('.gz') else open
    with opener(r1_path, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip('\n')
            fh.readline()  # +
            fh.readline()  # qual
            read_name = header[1:].split()[0]
            cb = seq[:cb_length]
            umi = seq[cb_length:cb_length + umi_length]
            mapping[read_name] = (cb, umi)
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--r1', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--cb-length', type=int, required=True)
    ap.add_argument('--umi-length', type=int, required=True)
    args = ap.parse_args()

    print('Parsing R1 barcodes...', file=sys.stderr)
    r1_map = parse_r1_barcodes(args.r1, args.cb_length, args.umi_length)
    print(f'  {len(r1_map)} read names loaded from R1', file=sys.stderr)

    tagged = 0
    untagged = 0
    with pysam.AlignmentFile(args.bam, 'rb') as inbam:
        with pysam.AlignmentFile(args.output, 'wb', header=inbam.header) as outbam:
            for read in inbam:
                key = read.query_name
                if key in r1_map:
                    cb, umi = r1_map[key]
                    read.set_tag('CB', cb, value_type='Z')
                    read.set_tag('UB', umi, value_type='Z')
                    tagged += 1
                else:
                    untagged += 1
                outbam.write(read)

    print(f'  Tagged {tagged} reads, {untagged} unmatched', file=sys.stderr)
    pysam.sort('-o', args.output + '.tmp', args.output)
    import shutil
    shutil.move(args.output + '.tmp', args.output)
    pysam.index(args.output)


if __name__ == '__main__':
    main()
