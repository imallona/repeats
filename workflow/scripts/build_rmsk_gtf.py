#!/usr/bin/env python3
"""
Convert a UCSC rmsk.txt.gz flat file (or per-chromosome MySQL query output)
to a TEtranscripts-compatible GTF.

UCSC rmsk.txt.gz columns (0-based index):
  0  bin
  1  swScore
  2  milliDiv
  3  milliDel
  4  milliIns
  5  genoName   (chromosome, e.g. chr1)
  6  genoStart  (0-based)
  7  genoEnd
  8  genoLeft
  9  strand     (+ or -)
  10 repName    (specific element name, e.g. L1PA2)
  11 repClass   (LINE, SINE, LTR, DNA, ...)
  12 repFamily  (L1, Alu, ERV1, ...)
  13 repStart
  14 repEnd
  15 repLeft
  16 id

Output GTF attributes match TEtranscripts format:
  gene_id "repName"; transcript_id "repName_dupN"; family_id "repFamily"; class_id "repClass";

The file is processed line by line so only one line is in memory at a time.
"""

import argparse
import gzip
import sys
from collections import defaultdict


def make_gtf_attributes(rep_name, locus_index, rep_family, rep_class):
    transcript_id = f'{rep_name}_dup{locus_index}'
    return (f'gene_id "{rep_name}"; transcript_id "{transcript_id}"; '
            f'family_id "{rep_family}"; class_id "{rep_class}";')


def convert_rmsk_to_gtf(input_path, output_path, allowed_chroms=None, min_length=50):
    opener = gzip.open if input_path.endswith('.gz') else open
    out_opener = gzip.open if output_path.endswith('.gz') else open

    locus_counts = defaultdict(int)
    written = 0
    skipped_chrom = 0
    skipped_short = 0

    with opener(input_path, 'rt') as in_fh, out_opener(output_path, 'wt') as out_fh:
        for line in in_fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 13:
                continue

            chrom = fields[5]
            start_0 = int(fields[6])
            end = int(fields[7])
            strand = fields[9]
            rep_name = fields[10]
            rep_class = fields[11]
            rep_family = fields[12]

            if allowed_chroms and chrom not in allowed_chroms:
                skipped_chrom += 1
                continue

            if end - start_0 < min_length:
                skipped_short += 1
                continue

            locus_counts[rep_name] += 1
            attrs = make_gtf_attributes(rep_name, locus_counts[rep_name], rep_family, rep_class)
            gtf_start = start_0 + 1

            out_fh.write(
                f'{chrom}\trmsk\texon\t{gtf_start}\t{end}\t.\t{strand}\t.\t{attrs}\n'
            )
            written += 1

    print(f'Written {written} GTF records', file=sys.stderr)
    if skipped_chrom:
        print(f'Skipped {skipped_chrom} records not in allowed chromosomes', file=sys.stderr)
    if skipped_short:
        print(f'Skipped {skipped_short} records shorter than {min_length} bp', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description='Convert UCSC rmsk.txt.gz to TEtranscripts GTF')
    ap.add_argument('--input', required=True, help='rmsk.txt.gz from UCSC or MySQL TSV output')
    ap.add_argument('--output', required=True, help='Output GTF (use .gtf.gz to compress)')
    ap.add_argument('--chromosomes', nargs='+', default=None,
                    help='Keep only these chromosomes (e.g. chr1 chr2)')
    ap.add_argument('--min-length', type=int, default=50)
    args = ap.parse_args()

    allowed = set(args.chromosomes) if args.chromosomes else None
    print(f'Converting {args.input}', file=sys.stderr)
    convert_rmsk_to_gtf(args.input, args.output, allowed, args.min_length)
    print(f'Output written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
