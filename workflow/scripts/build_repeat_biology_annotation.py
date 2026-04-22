#!/usr/bin/env python3
"""
Build per-subfamily biology annotation TSV for repeat loci.

For each RepeatMasker subfamily (gene_id == repName), aggregate across its loci:
  - n_loci, median milliDiv
  - family_age from subfamily name (classify_family_age)
  - family_age from median milliDiv (classify_family_age_from_divergence)
  - strand composition against overlapping gene annotation:
    frac_sense_gene, frac_antisense_gene, frac_both, frac_intergenic
    plus a dominant strand_vs_gene label ("sense" / "antisense" / "both" /
    "intergenic") taken from the most frequent per-locus class
  - frac_downstream_a_run: fraction of loci whose downstream genomic region
    (in mRNA-sense orientation) contains an A-run at least a_run_min_length long

Inputs:
  --rmsk           UCSC rmsk.txt.gz flatfile (needed for milliDiv)
  --genes-gtf      Decompressed genes GTF, chr-stripped (Ensembl layout)
  --genome-fasta   Decompressed genome FASTA, chr-stripped (Ensembl layout)
  --output         TSV output path

Pure classification logic lives in biology_annotations.py; this script only
does I/O orchestration (FASTA random access, merged gene-interval lookup,
per-subfamily aggregation). All helpers are stdlib-only so no new conda
dependency is introduced.
"""

import argparse
import bisect
import gzip
import statistics
import sys
from collections import defaultdict

from biology_annotations import (
    classify_family_age,
    classify_family_age_from_divergence,
    strand_vs_gene,
    longest_a_run,
)


def normalize_chrom(name):
    return name[3:] if name.startswith('chr') else name


def build_fasta_index(path):
    """Return dict chrom -> (length, seq_start_byte, line_bases, line_width).

    Assumes each sequence has a consistent line width (except the final line),
    matching standard Ensembl/UCSC FASTA conventions.
    """
    index = {}
    current = None
    seq_start = 0
    seq_length = 0
    line_bases = 0
    line_width = 0
    offset = 0
    with open(path, 'rb') as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            line_start = offset
            offset += len(line)
            if line.startswith(b'>'):
                if current is not None:
                    index[current] = (seq_length, seq_start, line_bases, line_width)
                current = line[1:].split()[0].decode()
                seq_start = offset
                seq_length = 0
                line_bases = 0
                line_width = 0
            else:
                stripped = line.rstrip(b'\r\n')
                if line_bases == 0:
                    line_bases = len(stripped)
                    line_width = len(line)
                seq_length += len(stripped)
    if current is not None:
        index[current] = (seq_length, seq_start, line_bases, line_width)
    return index


class FastaRandomAccess:
    """Random access into a plain FASTA using an in-memory index.

    Coordinates are 0-based half-open. Out-of-range queries return '' rather
    than raising, so the caller can treat unknown chromosomes as no-hit.
    """

    def __init__(self, path, index=None):
        self.path = path
        self.index = index if index is not None else build_fasta_index(path)
        self._fh = open(path, 'rb')

    def close(self):
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def fetch(self, chrom, start, end):
        meta = self.index.get(chrom)
        if meta is None:
            return ''
        length, seq_start, line_bases, line_width = meta
        if start < 0:
            start = 0
        if end > length:
            end = length
        if start >= end:
            return ''
        start_byte = seq_start + (start // line_bases) * line_width + (start % line_bases)
        end_byte = seq_start + (end // line_bases) * line_width + (end % line_bases)
        self._fh.seek(start_byte)
        raw = self._fh.read(end_byte - start_byte)
        return raw.replace(b'\n', b'').replace(b'\r', b'').decode('ascii').upper()


_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def downstream_sequence(fasta, chrom, start, end, strand, length):
    """Return the genomic region downstream of a repeat in mRNA-sense orientation.

    Plus-strand repeat: [end, end + length).
    Minus-strand repeat: reverse complement of [start - length, start).
    """
    if strand == '+':
        return fasta.fetch(chrom, end, end + length)
    flank = fasta.fetch(chrom, max(start - length, 0), start)
    return reverse_complement(flank)


class GeneIntervals:
    """Merged sorted gene intervals keyed by (chrom, strand) for overlap queries.

    Feed with add(chrom, start, end, strand), then call finalize() once to
    merge and sort per key. any_overlap(chrom, start, end, strand) returns a
    bool in O(log n).
    """

    def __init__(self):
        self._raw = defaultdict(list)
        self._starts = {}
        self._ends = {}

    def add(self, chrom, start, end, strand):
        self._raw[(chrom, strand)].append((start, end))

    def finalize(self):
        for key, ivs in self._raw.items():
            ivs.sort()
            merged = []
            cur_s, cur_e = ivs[0]
            for s, e in ivs[1:]:
                if s <= cur_e:
                    if e > cur_e:
                        cur_e = e
                else:
                    merged.append((cur_s, cur_e))
                    cur_s, cur_e = s, e
            merged.append((cur_s, cur_e))
            self._starts[key] = [s for s, _ in merged]
            self._ends[key] = [e for _, e in merged]
        self._raw = None

    def any_overlap(self, chrom, start, end, strand):
        key = (chrom, strand)
        starts = self._starts.get(key)
        if not starts:
            return False
        idx = bisect.bisect_right(starts, end - 1)
        if idx == 0:
            return False
        return self._ends[key][idx - 1] > start


def load_gene_intervals(genes_gtf_path):
    intervals = GeneIntervals()
    opener = gzip.open if genes_gtf_path.endswith('.gz') else open
    with opener(genes_gtf_path, 'rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8 or fields[2] != 'gene':
                continue
            chrom = normalize_chrom(fields[0])
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            if strand not in ('+', '-'):
                continue
            intervals.add(chrom, start, end, strand)
    intervals.finalize()
    return intervals


def iter_rmsk(path):
    """Yield (chrom, start, end, strand, milli_div, rep_name, rep_class, rep_family)."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 13:
                continue
            try:
                milli_div = int(fields[2])
            except ValueError:
                continue
            yield (
                normalize_chrom(fields[5]),
                int(fields[6]),
                int(fields[7]),
                fields[9],
                milli_div,
                fields[10],
                fields[11],
                fields[12],
            )


OUTPUT_COLUMNS = [
    'gene_id', 'family_id', 'class_id', 'n_loci',
    'median_millidiv', 'family_age', 'family_age_from_millidiv',
    'frac_sense_gene', 'frac_antisense_gene', 'frac_both', 'frac_intergenic',
    'strand_vs_gene', 'frac_downstream_a_run',
]


def aggregate(rmsk_path, gene_intervals, fasta, downstream_bp, a_run_min_length,
              allowed_chroms=None, min_length=50, log_fh=None):
    """Stream rmsk, compute per-locus biology, return per-subfamily aggregates.

    Returns a dict keyed by rep_name with accumulated counts ready to be
    divided into fractions at write time. Kept separate from write_annotation
    so unit tests can assert on the aggregate without touching disk.
    """
    per_subfamily = {}
    skipped_chrom = 0
    skipped_short = 0
    processed = 0

    for rec in iter_rmsk(rmsk_path):
        chrom, start, end, strand, milli_div, rep_name, rep_class, rep_family = rec
        if allowed_chroms and chrom not in allowed_chroms:
            skipped_chrom += 1
            continue
        if end - start < min_length:
            skipped_short += 1
            continue
        if chrom not in fasta.index:
            skipped_chrom += 1
            continue

        overlap_strands = []
        if gene_intervals.any_overlap(chrom, start, end, '+'):
            overlap_strands.append('+')
        if gene_intervals.any_overlap(chrom, start, end, '-'):
            overlap_strands.append('-')
        svg = strand_vs_gene(strand, overlap_strands)

        ds = downstream_sequence(fasta, chrom, start, end, strand, downstream_bp)
        a_hit = longest_a_run(ds) >= a_run_min_length

        entry = per_subfamily.get(rep_name)
        if entry is None:
            entry = {
                'family_id': rep_family,
                'class_id': rep_class,
                'n_loci': 0,
                'milli': [],
                'n_sense': 0,
                'n_antisense': 0,
                'n_both': 0,
                'n_intergenic': 0,
                'n_a_run_hit': 0,
            }
            per_subfamily[rep_name] = entry
        entry['n_loci'] += 1
        entry['milli'].append(milli_div)
        if svg == 'sense':
            entry['n_sense'] += 1
        elif svg == 'antisense':
            entry['n_antisense'] += 1
        elif svg == 'both':
            entry['n_both'] += 1
        else:
            entry['n_intergenic'] += 1
        if a_hit:
            entry['n_a_run_hit'] += 1
        processed += 1

    if log_fh is not None:
        log_fh.write(
            f'Processed {processed} loci; skipped {skipped_chrom} by chrom, '
            f'{skipped_short} by length\n'
        )
    return per_subfamily


def write_annotation(per_subfamily, output_path):
    with open(output_path, 'w') as out:
        out.write('\t'.join(OUTPUT_COLUMNS) + '\n')
        for rep_name in sorted(per_subfamily):
            e = per_subfamily[rep_name]
            n = e['n_loci']
            median_milli = statistics.median(e['milli']) if e['milli'] else None
            age_name = classify_family_age(rep_name)
            age_milli = (classify_family_age_from_divergence(median_milli)
                         if median_milli is not None else 'unknown')
            counts = {
                'sense': e['n_sense'],
                'antisense': e['n_antisense'],
                'both': e['n_both'],
                'intergenic': e['n_intergenic'],
            }
            dominant = max(counts, key=counts.get)
            out.write('\t'.join([
                rep_name,
                e['family_id'] or '',
                e['class_id'] or '',
                str(n),
                f'{median_milli:.3f}' if median_milli is not None else '',
                age_name,
                age_milli,
                f'{e["n_sense"] / n:.4f}',
                f'{e["n_antisense"] / n:.4f}',
                f'{e["n_both"] / n:.4f}',
                f'{e["n_intergenic"] / n:.4f}',
                dominant,
                f'{e["n_a_run_hit"] / n:.4f}',
            ]) + '\n')


def main():
    ap = argparse.ArgumentParser(description='Build per-subfamily repeat biology annotation TSV')
    ap.add_argument('--rmsk', required=True)
    ap.add_argument('--genes-gtf', required=True)
    ap.add_argument('--genome-fasta', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--downstream-bp', type=int, default=50,
                    help='Length of the downstream region scanned for internal A-runs')
    ap.add_argument('--a-run-min-length', type=int, default=12,
                    help='Minimum contiguous A count that counts as an A-run hit')
    ap.add_argument('--chromosomes', nargs='+', default=None,
                    help='Optional chromosome allow-list (chr prefix optional)')
    ap.add_argument('--min-length', type=int, default=50,
                    help='Skip rmsk records shorter than this many bases')
    args = ap.parse_args()

    allowed = None
    if args.chromosomes:
        allowed = {normalize_chrom(c) for c in args.chromosomes}

    sys.stderr.write(f'Indexing {args.genome_fasta}\n')
    with FastaRandomAccess(args.genome_fasta) as fasta:
        sys.stderr.write(f'Loading gene intervals from {args.genes_gtf}\n')
        gene_intervals = load_gene_intervals(args.genes_gtf)
        sys.stderr.write(f'Streaming {args.rmsk}\n')
        per_subfamily = aggregate(
            rmsk_path=args.rmsk,
            gene_intervals=gene_intervals,
            fasta=fasta,
            downstream_bp=args.downstream_bp,
            a_run_min_length=args.a_run_min_length,
            allowed_chroms=allowed,
            min_length=args.min_length,
            log_fh=sys.stderr,
        )
    sys.stderr.write(f'Writing {args.output}\n')
    write_annotation(per_subfamily, args.output)


if __name__ == '__main__':
    main()
