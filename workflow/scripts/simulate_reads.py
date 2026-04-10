#!/usr/bin/env python3
"""
Simulate scRNA-seq reads from repeat elements for aligner benchmarking.

Modes:
  smartseq2 - one gzipped FASTQ per cell plus a STARsolo manifest
  chromium  - paired R1 (barcode+UMI) and R2 (cDNA) for all cells combined

Ground truth is written as a long-format TSV with columns:
  cell_id (or cell_barcode), locus_id, repeat_id, family_id, class_id, true_count

  locus_id  = transcript_id (finest granularity: bare transcript_id, e.g. AluSz6_dup1)
  repeat_id = gene_id (coarser: multiple loci may share the same repeat_id)

Memory design: GTF coordinates are loaded per-chromosome into memory as compact
tuples (no sequences). The FASTA is streamed one chromosome at a time and only
chromosomes containing sampled intervals are loaded.
"""

import argparse
import gzip
import math
import os
import random
import sys
from collections import defaultdict


def parse_gtf_repeats_by_chrom(gtf_path, allowed_chroms=None, max_per_chrom=None):
    intervals_by_chrom = defaultdict(list)
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            if allowed_chroms and chrom not in allowed_chroms:
                continue
            start_0 = int(start) - 1
            end_0 = int(end)
            if end_0 - start_0 < 50:
                continue
            gene_id = parse_gtf_attribute(attrs, 'gene_id') or f'{chrom}_{start_0}_{end_0}'
            # transcript_id is the canonical per-locus identifier used in the locus_map
            # and pseudo-genome FASTA headers (e.g. AluSz6_dup1)
            transcript_id = (parse_gtf_attribute(attrs, 'transcript_id')
                             or f'{gene_id}_{chrom}_{start_0}')
            family_id = parse_gtf_attribute(attrs, 'family_id') or 'unknown'
            class_id = parse_gtf_attribute(attrs, 'class_id') or 'unknown'
            locus_id = transcript_id
            intervals_by_chrom[chrom].append(
                (start_0, end_0, locus_id, gene_id, family_id, class_id, strand))

    if max_per_chrom:
        rng = random.Random(0)
        for chrom in intervals_by_chrom:
            if len(intervals_by_chrom[chrom]) > max_per_chrom:
                intervals_by_chrom[chrom] = rng.sample(intervals_by_chrom[chrom], max_per_chrom)

    return intervals_by_chrom


def parse_gtf_attribute(attrs, key):
    for part in attrs.split(';'):
        part = part.strip()
        if not part.startswith(key):
            continue
        remainder = part[len(key):].strip()
        if remainder.startswith('"'):
            return remainder.strip('"')
        tokens = remainder.split()
        if tokens:
            return tokens[0].strip('"')
    return None


def stream_fasta_by_chrom(fasta_path, needed_chroms):
    opener = gzip.open if fasta_path.endswith('.gz') else open
    current_chrom = None
    seq_parts = []
    loading = False
    with opener(fasta_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if loading:
                    yield current_chrom, ''.join(seq_parts)
                current_chrom = line[1:].split()[0]
                loading = current_chrom in needed_chroms
                seq_parts = []
            elif loading:
                seq_parts.append(line.upper())
    if loading:
        yield current_chrom, ''.join(seq_parts)


def extract_repeat_sequence(chrom_seq, start, end, strand):
    seq = chrom_seq[start:end]
    if len(seq) < 50:
        return None
    if seq.count('N') / len(seq) > 0.1:
        return None
    return reverse_complement(seq) if strand == '-' else seq


def reverse_complement(seq):
    table = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(table)[::-1]


def sample_subseq(seq, read_length, rng, mutation_rate=0.001):
    if len(seq) <= read_length:
        return (seq + 'N' * read_length)[:read_length]
    offset = rng.randint(0, len(seq) - read_length)
    read = seq[offset:offset + read_length]
    bases = list(read)
    for i in range(len(bases)):
        if rng.random() < mutation_rate:
            bases[i] = rng.choice('ACGT')
    return ''.join(bases)


def sample_count_geometric(rng, mean=5.0, max_count=50):
    p = 1.0 / (1.0 + mean)
    u = max(rng.random(), 1e-10)
    count = int(math.ceil(math.log(u) / math.log(1.0 - p)))
    return max(1, min(count, max_count))


def make_qual(length, char='F'):
    return char * length


def safe_id(locus_id):
    return locus_id.replace(' ', '_').replace('/', '_')


def build_cell_plan(intervals_by_chrom, n_cells, mean_expressed_per_cell, rng):
    all_intervals = [
        (chrom, start, end, locus_id, gene_id, family_id, class_id, strand)
        for chrom, items in intervals_by_chrom.items()
        for start, end, locus_id, gene_id, family_id, class_id, strand in items
    ]

    cell_plan = {}
    for cell_index in range(n_cells):
        cell_id = f'cell_{cell_index + 1:03d}'
        n_expressed = max(1, int(mean_expressed_per_cell * (0.7 + 0.6 * rng.random())))
        sampled = rng.sample(all_intervals, min(n_expressed, len(all_intervals)))
        # key = locus_id (unique per repeat instance); value includes gene_id for aggregation
        cell_plan[cell_id] = {
            locus_id: (count, gene_id, family_id, class_id, chrom, start, end, strand)
            for chrom, start, end, locus_id, gene_id, family_id, class_id, strand in sampled
            for count in [sample_count_geometric(rng)]
        }

    return cell_plan


def build_locus_to_cells(cell_plan):
    locus_to_cells = defaultdict(list)
    for cell_id, locus_map in cell_plan.items():
        for locus_id, (count, gene_id, family_id, class_id, chrom, start, end, strand) \
                in locus_map.items():
            locus_to_cells[locus_id].append((cell_id, count, gene_id, family_id, class_id))
    return locus_to_cells


def build_chrom_locus_coords(cell_plan):
    chrom_loci = defaultdict(dict)
    for locus_map in cell_plan.values():
        for locus_id, (count, gene_id, family_id, class_id, chrom, start, end, strand) \
                in locus_map.items():
            chrom_loci[chrom][locus_id] = (start, end, strand, gene_id)
    return chrom_loci


def simulate_smartseq2(fasta_path, cell_plan, read_length, outdir, rng, mutation_rate=0.001):
    os.makedirs(outdir, exist_ok=True)
    locus_to_cells = build_locus_to_cells(cell_plan)
    chrom_locus_coords = build_chrom_locus_coords(cell_plan)
    needed_chroms = set(chrom_locus_coords.keys())

    cell_fastq_paths = {cell_id: os.path.join(outdir, f'{cell_id}.fastq.gz') for cell_id in cell_plan}
    cell_read_counts = defaultdict(int)
    ground_truth = defaultdict(dict)

    cell_handles = {cell_id: gzip.open(path, 'wt') for cell_id, path in cell_fastq_paths.items()}

    try:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_path, needed_chroms):
            for locus_id, (start, end, strand, gene_id) in chrom_locus_coords[chrom].items():
                repeat_seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if repeat_seq is None:
                    continue
                for cell_id, count, gene_id, family_id, class_id in locus_to_cells[locus_id]:
                    fq = cell_handles[cell_id]
                    read_idx = cell_read_counts[cell_id]
                    for i in range(count):
                        read = sample_subseq(repeat_seq, read_length, rng, mutation_rate)
                        fq.write(f'@{cell_id}_r{read_idx + i}_{safe_id(locus_id)}\n'
                                 f'{read}\n+\n{make_qual(len(read))}\n')
                    cell_read_counts[cell_id] += count
                    ground_truth[cell_id][locus_id] = (count, gene_id, family_id, class_id)
    finally:
        for fh in cell_handles.values():
            fh.close()

    for cell_id, total_reads in cell_read_counts.items():
        print(f'  {cell_id}: {total_reads} reads, {len(ground_truth[cell_id])} loci', file=sys.stderr)

    return cell_fastq_paths, ground_truth


def simulate_chromium(fasta_path, cell_plan, read_length, barcode_length, umi_length, outdir, rng, mutation_rate=0.001):
    os.makedirs(outdir, exist_ok=True)
    locus_to_cells = build_locus_to_cells(cell_plan)
    chrom_locus_coords = build_chrom_locus_coords(cell_plan)
    needed_chroms = set(chrom_locus_coords.keys())

    cell_ids = sorted(set(cell_plan.keys()))
    cell_id_to_barcode = {}
    used_barcodes = set()
    for cid in cell_ids:
        bc = ''.join(rng.choices('ACGT', k=barcode_length))
        while bc in used_barcodes:
            bc = ''.join(rng.choices('ACGT', k=barcode_length))
        used_barcodes.add(bc)
        cell_id_to_barcode[cid] = bc

    barcodes_path = os.path.join(outdir, 'barcodes.txt')
    with open(barcodes_path, 'w') as fh:
        for cid in cell_ids:
            fh.write(cell_id_to_barcode[cid] + '\n')

    barcode_map_path = os.path.join(outdir, 'barcode_to_cell_id.tsv')
    with open(barcode_map_path, 'w') as fh:
        fh.write('barcode\tcell_id\n')
        for cid in cell_ids:
            fh.write(f'{cell_id_to_barcode[cid]}\t{cid}\n')

    r1_path = os.path.join(outdir, 'R1.fastq.gz')
    r2_path = os.path.join(outdir, 'R2.fastq.gz')
    ground_truth = defaultdict(dict)
    total_reads = 0

    with gzip.open(r1_path, 'wt') as r1f, gzip.open(r2_path, 'wt') as r2f:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_path, needed_chroms):
            for locus_id, (start, end, strand, gene_id) in chrom_locus_coords[chrom].items():
                repeat_seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if repeat_seq is None:
                    continue
                for cell_id, count, gene_id, family_id, class_id in locus_to_cells[locus_id]:
                    barcode = cell_id_to_barcode[cell_id]
                    used_umis = set()
                    for _ in range(count):
                        umi = ''.join(rng.choices('ACGT', k=umi_length))
                        while umi in used_umis:
                            umi = ''.join(rng.choices('ACGT', k=umi_length))
                        used_umis.add(umi)
                        read_id = f'r{total_reads}_{barcode[:8]}_{umi}_{safe_id(locus_id)}'
                        r1_seq = barcode + umi
                        r2_seq = sample_subseq(repeat_seq, read_length, rng, mutation_rate)
                        r1f.write(f'@{read_id}\n{r1_seq}\n+\n{make_qual(len(r1_seq))}\n')
                        r2f.write(f'@{read_id}\n{r2_seq}\n+\n{make_qual(len(r2_seq))}\n')
                        total_reads += 1
                    ground_truth[barcode][locus_id] = (count, gene_id, family_id, class_id)

    print(f'  {total_reads} read pairs across {len(cell_ids)} cells', file=sys.stderr)
    return (r1_path, r2_path), ground_truth


def write_ground_truth(ground_truth, output_path, cell_column):
    with open(output_path, 'w') as fh:
        fh.write(f'{cell_column}\tlocus_id\trepeat_id\tfamily_id\tclass_id\ttrue_count\n')
        for cell_id in sorted(ground_truth):
            for locus_id in sorted(ground_truth[cell_id]):
                count, gene_id, family_id, class_id = ground_truth[cell_id][locus_id]
                fh.write(f'{cell_id}\t{locus_id}\t{gene_id}\t{family_id}\t{class_id}\t{count}\n')


def write_smartseq2_manifest(cell_fastq_paths, manifest_path):
    with open(manifest_path, 'w') as fh:
        for cell_id in sorted(cell_fastq_paths):
            fq_path = cell_fastq_paths[cell_id]
            fh.write(f'{fq_path}\t-\t{cell_id}\n')


def write_repeat_fasta(cell_plan, fasta_path, chrom_locus_coords, fasta_source_path):
    needed_chroms = set(chrom_locus_coords.keys())
    all_locus_coords = {
        locus_id: (chrom, start, end, strand)
        for chrom, loci in chrom_locus_coords.items()
        for locus_id, (start, end, strand, gene_id) in loci.items()
    }
    written = set()
    with open(fasta_path, 'w') as out_fh:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_source_path, needed_chroms):
            for locus_id, (l_chrom, start, end, strand) in all_locus_coords.items():
                if l_chrom != chrom or locus_id in written:
                    continue
                seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if seq is not None:
                    out_fh.write(f'>{safe_id(locus_id)}\n{seq}\n')
                    written.add(locus_id)


def main():
    ap = argparse.ArgumentParser(description='Simulate scRNA-seq reads from repeat elements')
    ap.add_argument('--mode', choices=['single_end', 'chromium'], required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--n-cells', type=int, default=20)
    ap.add_argument('--n-expressed', type=int, default=100)
    ap.add_argument('--read-length', type=int, default=90)
    ap.add_argument('--cb-length', type=int, default=16)
    ap.add_argument('--umi-length', type=int, default=12)
    ap.add_argument('--chromosomes', nargs='+', default=None)
    ap.add_argument('--max-repeats-per-chrom', type=int, default=None,
                    help='Cap intervals per chromosome to limit memory use on full genomes')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--mutation-rate', type=float, default=0.001,
                    help='Per-base substitution probability for simulated reads (default 0.001)')
    args = ap.parse_args()

    rng = random.Random(args.seed)
    # Strip chr prefix so config values like "chr10" match bare-number GTF/FASTA chromosome names
    allowed_chroms = (
        {c.lstrip("chr") if c.startswith("chr") else c for c in args.chromosomes}
        if args.chromosomes else None
    )

    print(f'Parsing {args.gtf}', file=sys.stderr)
    intervals_by_chrom = parse_gtf_repeats_by_chrom(
        args.gtf, allowed_chroms, args.max_repeats_per_chrom)
    total_intervals = sum(len(v) for v in intervals_by_chrom.values())
    print(f'  {total_intervals} repeat intervals across {len(intervals_by_chrom)} chromosomes',
          file=sys.stderr)

    if total_intervals == 0:
        sys.exit('No repeat intervals found. Check --chromosomes and input paths.')

    print(f'Planning {args.n_cells} cell assignments', file=sys.stderr)
    cell_plan = build_cell_plan(intervals_by_chrom, args.n_cells, args.n_expressed, rng)

    os.makedirs(args.outdir, exist_ok=True)
    ground_truth_path = os.path.join(args.outdir, 'ground_truth.tsv')
    repeat_fasta_path = os.path.join(args.outdir, 'simulated_repeats.fa')
    chrom_locus_coords = build_chrom_locus_coords(cell_plan)

    if args.mode == 'single_end':
        print(f'Simulating single-end reads (streaming FASTA by chrom)', file=sys.stderr)
        cell_fastq_paths, ground_truth = simulate_smartseq2(
            args.fasta, cell_plan, args.read_length, args.outdir, rng, mutation_rate=args.mutation_rate)
        write_ground_truth(ground_truth, ground_truth_path, cell_column='cell_id')
        manifest_path = os.path.join(args.outdir, 'manifest.tsv')
        write_smartseq2_manifest(cell_fastq_paths, manifest_path)
        print(f'Manifest written to {manifest_path}', file=sys.stderr)

    else:
        print(f'Simulating Chromium reads (streaming FASTA by chrom)', file=sys.stderr)
        (r1, r2), ground_truth = simulate_chromium(
            args.fasta, cell_plan, args.read_length,
            args.cb_length, args.umi_length, args.outdir, rng, mutation_rate=args.mutation_rate)
        write_ground_truth(ground_truth, ground_truth_path, cell_column='cell_barcode')
        print(f'R1: {r1}', file=sys.stderr)
        print(f'R2: {r2}', file=sys.stderr)

    print(f'Writing repeat FASTA to {repeat_fasta_path}', file=sys.stderr)
    write_repeat_fasta(cell_plan, repeat_fasta_path, chrom_locus_coords, args.fasta)
    print(f'Ground truth written to {ground_truth_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
