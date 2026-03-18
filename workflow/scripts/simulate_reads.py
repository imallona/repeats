#!/usr/bin/env python3
"""
Simulate scRNA-seq reads from repeat elements for aligner benchmarking.

Modes:
  smartseq2 - one gzipped FASTQ per cell plus a STARsolo manifest
  chromium  - paired R1 (barcode+UMI) and R2 (cDNA) for all cells combined

Ground truth is written as a long-format TSV with columns:
  cell_id (or cell_barcode), repeat_id, family_id, class_id, true_count

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
            family_id = parse_gtf_attribute(attrs, 'family_id') or 'unknown'
            class_id = parse_gtf_attribute(attrs, 'class_id') or 'unknown'
            intervals_by_chrom[chrom].append((start_0, end_0, gene_id, family_id, class_id, strand))

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


def sample_subseq(seq, read_length, rng):
    if len(seq) <= read_length:
        return (seq + 'N' * read_length)[:read_length]
    offset = rng.randint(0, len(seq) - read_length)
    read = seq[offset:offset + read_length]
    bases = list(read)
    for i in range(len(bases)):
        if rng.random() < 0.001:
            bases[i] = rng.choice('ACGT')
    return ''.join(bases)


def sample_count_geometric(rng, mean=5.0, max_count=50):
    p = 1.0 / (1.0 + mean)
    u = max(rng.random(), 1e-10)
    count = int(math.ceil(math.log(u) / math.log(1.0 - p)))
    return max(1, min(count, max_count))


def make_qual(length, char='F'):
    return char * length


def safe_id(gene_id):
    return gene_id.replace(' ', '_').replace('/', '_')


def build_cell_plan(intervals_by_chrom, n_cells, mean_expressed_per_cell, rng):
    all_intervals = [
        (chrom, start, end, gene_id, family_id, class_id, strand)
        for chrom, items in intervals_by_chrom.items()
        for start, end, gene_id, family_id, class_id, strand in items
    ]

    cell_plan = {}
    for cell_index in range(n_cells):
        cell_id = f'cell_{cell_index + 1:03d}'
        n_expressed = max(1, int(mean_expressed_per_cell * (0.7 + 0.6 * rng.random())))
        sampled = rng.sample(all_intervals, min(n_expressed, len(all_intervals)))
        cell_plan[cell_id] = {
            gene_id: (count, family_id, class_id, chrom, start, end, strand)
            for chrom, start, end, gene_id, family_id, class_id, strand in sampled
            for count in [sample_count_geometric(rng)]
        }

    return cell_plan


def build_gene_to_cells(cell_plan):
    gene_to_cells = defaultdict(list)
    for cell_id, gene_map in cell_plan.items():
        for gene_id, (count, family_id, class_id, chrom, start, end, strand) in gene_map.items():
            gene_to_cells[gene_id].append((cell_id, count, family_id, class_id))
    return gene_to_cells


def build_chrom_gene_coords(cell_plan):
    chrom_genes = defaultdict(dict)
    for gene_map in cell_plan.values():
        for gene_id, (count, family_id, class_id, chrom, start, end, strand) in gene_map.items():
            chrom_genes[chrom][gene_id] = (start, end, strand)
    return chrom_genes


def simulate_smartseq2(fasta_path, cell_plan, read_length, outdir, rng):
    os.makedirs(outdir, exist_ok=True)
    gene_to_cells = build_gene_to_cells(cell_plan)
    chrom_gene_coords = build_chrom_gene_coords(cell_plan)
    needed_chroms = set(chrom_gene_coords.keys())

    cell_fastq_paths = {cell_id: os.path.join(outdir, f'{cell_id}.fastq.gz') for cell_id in cell_plan}
    cell_read_counts = defaultdict(int)
    ground_truth = defaultdict(dict)

    cell_handles = {cell_id: gzip.open(path, 'wt') for cell_id, path in cell_fastq_paths.items()}

    try:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_path, needed_chroms):
            for gene_id, (start, end, strand) in chrom_gene_coords[chrom].items():
                repeat_seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if repeat_seq is None:
                    continue
                for cell_id, count, family_id, class_id in gene_to_cells[gene_id]:
                    fq = cell_handles[cell_id]
                    read_idx = cell_read_counts[cell_id]
                    for i in range(count):
                        read = sample_subseq(repeat_seq, read_length, rng)
                        fq.write(f'@{cell_id}_r{read_idx + i}_{safe_id(gene_id)}\n'
                                 f'{read}\n+\n{make_qual(len(read))}\n')
                    cell_read_counts[cell_id] += count
                    ground_truth[cell_id][gene_id] = (count, family_id, class_id)
    finally:
        for fh in cell_handles.values():
            fh.close()

    for cell_id, total_reads in cell_read_counts.items():
        print(f'  {cell_id}: {total_reads} reads, {len(ground_truth[cell_id])} repeats', file=sys.stderr)

    return cell_fastq_paths, ground_truth


def simulate_chromium(fasta_path, cell_plan, read_length, barcode_length, umi_length, outdir, rng):
    os.makedirs(outdir, exist_ok=True)
    gene_to_cells = build_gene_to_cells(cell_plan)
    chrom_gene_coords = build_chrom_gene_coords(cell_plan)
    needed_chroms = set(chrom_gene_coords.keys())

    cell_barcodes = sorted(set(cell_plan.keys()))
    barcodes_path = os.path.join(outdir, 'barcodes.txt')
    with open(barcodes_path, 'w') as fh:
        fh.write('\n'.join(cell_barcodes) + '\n')

    r1_path = os.path.join(outdir, 'R1.fastq.gz')
    r2_path = os.path.join(outdir, 'R2.fastq.gz')
    ground_truth = defaultdict(dict)
    total_reads = 0

    with gzip.open(r1_path, 'wt') as r1f, gzip.open(r2_path, 'wt') as r2f:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_path, needed_chroms):
            for gene_id, (start, end, strand) in chrom_gene_coords[chrom].items():
                repeat_seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if repeat_seq is None:
                    continue
                for cell_id, count, family_id, class_id in gene_to_cells[gene_id]:
                    used_umis = set()
                    for _ in range(count):
                        umi = ''.join(rng.choices('ACGT', k=umi_length))
                        while umi in used_umis:
                            umi = ''.join(rng.choices('ACGT', k=umi_length))
                        used_umis.add(umi)
                        read_id = f'r{total_reads}_{cell_id[:6]}_{umi}_{safe_id(gene_id)}'
                        r1_seq = cell_id + umi
                        r2_seq = sample_subseq(repeat_seq, read_length, rng)
                        r1f.write(f'@{read_id}\n{r1_seq}\n+\n{make_qual(len(r1_seq))}\n')
                        r2f.write(f'@{read_id}\n{r2_seq}\n+\n{make_qual(len(r2_seq))}\n')
                        total_reads += 1
                    ground_truth[cell_id][gene_id] = (count, family_id, class_id)

    print(f'  {total_reads} read pairs across {len(cell_barcodes)} cells', file=sys.stderr)
    return (r1_path, r2_path), ground_truth


def write_ground_truth(ground_truth, output_path, cell_column):
    with open(output_path, 'w') as fh:
        fh.write(f'{cell_column}\trepeat_id\tfamily_id\tclass_id\ttrue_count\n')
        for cell_id in sorted(ground_truth):
            for repeat_id in sorted(ground_truth[cell_id]):
                count, family_id, class_id = ground_truth[cell_id][repeat_id]
                fh.write(f'{cell_id}\t{repeat_id}\t{family_id}\t{class_id}\t{count}\n')


def write_smartseq2_manifest(cell_fastq_paths, manifest_path):
    with open(manifest_path, 'w') as fh:
        for cell_id in sorted(cell_fastq_paths):
            fq_path = cell_fastq_paths[cell_id]
            fh.write(f'{fq_path}\t-\t{cell_id}\n')


def write_repeat_fasta(cell_plan, fasta_path, chrom_gene_coords, fasta_source_path):
    needed_chroms = set(chrom_gene_coords.keys())
    all_gene_coords = {
        gene_id: (chrom, start, end, strand)
        for chrom, genes in chrom_gene_coords.items()
        for gene_id, (start, end, strand) in genes.items()
    }
    written = set()
    with open(fasta_path, 'w') as out_fh:
        for chrom, chrom_seq in stream_fasta_by_chrom(fasta_source_path, needed_chroms):
            for gene_id, (g_chrom, start, end, strand) in all_gene_coords.items():
                if g_chrom != chrom or gene_id in written:
                    continue
                seq = extract_repeat_sequence(chrom_seq, start, end, strand)
                if seq is not None:
                    out_fh.write(f'>{safe_id(gene_id)}\n{seq}\n')
                    written.add(gene_id)


def main():
    ap = argparse.ArgumentParser(description='Simulate scRNA-seq reads from repeat elements')
    ap.add_argument('--mode', choices=['smartseq2', 'chromium'], required=True)
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
    args = ap.parse_args()

    rng = random.Random(args.seed)
    allowed_chroms = set(args.chromosomes) if args.chromosomes else None

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
    chrom_gene_coords = build_chrom_gene_coords(cell_plan)

    if args.mode == 'smartseq2':
        print(f'Simulating SmartSeq2 reads (streaming FASTA by chrom)', file=sys.stderr)
        cell_fastq_paths, ground_truth = simulate_smartseq2(
            args.fasta, cell_plan, args.read_length, args.outdir, rng)
        write_ground_truth(ground_truth, ground_truth_path, cell_column='cell_id')
        manifest_path = os.path.join(args.outdir, 'manifest.tsv')
        write_smartseq2_manifest(cell_fastq_paths, manifest_path)
        print(f'Manifest written to {manifest_path}', file=sys.stderr)

    else:
        print(f'Simulating Chromium reads (streaming FASTA by chrom)', file=sys.stderr)
        (r1, r2), ground_truth = simulate_chromium(
            args.fasta, cell_plan, args.read_length,
            args.cb_length, args.umi_length, args.outdir, rng)
        write_ground_truth(ground_truth, ground_truth_path, cell_column='cell_barcode')
        print(f'R1: {r1}', file=sys.stderr)
        print(f'R2: {r2}', file=sys.stderr)

    print(f'Writing repeat FASTA to {repeat_fasta_path}', file=sys.stderr)
    write_repeat_fasta(cell_plan, repeat_fasta_path, chrom_gene_coords, args.fasta)
    print(f'Ground truth written to {ground_truth_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
