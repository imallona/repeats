"""
Unit tests for workflow/scripts/simulate_reads.py

Tests pure functions that do not require real FASTA/GTF files.
"""
import gzip
import io
import os
import random
import textwrap

import pytest

import simulate_reads as sr


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

def test_reverse_complement_known():
    assert sr.reverse_complement('ATCG') == 'CGAT'


def test_reverse_complement_all_bases():
    assert sr.reverse_complement('AACCGGTT') == 'AACCGGTT'


def test_reverse_complement_n_preserved():
    assert sr.reverse_complement('NNAANN') == 'NNTTN N'.replace(' ', '')
    assert sr.reverse_complement('NNN') == 'NNN'


def test_reverse_complement_single_base():
    assert sr.reverse_complement('A') == 'T'
    assert sr.reverse_complement('T') == 'A'
    assert sr.reverse_complement('C') == 'G'
    assert sr.reverse_complement('G') == 'C'


def test_reverse_complement_involution():
    seq = 'ATCGATCGNNATCG'
    assert sr.reverse_complement(sr.reverse_complement(seq)) == seq


# ---------------------------------------------------------------------------
# parse_gtf_attribute
# ---------------------------------------------------------------------------

def test_parse_gtf_attribute_basic():
    attrs = 'gene_id "AluSz6"; transcript_id "AluSz6_dup1"; family_id "Alu";'
    assert sr.parse_gtf_attribute(attrs, 'gene_id') == 'AluSz6'
    assert sr.parse_gtf_attribute(attrs, 'transcript_id') == 'AluSz6_dup1'
    assert sr.parse_gtf_attribute(attrs, 'family_id') == 'Alu'


def test_parse_gtf_attribute_missing_returns_none():
    attrs = 'gene_id "AluSz6";'
    assert sr.parse_gtf_attribute(attrs, 'class_id') is None


def test_parse_gtf_attribute_no_quotes():
    attrs = 'gene_id AluSz6; transcript_id AluSz6_dup1;'
    assert sr.parse_gtf_attribute(attrs, 'gene_id') == 'AluSz6'


# ---------------------------------------------------------------------------
# extract_repeat_sequence
# ---------------------------------------------------------------------------

def test_extract_repeat_sequence_forward_strand():
    chrom = 'A' * 100 + 'GATTACA' * 10 + 'T' * 100
    seq = sr.extract_repeat_sequence(chrom, 100, 170, '+')
    assert seq is not None
    assert seq == ('GATTACA' * 10)[:70]


def test_extract_repeat_sequence_minus_strand():
    seq_plus = 'ACGT' * 20
    chrom = 'N' * 10 + seq_plus + 'N' * 10
    seq = sr.extract_repeat_sequence(chrom, 10, 90, '-')
    assert seq == sr.reverse_complement(seq_plus)


def test_extract_repeat_sequence_too_many_n():
    chrom = 'N' * 200
    assert sr.extract_repeat_sequence(chrom, 0, 100, '+') is None


def test_extract_repeat_sequence_too_short():
    chrom = 'ACGT' * 100
    assert sr.extract_repeat_sequence(chrom, 0, 10, '+') is None  # < 50 bp


def test_extract_repeat_sequence_exactly_50bp_ok():
    chrom = 'ACGT' * 100
    seq = sr.extract_repeat_sequence(chrom, 0, 50, '+')
    assert seq is not None
    assert len(seq) == 50


def test_extract_repeat_sequence_n_just_below_threshold():
    # 9% N -> should pass (threshold = 10%)
    body = 'A' * 91 + 'N' * 9
    chrom = body + 'A' * 100
    seq = sr.extract_repeat_sequence(chrom, 0, 100, '+')
    assert seq is not None


def test_extract_repeat_sequence_n_just_above_threshold():
    # 11% N -> should fail (threshold is N/len > 0.1)
    body = 'A' * 89 + 'N' * 11
    chrom = body + 'A' * 100
    assert sr.extract_repeat_sequence(chrom, 0, 100, '+') is None


# ---------------------------------------------------------------------------
# sample_subseq
# ---------------------------------------------------------------------------

def test_sample_subseq_correct_length():
    rng = random.Random(42)
    seq = 'ACGT' * 100
    for _ in range(10):
        sub = sr.sample_subseq(seq, 90, rng)
        assert len(sub) == 90


def test_sample_subseq_short_sequence_padded():
    rng = random.Random(42)
    seq = 'ACGT' * 5  # 20 bp
    sub = sr.sample_subseq(seq, 90, rng)
    assert len(sub) == 90
    assert 'N' in sub  # padded with N


def test_sample_subseq_only_valid_bases():
    rng = random.Random(99)
    seq = 'ACGT' * 100
    sub = sr.sample_subseq(seq, 50, rng)
    assert all(b in 'ACGTN' for b in sub)


# ---------------------------------------------------------------------------
# sample_count_geometric
# ---------------------------------------------------------------------------

def test_sample_count_geometric_always_positive():
    rng = random.Random(0)
    for _ in range(1000):
        c = sr.sample_count_geometric(rng)
        assert c >= 1


def test_sample_count_geometric_respects_max():
    rng = random.Random(0)
    for _ in range(1000):
        c = sr.sample_count_geometric(rng, max_count=10)
        assert c <= 10


def test_sample_count_geometric_mean_approx():
    rng = random.Random(42)
    counts = [sr.sample_count_geometric(rng, mean=5.0) for _ in range(5000)]
    mean = sum(counts) / len(counts)
    # geometric with mean=5 should be around 5, allow ±1.5
    assert 3.5 <= mean <= 6.5


# ---------------------------------------------------------------------------
# build_cell_plan
# ---------------------------------------------------------------------------

def test_build_cell_plan_structure():
    intervals = {
        'chr1': [
            (0, 200, 'AluSz6_dup1', 'AluSz6', 'Alu', 'SINE', '+'),
            (500, 800, 'L1PA2_dup1', 'L1PA2', 'L1', 'LINE', '-'),
        ]
    }
    rng = random.Random(42)
    plan = sr.build_cell_plan(intervals, n_cells=3, mean_expressed_per_cell=2, rng=rng)
    assert len(plan) == 3
    for cell_id, locus_map in plan.items():
        assert cell_id.startswith('cell_')
        for locus_id, info in locus_map.items():
            count, gene_id, family_id, class_id, chrom, start, end, strand = info
            assert count >= 1


# ---------------------------------------------------------------------------
# build_locus_to_cells
# ---------------------------------------------------------------------------

def test_build_locus_to_cells_groups_correctly():
    plan = {
        'cell_001': {'locus_A': (3, 'gA', 'fA', 'cA', 'chr1', 0, 100, '+')},
        'cell_002': {'locus_A': (2, 'gA', 'fA', 'cA', 'chr1', 0, 100, '+'),
                     'locus_B': (5, 'gB', 'fB', 'cB', 'chr1', 200, 400, '+')},
    }
    l2c = sr.build_locus_to_cells(plan)
    assert len(l2c['locus_A']) == 2
    assert len(l2c['locus_B']) == 1
    cell_ids_for_A = [x[0] for x in l2c['locus_A']]
    assert 'cell_001' in cell_ids_for_A
    assert 'cell_002' in cell_ids_for_A


# ---------------------------------------------------------------------------
# write_ground_truth
# ---------------------------------------------------------------------------

def test_write_ground_truth_format(tmp_path):
    ground_truth = {
        'cell_001': {'locus_A': (5, 'gene_A', 'fam_A', 'cls_A')},
        'cell_002': {'locus_A': (2, 'gene_A', 'fam_A', 'cls_A'),
                     'locus_B': (7, 'gene_B', 'fam_B', 'cls_B')},
    }
    out = tmp_path / 'gt.tsv'
    sr.write_ground_truth(ground_truth, str(out), cell_column='cell_id')

    with open(out) as fh:
        lines = fh.readlines()

    # header
    header = lines[0].rstrip('\n').split('\t')
    assert header[0] == 'cell_id'
    assert 'locus_id' in header
    assert 'true_count' in header

    # data rows sorted by cell_id then locus_id
    rows = [l.rstrip('\n').split('\t') for l in lines[1:]]
    assert len(rows) == 3
    # first row: cell_001 locus_A
    assert rows[0][0] == 'cell_001'
    assert rows[0][5] == '5'


def test_write_ground_truth_total_rows(tmp_path):
    gt = {'c1': {'L1': (1, 'L1', 'L1', 'LINE'), 'L2': (2, 'L2', 'L2', 'LINE')},
          'c2': {'L1': (3, 'L1', 'L1', 'LINE')}}
    out = tmp_path / 'gt.tsv'
    sr.write_ground_truth(gt, str(out), 'cell_id')
    with open(out) as fh:
        lines = fh.readlines()
    assert len(lines) == 4  # header + 3 data rows


# ---------------------------------------------------------------------------
# parse_gtf_repeats_by_chrom  (uses in-memory GTF string)
# ---------------------------------------------------------------------------

def test_parse_gtf_repeats_by_chrom_minimal(tmp_path):
    gtf_content = (
        'chr1\trmsk\texon\t101\t500\t.\t+\t.\t'
        'gene_id "AluSz6"; transcript_id "AluSz6_dup1"; '
        'family_id "Alu"; class_id "SINE";\n'
        'chr2\trmsk\texon\t1001\t1600\t.\t-\t.\t'
        'gene_id "L1PA2"; transcript_id "L1PA2_dup1"; '
        'family_id "L1"; class_id "LINE";\n'
    )
    gtf = tmp_path / 'repeats.gtf'
    gtf.write_text(gtf_content)
    intervals = sr.parse_gtf_repeats_by_chrom(str(gtf))
    assert 'chr1' in intervals
    assert 'chr2' in intervals
    assert len(intervals['chr1']) == 1
    locus = intervals['chr1'][0]
    # (start_0, end_0, locus_id, gene_id, family_id, class_id, strand)
    assert locus[2] == 'AluSz6_dup1'
    assert locus[3] == 'AluSz6'


def test_parse_gtf_repeats_by_chrom_filters_short(tmp_path):
    gtf_content = (
        'chr1\trmsk\texon\t101\t140\t.\t+\t.\t'  # only 39 bp < 50 -> filtered
        'gene_id "Tiny"; transcript_id "Tiny_dup1"; '
        'family_id "Tiny"; class_id "DNA";\n'
    )
    gtf = tmp_path / 'repeats.gtf'
    gtf.write_text(gtf_content)
    intervals = sr.parse_gtf_repeats_by_chrom(str(gtf))
    assert len(intervals) == 0


def test_parse_gtf_repeats_by_chrom_allowed_chroms(tmp_path):
    gtf_content = (
        'chr1\trmsk\texon\t101\t500\t.\t+\t.\t'
        'gene_id "AluSz6"; transcript_id "AluSz6_dup1"; '
        'family_id "Alu"; class_id "SINE";\n'
        'chr2\trmsk\texon\t1001\t1600\t.\t-\t.\t'
        'gene_id "L1PA2"; transcript_id "L1PA2_dup1"; '
        'family_id "L1"; class_id "LINE";\n'
    )
    gtf = tmp_path / 'repeats.gtf'
    gtf.write_text(gtf_content)
    intervals = sr.parse_gtf_repeats_by_chrom(str(gtf), allowed_chroms={'chr1'})
    assert 'chr1' in intervals
    assert 'chr2' not in intervals
