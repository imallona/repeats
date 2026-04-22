"""
Unit tests for workflow/scripts/build_repeat_biology_annotation.py

The pure classification helpers are covered by test_biology_annotations.py;
here we exercise the I/O-orchestration pieces: FASTA random access by .fai
arithmetic, strand-aware merged-interval overlap, and the end-to-end
aggregation from a tiny rmsk fixture.
"""
import gzip
import os
import textwrap

import pytest

import build_repeat_biology_annotation as bra


def write_fasta(tmp_path, sequences, line_width=10):
    """Write a multi-record FASTA with a fixed per-chrom line width."""
    path = os.path.join(tmp_path, 'genome.fa')
    with open(path, 'w') as fh:
        for name, seq in sequences:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + '\n')
    return path


def test_normalize_chrom():
    assert bra.normalize_chrom('chr1') == '1'
    assert bra.normalize_chrom('1') == '1'
    assert bra.normalize_chrom('chrX') == 'X'
    assert bra.normalize_chrom('chrMT') == 'MT'


def test_reverse_complement():
    assert bra.reverse_complement('ACGT') == 'ACGT'
    assert bra.reverse_complement('AAAA') == 'TTTT'
    assert bra.reverse_complement('ACGTN') == 'NACGT'
    assert bra.reverse_complement('') == ''


def test_fasta_random_access_spans_lines(tmp_path):
    seq = 'ACGTACGTAC' + 'GTACGTACGT' + 'ACGT'
    path = write_fasta(str(tmp_path), [('1', seq)], line_width=10)
    with bra.FastaRandomAccess(path) as fa:
        assert fa.index['1'][0] == len(seq)
        assert fa.fetch('1', 0, 4) == 'ACGT'
        assert fa.fetch('1', 8, 14) == 'ACGTAC'
        assert fa.fetch('1', 0, len(seq)) == seq
        assert fa.fetch('1', 10, 14) == 'GTAC'


def test_fasta_random_access_unknown_chrom_and_bounds(tmp_path):
    path = write_fasta(str(tmp_path), [('1', 'ACGTACGT')], line_width=10)
    with bra.FastaRandomAccess(path) as fa:
        assert fa.fetch('missing', 0, 5) == ''
        assert fa.fetch('1', 5, 100) == 'CGT'
        assert fa.fetch('1', 10, 20) == ''
        assert fa.fetch('1', -5, 4) == 'ACGT'


def test_fasta_random_access_multiple_records(tmp_path):
    path = write_fasta(str(tmp_path), [
        ('1', 'AAAAAAAAAA'),
        ('2', 'CCCCCCCCCC'),
    ], line_width=10)
    with bra.FastaRandomAccess(path) as fa:
        assert fa.fetch('1', 0, 10) == 'A' * 10
        assert fa.fetch('2', 0, 10) == 'C' * 10
        assert fa.fetch('1', 3, 7) == 'AAAA'
        assert fa.fetch('2', 3, 7) == 'CCCC'


def test_downstream_sequence_plus_strand(tmp_path):
    path = write_fasta(str(tmp_path), [('1', 'CGCGAAAAAAAA')], line_width=20)
    with bra.FastaRandomAccess(path) as fa:
        assert bra.downstream_sequence(fa, '1', 0, 4, '+', 8) == 'AAAAAAAA'


def test_downstream_sequence_minus_strand(tmp_path):
    seq = 'TTTTTTTTCGCG'
    path = write_fasta(str(tmp_path), [('1', seq)], line_width=20)
    with bra.FastaRandomAccess(path) as fa:
        ds = bra.downstream_sequence(fa, '1', 8, 12, '-', 8)
        assert ds == 'AAAAAAAA'


def test_gene_intervals_any_overlap():
    g = bra.GeneIntervals()
    g.add('1', 100, 200, '+')
    g.add('1', 150, 250, '+')
    g.add('1', 400, 500, '+')
    g.add('1', 100, 200, '-')
    g.finalize()

    assert g.any_overlap('1', 120, 130, '+') is True
    assert g.any_overlap('1', 249, 260, '+') is True
    assert g.any_overlap('1', 250, 260, '+') is False
    assert g.any_overlap('1', 250, 400, '+') is False
    assert g.any_overlap('1', 499, 501, '+') is True
    assert g.any_overlap('1', 500, 600, '+') is False
    assert g.any_overlap('1', 120, 130, '-') is True
    assert g.any_overlap('2', 120, 130, '+') is False


def test_gene_intervals_merge_and_boundary():
    g = bra.GeneIntervals()
    g.add('1', 0, 100, '+')
    g.add('1', 100, 200, '+')
    g.add('1', 250, 300, '+')
    g.finalize()
    assert g.any_overlap('1', 50, 150, '+') is True
    assert g.any_overlap('1', 199, 201, '+') is True
    assert g.any_overlap('1', 200, 250, '+') is False
    assert g.any_overlap('1', 299, 301, '+') is True


def test_load_gene_intervals_from_gtf(tmp_path):
    gtf = os.path.join(str(tmp_path), 'genes.gtf')
    with open(gtf, 'w') as fh:
        fh.write('#comment line\n')
        fh.write('1\thavana\tgene\t101\t200\t.\t+\t.\tgene_id "g1";\n')
        fh.write('1\thavana\texon\t101\t150\t.\t+\t.\tgene_id "g1";\n')
        fh.write('1\thavana\tgene\t301\t400\t.\t-\t.\tgene_id "g2";\n')
        fh.write('chr2\thavana\tgene\t1\t50\t.\t+\t.\tgene_id "g3";\n')
    ivs = bra.load_gene_intervals(gtf)
    assert ivs.any_overlap('1', 120, 130, '+') is True
    assert ivs.any_overlap('1', 120, 130, '-') is False
    assert ivs.any_overlap('1', 320, 330, '-') is True
    assert ivs.any_overlap('2', 10, 20, '+') is True


def _write_rmsk(tmp_path, rows, gzipped=True):
    path = os.path.join(str(tmp_path),
                        'rmsk.txt.gz' if gzipped else 'rmsk.txt')
    opener = gzip.open if gzipped else open
    with opener(path, 'wt') as fh:
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')
    return path


def _rmsk_row(milli, chrom, start, end, strand, name, cls, family):
    return [0, 1000, milli, 0, 0, chrom, start, end, -100,
            strand, name, cls, family, 0, end - start, 0, 1]


def test_iter_rmsk_parses_core_fields(tmp_path):
    path = _write_rmsk(tmp_path, [
        _rmsk_row(30, 'chr1', 100, 400, '+', 'AluY', 'SINE', 'Alu'),
        _rmsk_row(250, 'chr2', 0, 60, '-', 'L1MA1', 'LINE', 'L1'),
    ])
    records = list(bra.iter_rmsk(path))
    assert len(records) == 2
    chrom, start, end, strand, milli, name, cls, family = records[0]
    assert chrom == '1'
    assert (start, end, strand, milli, name, cls, family) == (
        100, 400, '+', 30, 'AluY', 'SINE', 'Alu')
    assert records[1][0] == '2'
    assert records[1][4] == 250


def test_aggregate_end_to_end(tmp_path):
    # Genome: two short chromosomes, each with an A-run downstream of a plus
    # strand locus and a T-run upstream of a minus strand locus (which becomes
    # an A-run in mRNA-sense orientation).
    chrom1 = 'ACGTACGTAC' + 'AAAAAAAAAAAA' + 'GGGGGGGG'
    chrom2 = 'GGGGGGGG' + 'TTTTTTTTTTTT' + 'CGCGCGCGCG'
    fasta = write_fasta(str(tmp_path), [('1', chrom1), ('2', chrom2)],
                        line_width=20)

    gtf = os.path.join(str(tmp_path), 'genes.gtf')
    with open(gtf, 'w') as fh:
        fh.write('1\thav\tgene\t5\t50\t.\t+\t.\tgene_id "g1";\n')

    rmsk_rows = [
        _rmsk_row(20, 'chr1', 0, 10, '+', 'AluY', 'SINE', 'Alu'),
        _rmsk_row(300, 'chr1', 0, 10, '+', 'AluJb', 'SINE', 'Alu'),
        _rmsk_row(40, 'chr2', 20, 30, '-', 'L1HS', 'LINE', 'L1'),
    ]
    rmsk = _write_rmsk(tmp_path, rmsk_rows)

    with bra.FastaRandomAccess(fasta) as fa:
        ivs = bra.load_gene_intervals(gtf)
        agg = bra.aggregate(
            rmsk_path=rmsk,
            gene_intervals=ivs,
            fasta=fa,
            downstream_bp=12,
            a_run_min_length=8,
            min_length=5,
        )

    assert set(agg.keys()) == {'AluY', 'AluJb', 'L1HS'}

    aluy = agg['AluY']
    assert aluy['n_loci'] == 1
    assert aluy['n_sense'] == 1
    assert aluy['n_intergenic'] == 0
    assert aluy['n_a_run_hit'] == 1

    alujb = agg['AluJb']
    assert alujb['n_sense'] == 1
    assert alujb['milli'] == [300]

    l1hs = agg['L1HS']
    assert l1hs['n_intergenic'] == 1
    assert l1hs['n_a_run_hit'] == 1


def test_write_annotation_columns_and_values(tmp_path):
    per_subfamily = {
        'AluY': {
            'family_id': 'Alu', 'class_id': 'SINE',
            'n_loci': 4, 'milli': [10, 20, 30, 40],
            'n_sense': 2, 'n_antisense': 1, 'n_both': 0, 'n_intergenic': 1,
            'n_a_run_hit': 3,
        },
    }
    out_path = os.path.join(str(tmp_path), 'out.tsv')
    bra.write_annotation(per_subfamily, out_path)
    with open(out_path) as fh:
        header = fh.readline().rstrip().split('\t')
        row = fh.readline().rstrip().split('\t')
    assert header == bra.OUTPUT_COLUMNS
    fields = dict(zip(header, row))
    assert fields['gene_id'] == 'AluY'
    assert fields['family_id'] == 'Alu'
    assert fields['class_id'] == 'SINE'
    assert fields['n_loci'] == '4'
    assert fields['family_age'] == 'young'
    assert fields['family_age_from_millidiv'] == 'young'
    assert fields['frac_sense_gene'] == '0.5000'
    assert fields['frac_antisense_gene'] == '0.2500'
    assert fields['frac_both'] == '0.0000'
    assert fields['frac_intergenic'] == '0.2500'
    assert fields['strand_vs_gene'] == 'sense'
    assert fields['frac_downstream_a_run'] == '0.7500'


def test_write_annotation_empty_milli_branch(tmp_path):
    # Defensive: write_annotation should cope with a subfamily that has no
    # milliDiv observations (e.g. all rows hit a parse failure upstream).
    per_subfamily = {
        'X': {
            'family_id': '', 'class_id': '',
            'n_loci': 1, 'milli': [],
            'n_sense': 0, 'n_antisense': 0, 'n_both': 0, 'n_intergenic': 1,
            'n_a_run_hit': 0,
        },
    }
    out_path = os.path.join(str(tmp_path), 'out.tsv')
    bra.write_annotation(per_subfamily, out_path)
    with open(out_path) as fh:
        header = fh.readline()
        row = fh.readline().rstrip().split('\t')
    fields = dict(zip(bra.OUTPUT_COLUMNS, row))
    assert fields['median_millidiv'] == ''
    assert fields['family_age_from_millidiv'] == 'unknown'
    assert fields['strand_vs_gene'] == 'intergenic'
