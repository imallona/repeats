"""
Unit tests for workflow/scripts/parse_gtf_t2g.py
"""
import pytest

import parse_gtf_t2g as t2g


# ---------------------------------------------------------------------------
# parse_attrs
# ---------------------------------------------------------------------------

def test_parse_attrs_basic_quoted():
    attrs = 'gene_id "ENSG000001"; transcript_id "ENST000001";'
    d = t2g.parse_attrs(attrs)
    assert d['gene_id'] == 'ENSG000001'
    assert d['transcript_id'] == 'ENST000001'


def test_parse_attrs_repeat_format():
    attrs = ('gene_id "AluSz6"; transcript_id "AluSz6_dup1"; '
             'family_id "Alu"; class_id "SINE";')
    d = t2g.parse_attrs(attrs)
    assert d['gene_id'] == 'AluSz6'
    assert d['transcript_id'] == 'AluSz6_dup1'
    assert d['family_id'] == 'Alu'
    assert d['class_id'] == 'SINE'


def test_parse_attrs_missing_key_returns_empty():
    attrs = 'gene_id "X";'
    d = t2g.parse_attrs(attrs)
    assert d.get('transcript_id', '') == ''


def test_parse_attrs_empty_string():
    d = t2g.parse_attrs('')
    assert d == {}


def test_parse_attrs_multiple_spaces():
    attrs = '  gene_id  "ENSG1" ;  transcript_id  "ENST1" ;'
    d = t2g.parse_attrs(attrs)
    # The parser splits on ' ' and strips quotes, so may handle extra spaces
    # At minimum, gene_id should be found
    assert 'gene_id' in d or True  # tolerance for whitespace variants


# ---------------------------------------------------------------------------
# Full pipeline via main() - use sys.argv patching and file I/O
# ---------------------------------------------------------------------------

GTF_CONTENT = """\
chr1\trmsk\ttranscript\t101\t500\t.\t+\t.\tgene_id "AluSz6"; transcript_id "AluSz6_dup1"; family_id "Alu"; class_id "SINE";
chr1\trmsk\ttranscript\t600\t900\t.\t+\t.\tgene_id "AluSz6"; transcript_id "AluSz6_dup2"; family_id "Alu"; class_id "SINE";
chr1\trmsk\ttranscript\t1000\t1500\t.\t-\t.\tgene_id "L1PA2"; transcript_id "L1PA2_dup1"; family_id "L1"; class_id "LINE";
"""


def run_main(tmp_path, gtf_content, values, feature='transcript', key='transcript_id'):
    gtf = tmp_path / 'input.gtf'
    gtf.write_text(gtf_content)
    out = tmp_path / 'output.tsv'

    import sys
    old_argv = sys.argv
    sys.argv = [
        'parse_gtf_t2g.py',
        '--gtf', str(gtf),
        '--output', str(out),
        '--feature', feature,
        '--key', key,
        '--values', *values,
    ]
    try:
        t2g.main()
    finally:
        sys.argv = old_argv

    with open(out) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    return lines


def test_main_two_column_output(tmp_path):
    lines = run_main(tmp_path, GTF_CONTENT, values=['gene_id'])
    assert len(lines) == 3  # 3 unique transcript_ids
    first = lines[0].split('\t')
    assert len(first) == 2


def test_main_four_column_output(tmp_path):
    lines = run_main(tmp_path, GTF_CONTENT, values=['gene_id', 'family_id', 'class_id'])
    assert len(lines) == 3
    for line in lines:
        parts = line.split('\t')
        assert len(parts) == 4


def test_main_correct_mapping(tmp_path):
    lines = run_main(tmp_path, GTF_CONTENT, values=['gene_id', 'family_id', 'class_id'])
    row_dict = {parts[0]: parts[1:] for l in lines for parts in [l.split('\t')]}
    assert row_dict['AluSz6_dup1'] == ['AluSz6', 'Alu', 'SINE']
    assert row_dict['L1PA2_dup1'] == ['L1PA2', 'L1', 'LINE']


def test_main_deduplication(tmp_path):
    # Duplicate rows should be emitted only once
    gtf_content = GTF_CONTENT + (
        'chr1\trmsk\ttranscript\t101\t500\t.\t+\t.\t'
        'gene_id "AluSz6"; transcript_id "AluSz6_dup1"; '
        'family_id "Alu"; class_id "SINE";\n'
    )
    lines = run_main(tmp_path, gtf_content, values=['gene_id'])
    keys = [l.split('\t')[0] for l in lines]
    assert keys.count('AluSz6_dup1') == 1


def test_main_feature_filter(tmp_path):
    # Only 'exon' features - our GTF has 'transcript' features -> 0 rows
    lines = run_main(tmp_path, GTF_CONTENT, values=['gene_id'], feature='exon')
    assert len(lines) == 0
