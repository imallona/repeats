"""
Unit tests for workflow/scripts/build_rmsk_gtf.py
"""
import gzip
import os

import pytest

import build_rmsk_gtf as rmsk


# ---------------------------------------------------------------------------
# make_gtf_attributes
# ---------------------------------------------------------------------------

def test_make_gtf_attributes_format():
    attr = rmsk.make_gtf_attributes('AluSz6', 3, 'Alu', 'SINE')
    assert 'gene_id "AluSz6"' in attr
    assert 'transcript_id "AluSz6_dup3"' in attr
    assert 'family_id "Alu"' in attr
    assert 'class_id "SINE"' in attr


def test_make_gtf_attributes_locus_index():
    a1 = rmsk.make_gtf_attributes('L1PA2', 1, 'L1', 'LINE')
    a5 = rmsk.make_gtf_attributes('L1PA2', 5, 'L1', 'LINE')
    assert 'transcript_id "L1PA2_dup1"' in a1
    assert 'transcript_id "L1PA2_dup5"' in a5


# ---------------------------------------------------------------------------
# convert_rmsk_to_gtf
# ---------------------------------------------------------------------------


def rmsk_row(chrom='chr1', start=100, end=400, strand='+',
             name='AluSz6', rep_class='SINE', family='Alu'):
    return (f'0\t1000\t100\t0\t0\t{chrom}\t{start}\t{end}\t-100\t'
            f'{strand}\t{name}\t{rep_class}\t{family}\t0\t300\t0\t1\n')


def test_convert_rmsk_to_gtf_basic(tmp_path):
    infile = tmp_path / 'rmsk.txt'
    infile.write_text(rmsk_row())
    outfile = tmp_path / 'out.gtf'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile))
    lines = outfile.read_text().strip().split('\n')
    assert len(lines) == 1
    fields = lines[0].split('\t')
    assert fields[0] == 'chr1'
    assert fields[2] == 'exon'
    assert 'gene_id "AluSz6"' in fields[8]
    assert 'transcript_id "AluSz6_dup1"' in fields[8]


def test_convert_rmsk_to_gtf_increments_dup_index(tmp_path):
    infile = tmp_path / 'rmsk.txt'
    infile.write_text(
        rmsk_row(name='AluSz6', start=100, end=400) +
        rmsk_row(name='AluSz6', start=500, end=900)
    )
    outfile = tmp_path / 'out.gtf'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile))
    lines = outfile.read_text().strip().split('\n')
    assert len(lines) == 2
    assert 'transcript_id "AluSz6_dup1"' in lines[0]
    assert 'transcript_id "AluSz6_dup2"' in lines[1]


def test_convert_rmsk_to_gtf_filters_short(tmp_path):
    infile = tmp_path / 'rmsk.txt'
    # Only 30 bp wide -> below default min_length=50
    infile.write_text(rmsk_row(start=100, end=130))
    outfile = tmp_path / 'out.gtf'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile))
    assert outfile.read_text().strip() == ''


def test_convert_rmsk_to_gtf_chromosome_filter(tmp_path):
    infile = tmp_path / 'rmsk.txt'
    infile.write_text(
        rmsk_row(chrom='chr1') +
        rmsk_row(chrom='chr2')
    )
    outfile = tmp_path / 'out.gtf'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile), allowed_chroms={'chr1'})
    lines = [l for l in outfile.read_text().strip().split('\n') if l]
    assert len(lines) == 1
    assert lines[0].startswith('chr1')


def test_convert_rmsk_to_gtf_gtf_start_is_one_based(tmp_path):
    # RMSK uses 0-based start; GTF should be 1-based
    infile = tmp_path / 'rmsk.txt'
    infile.write_text(rmsk_row(start=100, end=400))
    outfile = tmp_path / 'out.gtf'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile))
    fields = outfile.read_text().strip().split('\t')
    assert fields[3] == '101'   # 100 + 1
    assert fields[4] == '400'   # end unchanged


def test_convert_rmsk_to_gtf_compressed_output(tmp_path):
    infile = tmp_path / 'rmsk.txt'
    infile.write_text(rmsk_row())
    outfile = tmp_path / 'out.gtf.gz'
    rmsk.convert_rmsk_to_gtf(str(infile), str(outfile))
    with gzip.open(outfile, 'rt') as fh:
        content = fh.read()
    assert 'AluSz6' in content
