"""
Unit tests for count_starsolo_locus.py pure functions and mocked subprocess paths.
"""
import os
import sys
import textwrap
from unittest.mock import MagicMock, patch

import pytest

import count_starsolo_locus as csl


# ---------------------------------------------------------------------------
# parse_locus_id
# ---------------------------------------------------------------------------

def test_parse_locus_id_plus_strand():
    result = csl.parse_locus_id("AluSz_dup1::chr10:10000-10435(+)")
    assert result == ("chr10", 10000, 10435, "+")


def test_parse_locus_id_minus_strand():
    result = csl.parse_locus_id("L1M4_dup3::10:65000-65900(-)")
    assert result == ("10", 65000, 65900, "-")


def test_parse_locus_id_malformed_returns_none():
    assert csl.parse_locus_id("no_double_colon") is None
    assert csl.parse_locus_id("bad::coords") is None


# ---------------------------------------------------------------------------
# load_intervals
# ---------------------------------------------------------------------------

def test_load_intervals_basic(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(
        ">AluSz_dup1::10:10000-10435(+)\nACGT\n"
        ">L1M4_dup3::10:65000-65900(-)\nACGT\n"
    )
    chrom_starts, chrom_intervals = csl.load_intervals(str(fa))
    assert "10" in chrom_intervals
    assert len(chrom_intervals["10"]) == 2
    # sorted by start
    assert chrom_intervals["10"][0][0] == 10000
    assert chrom_intervals["10"][1][0] == 65000
    assert chrom_starts["10"] == [10000, 65000]


def test_load_intervals_skips_non_header_lines(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">AluSz_dup1::10:10000-10200(+)\nACGTACGT\n")
    _, intervals = csl.load_intervals(str(fa))
    assert len(intervals["10"]) == 1
    assert intervals["10"][0][2] == "AluSz_dup1"


def test_load_intervals_deduplicates(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(
        ">dup1::10:100-200(+)\nACGT\n"
        ">dup1::10:100-200(+)\nACGT\n"
    )
    _, intervals = csl.load_intervals(str(fa))
    assert len(intervals["10"]) == 1


# ---------------------------------------------------------------------------
# find_locus
# ---------------------------------------------------------------------------

def make_intervals():
    chrom_intervals = {"10": [(10000, 10200, "AluSz_dup1"), (20000, 20500, "L1_dup2")]}
    chrom_starts = {"10": [10000, 20000]}
    return chrom_starts, chrom_intervals


def test_find_locus_hit_first():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 10100, cs, ci) == "AluSz_dup1"


def test_find_locus_hit_second():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 20300, cs, ci) == "L1_dup2"


def test_find_locus_boundary_start():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 10000, cs, ci) == "AluSz_dup1"


def test_find_locus_boundary_end_exclusive():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 10200, cs, ci) is None


def test_find_locus_before_all():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 5000, cs, ci) is None


def test_find_locus_between_intervals():
    cs, ci = make_intervals()
    assert csl.find_locus("10", 15000, cs, ci) is None


def test_find_locus_unknown_chrom():
    cs, ci = make_intervals()
    assert csl.find_locus("chrX", 10100, cs, ci) is None


# ---------------------------------------------------------------------------
# process_smartseq2 with mocked subprocess
# ---------------------------------------------------------------------------

def make_fasta(tmp_path):
    fa = tmp_path / "rep.fa"
    fa.write_text(">AluSz_dup1::10:10000-10200(+)\nACGT\n>L1_dup2::10:20000-20500(+)\nACGT\n")
    return str(fa)


def _mock_proc(lines):
    proc = MagicMock()
    proc.stdout = iter([l.encode() for l in lines])
    return proc


def test_process_smartseq2_assigns_reads(tmp_path):
    chrom_starts, chrom_intervals = csl.load_intervals(make_fasta(tmp_path))
    sam = [
        "cell_001_r1_AluSz_dup1\t0\t10\t10050\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\n",
        "cell_001_r2_AluSz_dup1\t0\t10\t10080\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\n",
        "cell_002_r3_L1_dup2\t0\t10\t20100\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\n",
    ]
    with patch("count_starsolo_locus.subprocess.Popen", return_value=_mock_proc(sam)):
        counts, cbs, loci = csl.process_smartseq2("/fake.bam", chrom_starts, chrom_intervals, "unique")
    assert "cell_001" in cbs
    assert "cell_002" in cbs
    assert counts["cell_001"]["AluSz_dup1"] == 2
    assert counts["cell_002"]["L1_dup2"] == 1


def test_process_smartseq2_unique_filter(tmp_path):
    chrom_starts, chrom_intervals = csl.load_intervals(make_fasta(tmp_path))
    # NH:i:3 = multimapper, should be dropped in unique mode
    sam = [
        "cell_001_r1_AluSz_dup1\t0\t10\t10050\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:3\n",
    ]
    with patch("count_starsolo_locus.subprocess.Popen", return_value=_mock_proc(sam)):
        counts, cbs, loci = csl.process_smartseq2("/fake.bam", chrom_starts, chrom_intervals, "unique")
    assert len(cbs) == 0


def test_process_smartseq2_multi_mode_keeps_multimappers(tmp_path):
    chrom_starts, chrom_intervals = csl.load_intervals(make_fasta(tmp_path))
    sam = [
        "cell_001_r1_AluSz_dup1\t0\t10\t10050\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:3\n",
    ]
    with patch("count_starsolo_locus.subprocess.Popen", return_value=_mock_proc(sam)):
        counts, cbs, loci = csl.process_smartseq2("/fake.bam", chrom_starts, chrom_intervals, "multi")
    assert "cell_001" in cbs


# ---------------------------------------------------------------------------
# process_chromium with mocked subprocess
# ---------------------------------------------------------------------------

def test_process_chromium_assigns_reads(tmp_path):
    chrom_starts, chrom_intervals = csl.load_intervals(make_fasta(tmp_path))
    sam = [
        "r1\t0\t10\t10050\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\tCB:Z:ACGAACGAAGTGAGCT\tUB:Z:GGATGACGAAGG\n",
        "r2\t0\t10\t10080\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\tCB:Z:ACGAACGAAGTGAGCT\tUB:Z:TTTTTTTTTAAA\n",
        "r3\t0\t10\t20100\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\tCB:Z:CTAAGCCACCCTGAAG\tUB:Z:AAAAAAAAAAAA\n",
    ]
    with patch("count_starsolo_locus.sort_bam_by_cb", return_value="/tmp/fake_sorted.bam"), \
         patch("count_starsolo_locus.subprocess.Popen", return_value=_mock_proc(sam)), \
         patch("os.unlink"):
        counts, cbs, loci = csl.process_chromium(
            "/fake.bam", chrom_starts, chrom_intervals, "unique", 1)
    assert "ACGAACGAAGTGAGCT" in cbs
    assert "CTAAGCCACCCTGAAG" in cbs
    # UMI deduplication: 2 distinct UMIs for cell1/AluSz_dup1
    assert counts["ACGAACGAAGTGAGCT"]["AluSz_dup1"] == 2
    assert counts["CTAAGCCACCCTGAAG"]["L1_dup2"] == 1


def test_process_chromium_drops_no_cb(tmp_path):
    chrom_starts, chrom_intervals = csl.load_intervals(make_fasta(tmp_path))
    sam = [
        "r1\t0\t10\t10050\t255\t90M\t*\t0\t0\tACGT\tFFFF\tNH:i:1\n",  # no CB tag
    ]
    with patch("count_starsolo_locus.sort_bam_by_cb", return_value="/tmp/fake.bam"), \
         patch("count_starsolo_locus.subprocess.Popen", return_value=_mock_proc(sam)), \
         patch("os.unlink"):
        counts, cbs, loci = csl.process_chromium(
            "/fake.bam", chrom_starts, chrom_intervals, "unique", 1)
    assert len(cbs) == 0
