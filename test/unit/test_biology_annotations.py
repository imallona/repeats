"""
Unit tests for workflow/scripts/biology_annotations.py
"""
import math

import pytest

import biology_annotations as ba


@pytest.mark.parametrize('rep_name, expected', [
    ('L1HS', 'young'),
    ('L1PA2', 'young'),
    ('L1PA7', 'young'),
    ('L1PA8', 'medium'),
    ('L1PA10', 'medium'),
    ('L1PA15', 'medium'),
    ('L1PA16', 'medium'),
    ('L1PA17', 'old'),
    ('L1PA25', 'old'),
    ('L1PB1', 'old'),
    ('L1MA1', 'old'),
    ('L1MA4', 'old'),
    ('L1ME2', 'old'),
    ('AluY', 'young'),
    ('AluYa5', 'young'),
    ('AluSx1', 'medium'),
    ('AluSg', 'medium'),
    ('AluJb', 'old'),
    ('FLAM_A', 'old'),
    ('FRAM', 'old'),
    ('SVA_A', 'young'),
    ('SVA_F', 'young'),
    ('HERVK-int', 'young'),
    ('HERVH-int', 'young'),
    ('LTR5_Hs', 'young'),
    ('MER4A', 'unknown'),
    ('MLT1A', 'unknown'),
    ('(CA)n', 'unknown'),
    ('', 'unknown'),
])
def test_classify_family_age(rep_name, expected):
    assert ba.classify_family_age(rep_name) == expected


def test_classify_family_age_none_input():
    assert ba.classify_family_age(None) == 'unknown'


@pytest.mark.parametrize('milli, expected', [
    (0, 'young'),
    (49, 'young'),
    (50, 'young'),
    (51, 'medium'),
    (150, 'medium'),
    (200, 'medium'),
    (201, 'old'),
    (500, 'old'),
    (1000, 'old'),
    (-1, 'unknown'),
    (None, 'unknown'),
    ('bad', 'unknown'),
])
def test_classify_family_age_from_divergence(milli, expected):
    assert ba.classify_family_age_from_divergence(milli) == expected


def test_classify_family_age_from_divergence_nan():
    assert ba.classify_family_age_from_divergence(float('nan')) == 'unknown'


def test_classify_family_age_from_divergence_accepts_string_number():
    assert ba.classify_family_age_from_divergence('75') == 'medium'


def test_strand_vs_gene_intergenic():
    assert ba.strand_vs_gene('+', []) == 'intergenic'
    assert ba.strand_vs_gene('-', []) == 'intergenic'


def test_strand_vs_gene_sense_only():
    assert ba.strand_vs_gene('+', ['+']) == 'sense'
    assert ba.strand_vs_gene('-', ['-', '-']) == 'sense'


def test_strand_vs_gene_antisense_only():
    assert ba.strand_vs_gene('+', ['-']) == 'antisense'
    assert ba.strand_vs_gene('-', ['+']) == 'antisense'


def test_strand_vs_gene_both():
    assert ba.strand_vs_gene('+', ['+', '-']) == 'both'
    assert ba.strand_vs_gene('-', ['-', '+']) == 'both'


def test_strand_vs_gene_dedupes_overlaps():
    # Many overlapping same-strand genes should still give "sense".
    assert ba.strand_vs_gene('+', ['+', '+', '+']) == 'sense'


@pytest.mark.parametrize('seq, expected', [
    ('', 0),
    ('CGCGCG', 0),
    ('A', 1),
    ('AAA', 3),
    ('AAACCCAAAA', 4),
    ('AACCAAAA', 4),
    ('aaaa', 4),
    ('AANAAAAA', 5),
    ('CCAAACCAAAACCCAA', 4),
    ('TTTT', 0),
])
def test_longest_a_run(seq, expected):
    assert ba.longest_a_run(seq) == expected


def test_longest_a_run_none_input():
    assert ba.longest_a_run(None) == 0


def test_longest_a_run_ignores_case():
    assert ba.longest_a_run('aAaAaA') == 6
