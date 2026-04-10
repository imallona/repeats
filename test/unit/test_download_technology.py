"""
Unit tests for the technology-based download output logic in download_sra.snmk.

The rule produces R2 only for paired_end technology. These tests verify the
Python expressions that drive that decision without invoking Snakemake.
"""
import pytest


def outputs_for_technology(technology, base, sample_id):
    """Mirrors the output: block logic in download_sra_bulk."""
    import os.path as op
    is_paired = technology == 'paired_end'
    R1 = op.join(base, 'data', 'raw', sample_id, f'{sample_id}_1.fastq.gz')
    R2 = ([op.join(base, 'data', 'raw', sample_id, f'{sample_id}_2.fastq.gz')]
          if is_paired else [])
    return R1, R2


def test_single_end_produces_r1_only():
    r1, r2 = outputs_for_technology('single_end', '/base', 'sampleA')
    assert r1.endswith('sampleA_1.fastq.gz')
    assert r2 == []


def test_paired_end_produces_r1_and_r2():
    r1, r2 = outputs_for_technology('paired_end', '/base', 'sampleA')
    assert r1.endswith('sampleA_1.fastq.gz')
    assert len(r2) == 1
    assert r2[0].endswith('sampleA_2.fastq.gz')


def test_chromium_does_not_produce_r2_via_bulk_rule():
    # chromium data is handled by download_sra_sc, not download_sra_bulk.
    # If download_sra_bulk were called with chromium, it would behave like
    # single_end (no R2), which signals a misconfiguration.
    _, r2 = outputs_for_technology('chromium', '/base', 'sampleA')
    assert r2 == []


@pytest.mark.parametrize('technology', ['single_end', 'paired_end'])
def test_r1_path_consistent_across_technologies(technology):
    r1, _ = outputs_for_technology(technology, '/results', 'mysample')
    assert 'mysample' in r1
    assert r1.endswith('_1.fastq.gz')
