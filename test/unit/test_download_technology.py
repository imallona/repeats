"""
Unit tests for the library_layout-based download output logic in download_sra.snmk.

The rule produces R2 only for library_layout: paired. These tests verify the
Python expressions that drive that decision without invoking Snakemake.
"""
import pytest


def outputs_for_layout(library_layout, base, sample_id):
    """Mirrors the output: block logic in download_sra_bulk."""
    import os.path as op
    is_paired = library_layout == 'paired'
    R1 = op.join(base, 'data', 'raw', sample_id, f'{sample_id}_1.fastq.gz')
    R2 = ([op.join(base, 'data', 'raw', sample_id, f'{sample_id}_2.fastq.gz')]
          if is_paired else [])
    return R1, R2


def test_single_layout_produces_r1_only():
    r1, r2 = outputs_for_layout('single', '/base', 'sampleA')
    assert r1.endswith('sampleA_1.fastq.gz')
    assert r2 == []


def test_paired_layout_produces_r1_and_r2():
    r1, r2 = outputs_for_layout('paired', '/base', 'sampleA')
    assert r1.endswith('sampleA_1.fastq.gz')
    assert len(r2) == 1
    assert r2[0].endswith('sampleA_2.fastq.gz')


def test_unknown_layout_defaults_to_no_r2():
    _, r2 = outputs_for_layout('unknown', '/base', 'sampleA')
    assert r2 == []


@pytest.mark.parametrize('layout', ['single', 'paired'])
def test_r1_path_consistent_across_layouts(layout):
    r1, _ = outputs_for_layout(layout, '/results', 'mysample')
    assert 'mysample' in r1
    assert r1.endswith('_1.fastq.gz')
