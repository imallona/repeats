"""
Snakemake dry-run tests.

These tests run `snakemake --dry-run` on both the main Snakefile and the
test Snakefile to verify that the DAG can be constructed without errors.
They do not execute any actual rules and require only snakemake to be installed.

Marked as `workflow` so they can be excluded with `-m "not workflow"` when
snakemake is not in the PATH.
"""
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
WORKFLOW_DIR = os.path.join(REPO_ROOT, 'workflow')
TEST_WORKFLOW = os.path.join(REPO_ROOT, 'test', 'workflow', 'Snakefile_test')
SMARTSEQ2_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'simulation_smartseq2.yaml')
CHROMIUM_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'simulation_chromium.yaml')
NEG_CTRL_CFG = os.path.join(REPO_ROOT, 'test', 'workflow', 'configs',
                             'test_negative_control.yaml')


def snakemake_available():
    try:
        subprocess.run(['snakemake', '--version'], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


SKIP_IF_NO_SNAKEMAKE = pytest.mark.skipif(
    not snakemake_available(),
    reason='snakemake not found in PATH'
)


def run_dryrun(snakefile, configfile, workdir=WORKFLOW_DIR):
    cmd = [
        'snakemake',
        '-s', snakefile,
        '--configfile', configfile,
        '--dry-run',
        '--cores', '1',
        '--quiet',
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_smartseq2():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), SMARTSEQ2_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for simulation_smartseq2.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_chromium():
    if not os.path.exists(CHROMIUM_CFG):
        pytest.skip('simulation_chromium.yaml not found')
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), CHROMIUM_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for simulation_chromium.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_test_snakefile_dryrun_negative_control():
    r = run_dryrun(TEST_WORKFLOW, NEG_CTRL_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for test negative control config:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_lint():
    """snakemake --lint catches common Snakemake style/correctness issues."""
    cmd = [
        'snakemake',
        '--lint',
        '--configfile', SMARTSEQ2_CFG,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=WORKFLOW_DIR)
    # lint may return non-zero for warnings; we only fail on errors
    # (lint output on stderr distinguishes errors from warnings)
    assert 'Error' not in r.stdout, f'Snakemake lint errors:\n{r.stdout}'
