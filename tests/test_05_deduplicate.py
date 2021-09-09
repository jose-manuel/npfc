"""
Module test_04_duplicate
====================
Tests for the duplicate module.
"""
# standard
from pathlib import Path
# tests
import pytest
from npfc import load
from npfc import deduplicate
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def input_files_dupl():
    p = Path('tests/tmp')
    return [str(f) for f in list(p.glob('test_save_dupl_00[1-4].csv.gz'))]


@pytest.fixture
def ref_file_csv():
    return 'tests/tmp/test_dupl_ref.csv'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_remove_dupl_wo_ref_file(input_files_dupl):
    """Remove duplicates found in chunks individually without using a reference file."""

    # without a ref file: checks perfomed only intra chunks
    passed = 0
    filtered = 0
    total = 0
    for f in input_files_dupl:
        df = load.file(f)
        total_curr = len(df.index)
        df = deduplicate.filter_duplicates(df, ref_file=None)
        passed_curr = len(df.index)
        passed += passed_curr
        filtered += total_curr - passed_curr
        total += total_curr

    assert passed == 6 and filtered == 1


def test_remove_dupl_w_ref_file_csv(input_files_dupl, ref_file_csv):

    # make sure there is no previous ref file existing
    p = Path(ref_file_csv)
    if p.exists():
        p.unlink()
    assert p.exists() is False
    # without a ref file: checks perfomed only intra chunks
    passed = 0
    filtered = 0
    total = 0
    for f in input_files_dupl:
        df = load.file(f)
        total_curr = len(df.index)
        df = deduplicate.filter_duplicates(df, ref_file=ref_file_csv)
        passed_curr = len(df.index)
        passed += passed_curr
        filtered += total_curr - passed_curr
        total += total_curr

    assert Path(ref_file_csv).exists()
    assert passed == 4 and filtered == 3
