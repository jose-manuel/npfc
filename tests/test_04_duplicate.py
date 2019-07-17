"""
Module test_04_duplicate
====================
Tests for the duplicate module.
"""
# standard
from pathlib import Path
import warnings
# data science
import pandas as pd
# chemoinformatics
from rdkit import RDLogger
# tests
import pytest
from npfc.duplicate import DuplicateFilter
from npfc import load
# configure logging
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
# debug
# import logging
# logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def input_files_dupl():
    p = Path('tests/tmp')
    return [str(f) for f in list(p.glob('test_save_dupl_00[1-4].csv.gz'))]


@pytest.fixture
def ref_file():
    return 'tests/tmp/test_dupl_ref.hdf'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init_ref(ref_file):
    """Make sure ref file is computed during this test."""
    p = Path(ref_file)
    if p.exists():
        p.unlink()
    assert p.exists() is False


def test_remove_dupl(input_files_dupl, ref_file):
    """Remove duplicates accross chunks using a syn file"""

    # without a ref file: checks perfomed only intra chunks
    passed = 0
    filtered = 0
    duplicate_filter = DuplicateFilter(ref_file=None)
    for f in input_files_dupl:
        df = load.file(f)
        df = duplicate_filter.mark_dupl(df)
        passed += len(df[df["status"] == "passed"])
        filtered += len(df[df["status"] == "filtered"])

    assert passed == 6 and filtered == 1

    # with a ref file: checks performed intra and inter chunks
    passed = 0
    filtered = 0
    duplicate_filter = DuplicateFilter(ref_file=ref_file)
    for f in input_files_dupl:
        df = load.file(f)
        df = duplicate_filter.mark_dupl(df)
        passed += len(df[df["status"] == "passed"])
        filtered += len(df[df["status"] == "filtered"])

    assert passed == 4 and filtered == 3
