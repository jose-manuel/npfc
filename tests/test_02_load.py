"""
Module test_02_load
====================
Tests for the load module.

Uses the files generated by the save test suite.
"""

# standard
from pathlib import Path
# chemoinformatics
from rdkit import Chem
from rdkit import RDLogger
# tests
import pytest
from npfc import load

# disable RDKit logs on the console
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


@pytest.fixture
def input_hdf():
    return "tests/tmp/test_save_simple.hdf"


@pytest.fixture
def input_csv():
    return "tests/tmp/test_save_simple.csv"


@pytest.fixture
def input_sdf():
    return "tests/tmp/test_save_simple.sdf"



def test_from_sdf(input_sdf):
    """Read a hdf file with encoded molecules."""
    # without properties
    df = load.from_sdf(input_sdf, keep_props=False)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm']
    assert isinstance(df.iloc[0]['mol'], Chem.Mol)
    # with properties
    df = load.from_sdf(input_sdf, keep_props=True)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm', '_Name']
    assert isinstance(df.iloc[0]['mol'], Chem.Mol)

def test_from_hdf(input_hdf):
    """Read a hdf file with encoded molecules."""
    df = load.from_hdf(input_hdf)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm']
    assert isinstance(df.iloc[0]['mol'], Chem.Mol)


def test_from_csv(input_csv):
    """Read a csv file with encoded molecules."""
    df = load.from_csv(input_csv)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm']
    assert isinstance(df.iloc[0]['mol'], Chem.Mol)
