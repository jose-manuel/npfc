"""
Module test_02_load
====================
Tests for the load module.

Uses the files generated by the save test suite.
"""
# standard
from shutil import copy
from itertools import product
# chemoinformatics
from rdkit.Chem import Mol
# tests
import pytest
from npfc import load
# debug
import logging
logging.basicConfig(level=logging.INFO)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def input_csv():
    return "tests/tmp/test_save.csv"


@pytest.fixture
def input_sdf():
    return "tests/tmp/test_save.sdf"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_file_sdf_without_props(input_sdf):
    """Read a sdf file without properties."""
    # without properties
    df = load.file(input_sdf, keep_props=False, col_idm='_Name')
    assert len(df.index) == 5
    assert list(df.columns.values) == ['idm', 'mol']
    assert isinstance(df.iloc[0]['mol'], Mol)

def test_file_sdf_with_props_1(input_sdf):
    """Read a sdf file with properties including idm and using it for identification."""
    df = load.file(input_sdf, keep_props=True, col_idm='idm')
    assert len(df.index) == 5
    cols = sorted(list(df.columns.values))
    assert cols == ['_Name', 'idm', 'mol', 'prop']
    assert isinstance(df.iloc[0]['mol'], Mol)

def test_file_sdf_with_props_2(input_sdf):
    """Read a sdf file with properties including idm but not using it for identification."""
    df = load.file(input_sdf, keep_props=True, col_idm='_Name')
    assert len(df.index) == 5
    cols = sorted(list(df.columns.values))
    assert cols == ['_Name', 'idm', 'idm.1', 'mol', 'prop']
    assert isinstance(df.iloc[0]['mol'], Mol)

def test_file_sdf_gz(input_sdf):
    """Read a compressed sdf file."""
    # without properties
    df = load.file(f"{input_sdf}.gz", keep_props=False, col_idm='_Name')
    assert len(df.index) == 5
    assert list(df.columns.values) == ['idm', 'mol']
    assert isinstance(df.iloc[0]['mol'], Mol)


def test_file_csv(input_csv):
    """Read a csv file with encoded molecules."""
    # without properties
    df = load.file(input_csv, keep_props=False)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm']
    assert isinstance(df.iloc[0]['mol'], Mol)
    # with properties
    df = load.file(input_csv, keep_props=True)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm', 'prop']
    assert isinstance(df.iloc[0]['mol'], Mol)


def test_file_csv_gz(input_csv):
    """Read a csv file with encoded molecules."""
    df = load.file(input_csv + '.gz', keep_props=False)
    assert len(df.index) == 5
    assert list(df.columns.values) == ['mol', 'idm']
    assert isinstance(df.iloc[0]['mol'], Mol)


def test_count_sdf(input_sdf):
    """Count molecules in a SDF file"""
    count = load.count_mols(input_sdf)
    assert count == 5


def test_count_sdf_gz(input_sdf):
    """Count molecules in a compressed SDF file"""
    count = load.count_mols(f"{input_sdf}.gz", keep_uncompressed=True)  # keep uncompressed file in case of repeated tests
    assert count == 5


def test_count_csv(input_csv):
    """Count molecules in a CSV file"""
    count = load.count_mols(input_csv)  # headers not counted
    assert count == 5


def test_count_csv_gz(input_csv):
    """Count molecules in a compressed CSV file"""
    count = load.count_mols(f"{input_csv}.gz", keep_uncompressed=True)  # headers not counted
    assert count == 5
