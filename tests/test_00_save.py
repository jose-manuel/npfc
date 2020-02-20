"""
Module test_00_save
====================
Tests for the save module.

The whole reason why tests are numeroted is because I want the save tests to be
executed before the load tests, so that temporary files can be used before they
are erased.
"""
# standard
from pathlib import Path
import shutil
# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import save
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def output_file_prefix():
    return 'tests/tmp/test_save'


@pytest.fixture
def df_mols():
    """An example of a DataFrame with molecules."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'FC1CCCCC1', 'OC1CCCCC1', 'NC1CCCCC1', '[Cl]C1CCCCC1'],
                       'idm': ['mol1', 'mol2', 'mol3', 'mol4', 'mol5'],
                       })
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    df['prop'] = [f"MOL{i}" for i in range(len(df.index))]  # triggers pandas warnings
    return df


@pytest.fixture
def df_mols_dupl():
    """Example of a DataFrame with some duplicate molecules. Supposed to be exported as chunks."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'C1CCCCC1', 'C1CCCCC1', 'NC1CCCCC1', 'FC1CCCCC1', '[Cl]C1CCCCC1', 'C1CCCCC1'],   # 4/7 are dupl
                       })
    df['idm'] = [f"mol_{i}" for i in range(len(df.index))]

    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    return df


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init(output_file_prefix):
    """Initialize the tests/temp folder and check if the saver object is initialized with expected attributes."""
    # init temp folder
    WD = Path(output_file_prefix).resolve().parent
    # delete previous runs
    shutil.rmtree(str(WD))
    # recreate an empty WD
    WD.mkdir(parents=True)


def test_save_compressed(df_mols, output_file_prefix):
    """Save outputs into gzip format. Has to be first as save method removes byproducts."""
    # test simple export but to an archive (for csv only)
    outputs_csv_compressed = save.file(df_mols, output_file_prefix + '.csv.gz')
    assert Path(outputs_csv_compressed[0][0]).is_file() is True and outputs_csv_compressed[0][1] == 5
    outputs_csv_compressed = save.file(df_mols, output_file_prefix + '.sdf.gz')
    assert Path(outputs_csv_compressed[0][0]).is_file() is True and outputs_csv_compressed[0][1] == 5


def test_save_simple(df_mols, output_file_prefix):
    """Test several cases (shuffle, encode) for saving molecules into different formats (hdf, csv)."""
    # hdf
    outputs_hdf = save.file(df_mols, output_file_prefix + '.hdf')
    assert Path(outputs_hdf[0][0]).is_file() is True and outputs_hdf[0][1] == 5
    # csv
    outputs_csv = save.file(df_mols, output_file_prefix + '.csv')
    assert Path(outputs_csv[0][0]).is_file() is True and outputs_csv[0][1] == 5
    # sdf
    outputs_sdf = save.file(df_mols, output_file_prefix + '.sdf')
    assert Path(outputs_sdf[0][0]).is_file() is True and outputs_sdf[0][1] == 5


def test_save_chunks(df_mols, output_file_prefix):

    # test create chunks
    chunk_size = 2
    # hdf
    outputs_hdf = save.file(df_mols, output_file_prefix + '_chunks.hdf', chunk_size=chunk_size)
    assert Path(outputs_hdf[0][0]).is_file() and outputs_hdf[0][1] == 2
    assert Path(outputs_hdf[1][0]).is_file() and outputs_hdf[1][1] == 2
    assert Path(outputs_hdf[2][0]).is_file() and outputs_hdf[2][1] == 1
    # csv
    outputs_csv = save.file(df_mols, output_file_prefix + '_chunks.csv', chunk_size=chunk_size)
    assert Path(outputs_csv[0][0]).is_file() and outputs_csv[0][1] == 2
    assert Path(outputs_csv[1][0]).is_file() and outputs_csv[1][1] == 2
    assert Path(outputs_csv[2][0]).is_file() and outputs_csv[2][1] == 1
    # sdf
    outputs_sdf = save.file(df_mols, output_file_prefix + '_chunks.sdf', chunk_size=chunk_size)
    assert Path(outputs_sdf[0][0]).is_file() and outputs_sdf[0][1] == 2
    assert Path(outputs_sdf[1][0]).is_file() and outputs_sdf[1][1] == 2
    assert Path(outputs_sdf[2][0]).is_file() and outputs_sdf[2][1] == 1


def test_save_func_dupl(df_mols_dupl, output_file_prefix):
    """Save molecules using the save function instead of using a Saver object.
    Produced chunks are used for testing the removal of duplicates."""
    outputs_csv = save.file(df_mols_dupl, output_file_prefix + '_dupl.csv.gz', chunk_size=2)
    assert Path(outputs_csv[0][0]).is_file() and outputs_csv[0][1] == 2
    assert Path(outputs_csv[1][0]).is_file() and outputs_csv[1][1] == 2
    assert Path(outputs_csv[2][0]).is_file() and outputs_csv[2][1] == 2
    assert Path(outputs_csv[3][0]).is_file() and outputs_csv[3][1] == 1
