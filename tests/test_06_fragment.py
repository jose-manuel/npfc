"""
tests.test_fcc
~~~~~~~~~~
Tests for the npfc.fcc module.
"""

# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import save
from npfc import fragment
from npfc import utils
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

@pytest.fixture
def df_frags():
    """Two fragments used for most of substructue searches."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles(x) for x in ['C1NCCC1', 'C1CCCOC1']]}, index=['QA', 'QB'])


@pytest.fixture
def df_mols_match():
    """One molecule with multiple matches."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2N1)C1CCCOC1')]}, index=['MOL1'])

@pytest.fixture
def df_mols_no_match():
    """One molecule without any match"""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2S1)C1CCCSC1')]}, index=['MOL2'])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_find_fragments_hit(df_mols_match, df_frags):
    """Find fragments in a molecule that contain them."""
    df_aidxf = fragment.find(df_mols_match, df_frags)
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert list(df_aidxf.index) == [0, 1]
    assert list(df_aidxf['idm']) == ['MOL1'] * 2
    assert list(df_aidxf['idf']) == ['QA', 'QB']
    assert list(df_aidxf['_aidxf']) == [{1, 2 ,3 ,8 , 9}, {10, 11, 12, 13 ,14 ,15}]
    assert list(df_aidxf['mol_perc']) == [31.0, 38.0]


def test_find_fragments_no_hit(df_mols_no_match, df_frags):
    """Do not find fragments if the molecule does not contain them."""
    df_aidxf = fragment.find(df_mols_no_match, df_frags)
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf.index) == 0
