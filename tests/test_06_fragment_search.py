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
from npfc import fragment_search
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def df_frags():
    """Two fragments used for most of substructue searches."""
    df = pd.DataFrame([['QA', Chem.MolFromSmiles('C1NCCC1'), {1: '1', 0: '2a', 2: '2b', 3: '3a', 4: '3b'}],
                       ['QB', Chem.MolFromSmiles('C1CCCOC1'), {1: '1', 0: '2a', 2: '2b', 4: '3', 3: '4a', 5: '4b'}],
                       ], columns=['idf', 'mol', '_fcp_labels'])
    df.index = [x for x in df['idf'].astype(str)]
    return df


@pytest.fixture
def df_mols_match():
    """One molecule with multiple matches."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2N1)C1CCCOC1')], 'inchikey': ['UUKYKIBPEZXAON-UHFFFAOYSA-N']}, index=['MOL1'])


@pytest.fixture
def df_mols_no_match():
    """One molecule without any match"""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2S1)C1CCCSC1')], 'inchikey': ['BMTGHHFFWCOFTB-UHFFFAOYSA-N']}, index=['MOL2'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_find_fragments_hit(df_mols_match, df_frags):
    """Find fragments in a molecule that contain them."""
    df_aidxf = fragment_search.get_fragment_hits(df_mols_match, df_frags)
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert list(df_aidxf.index) == [0, 1]
    assert list(df_aidxf['idm']) == ['MOL1'] * 2
    assert list(df_aidxf['idf']) == ['QA', 'QB']
    assert list(df_aidxf['_aidxf']) == [(1, 9, 8, 3, 2), (10, 11, 12, 13, 14, 15)]
    assert list(df_aidxf['mol_perc']) == [31.0, 38.0]


def test_find_fragments_no_hit(df_mols_no_match, df_frags):
    """Do not find fragments if the molecule does not contain them."""
    df_aidxf = fragment_search.get_fragment_hits(df_mols_no_match, df_frags)
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf.index) == 0
