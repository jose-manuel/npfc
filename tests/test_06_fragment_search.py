"""
tests.test_fcc
~~~~~~~~~~~~~~~
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
logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def df_frags():
    """Two fragments used for most of substructue searches."""
    df = pd.DataFrame([['QA', Chem.MolFromSmiles('C1NCCC1'), {1: '1', 0: '2a', 2: '2b', 3: '3a', 4: '3b'}],
                       ['QB', Chem.MolFromSmiles('C1CCCOC1'), {1: '1', 0: '2a', 2: '2b', 4: '3', 3: '4a', 5: '4b'}],
                       ], columns=['idf', 'mol', '_fcp_labels'], index=['QA', 'QB'])
    return df


@pytest.fixture
def df_mols_match():
    """One molecule with multiple matches."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2N1)C1CCCOC1')], 'inchikey': ['UUKYKIBPEZXAON-UHFFFAOYSA-N']}, index=['MOL1'])


@pytest.fixture
def df_mols_no_match():
    """One molecule without any match"""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(C1CC2CCCCC2S1)C1CCCSC1')], 'inchikey': ['BMTGHHFFWCOFTB-UHFFFAOYSA-N']}, index=['MOL2'])


@pytest.fixture
def df_mols_tautomers_1():
    """Two tautomers of a molecule which result in different fragment hits"""
    return pd.DataFrame([[Chem.MolFromSmiles('O=C1CCCCC1C'), ''],
                         [Chem.MolFromSmiles('OC1=CCCCC1C'), ''],
                         ], columns=['mol', 'inchikey'], index=['TAUT_1A', 'TAUT_1B'])


@pytest.fixture
def df_frags_tautomers_1():
    """Two fragments used for most of substructue searches."""
    df = pd.DataFrame([['QC', Chem.MolFromSmiles('C1CCCCC1=O'), {0: '1a', 1: '2a', 2: '3', 3: '2b', 4: '1b', 5: '4', 6: '5'}],
                       ], columns=['idf', 'mol', '_fcp_labels'])
    df.index = [x for x in df['idf'].astype(str)]
    return df


@pytest.fixture
def df_mols_tautomers_2():
    """Two tautomers of a molecule which result in different fragment hits"""
    return pd.DataFrame([[Chem.MolFromSmiles('CC1(c2cccc(NC(=O)c3ccc(Br)cn3)c2)COC(C)(C(F)(F)F)C(N)=N1'), ''],
                         [Chem.MolFromSmiles('CC1(c2cccc(NC(=O)c3ccc(Br)cn3)c2)COC(C)(C(F)(F)F)C(=N)N1'), ''],
                         ], columns=['mol', 'inchikey'], index=['TAUT_2A', 'TAUT_2B'])


@pytest.fixture
def df_frags_tautomers_2():
    """Two fragments used for most of substructue searches."""
    df = pd.DataFrame([['QD', Chem.MolFromSmiles('C1COCCN1'), {0: '1a', 1: '2a', 2: '3', 3: '2b', 4: '1b', 5: '4'}],
                       ], columns=['idf', 'mol', '_fcp_labels'])
    df.index = [x for x in df['idf'].astype(str)]
    return df


@pytest.fixture
def df_mols_tautomers_3():
    """Two tautomers of a molecule which result in different fragment hits"""
    return pd.DataFrame([[Chem.MolFromSmiles('CCN(CC)c1ccc2nc3ccc(N(CC)CC)cc3[o+]c2c1'), ''],
                         [Chem.MolFromSmiles('CCN(CC)c1ccc2nc3ccc(=[N+](CC)CC)cc-3oc2c1'), ''],
                         ], columns=['mol', 'inchikey'], index=['TAUT_3A', 'TAUT_3B'])


@pytest.fixture
def df_frags_tautomers_3():
    """Two fragments used for most of substructue searches."""
    df = pd.DataFrame([['QE', Chem.MolFromSmiles('c1ccccc1'), {0: '1a', 1: '1b', 2: '1c', 3: '1d', 4: '1e', 5: '1f'}],
                       ], columns=['idf', 'mol', '_fcp_labels'])
    df.index = [x for x in df['idf'].astype(str)]
    return df


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_find_fragments_hit(df_mols_match, df_frags):
    """Find fragments in a molecule that contain them."""
    df_aidxf = fragment_search.get_fragment_hits(df_mols_match, df_frags, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert list(df_aidxf.index) == [0, 1]
    assert list(df_aidxf['idm']) == ['MOL1'] * 2
    assert list(df_aidxf['idf']) == ['QA', 'QB']
    assert list(df_aidxf['_aidxf']) == [(1, 9, 8, 3, 2), (10, 11, 12, 13, 14, 15)]
    assert list(df_aidxf['mol_perc']) == [31.0, 38.0]


def test_find_fragments_no_hit(df_mols_no_match, df_frags):
    """Do not find fragments if the molecule does not contain them."""
    df_aidxf = fragment_search.get_fragment_hits(df_mols_no_match, df_frags, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf.index) == 0


def test_find_fragments_hit_tautomers_1(df_mols_tautomers_1, df_frags_tautomers_1):
    """Check that tautomer-dependant search fails to identify both hits but tautomer-independant searchs succeeds (#1)"""
    # tautomer-dependant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_1, df_frags_tautomers_1, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 1
    assert list(df_aidxf['idm']) == ['TAUT_1A']
    # tautomer-independant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_1, df_frags_tautomers_1, tautomer=True, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 2
    assert list(df_aidxf['idm']) == ['TAUT_1A', 'TAUT_1B']


def test_find_fragments_hit_tautomers_2(df_mols_tautomers_2, df_frags_tautomers_2, col_to_index_mols='', col_to_index_frags=''):
    """Check that tautomer-dependant search fails to identify both hits but tautomer-independant searchs succeeds (#2)"""
    # tautomer-dependant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_2, df_frags_tautomers_2, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 1
    assert list(df_aidxf['idm']) == ['TAUT_2B']
    # tautomer-independant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_2, df_frags_tautomers_2, tautomer=True, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 2
    assert list(df_aidxf['idm']) == ['TAUT_2A', 'TAUT_2B']


def test_find_fragments_hit_tautomers_3(df_mols_tautomers_3, df_frags_tautomers_3, col_to_index_mols='', col_to_index_frags=''):
    """Check that tautomer-independant does not outperform tautomer-dependant on this particular case"""
    # tautomer-dependant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_3, df_frags_tautomers_3, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 3
    assert list(df_aidxf['idm']) == ['TAUT_3A', 'TAUT_3A', 'TAUT_3B']
    # tautomer-independant
    df_aidxf = fragment_search.get_fragment_hits(df_mols_tautomers_3, df_frags_tautomers_3, tautomer=True, col_to_index_mols='', col_to_index_frags='')
    logging.debug(f"Fragment Hits:\n{df_aidxf}\n")
    assert len(df_aidxf) == 3
    assert list(df_aidxf['idm']) == ['TAUT_3A', 'TAUT_3A', 'TAUT_3B']
