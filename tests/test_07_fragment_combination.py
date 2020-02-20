"""
tests.test_fragment_combination
===============================
Tests for the npfc.fragment_combination module.
"""

# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import save
from npfc import fragment_search
from npfc import fragment_combination
from npfc import utils
# debug
import logging
logging.basicConfig(level=logging.INFO)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def df_aidxf_fs():
    """Example df_aidxf of a molecule with the fusion spiro fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC2(C1)CCNC2')
    return pd.DataFrame([
                         ['mol_fs', 'QA', 0, frozenset([0, 1, 2, 3, 4]), 50.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_fs', 'QB', 0, frozenset([2, 5, 6, 7, 8, 9]), 60.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fe():
    """Example df_aidxf of a molecule with the fusion edge fragment combination."""
    mol = Chem.MolFromSmiles('C1CC2CCOCC2N1')
    return pd.DataFrame([
                         ['mol_fe', 'QA', 0, frozenset([0, 1, 2, 7, 8]), 56.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_fe', 'QB', 0, frozenset([2, 3, 4, 5, 6, 7]), 67.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fb():
    """Example df_aidxf of a molecule with the fusion bridge fragment combination."""
    mol = Chem.MolFromSmiles('C1NC2COCC1C2')
    return pd.DataFrame([
                         ['mol_fb', 'QA', 0, frozenset([0, 1, 2, 3, 4]), 56.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_fb', 'QB', 0, frozenset([2, 3, 4, 5, 6, 7]), 67.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fl():
    """Example df_aidxf of a molecule with the fusion linker fragment combination."""
    mol = Chem.MolFromSmiles('C1CCCC2COCC(CC1)CN2')
    return pd.DataFrame([
                         ['mol_fl', 'QC', 0, frozenset([0, 1, 2, 3, 7, 8, 9, 10, 11, 12]), 77.0, mol, Chem.MolFromSmiles('C1CCCCNCCCC1')],
                         ['mol_fl', 'QD', 0, frozenset([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]), 85.0, mol, Chem.MolFromSmiles('C1CCCCCOCCCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ca():
    """Example df_aidxf of a molecule with a annulated connection fragment combination."""
    mol = Chem.MolFromSmiles('C1CC2C(CCC3C2CCC2CCOCC32)CN1')
    return pd.DataFrame([
                         ['mol_ca', 'QA', 0, frozenset([10, 11, 12, 13, 14, 15]), 33.0, mol, Chem.MolFromSmiles('C1CCNC1'), 'RRJUGKFDAVWIGB-UHFFFAOYNA-N'],
                         ['mol_ca', 'QB', 0, frozenset([2, 1, 0, 17, 16, 3]), 33.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'RRJUGKFDAVWIGB-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_cm1():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC(C2CCNC2)C1')
    return pd.DataFrame([
                         ['mol_cm', 'QA', 0, frozenset([5, 6, 7, 8, 9]), 45.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cm', 'QB', 0, frozenset([0, 1, 2, 3, 4, 10]), 55.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cm2():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination. One of the fragments has a ring attached to it in between fragments."""
    mol = Chem.MolFromSmiles('C1CC(CCN1)C1CCC2CCOCC2C1')
    return pd.DataFrame([
                         ['mol_cm', 'QA', 0, frozenset([0, 1, 2, 3, 4, 5]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1'), 'WXNLMLWVMZMRGS-UHFFFAOYNA-N'],
                         ['mol_cm', 'QB', 0, frozenset([9, 10, 11, 12, 13, 14]), 38.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'WXNLMLWVMZMRGS-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_cm3():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination. The shortest path goes through another ring but does not form a ring system."""
    mol = Chem.MolFromSmiles('C1C(CC1C1CCCOC1)C1CCNCC1')
    return pd.DataFrame([
                         ['mol_cm', 'QA', 0, frozenset([4, 5, 6, 7, 8, 9]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1'), 'FRFCBPABRGKJMJ-UHFFFAOYNA-N'],
                         ['mol_cm', 'QB', 0, frozenset([10, 11, 12, 13, 14, 15]), 38.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'FRFCBPABRGKJMJ-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_cbs():
    """Example df_aidxf of a molecule with a bipodal connection spiro fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC2(C1)CCC1CNCC12')
    return pd.DataFrame([
                         ['mol_cbs', 'QA', 0, frozenset([8, 9, 10, 11, 12]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cbs', 'QB', 0, frozenset([0, 1, 2, 3, 4, 5]), 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbe():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal edge."""
    mol = Chem.MolFromSmiles('C1NCC2C1CCC1OCCCC21')
    return pd.DataFrame([
                         ['mol_cbe', 'QA', 0, frozenset([0, 1, 2, 3, 4]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cbe', 'QB', 0, frozenset([7, 8, 9, 10, 10, 12]), 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbb():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal bridged."""
    mol = Chem.MolFromSmiles('C1COC2C3CC(CN3)C2C1')
    return pd.DataFrame([
                         ['mol_cbb', 'QA', 0, frozenset([4, 5, 6, 7, 8]), 45.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cbb', 'QB', 0, frozenset([0, 1, 2, 3, 9, 10]), 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbl():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal linker."""
    mol = Chem.MolFromSmiles('C1CCCC2NCC(CC1)C1CCCOC21')
    return pd.DataFrame([
                         ['mol_cbl', 'QC', 0, frozenset([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), 62.0, mol, Chem.MolFromSmiles('C1CCCCNCCCC1')],
                         ['mol_cbl', 'QB', 0, frozenset([10, 11, 12, 13, 14, 15]), 38.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cts():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal spiro."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NCC4CC2(CO1)C43')
    return pd.DataFrame([
                         ['mol_cts', 'QA', 0, frozenset([5, 6, 7, 8, 13]), 36.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cts', 'QB', 0, frozenset([0, 1, 2, 10, 11, 12]), 43.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cte():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal edge."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NCC4CCC(O1)C2C43')
    return pd.DataFrame([
                         ['mol_cte', 'QA', 0, frozenset([5, 6, 7, 8, 14]), 33.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cte', 'QB', 0, frozenset([0, 1, 2, 11, 12, 13]), 40.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ctb():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal bridged."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NC4CC3C2C4O1')
    return pd.DataFrame([
                         ['mol_ctb', 'QA', 0, frozenset([5, 6, 7, 8, 9]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_ctb', 'QB', 0, frozenset([0, 1, 2, 10, 11, 12]), 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ctl():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal linker."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NC4CCCCCCC3C2C4O1')
    return pd.DataFrame([
                         ['mol_ctl', 'QC', 0, frozenset([5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), 56.0, mol, Chem.MolFromSmiles('C1CCCCNCCCC1')],
                         ['mol_ctl', 'QB', 0, frozenset([0, 1, 2, 15, 16, 17]), 33.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cos():
    """Example df_aidxf of a molecule with a fragment combination of type connection other spiro."""
    mol = Chem.MolFromSmiles('C1OC2CCC3NC4CC1CC21CCC4C31')
    return pd.DataFrame([
                         ['mol_cos', 'QA', 0, frozenset([5, 6, 7, 14, 15]), 38.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_cos', 'QB', 0, frozenset([0, 1, 2, 9, 10, 11]), 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_coe():
    """Example df_aidxf of a molecule with a fragment combination of type connection other edge."""
    mol = Chem.MolFromSmiles('C1OC2CCC3NC4CC1C1CCC4C3C21')
    return pd.DataFrame([
                         ['mol_coe', 'QA', 0, frozenset([5, 6, 7, 13, 14]), 31.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_coe', 'QB', 0, frozenset([0, 1, 2, 9, 10, 15]), 38.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cob():
    """Example df_aidxf of a molecule with a fragment combination of type connection other bridged."""
    mol = Chem.MolFromSmiles('C1OC2CCC3OC4CC1C1CCC(C4)C3C21')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_cob', 'QB', 0, frozenset([5, 6, 7, 13, 14, 15]), 35.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                         ['mol_cob', 'QB', 1, frozenset([0, 1, 2, 9, 10, 16]), 35.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_col():   # TODO: replace one QB with QC
    """Example df_aidxf of a molecule with a fragment combination of type connection other linker."""
    mol = Chem.MolFromSmiles('C1CCCC2COC3CCC4OC(CC1)CC1CCC2C3C14')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_cob', 'QB', 0, frozenset([16, 15, 12, 11, 10, 21]), 35.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                         ['mol_cob', 'QB', 1, frozenset([19, 4, 5, 6, 7, 20]), 35.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ffs():
    """Example df_aidxf of a molecule with a fragment combination of type fusion false positive substructure."""
    mol = Chem.MolFromSmiles('CCC1CCNC1')
    return pd.DataFrame([
                         ['mol_ffs', 'QA', 0, frozenset([2, 3, 4, 5, 6]), 71.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_ffs', 'QE', 1, frozenset([1, 2, 3, 4, 5, 6]), 86.0, mol, Chem.MolFromSmiles('CC1CCNC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cfc():
    """Example df_aidxf of a molecule with a fragment combination of type connection false positive cutoff."""
    mol = Chem.MolFromSmiles('C(CCC1CCCOC1)CC1CCNC1')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_ffs', 'QA', 0, frozenset([10, 11, 12, 13, 14]), 71.0, mol, Chem.MolFromSmiles('C1CCNC1')],
                         ['mol_ffs', 'QB', 1, frozenset([3, 4, 5, 6, 7, 8]), 86.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ffo():
    """Example df_aidxf of a molecule with a fragment combination of type fusion false positive overlap."""
    mol = Chem.MolFromSmiles('C1NCC2CC3COCCC3CC12')
    return pd.DataFrame([
                         ['mol_ffo', 'QG', 0, frozenset([0, 1, 2, 3, 4, 5, 10, 11, 12]), 69.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12')],
                         ['mol_ffo', 'QH', 1, frozenset([3, 4, 5, 6, 7, 8, 9, 10, 11, 12]), 77.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1')],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_macro_cbs():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbs. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CCC2CCCC3(CCCOC3)CCC3NCCCC3CCCCC2C1')
    return pd.DataFrame([
                         ['macro_cbs', 'QA', 0, frozenset([15, 16, 17, 18, 19, 20]), 22.0, mol, Chem.MolFromSmiles('C1CCNCC1'), 'LGLPBIYRGSLKKJ-UHFFFAOYNA-N'],
                         ['macro_cbs', 'QB', 1, frozenset([7, 8, 9, 10, 11, 12]), 22.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'LGLPBIYRGSLKKJ-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_macro_cbe():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbe. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CCC2CCCC3CCOCC3CCC3NCCCC3CCCCC2C1')
    return pd.DataFrame([
                         ['macro_cbe', 'QA', 0, frozenset([19, 18, 17, 16, 15, 20]), 22.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12'), 'FQSALWOQPKLPEO-UHFFFAOYNA-N'],
                         ['macro_cbe', 'QB', 1, frozenset([7, 8, 9, 10, 11, 12]), 22.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1'), 'FQSALWOQPKLPEO-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_macro_cbb():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbb. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CC2CC(CCC3CC(CCO3)CCCCC3CCC(CC3)CCCCCC2)N1')
    return pd.DataFrame([
                         ['macro_cbb', 'QA', 0, frozenset([2, 1, 0, 29, 4, 3]), 20.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12'), 'DYGMURSGZXWMST-UHFFFAOYNA-N'],
                         ['macro_cbb', 'QB', 1, frozenset([9, 8, 7, 12, 11, 10]), 20.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1'), 'DYGMURSGZXWMST-UHFFFAOYNA-N'],
                        ], columns=['idm', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_classify_fs(df_aidxf_fs):
    """Check if fusion spiro fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fs)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'spiro' and result['subtype'] == '' and result['abbrev'] == 'fs'


def test_classify_fe(df_aidxf_fe):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fe)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'edge' and result['subtype'] == '' and result['abbrev'] == 'fe'


def test_classify_fb(df_aidxf_fb):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fb)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'bridged' and result['subtype'] == '' and result['abbrev'] == 'fb'


def test_classify_fl(df_aidxf_fl):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fl)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'linker' and result['subtype'] == '' and result['abbrev'] == 'fl'


def test_classify_cm1(df_aidxf_cm1):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm1)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cm'


def test_classify_cm2(df_aidxf_cm2):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm2)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cm'


def test_classify_cm3(df_aidxf_cm3):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm3)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cm'


def test_classify_cbs(df_aidxf_cbs):
    """Check if fragment combinations of type connection bipodal spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbs)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'


def test_classify_cbe(df_aidxf_cbe):
    """Check if fragment combinations of type connection bipodal edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbe)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cbe'


def test_classify_cbb(df_aidxf_cbb):
    """Check if fragment combinations of type connection bipodal bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbb)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'


@pytest.mark.skip
def test_classify_cbl(df_aidxf_cbl):
    """Check if fragment combinations of type connection bipodal linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbl)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'linker' and result['abbrev'] == 'cbl'


def test_classify_cts(df_aidxf_cts):
    """Check if fragment combinations of type connection tripodal spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cts)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cts'


def test_classify_cte(df_aidxf_cte):
    """Check if fragment combinations of type connection tripodal edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cte)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cte'


def test_classify_ctb(df_aidxf_ctb):
    """Check if fragment combinations of type connection tripodal bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ctb)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'ctb'


@pytest.mark.skip
def test_classify_ctl(df_aidxf_ctl):
    """Check if fragment combinations of type connection bipodal linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ctl)
    result = df_fcc.iloc[0]
    # currently found as ctr instead of ctl
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'linker' and result['abbrev'] == 'ctl'


def test_classify_cos(df_aidxf_cos):
    """Check if fragment combinations of type connection other spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cos)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'spiro' and result['abbrev'] == 'cos'


def test_classify_coe(df_aidxf_coe):
    """Check if fragment combinations of type connection other edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_coe)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'edge' and result['abbrev'] == 'coe'


def test_classify_cob(df_aidxf_cob):
    """Check if fragment combinations of type connection other bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cob)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'bridged' and result['abbrev'] == 'cob'


@pytest.mark.skip
def test_classify_col(df_aidxf_col):
    """Check if fragment combinations of type connection other linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_col)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'linker' and result['abbrev'] == 'col'


def test_classify_ffs(df_aidxf_ffs):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ffs)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'substructure' and result['abbrev'] == 'ffs'


def test_classify_ffo(df_aidxf_ffo):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ffo)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'overlap' and result['abbrev'] == 'ffo'


def test_classify_cfc(df_aidxf_cfc):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cfc, clear_cfc=False)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'false_positive' and result['subtype'] == 'cutoff' and result['abbrev'] == 'cfc'


def test_classify_cbs_in_macro(df_aidxf_macro_cbs):
    """Check if cbs fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbs)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'


def test_classify_cbe_in_macro(df_aidxf_macro_cbe):
    """Check if cbe fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbe)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cbe'


def test_classify_cbb_in_macro(df_aidxf_macro_cbb):
    """Check if cbb fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbb)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'
