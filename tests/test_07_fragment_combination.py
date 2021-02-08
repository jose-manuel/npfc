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
from npfc import fragment_combination
# debug
import logging
logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def df_aidxf_fs():
    """Example df_aidxf of a molecule with the fusion spiro fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC2(C1)CCNCC2')
    return pd.DataFrame([
                         ['mol_fs', 'XXX', 'QA', 0, [6, 4, 10, 9, 8, 7], 50.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_fs', 'XXX', 'QB', 0, [5, 4, 3, 2, 1, 0], 60.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fe():
    """Example df_aidxf of a molecule with the fusion edge fragment combination."""
    mol = Chem.MolFromSmiles('C1CC2CCOCC2NC1')
    return pd.DataFrame([
                         ['mol_fe', 'XXX', 'QA', 0, [0, 1, 2, 7, 8, 9], 56.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_fe', 'XXX', 'QB', 0, [2, 3, 4, 5, 6, 7], 67.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fb():
    """Example df_aidxf of a molecule with the fusion bridge fragment combination."""
    mol = Chem.MolFromSmiles('C1CNC2COCC1C2')
    return pd.DataFrame([
                         ['mol_fb', 'XXX', 'QA', 0, [0, 7, 8, 3, 2, 1], 56.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_fb', 'XXX', 'QB', 0, [8, 7, 6, 5, 4, 3], 67.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_fl():
    """Example df_aidxf of a molecule with the fusion linker fragment combination."""
    mol = Chem.MolFromSmiles('C1CCCC2COCC(CC1)CN2')
    return pd.DataFrame([
                         ['mol_fl', 'XXX', 'QC', 0, [0, 1, 2, 3, 4, 12, 11, 8, 9, 10], 77.0, mol, Chem.MolFromSmiles('C1CCCCNCCCC1')],
                         ['mol_fl', 'XXX', 'QD', 0, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 85.0, mol, Chem.MolFromSmiles('C1CCCCCOCCCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ca():
    """Example df_aidxf of a molecule with a annulated connection fragment combination."""
    mol = Chem.MolFromSmiles('C1CNCC2C1CCC1C2CCC2CCOCC12')
    return pd.DataFrame([
                         ['mol_ca', 'XXX', 'QA', 0, [0, 5, 4, 3, 2, 1], 33.0, mol, Chem.MolFromSmiles('C1CCCNC1'), 'RRJUGKFDAVWIGB-UHFFFAOYNA-N'],
                         ['mol_ca', 'XXX', 'QB', 0, [12, 13, 14, 15, 16, 17], 33.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'RRJUGKFDAVWIGB-UHFFFAOYNA-N'],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_cm1():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC(C2CCCNC2)C1')
    return pd.DataFrame([
                         ['mol_cm', 'XXX', 'QA', 0, [5, 6, 7, 8, 9, 10], 45.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cm', 'XXX', 'QB', 0, [11, 4, 3, 2, 1, 0], 55.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cm2():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination. One of the fragments has a ring attached to it in between fragments."""
    mol = Chem.MolFromSmiles('C1CC(CCN1)C1CCC2CCOCC2C1')
    return pd.DataFrame([
                         ['mol_cm', 'WXNLMLWVMZMRGS-UHFFFAOYNA-N', 'QA', 0, [1, 2, 3, 4, 5, 0], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cm', 'WXNLMLWVMZMRGS-UHFFFAOYNA-N', 'QB', 0, [9, 10, 11, 12, 13, 14], 38.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cm3():
    """Example df_aidxf of a molecule with a monopodal connection fragment combination. The shortest path goes through another ring but does not form a ring system."""
    mol = Chem.MolFromSmiles('C1C(CC1C1CCCOC1)C1CCNCC1')
    return pd.DataFrame([
                         ['mol_cm', 'FRFCBPABRGKJMJ-UHFFFAOYNA-N', 'QA', 0, [11, 10, 15, 14, 13, 12], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cm', 'FRFCBPABRGKJMJ-UHFFFAOYNA-N', 'QB', 0, [5, 4, 9, 8, 7, 6], 38.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbs():
    """Example df_aidxf of a molecule with a bipodal connection spiro fragment combination."""
    mol = Chem.MolFromSmiles('C1COCC2(C1)CCC1CCNCC12')
    return pd.DataFrame([
                         ['mol_cbs', 'XXX', 'QA', 0, [9, 8, 13, 12, 11, 10], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cbs', 'XXX', 'QB', 0, [5, 4, 3, 2, 1, 0], 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbe():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal edge."""
    mol = Chem.MolFromSmiles('C1CNCC2C1CCC1OCCCC21')
    return pd.DataFrame([
                         ['mol_cbe', 'XXX', 'QA', 0, [0, 5, 4, 3, 2, 1], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cbe', 'XXX', 'QB', 0, [12, 11, 10, 9, 8, 13], 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbb():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal bridged."""
    mol = Chem.MolFromSmiles('C1COC2C3CC(CCN3)C2C1')
    return pd.DataFrame([
                         ['mol_cbb', 'XXX', 'QA', 0, [5, 6, 7, 8, 9, 4], 45.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cbb', 'XXX', 'QB', 0, [11, 10, 3, 2, 1, 0], 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cbl():
    """Example df_aidxf of a molecule with a fragment combination of type connection bipodal linker."""
    mol = Chem.MolFromSmiles('C1CC2CC1CCCCC1CCC(CCCC3CCC2C3)NC1')
    return pd.DataFrame([
                         ['mol_cbl', 'XXX', 'QA', 0, [10, 9, 22, 21, 12, 11], 62.0, mol, Chem.MolFromSmiles('C1CCNCC1')],
                         ['mol_cbl', 'XXX', 'QI', 0, [0, 4, 3, 2, 1, 19, 18, 17, 16, 20], 38.0, mol, Chem.MolFromSmiles('C1CCC(C1)C1CCCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cts():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal spiro."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NCCC4CC2(CO1)C43')
    return pd.DataFrame([
                         ['mol_cts', 'XXX', 'QA', 0, [8, 9, 14, 5, 6, 7], 36.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cts', 'XXX', 'QB', 0, [2, 1, 0, 13, 12, 11], 43.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cte():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal edge."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NCCC4CCC(O1)C2C43')
    return pd.DataFrame([
                         ['mol_cte', 'XXX', 'QA', 0, [8, 9, 15, 5, 6, 7], 33.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cte', 'XXX', 'QB', 0, [2, 1, 0, 13, 12, 14], 40.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ctb():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal bridged."""
    mol = Chem.MolFromSmiles('C1CC2CCC3NCC4CC3C2C4O1')
    return pd.DataFrame([
                         ['mol_ctb', 'XXX', 'QA', 0, [8, 9, 10, 5, 6, 7], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_ctb', 'XXX', 'QB', 0, [2, 1, 0, 13, 12, 11], 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ctl():
    """Example df_aidxf of a molecule with a fragment combination of type connection tripodal linker."""
    mol = Chem.MolFromSmiles('C1CC2CC1CCCCC1CCNC3CCCCCCC4CC2CC4CCCC13')
    return pd.DataFrame([
                         ['mol_ctl', 'XXX', 'QA', 0, [9, 10, 11, 12, 13, 28], 56.0, mol, Chem.MolFromSmiles('C1CCNCC1')],
                         ['mol_ctl', 'XXX', 'QI', 0, [0, 4, 3, 2, 1, 22, 21, 20, 24, 23], 33.0, mol, Chem.MolFromSmiles('C1CCC(C1)C1CCCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cos():
    """Example df_aidxf of a molecule with a fragment combination of type connection other spiro."""
    mol = Chem.MolFromSmiles('C1OC2CCC3NCC4CC1CC21CCC4C31')
    return pd.DataFrame([
                         ['mol_cos', 'XXX', 'QA', 0, [8, 15, 16, 5, 6, 7], 38.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cos', 'XXX', 'QB', 0, [11, 10, 0, 1, 2, 12], 46.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_coe():
    """Example df_aidxf of a molecule with a fragment combination of type connection other edge."""
    mol = Chem.MolFromSmiles('C1OC2CCC3NCC4CC1C1CCC4C3C21')
    return pd.DataFrame([
                         ['mol_coe', 'XXX', 'QA', 0, [8, 14, 15, 5, 6, 7], 31.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_coe', 'XXX', 'QB', 0, [11, 10, 0, 1, 2, 16], 38.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cob():
    """Example df_aidxf of a molecule with a fragment combination of type connection other bridged."""
    mol = Chem.MolFromSmiles('C1CC2C3CNC4CCC5OC(CC1C5C24)C3')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_cob', 'XXX', 'QA', 0, [3, 2, 15, 6, 5, 4], 35.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cob', 'XXX', 'QB', 0, [13, 12, 11, 10, 9, 14], 35.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_col():   # TODO: replace one QB with QC
    """Example df_aidxf of a molecule with a fragment combination of type connection other linker."""
    mol = Chem.MolFromSmiles('C1CC2C3CC4CCCCCCC5NCC6CCCCCCC2C1CCCCC6C5CCCC4C3')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_cob', 'XXX', 'QA', 0, [28, 29, 12, 13, 14, 15], 35.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_cob', 'XXX', 'QI', 0, [0, 23, 22, 2, 1, 3, 4, 5, 33, 34], 35.0, mol, Chem.MolFromSmiles('C1CCC(C1)C1CCCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ffs():
    """Example df_aidxf of a molecule with a fragment combination of type fusion false positive substructure."""
    mol = Chem.MolFromSmiles('CCC1CCCNC1')
    return pd.DataFrame([
                         ['mol_ffs', 'XXX', 'QA', 0, [2, 3, 4, 5, 6, 7], 71.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_ffs', 'XXX', 'QE', 0, [1, 2, 3, 4, 5, 6, 7], 86.0, mol, Chem.MolFromSmiles('CC1CCCNC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_cfc():
    """Example df_aidxf of a molecule with a fragment combination of type connection false positive cutoff."""
    mol = Chem.MolFromSmiles('C(CCC1CCCOC1)CC1CCCNC1')  # same QB twice because QA is too small for this case
    return pd.DataFrame([
                         ['mol_ffs', 'XXX', 'QA', 0, [10, 11, 12, 13, 14, 15], 71.0, mol, Chem.MolFromSmiles('C1CCCNC1')],
                         ['mol_ffs', 'XXX', 'QB', 0, [4, 3, 8, 7, 6, 5], 86.0, mol, Chem.MolFromSmiles('C1CCOCC1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_ffo():
    """Example df_aidxf of a molecule with a fragment combination of type fusion false positive overlap."""
    mol = Chem.MolFromSmiles('C1NCCC2CC3COCCC3CC12')
    return pd.DataFrame([
                         ['mol_ffo', 'XXX', 'QG', 0, [0, 1, 2, 3, 4, 5, 6, 11, 12, 13], 69.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12')],
                         ['mol_ffo', 'XXX', 'QH', 0, [13, 4, 5, 6, 7, 8, 9, 10, 11, 12], 77.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1')],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag']
                        )


@pytest.fixture
def df_aidxf_macro_cbs():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbs. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CCC2CCCC3(CCCOC3)CCC3NCCCC3CCCCC2C1')
    return pd.DataFrame([
                         ['macro_cbs', 'XXX', 'QA', 0, [19, 18, 17, 16, 15, 20], 22.0, mol, Chem.MolFromSmiles('C1CCNCC1'), 'LGLPBIYRGSLKKJ-UHFFFAOYNA-N'],
                         ['macro_cbs', 'XXX', 'QB', 0, [8, 7, 12, 11, 10, 9], 22.0, mol, Chem.MolFromSmiles('C1CCOCC1'), 'LGLPBIYRGSLKKJ-UHFFFAOYNA-N'],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_macro_cbe():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbe. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CCC2CCCC3CCOCC3CCC3NCCCC3CCCCC2C1')
    return pd.DataFrame([
                         ['macro_cbe', 'XXX', 'QA', 0, [19, 18, 17, 16, 15, 20], 22.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12'), 'FQSALWOQPKLPEO-UHFFFAOYNA-N'],
                         ['macro_cbe', 'XXX', 'QB', 0, [7, 8, 9, 10, 11, 12], 22.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1'), 'FQSALWOQPKLPEO-UHFFFAOYNA-N'],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


@pytest.fixture
def df_aidxf_macro_cbb():
    """Example df_aidxf of a macrocycle with a fragment combination of type cbb. This tests that redundant intermediary rings are filtered out."""
    mol = Chem.MolFromSmiles('C1CC2CC(CCC3CC(CCO3)CCCCC3CCC(CC3)CCCCCC2)N1')
    return pd.DataFrame([
                         ['macro_cbb', 'XXX', 'QA', 0, [2, 1, 0, 29, 4, 3], 20.0, mol, Chem.MolFromSmiles('C1NCC2CCCCC12'), 'DYGMURSGZXWMST-UHFFFAOYNA-N'],
                         ['macro_cbb', 'XXX', 'QB', 0, [9, 8, 7, 12, 11, 10], 20.0, mol, Chem.MolFromSmiles('C1CCC2COCCC2C1'), 'DYGMURSGZXWMST-UHFFFAOYNA-N'],
                        ], columns=['idm', 'inchikey', 'idf', 'idf_idx', '_aidxf', 'mol_perc', 'mol', 'mol_frag', 'inchikey']
                        )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_classify_fs(df_aidxf_fs):
    """Check if fusion spiro fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fs)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1[fs]QB:0@1'


def test_classify_fe(df_aidxf_fe):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fe)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@2,3[fe]QB:0@0,5'


def test_classify_fb(df_aidxf_fb):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fb)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2,3[fb]QB:0@0,1,5'


def test_classify_fl(df_aidxf_fl):
    """Check if fusion edge fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_fl)
    assert df_fcc.iloc[0]['fc'] == 'QC:0@0,1,2,3,4,7,8,9[fl]QD:0@0,1,2,3,4,8,9,10'


def test_classify_ca(df_aidxf_ca):
    """Check if connection anulated fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ca)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2[ca]QB:0@0,5'


def test_classify_cm1(df_aidxf_cm1):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm1)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0[cm]QB:0@1'


def test_classify_cm2(df_aidxf_cm2):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm2)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1[cm]QB:0@5'


def test_classify_cm3(df_aidxf_cm3):
    """Check if connection monopodal fragment combinations are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cm3)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1[cm]QB:0@1'


def test_classify_cbs(df_aidxf_cbs):
    """Check if fragment combinations of type connection bipodal spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbs)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2[cbs]QB:0@1'


def test_classify_cbe(df_aidxf_cbe):
    """Check if fragment combinations of type connection bipodal edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbe)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2[cbe]QB:0@4,5'


def test_classify_cbb(df_aidxf_cbb):
    """Check if fragment combinations of type connection bipodal bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbb)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,5[cbb]QB:0@1,2'



def test_classify_cbl(df_aidxf_cbl):
    """Check if fragment combinations of type connection bipodal linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cbl)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,4,5[cbl]QI:0@1,2,3,5,8,9'


def test_classify_cts(df_aidxf_cts):
    """Check if fragment combinations of type connection tripodal spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cts)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2,3[cts]QB:0@0,5'


def test_classify_cte(df_aidxf_cte):
    """Check if fragment combinations of type connection tripodal edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cte)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@1,2,3[cte]QB:0@0,4,5'


def test_classify_ctb(df_aidxf_ctb):
    """Check if fragment combinations of type connection tripodal bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ctb)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,3[ctb]QB:0@0,4,5'


def test_classify_ctl(df_aidxf_ctl):
    """Check if fragment combinations of type connection bipodal linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ctl)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,4,5[ctl]QI:0@1,2,3,5,7,8,9'


def test_classify_cos(df_aidxf_cos):
    """Check if fragment combinations of type connection other spiro are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cos)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,3[cos]QB:0@0,1,4,5'


def test_classify_coe(df_aidxf_coe):
    """Check if fragment combinations of type connection other edge are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_coe)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,3[coe]QB:0@0,1,4,5'


def test_classify_cob(df_aidxf_cob):
    """Check if fragment combinations of type connection other bridged are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cob)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,3[cob]QB:0@0,1,2,4,5'


def test_classify_col(df_aidxf_col):  # identified as cob instead of col
    """Check if fragment combinations of type connection other linker are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_col)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,5[col]QI:0@1,2,3,5,7,8,9'


def test_classify_ffs(df_aidxf_ffs):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ffs)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,1,2,3,4,5[ffs]QE:0@1,2,3,4,5,6'


def test_classify_ffo(df_aidxf_ffo):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_ffo)
    assert df_fcc.iloc[0]['fc'] == 'QG:0@4,5,6,7,8,9[ffo]QH:0@0,1,2,3,8,9'


def test_classify_cfc(df_aidxf_cfc):
    """Check if fragment combinations of type fusion false positive substructure are identified."""
    df_fcc = fragment_combination.classify_df(df_aidxf_cfc, clear_cfc=False)
    assert df_fcc.iloc[0]['fc'] == 'QA:0[cfc]QB:0'


def test_classify_cbs_in_macro(df_aidxf_macro_cbs):
    """Check if cbs fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbs)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@4,5[cbs]QB:0@1'


def test_classify_cbe_in_macro(df_aidxf_macro_cbe):
    """Check if cbe fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbe)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@4,5[cbe]QB:0@0,5'


def test_classify_cbb_in_macro(df_aidxf_macro_cbb):
    """Check if cbb fragment combinations are correctly identified in macrocycles."""
    df_fcc = fragment_combination.classify_df(df_aidxf_macro_cbb)
    assert df_fcc.iloc[0]['fc'] == 'QA:0@0,4,5[cbb]QB:0@0,1,2'
