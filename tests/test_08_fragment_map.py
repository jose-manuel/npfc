"""
tests.test_fragment_map
=======================
Tests for the npfc.fragment_map module.
"""

# data handling
from pandas import DataFrame
# chemoinformatics
from rdkit import Chem
# tests
import pytest
# dev
from npfc import fragment_map
# debug
import logging
logging.basicConfig(level=logging.INFO)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def df_fc_simple():
    """Simplest scenario case with a map of 3 fragments."""
    mol = Chem.MolFromSmiles('C(C1CC(CN1)C1CCCOC1)C1CCSCC1')
    return DataFrame([
                      ['mol_cmo', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4]), frozenset([6, 7, 8, 9, 10, 11]), 11, mol],
                      ['mol_cmo', 'QA', 0, 'QA:0', 'QC', 0, 'QC:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4]), frozenset([12, 13, 14, 15, 16, 17]), 11, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_redundant():
    """Scenario case with a map of 2 redundant fragments ."""
    mol = Chem.MolFromSmiles('C1CNC(C1)C1CCNC1')
    return DataFrame([
                      ['mol_fc_redundant', 'QA', 0, 'QA:0', 'QA', 1, 'QA:1', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4]), frozenset([5, 6, 7, 8, 9]), 10, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_circular():
    """Scenario case with a circular map."""
    mol = Chem.MolFromSmiles('C1NC2CC1CC1COCC(C1)CC1CC(CCS1)C2')
    return DataFrame([
                      ['mol_fc_circular', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4]), frozenset([6, 7, 8, 9, 10, 11]), 20, mol],
                      ['mol_fc_circular', 'QA', 0, 'QA:0', 'QC', 0, 'QC:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4]), frozenset([13, 14, 15, 16, 17, 18]), 20, mol],
                      ['mol_fc_circular', 'QB', 0, 'QA:0', 'QC', 0, 'QC:0', 'cmo', 'connection', 'monopodal', '', frozenset([6, 7, 8, 9, 10, 11]), frozenset([13, 14, 15, 16, 17, 18]), 20, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_independant():
    """Scenario case with 2 independant submaps."""
    mol = Chem.MolFromSmiles('C(CCC1CSC2CCCCC2C1)CC1COCC(CC2CCNC2)C1')
    return DataFrame([
                      ['mol_fc_independant', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cmo', 'connection', 'monopodal', '', frozenset([20, 21, 22, 23, 24]), frozenset([14, 15, 16, 17, 18, 25]), 26, mol],
                      ['mol_fc_independant', 'QC', 0, 'QC:0', 'QD', 0, 'QD:0', 'fed', 'fusion', 'edge', '', frozenset([3, 4, 5, 6, 11, 12]), frozenset([6, 7, 8, 9, 10, 11]), 26, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])

@pytest.fixture
def df_fc_overlap_1():
    """Scenario case with 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1CCCCC1O')
    return DataFrame([
                      ['mol_fc_overlap_1', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 6]), frozenset([1, 2, 3, 4, 5, 6, 7]), 8, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])

@pytest.fixture
def df_fc_overlap_2():
    """Scenario case with 2 fragments overlapping, bound to a common fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CCC1O)C1CCC1')
    return DataFrame([
                      ['mol_fc_overlap_2', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 6]), frozenset([1, 2, 3, 4, 5, 6, 7]), 12, mol],
                      ['mol_fc_overlap_2', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4, 5, 6]), frozenset([8, 9, 10, 11]), 12, mol],
                      ['mol_fc_overlap_2', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 6, 7]), frozenset([8, 9, 10, 11]), 12, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_overlap_3():
    """Scenario case with 2 fragments overlapping, bound to a common fragment, itself bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CCC1O)C1CC(CC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_3', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 6]), frozenset([1, 2, 3, 4, 5, 6, 7]), 16, mol],
                      ['mol_fc_overlap_3', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4, 5, 6]), frozenset([8, 9, 10, 15]), 16, mol],
                      ['mol_fc_overlap_3', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 6, 7]), frozenset([8, 9, 10, 15]), 16, mol],
                      ['mol_fc_overlap_3', 'O3', 0, 'O3:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([8, 9, 10, 15]), frozenset([12, 13, 14]), 16, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])

@pytest.fixture
def df_fc_overlap_4():
    """Scenario case with 2 independant sets of 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1CC(CCCCC2CCC(O)C(N)C2)CC1O')
    return DataFrame([
                      ['mol_fc_overlap_4', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([8, 9, 10, 11, 13, 14, 15]), frozenset([8, 9, 10, 11, 12, 13, 15]), 19, mol],
                      ['mol_fc_overlap_4', 'O4', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([1, 2, 3, 16, 17, 18]), frozenset([0, 1, 2, 3, 16, 17]), 19, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_overlap_5():
    """Scenario case with 2 sets of 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1C(O)CCC1C1CCC(O)C(N)C1')
    return DataFrame([
                      ['mol_fc_overlap_5', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([7, 8, 9, 10, 12, 13, 14]), frozenset([7, 8, 9, 10, 11, 12, 14]), 15, mol],
                      ['mol_fc_overlap_5', 'O1', 0, 'O1:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 12, 13, 14]), frozenset([1, 2, 3, 4, 5, 6]), 15, mol],
                      ['mol_fc_overlap_5', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 12, 13, 14]), frozenset([0, 1, 2, 4, 5, 6]), 15, mol],
                      ['mol_fc_overlap_5', 'O2', 0, 'O2:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 11, 12, 14]), frozenset([1, 2, 3, 4, 5, 6]), 15, mol],
                      ['mol_fc_overlap_5', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 11, 12, 14]), frozenset([0, 1, 2, 4, 5, 6]), 15, mol],
                      ['mol_fc_overlap_5', 'O4', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([1, 2, 3, 4, 5, 6]), frozenset([0, 1, 2, 4, 5, 6]), 15, mol],
                      ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_overlap_6():
    """Scenario case with 2 sets of 2 overlapping fragments, 1 set being bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1C(O)CCC1C1CC(N)C(O)CC1CCC1CC(CCC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_6', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([7, 8, 9, 10, 11, 13, 14]), frozenset([7, 8, 9, 11, 12, 13, 14]), 26, mol],
                      ['mol_fc_overlap_6', 'O1', 0, 'O1:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 11, 13, 14]), frozenset([1, 2, 3, 4, 5, 6]), 26, mol],
                      ['mol_fc_overlap_6', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 11, 13, 14]), frozenset([0, 1, 2, 4, 5, 6]), 26, mol],
                      ['mol_fc_overlap_6', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 10, 11, 13, 14]), frozenset([17, 18, 19, 25]), 26, mol],
                      ['mol_fc_overlap_6', 'O2', 0, 'O2:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 11, 12, 13, 14]), frozenset([1, 2, 3, 4, 5, 6]), 26, mol],
                      ['mol_fc_overlap_6', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([7, 8, 9, 11, 12, 13, 14]), frozenset([0, 1, 2, 4, 5, 6]), 26, mol],
                      ['mol_fc_overlap_6', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'false_positive', '', frozenset([7, 8, 9, 11, 12, 13, 14]), frozenset([17, 18, 19, 25]), 26, mol],
                      ['mol_fc_overlap_6', 'O4', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'false_positive', frozenset([1, 2, 3, 4, 5, 6]), frozenset([0, 1, 2, 4, 5, 6]), 26, mol],
                      ['mol_fc_overlap_6', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cmo', 'connection', 'monopodal', '', frozenset([17, 18, 19, 25]), frozenset([24, 22, 23]), 26, mol],
                  ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])



@pytest.fixture
def df_fc_overlap_7():
    """Scenario case with 3 overlapping fragments, bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CCC1')
    return DataFrame([
                      ['mol_fc_overlap_7', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([1, 2, 3, 4, 5, 7, 8]), 13, mol],
                      ['mol_fc_overlap_7', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([1, 2, 3, 4, 5, 6, 7]), 13, mol],
                      ['mol_fc_overlap_7', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([9, 10, 11, 12]), 13, mol],
                      ['mol_fc_overlap_7', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([1, 2, 3, 4, 5, 7, 8]), frozenset([1, 2, 3, 4, 5, 6, 7]), 13, mol],
                      ['mol_fc_overlap_7', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 7, 8]), frozenset([9, 10, 11, 12]), 13, mol],
                      ['mol_fc_overlap_7', 'O3', 0, 'O3:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 6, 7]), frozenset([9, 10, 11, 12]), 13, mol],
                  ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_overlap_8():
    """Scenario case with 3 overlapping fragments, bound to a common combination of 2 fragments."""
    mol = Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CC(CCC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_8', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([1, 2, 3, 4, 5, 7, 8]), 18, mol],
                      ['mol_fc_overlap_8', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([1, 2, 3, 4, 5, 6, 7]), 18, mol],
                      ['mol_fc_overlap_8', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 3, 4, 5, 7]), frozenset([9, 10, 11, 17]), 18, mol],
                      ['mol_fc_overlap_8', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([1, 2, 3, 4, 5, 7, 8]), frozenset([1, 2, 3, 4, 5, 6, 7]), 18, mol],
                      ['mol_fc_overlap_8', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 7, 8]), frozenset([9, 10, 11, 17]), 18, mol],
                      ['mol_fc_overlap_8', 'O3', 0, 'O3:0', 'O6', 0, 'O6:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 5, 6, 7]), frozenset([9, 10, 11, 17]), 18, mol],
                      ['mol_fc_overlap_8', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cmo', 'connection', 'monopodal', '', frozenset([9, 10, 11, 17]), frozenset([16, 14, 15]), 18, mol],
                  ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


@pytest.fixture
def df_fc_overlap_9():
    """Scenario case with 3 sets of 2 overlapping fragments, bound to a common combination of 2 redundqnt fragments."""
    mol = Chem.MolFromSmiles('NC1C(O)C(CCCC2CC2CCC2CC2)C1CCC1CC(C(N)C1O)C1CCC(O)C(N)C1')
    return DataFrame([
                      ['mol_fc_overlap_9', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([32, 33, 26, 27, 28, 29, 31]), frozenset([33, 26, 27, 28, 29, 30, 31]), 34, mol],
                      ['mol_fc_overlap_9', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([32, 33, 26, 27, 28, 29, 31]), frozenset([19, 20, 21, 22, 23, 24]), 34, mol],
                      ['mol_fc_overlap_9', 'O1', 0, 'O1:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([32, 33, 26, 27, 28, 29, 31]), frozenset([19, 20, 21, 22, 24, 25]), 34, mol],
                      ['mol_fc_overlap_9', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([33, 26, 27, 28, 29, 30, 31]), frozenset([19, 20, 21, 22, 23, 24]), 34, mol],
                      ['mol_fc_overlap_9', 'O2', 0, 'O2:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([33, 26, 27, 28, 29, 30, 31]), frozenset([19, 20, 21, 22, 24, 25]), 34, mol],
                      ['mol_fc_overlap_9', 'O7', 0, 'O7:0', 'O7', 1, 'O7:1', 'cmo', 'connection', 'monopodal', '', frozenset([8, 9, 10]), frozenset([13, 14, 15]), 34, mol],
                      ['mol_fc_overlap_9', 'O7', 0, 'O7:0', 'O8', 0, 'O8:0', 'cmo', 'connection', 'monopodal', '', frozenset([8, 9, 10]), frozenset([0, 1, 2, 4, 16]), 34, mol],
                      ['mol_fc_overlap_9', 'O7', 0, 'O7:0', 'O9', 0, 'O9:0', 'cmo', 'connection', 'monopodal', '', frozenset([8, 9, 10]), frozenset([1, 2, 3, 4, 16]), 34, mol],
                      ['mol_fc_overlap_9', 'O8', 0, 'O8:0', 'O9', 0, 'O9:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([0, 1, 2, 4, 16]), frozenset([1, 2, 3, 4, 16]), 34, mol],
                      ['mol_fc_overlap_9', 'O8', 0, 'O8:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 4, 16]), frozenset([19, 20, 21, 22, 23, 24]), 34, mol],
                      ['mol_fc_overlap_9', 'O8', 0, 'O8:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([0, 1, 2, 4, 16]), frozenset([19, 20, 21, 22, 24, 25]), 34, mol],
                      ['mol_fc_overlap_9', 'O9', 0, 'O9:0', 'O5', 0, 'O5:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 16]), frozenset([19, 20, 21, 22, 23, 24]), 34, mol],
                      ['mol_fc_overlap_9', 'O9', 0, 'O9:0', 'O4', 0, 'O4:0', 'cmo', 'connection', 'monopodal', '', frozenset([1, 2, 3, 4, 16]), frozenset([19, 20, 21, 22, 24, 25]), 34, mol],
                      ['mol_fc_overlap_9', 'O5', 0, 'O5:0', 'O4', 0, 'O4:0', 'ffo', 'fusion', 'false_positive', 'overlap', frozenset([19, 20, 21, 22, 23, 24]), frozenset([19, 20, 21, 22, 24, 25]), 34, mol],
              ], columns=['idm', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'abbrev', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_map_fc_simple(df_fc_simple):
    """Check if """
    # default: 3 <= n <= 5
    df_map = fragment_map.generate(df_fc_simple)
    result = df_map.iloc[0]
    assert result['fmap_str'] == 'QA:0[cmo]QB:0-QA:0[cmo]QC:0' and result['nfrags_u'] == 3
    # when n > max (2)
    df_map = fragment_map.generate(df_fc_simple, max_frags=2)
    assert len(df_map.index) == 0
    # when n < min (4)
    df_map = fragment_map.generate(df_fc_simple, min_frags=4)
    assert len(df_map.index) == 0


def test_map_fc_redundant(df_fc_redundant):
    """"""
    df_map = fragment_map.generate(df_fc_redundant)
    result = df_map.iloc[0]
    assert result['fmap_str'] == 'QA:0[cmo]QA:1' and result['nfrags'] == 2 and result['nfrags_u'] == 1


def test_map_fc_circular(df_fc_circular):
    """"""
    df_map = fragment_map.generate(df_fc_circular)
    result = df_map.iloc[0]
    assert result['fmap_str'] == 'QA:0[cmo]QB:0-QA:0[cmo]QC:0-QA:0[cmo]QC:0' and result['nfrags'] == 3 and result['nfrags_u'] == 3


def test_map_fc_independant(df_fc_independant):
    """"""
    df_map = fragment_map.generate(df_fc_independant)
    assert len(df_map.index) == 2
    result1 = df_map.iloc[0]
    assert result1['fmap_str'] == 'QA:0[cmo]QB:0' and result1['nfrags'] == 2 and result1['nfrags_u'] == 2
    result2 = df_map.iloc[1]
    assert result2['fmap_str'] == 'QC:0[fed]QD:0' and result2['nfrags'] == 2 and result2['nfrags_u'] == 2


def test_map_fc_overlap_1(df_fc_overlap_1):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_1)
    assert len(df_map.index) == 0


def test_map_fc_overlap_2(df_fc_overlap_2):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_2)
    assert len(df_map.index) == 2
    assert sorted(list(df_map['fmap_str'].values)) == ['O1:0[cmo]O3:0',
                                                       'O2:0[cmo]O3:0'
                                                       ]

def test_map_fc_overlap_3(df_fc_overlap_3):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_3)
    assert len(df_map.index) == 2
    assert sorted(list(df_map['fmap_str'].values)) == ['O3:0[cmo]O4:0-O1:0[cmo]O3:0',
                                                       'O3:0[cmo]O4:0-O2:0[cmo]O3:0'
                                                       ]

def test_map_fc_overlap_4(df_fc_overlap_4):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_4)
    assert len(df_map.index) == 0


def test_map_fc_overlap_5(df_fc_overlap_5):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_5)
    assert len(df_map.index) == 4
    # test each of the fmaps
    # for i in range(len(df_map.index)):
    #     print(df_map.iloc[i]['fmap_str'])
    assert sorted(list(df_map['fmap_str'].values)) == ['O1:0[cmo]O4:0',
                                                       'O1:0[cmo]O5:0',
                                                       'O2:0[cmo]O4:0',
                                                       'O2:0[cmo]O5:0'
                                                       ]

def test_map_fc_overlap_6(df_fc_overlap_6):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_6)
    assert len(df_map.index) == 4
    assert sorted(list(df_map['fmap_str'].values)) == ['O6:0[cmo]O7:0-O1:0[cmo]O4:0-O1:0[cmo]O6:0',
                                                       'O6:0[cmo]O7:0-O1:0[cmo]O5:0-O1:0[cmo]O6:0',
                                                       'O6:0[cmo]O7:0-O2:0[cmo]O4:0-O2:0[cmo]O6:0',
                                                       'O6:0[cmo]O7:0-O2:0[cmo]O5:0-O2:0[cmo]O6:0'
                                                       ]


def test_map_fc_overlap_7(df_fc_overlap_7):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_7)
    assert len(df_map.index) == 3
    assert sorted(list(df_map['fmap_str'].values)) == ['O1:0[cmo]O6:0',
                                                       'O2:0[cmo]O6:0',
                                                       'O3:0[cmo]O6:0',
                                                       ]

def test_map_fc_overlap_8(df_fc_overlap_8):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_8)
    assert len(df_map.index) == 3
    assert sorted(list(df_map['fmap_str'].values)) == ['O6:0[cmo]O7:0-O1:0[cmo]O6:0',
                                                       'O6:0[cmo]O7:0-O2:0[cmo]O6:0',
                                                       'O6:0[cmo]O7:0-O3:0[cmo]O6:0',
                                                       ]

def test_map_fc_overlap_9(df_fc_overlap_9):
    """"""
    df_map = fragment_map.generate(df_fc_overlap_9)
    assert len(df_map.index) == 8

    assert sorted(list(df_map['fmap_str'].values)) == ['O7:0[cmo]O7:1-O1:0[cmo]O4:0-O7:0[cmo]O8:0-O8:0[cmo]O4:0',
                                                       'O7:0[cmo]O7:1-O1:0[cmo]O4:0-O7:0[cmo]O9:0-O9:0[cmo]O4:0',
                                                       'O7:0[cmo]O7:1-O1:0[cmo]O5:0-O7:0[cmo]O8:0-O8:0[cmo]O5:0',
                                                       'O7:0[cmo]O7:1-O1:0[cmo]O5:0-O7:0[cmo]O9:0-O9:0[cmo]O5:0',
                                                       'O7:0[cmo]O7:1-O2:0[cmo]O4:0-O7:0[cmo]O8:0-O8:0[cmo]O4:0',
                                                       'O7:0[cmo]O7:1-O2:0[cmo]O4:0-O7:0[cmo]O9:0-O9:0[cmo]O4:0',
                                                       'O7:0[cmo]O7:1-O2:0[cmo]O5:0-O7:0[cmo]O8:0-O8:0[cmo]O5:0',
                                                       'O7:0[cmo]O7:1-O2:0[cmo]O5:0-O7:0[cmo]O9:0-O9:0[cmo]O5:0',
                                                       ]
