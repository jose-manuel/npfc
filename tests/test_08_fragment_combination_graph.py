"""
tests.test_fragment_combination_graph
=======================================
Tests for the npfc.fragment_combination_graph module.
"""

# data handling
from pandas import DataFrame
# chemoinformatics
from rdkit import Chem
# tests
import pytest
# dev
from npfc import fragment_combination_graph
# debug
import logging
logging.basicConfig(level=logging.INFO)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# fragments for testing normal fgraphs
mol_qa = Chem.MolFromSmiles('C1CCNC1')
mol_qb = Chem.MolFromSmiles('C1CCOC1')
mol_qc = Chem.MolFromSmiles('C1CCSC1')
mol_qd = Chem.MolFromSmiles('C1CCCC1')
# fragments for testing overlaps
mol_o1 = Chem.MolFromSmiles('NC1CCCCC1')
mol_o2 = Chem.MolFromSmiles('OC1CCCCC1')
mol_o3 = Chem.MolFromSmiles('SC1CCCCC1')
mol_o4 = Chem.MolFromSmiles('OC1CCCC1')
mol_o5 = Chem.MolFromSmiles('NC1CCCC1')
mol_o6 = Chem.MolFromSmiles('C1CCC1')
mol_o7 = Chem.MolFromSmiles('C1CC1')
mol_o8 = Chem.MolFromSmiles('NC1CCC1')
mol_o9 = Chem.MolFromSmiles('OC1CCC1')


@pytest.fixture
def df_fc_simple():
    """Simplest scenario case with a fcg of 3 fragments."""
    mol = Chem.MolFromSmiles('C(C1CCSCC1)C1CCC(CN1)C1CCCOC1')
    return DataFrame([
                      ['mol_cm', 'XXX', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cm', 'connection', 'monopodal', '', [8, 9, 10, 11, 12, 7], [14, 13, 18, 17, 16, 15], 11, mol, mol_qa, mol_qb, 'QA:0@2[cm]QB:0@1'],
                      ['mol_cm', 'XXX', 'QA', 0, 'QA:0', 'QC', 0, 'QC:0', 'cm', 'connection', 'monopodal', '', [8, 9, 10, 11, 12, 7], [1, 2, 3, 4, 5, 6], 11, mol, mol_qa, mol_qc, 'QA:0@5[cm]QC:0@0'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_redundant():
    """Scenario case with a fcg of 2 redundant fragments ."""
    mol = Chem.MolFromSmiles('C1CCC(NC1)C1CCCNC1')
    return DataFrame([
                      ['mol_fc_redundant', 'XXX', 'QA', 0, 'QA:0', 'QA', 1, 'QA:1', 'cm', 'connection', 'monopodal', '', [0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11], 10, mol, mol_qa, mol_qa, 'QA:0@3[cm]QA:1@0'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_circular():
    """Scenario case with a circular fcg."""
    mol = Chem.MolFromSmiles('C1COCC(C1)C(C1CCNCC1)C1CCCSC1')
    return DataFrame([
                      ['mol_fc_circular', 'XXX', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cm', 'connection', 'monopodal', '', (8, 7, 12, 11, 10, 9), (5, 4, 3, 2, 1, 0), 20, mol, mol_qa, mol_qb, 'QA:0@1[cm]QB:0@1'],
                      ['mol_fc_circular', 'XXX', 'QA', 0, 'QA:0', 'QC', 0, 'QC:0', 'cm', 'connection', 'monopodal', '', (8, 7, 12, 11, 10, 9), (14, 13, 18, 17, 16, 15), 20, mol, mol_qa, mol_qc, 'QA:0@1[cm]QC:0@1'],
                      ['mol_fc_circular', 'XXX', 'QB', 0, 'QA:0', 'QC', 0, 'QC:0', 'cm', 'connection', 'monopodal', '', (5, 4, 3, 2, 1, 0), (14, 13, 18, 17, 16, 15), 20, mol, mol_qb, mol_qc, 'QB:0@1[cm]QC:0@1'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_independant():
    """Scenario case with 2 independant subfcgs."""
    mol = Chem.MolFromSmiles('C(CCC1CSC2CCCCC2C1)CC1COCC(CC2CCCNC2)C1')
    return DataFrame([
                      ['mol_fc_independant', 'XXX', 'QA', 0, 'QA:0', 'QB', 0, 'QB:0', 'cm', 'connection', 'monopodal', '', (20, 21, 22, 23, 24, 25), (26, 18, 17, 16, 15, 14), 26, mol, mol_qa, mol_qb, 'QA:0@0[cm]QB:0@1'],
                      ['mol_fc_independant', 'XXX', 'QC', 0, 'QC:0', 'QD', 0, 'QD:0', 'fe', 'fusion', 'edge', '', (12, 11, 6, 5, 4, 3), (6, 7, 8, 9, 10, 11), 26, mol, mol_qc, mol_qd, 'QC:0@1,2[fe]QD:0@0,5'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_1():
    """Scenario case with 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1CCCCC1O')
    return DataFrame([
                      ['mol_fc_overlap_1', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (7, 6, 5, 4, 3, 2, 1), (0, 1, 2, 3, 4, 5, 6), 8, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_2():
    """Scenario case with 2 fragments overlapping, bound to a common fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CCC1O)C1CCC1')
    return DataFrame([
                      ['mol_fc_overlap_2', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (7, 6, 5, 4, 3, 2, 1), (0, 1, 2, 3, 4, 5, 6), 12, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_2', 'XXX', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'cm', 'connection', 'monopodal', '', (7, 6, 5, 4, 3, 2, 1), (8, 9, 10, 11), 12, mol, mol_o1, mol_o3, 'O1:0@4[cm]O3:0@0'],
                      ['mol_fc_overlap_2', 'XXX', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'cm', 'connection', 'monopodal', '', (0, 1, 2, 3, 4, 5, 6), (8, 9, 10, 11), 12, mol, mol_o2, mol_o3, 'O2:0@3[cm]O3:0@0'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_3():
    """Scenario case with 2 fragments overlapping, bound to a common fragment, itself bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CCC1O)C1CC(CC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_3', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (7, 6, 5, 4, 3, 2, 1), (0, 1, 2, 3, 4, 5, 6), 16, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_3', 'XXX', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (7, 6, 5, 4, 3, 2, 1), (8, 9, 10, 15), 16, mol, mol_o1, mol_o3, 'O1:0@4[cm]O6:0@0'],
                      ['mol_fc_overlap_3', 'XXX', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (0, 1, 2, 3, 4, 5, 6), (8, 9, 10, 15), 16, mol, mol_o2, mol_o3, 'O2:0@3[cm]O6:0@0'],
                      ['mol_fc_overlap_3', 'XXX', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cm', 'connection', 'monopodal', '', (8, 9, 10, 15), (12, 13, 14), 16, mol, mol_o3, mol_o4, 'O6:0@2[cm]O7:0@0'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_4():
    """Scenario case with 2 independant sets of 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1CC(CCCCC2CCC(O)C(N)C2)CC1O')
    return DataFrame([
                      ['mol_fc_overlap_4', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (12, 11, 10, 9, 8, 15, 13), (14, 13, 11, 10, 9, 8, 15), 19, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_4', 'XXX', 'O5', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'overlap', (18, 17, 16, 3, 2, 1), (0, 1, 2, 3, 16, 17), 19, mol, mol_o4, mol_o5, 'O4:0@1,2,3,4,5[ffo]O5:0@1,2,3,4,5'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_5():
    """Scenario case with 2 sets of 2 fragments overlapping."""
    mol = Chem.MolFromSmiles('NC1C(O)CCC1C1CCC(O)C(N)C1')
    return DataFrame([
                      ['mol_fc_overlap_5', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (11, 10, 9, 8, 7, 14, 12), (13, 12, 10, 9, 8, 7, 14), 15, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_5', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (11, 10, 9, 8, 7, 14, 12), (3, 2, 1, 6, 5, 4), 15, mol, mol_o1, mol_o4, 'O1:0@4[cm]O4:0@3'],
                      ['mol_fc_overlap_5', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (11, 10, 9, 8, 7, 14, 12), (0, 1, 2, 4, 5, 6), 15, mol, mol_o1, mol_o5, 'O1:0@4[cm]O5:0@5'],
                      ['mol_fc_overlap_5', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (13, 12, 10, 9, 8, 7, 14), (3, 2, 1, 6, 5, 4), 15, mol, mol_o2, mol_o4, 'O2:0@5[cm]O4:0@3'],
                      ['mol_fc_overlap_5', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (13, 12, 10, 9, 8, 7, 14), (0, 1, 2, 4, 5, 6), 15, mol, mol_o2, mol_o5, 'O2:0@5[cm]O5:0@5'],
                      ['mol_fc_overlap_5', 'XXX', 'O5', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'overlap', (3, 2, 1, 6, 5, 4), (0, 1, 2, 4, 5, 6), 15, mol, mol_o4, mol_o5, 'O4:0@1,2,3,4,5[ffo]O5:0@1,2,3,4,5'],
                      ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_6():
    """Scenario case with 2 sets of 2 overlapping fragments, 1 set being bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1C(O)CCC1C1CC(N)C(O)CC1CCC1CC(CCC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_6', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (12, 11, 9, 8, 7, 14, 13), (10, 9, 8, 7, 14, 13, 11), 26, mol, mol_o1, mol_o2, 'O1:0@4[cm]O4:0@3'],
                      ['mol_fc_overlap_6', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (12, 11, 9, 8, 7, 14, 13), (3, 2, 1, 6, 5, 4), 26, mol, mol_o1, mol_o4, 'O1:0@4[cm]O5:0@5'],
                      ['mol_fc_overlap_6', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (12, 11, 9, 8, 7, 14, 13), (0, 1, 2, 4, 5, 6), 26, mol, mol_o1, mol_o5, 'O1:0@5[cm]O6:0@0'],
                      ['mol_fc_overlap_6', 'XXX', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (12, 11, 9, 8, 7, 14, 13), (17, 18, 19, 25), 26, mol, mol_o1, mol_o6, 'O2:0@3[cm]O4:0@3'],
                      ['mol_fc_overlap_6', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (10, 9, 8, 7, 14, 13, 11), (3, 2, 1, 6, 5, 4), 26, mol, mol_o2, mol_o4, 'O2:0@3[cm]O5:0@5'],
                      ['mol_fc_overlap_6', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (10, 9, 8, 7, 14, 13, 11), (0, 1, 2, 4, 5, 6), 26, mol, mol_o2, mol_o5, 'O2:0@4[cm]O6:0@0'],
                      ['mol_fc_overlap_6', 'XXX', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'false_positive', '', (10, 9, 8, 7, 14, 13, 11), (17, 18, 19, 25), 26, mol, mol_o2, mol_o6, 'O2:0@4[cm]O6:0@0'],
                      ['mol_fc_overlap_6', 'XXX', 'O5', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'false_positive', (3, 2, 1, 6, 5, 4), (0, 1, 2, 4, 5, 6), 26, mol, mol_o4, mol_o5, 'O4:0@1,2,3,4,5[ffo]O5:0@1,2,3,4,5'],
                      ['mol_fc_overlap_6', 'XXX', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cm', 'connection', 'monopodal', '', (17, 18, 19, 25), (22, 23, 24), 26, mol, mol_o6, mol_o7, 'O6:0@2[cm]O7:0@0'],
                  ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_7():
    """Scenario case with 3 overlapping fragments, bound to another fragment."""
    mol = Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CCC1')
    return DataFrame([
                      ['mol_fc_overlap_7', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (8, 7, 5, 4, 3, 2, 1), (0, 1, 2, 3, 4, 5, 7), 13, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_7', 'XXX', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', (8, 7, 5, 4, 3, 2, 1), (6, 5, 4, 3, 2, 1, 7), 13, mol, mol_o1, mol_o3, 'O1:0@1,2,3,4,5,6[ffo]O3:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_7', 'XXX', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (8, 7, 5, 4, 3, 2, 1), (9, 10, 11, 12), 13, mol, mol_o1, mol_o6, 'O1:0@4[cm]O6:0@0'],
                      ['mol_fc_overlap_7', 'XXX', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', (0, 1, 2, 3, 4, 5, 7), (6, 5, 4, 3, 2, 1, 7), 13, mol, mol_o3, mol_o3, 'O2:0@1,2,3,4,5,6[ffo]O3:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_7', 'XXX', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (0, 1, 2, 3, 4, 5, 7), (9, 10, 11, 12), 13, mol, mol_o2, mol_o6, 'O2:0@3[cm]O6:0@0'],
                      ['mol_fc_overlap_7', 'XXX', 'O3', 0, 'O3:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (6, 5, 4, 3, 2, 1, 7), (9, 10, 11, 12), 13, mol, mol_o3, mol_o6, 'O3:0@3[cm]O6:0@0'],
                  ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_8():
    """Scenario case with 3 overlapping fragments, bound to a common combination of 2 fragments."""
    mol = Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CC(CCC2CC2)C1')
    return DataFrame([
                      ['mol_fc_overlap_8', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (8, 7, 5, 4, 3, 2, 1), (0, 1, 2, 3, 4, 5, 7), 18, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_8', 'XXX', 'O1', 0, 'O1:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', (8, 7, 5, 4, 3, 2, 1), (6, 5, 4, 3, 2, 1, 7), 18, mol, mol_o1, mol_o3, 'O1:0@1,2,3,4,5,6[ffo]O3:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_8', 'XXX', 'O1', 0, 'O1:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (8, 7, 5, 4, 3, 2, 1), (9, 10, 11, 17), 18, mol, mol_o1, mol_o6, 'O1:0@4[cm]O6:0@0'],
                      ['mol_fc_overlap_8', 'XXX', 'O2', 0, 'O2:0', 'O3', 0, 'O3:0', 'ffo', 'fusion', 'false_positive', 'overlap', (0, 1, 2, 3, 4, 5, 7), (6, 5, 4, 3, 2, 1, 7), 18, mol, mol_o2, mol_o3, 'O2:0@1,2,3,4,5,6[ffo]O3:0@1,2,3,4,5,6'],
                      ['mol_fc_overlap_8', 'XXX', 'O2', 0, 'O2:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (0, 1, 2, 3, 4, 5, 7), (9, 10, 11, 17), 18, mol, mol_o2, mol_o6, 'O2:0@3[cm]O6:0@0'],
                      ['mol_fc_overlap_8', 'XXX', 'O3', 0, 'O3:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (6, 5, 4, 3, 2, 1, 7), (9, 10, 11, 17), 18, mol, mol_o3, mol_o6, 'O3:0@3[cm]O6:0@0'],
                      ['mol_fc_overlap_8', 'XXX', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cm', 'connection', 'monopodal', '', (9, 10, 11, 17), (14, 15, 16), 18, mol, mol_o6, mol_o7, 'O6:0@2[cm]O7:0@0'],
                  ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


@pytest.fixture
def df_fc_overlap_9():
    """Scenario case with 3 sets of 2 overlapping fragments, bound to a common combination of 2 redundant fragments."""
    mol = Chem.MolFromSmiles('NC1C(O)C(CCCC2CC2CCC2CC2)C1CCC1CC(C(N)C1O)C1CCC(O)C(N)C1')
    return DataFrame([
                    ['mol_fc_overlap_9', 'XXX', 'O1', 0, 'O1:0', 'O2', 0, 'O2:0', 'ffo', 'fusion', 'false_positive', 'overlap', (30, 29, 28, 27, 26, 33, 31), (32, 31, 29, 28, 27, 26, 33), 34, mol, mol_o1, mol_o2, 'O1:0@1,2,3,4,5,6[ffo]O2:0@1,2,3,4,5,6'],
                    ['mol_fc_overlap_9', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (30, 29, 28, 27, 26, 33, 31), (25, 24, 22, 21, 20, 19), 34, mol, mol_o1, mol_o4, 'O1:0@4[cm]O4:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O1', 0, 'O1:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (30, 29, 28, 27, 26, 33, 31), (23, 22, 21, 20, 19, 24), 34, mol, mol_o1, mol_o5, 'O1:0@4[cm]O5:0@2'],
                    ['mol_fc_overlap_9', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O4:0', 'cm', 'connection', 'monopodal', '', (32, 31, 29, 28, 27, 26, 33), (25, 24, 22, 21, 20, 19), 34, mol, mol_o2, mol_o4, 'O2:0@5[cm]O4:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O2', 0, 'O2:0', 'O5', 0, 'O5:0', 'cm', 'connection', 'monopodal', '', (32, 31, 29, 28, 27, 26, 33), (23, 22, 21, 20, 19, 24), 34, mol, mol_o2, mol_o5, 'O2:0@5[cm]O5:0@2'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O4:0', 'O5', 0, 'O5:0', 'ffo', 'fusion', 'false_positive', 'overlap', (25, 24, 22, 21, 20, 19), (23, 22, 21, 20, 19, 24), 34, mol, mol_o4, mol_o5, 'O4:0@1,2,3,4,5[ffo]O5:0@1,2,3,4,5'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O4:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (25, 24, 22, 21, 20, 19), (1, 2, 4, 16), 34, mol, mol_o4, mol_o6, 'O4:0@5[cm]O6:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O4:0', 'O8', 0, 'O8:0', 'cm', 'connection', 'monopodal', '', (25, 24, 22, 21, 20, 19), (0, 1, 2, 4, 16), 34, mol, mol_o4, mol_o8, 'O4:0@5[cm]O8:0@4'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O4:0', 'O9', 0, 'O9:0', 'cm', 'connection', 'monopodal', '', (25, 24, 22, 21, 20, 19), (3, 2, 1, 16, 4), 34, mol, mol_o4, mol_o9, 'O4:0@5[cm]O9:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O5:0', 'O6', 0, 'O6:0', 'cm', 'connection', 'monopodal', '', (23, 22, 21, 20, 19, 24), (1, 2, 4, 16), 34, mol, mol_o5, mol_o6, 'O5:0@4[cm]O6:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O5:0', 'O8', 0, 'O8:0', 'cm', 'connection', 'monopodal', '', (23, 22, 21, 20, 19, 24), (0, 1, 2, 4, 16), 34, mol, mol_o5, mol_o8, 'O5:0@4[cm]O8:0@4'],
                    ['mol_fc_overlap_9', 'XXX', 'O5', 0, 'O5:0', 'O9', 0, 'O9:0', 'cm', 'connection', 'monopodal', '', (23, 22, 21, 20, 19, 24), (3, 2, 1, 16, 4), 34, mol, mol_o5, mol_o9, 'O5:0@4[cm]O9:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O6', 0, 'O6:0', 'O7', 0, 'O7:0', 'cm', 'connection', 'monopodal', '', (1, 2, 4, 16), (8, 9, 10), 34, mol, mol_o6, mol_o7, 'O6:0@2[cm]O7:0@0'],
                    ['mol_fc_overlap_9', 'XXX', 'O6', 0, 'O6:0', 'O8', 0, 'O8:0', 'ffs', 'fusion', 'false_positive', 'substructure', (1, 2, 4, 16), (0, 1, 2, 4, 16), 34, mol, mol_o6, mol_o8, 'O6:0@0,1,2,3[ffs]O8:0@1,2,3,4'],
                    ['mol_fc_overlap_9', 'XXX', 'O6', 0, 'O6:0', 'O9', 0, 'O9:0', 'ffs', 'fusion', 'false_positive', 'substructure', (1, 2, 4, 16), (3, 2, 1, 16, 4), 34, mol, mol_o6, mol_o9, 'O6:0@0,1,2,3[ffs]O9:0@1,2,3,4'],
                    ['mol_fc_overlap_9', 'XXX', 'O7', 0, 'O7:0', 'O7', 1, 'O7:1', 'cm', 'connection', 'monopodal', '', (8, 9, 10), (13, 14, 15), 34, mol, mol_o7, mol_o7, 'O7:0@2[cm]O7:1@0'],
                    ['mol_fc_overlap_9', 'XXX', 'O7', 0, 'O7:0', 'O8', 0, 'O8:0', 'cm', 'connection', 'monopodal', '', (8, 9, 10), (0, 1, 2, 4, 16), 34, mol, mol_o7, mol_o8, 'O7:0@0[cm]O8:0@3'],
                    ['mol_fc_overlap_9', 'XXX', 'O7', 0, 'O7:0', 'O9', 0, 'O9:0', 'cm', 'connection', 'monopodal', '', (8, 9, 10), (3, 2, 1, 16, 4), 34, mol, mol_o7, mol_o9, 'O7:0@0[cm]O9:0@4'],
                    ['mol_fc_overlap_9', 'XXX', 'O8', 0, 'O8:0', 'O9', 0, 'O9:0', 'ffo', 'fusion', 'false_positive', 'overlap', (0, 1, 2, 4, 16), (3, 2, 1, 16, 4), 34, mol, mol_o8, mol_o9, 'O8:0@1,2,3,4[ffo]O9:0@1,2,3,4'],
              ], columns=['idm', 'inchikey', 'idf1', 'idxf1', 'fid1', 'idf2', 'idxf2', 'fid2', 'fcc', 'category', 'type', 'subtype', '_aidxf1', '_aidxf2', 'hac', 'mol', 'mol_frag_1', 'mol_frag_2', 'fc'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_fcg_fc_simple(df_fc_simple):
    """Compute a fragment graph from fragment combinations"""
    # default: 3 <= n <= 5
    df_fcg = fragment_combination_graph.generate(df_fc_simple)
    result = df_fcg.iloc[0]
    assert result['fcg_str'] == 'QA:0@2[cm]QB:0@1-QA:0@5[cm]QC:0@0' and result['nfrags_u'] == 3
    # when n > max (2)
    df_fcg = fragment_combination_graph.generate(df_fc_simple, max_frags=2)
    assert len(df_fcg.index) == 0
    # when n < min (4)
    df_fcg = fragment_combination_graph.generate(df_fc_simple, min_frags=4)
    assert len(df_fcg.index) == 0


def test_fcg_fc_redundant(df_fc_redundant):
    """Compute a fragment graph from fragment combinations with multiple occurrences from a same fragment"""
    df_fcg = fragment_combination_graph.generate(df_fc_redundant)
    result = df_fcg.iloc[0]
    assert result['fcg_str'] == 'QA:0@3[cm]QA:1@0' and result['nfrags'] == 2 and result['nfrags_u'] == 1


def test_fcg_fc_circular(df_fc_circular):
    """Compute a circular fragment graph from fragment combination"""
    df_fcg = fragment_combination_graph.generate(df_fc_circular)
    result = df_fcg.iloc[0]
    assert result['fcg_str'] == 'QA:0@1[cm]QB:0@1-QA:0@1[cm]QC:0@1-QB:0@1[cm]QC:0@1' and result['nfrags'] == 3 and result['nfrags_u'] == 3


def test_fcg_fc_independant(df_fc_independant):
    """Compute 2 fragment graphs from 2 independent fragment combinations """
    df_fcg = fragment_combination_graph.generate(df_fc_independant)
    assert len(df_fcg.index) == 2
    result1 = df_fcg.iloc[0]
    assert result1['fcg_str'] == 'QA:0@0[cm]QB:0@1' and result1['nfrags'] == 2 and result1['nfrags_u'] == 2
    result2 = df_fcg.iloc[1]
    assert result2['fcg_str'] == 'QC:0@1,2[fe]QD:0@0,5' and result2['nfrags'] == 2 and result2['nfrags_u'] == 2


def test_fcg_fc_overlap_1(df_fc_overlap_1):
    """Split 2 overlapping fragments that result in no fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_1)
    assert len(df_fcg.index) == 0


def test_fcg_fc_overlap_2(df_fc_overlap_2):
    """Split 2 overlapping fragments that result in 2 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_2)
    assert len(df_fcg.index) == 2
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O1:0@4[cm]O3:0@0',
                                                      'O2:0@3[cm]O3:0@0'
                                                      ]


def test_fcg_fc_overlap_3(df_fc_overlap_3):
    """Split 2 overlapping fragments that result in longer 2 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_3)
    assert len(df_fcg.index) == 2
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O6:0@2[cm]O7:0@0-O1:0@4[cm]O6:0@0',
                                                      'O6:0@2[cm]O7:0@0-O2:0@3[cm]O6:0@0'
                                                      ]


def test_fcg_fc_overlap_4(df_fc_overlap_4):
    """Split 2x2 overlapping fragments that result in 0 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_4)
    assert len(df_fcg.index) == 0


def test_fcg_fc_overlap_5(df_fc_overlap_5):
    """Split 2x2 overlapping fragments that result in 4 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_5)
    assert len(df_fcg.index) == 4
    # test each of the fgraphs
    # for i in range(len(df_fcg.index)):
    #     print(df_fcg.iloc[i]['fcg_str'])
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O1:0@4[cm]O4:0@3',
                                                      'O1:0@4[cm]O5:0@5',
                                                      'O2:0@5[cm]O4:0@3',
                                                      'O2:0@5[cm]O5:0@5',
                                                      ]


def test_fcg_fc_overlap_6(df_fc_overlap_6):
    """Split 2x2 overlapping fragments that result in 4 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_6)
    assert len(df_fcg.index) == 4
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O6:0@2[cm]O7:0@0-O1:0@4[cm]O5:0@5-O2:0@3[cm]O4:0@3',
                                                      'O6:0@2[cm]O7:0@0-O1:0@5[cm]O6:0@0-O2:0@3[cm]O4:0@3',
                                                      'O6:0@2[cm]O7:0@0-O2:0@3[cm]O5:0@5-O2:0@4[cm]O6:0@0',
                                                      'O6:0@2[cm]O7:0@0-O2:0@4[cm]O6:0@0-O2:0@4[cm]O6:0@0',
                                                      ]


def test_fcg_fc_overlap_7(df_fc_overlap_7):
    """Split 3x2 overlapping fragments that result in 3 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_7)
    assert len(df_fcg.index) == 3
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O1:0@4[cm]O6:0@0',
                                                      'O2:0@3[cm]O6:0@0',
                                                      'O3:0@3[cm]O6:0@0',
                                                      ]


def test_fcg_fc_overlap_8(df_fc_overlap_8):
    """Split 3x2 overlapping fragments that result in 3 longer fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_8)
    assert len(df_fcg.index) == 3
    assert sorted(list(df_fcg['fcg_str'].values)) == ['O6:0@2[cm]O7:0@0-O1:0@4[cm]O6:0@0',
                                                      'O6:0@2[cm]O7:0@0-O2:0@3[cm]O6:0@0',
                                                      'O6:0@2[cm]O7:0@0-O3:0@3[cm]O6:0@0',
                                                      ]


def test_fcg_fc_overlap_9(df_fc_overlap_9):
    """Split 2x2x2 overlapping fragments that result in 8 fragment graphs"""
    df_fcg = fragment_combination_graph.generate(df_fc_overlap_9)
    assert len(df_fcg.index) == 8

    assert sorted(list(df_fcg['fcg_str'].values)) == ['O7:0@2[cm]O7:1@0-O1:0@4[cm]O4:0@3-O4:0@5[cm]O8:0@4-O7:0@0[cm]O8:0@3',
                                                      'O7:0@2[cm]O7:1@0-O1:0@4[cm]O4:0@3-O4:0@5[cm]O9:0@3-O7:0@0[cm]O9:0@4',
                                                      'O7:0@2[cm]O7:1@0-O1:0@4[cm]O5:0@2-O5:0@4[cm]O8:0@4-O7:0@0[cm]O8:0@3',
                                                      'O7:0@2[cm]O7:1@0-O1:0@4[cm]O5:0@2-O5:0@4[cm]O9:0@3-O7:0@0[cm]O9:0@4',
                                                      'O7:0@2[cm]O7:1@0-O2:0@5[cm]O4:0@3-O4:0@5[cm]O8:0@4-O7:0@0[cm]O8:0@3',
                                                      'O7:0@2[cm]O7:1@0-O2:0@5[cm]O4:0@3-O4:0@5[cm]O9:0@3-O7:0@0[cm]O9:0@4',
                                                      'O7:0@2[cm]O7:1@0-O2:0@5[cm]O5:0@2-O5:0@4[cm]O8:0@4-O7:0@0[cm]O8:0@3',
                                                      'O7:0@2[cm]O7:1@0-O2:0@5[cm]O5:0@2-O5:0@4[cm]O9:0@3-O7:0@0[cm]O9:0@4',
                                                      ]
