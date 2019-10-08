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
from npfc import fragment_combination
from npfc import utils
# debug
import logging
logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



@pytest.fixture
def df_frags():
    """Two fragments used for most of substructue searches."""
    df_frags = pd.DataFrame({'smiles': ['C1NCCC1', 'C1CCCOC1']}, index=['QA', 'QB'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

#
# @pytest.fixture
# def df_case_repeated_frag():
#     """A molecule with a fragment repeated twice to test how colormap is computed."""
#     df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('OCCC(C1CCCC1)C1CCC2CCCCC2C1')]}, index=['REPEATEDFRAGS'])
#     df_frags = pd.DataFrame({'smiles': ['C1CCCCC1', 'C1CCCC1']}, index=['A', 'B'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return (df_mols, df_frags)
#
#
# @pytest.fixture
# def df_case_chembl_1():
#     """A molecule from real-life study with false positive combinations such as ffo and ffs for testing clean and map functions."""
#     df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('CS(=O)(=O)c1ccc2nc(O)c(c(O)c2c1)c3ccccc3')]}, index=['CHEMBL209576'])
#     df_frags = pd.DataFrame({'smiles': ['c1ccc2ncccc2c1', 'c1ccccc1', 'Oc1ccccn1', 'c1ccncc1']}, index=['2', '32', '320', '328'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return (df_mols, df_frags)
#
#
# @pytest.fixture
# def df_case_chembl_2():
#     """A molecule with wrong aidxfs in the final output."""
#     df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('COC1=CC=C2CC3C4C[C@]5(CCCCC6=CC=CC=C6)CO[C@@H]5C5OC1=C2C45CCN3C')]}, index=['CHEMBL10006'])
#     df_frags = pd.DataFrame({'smiles': ['C1OC2CCCCC12', 'C1OC2=C3C(CC4CC13CCN4)=CC=C2']}, index=['678', '1141'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return (df_mols, df_frags)
#
#
# @pytest.fixture
# def df_case_chembl_3():
#     """A molecule with 4 alternative fragment maps that should be only one."""
#     df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('COc1ccc2[nH]n(CC3(C)NC(=O)NC3=O)c(=O)c2c1')]}, index=['CHEMBL1078396'])
#     df_frags = pd.DataFrame({'smiles': ['O=c1cc[nH][nH]1', 'O=C1CNCN1', 'O=C1NCCN1', 'c1ccc2[nH]ncc2c1']}, index=['77', '358', '676', '1414'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return (df_mols, df_frags)
#
# @pytest.fixture
# def df_case_dnp_1():
#     """A molecule with mislabeld as cbr instead of ffo."""
#     df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('OCC(O)C1=CC(O)=C(O)N1')]}, index=['00242908-  -001'])
#     df_frags = pd.DataFrame({'smiles': ['OC1=CC=CN1', 'CC1=CC=CN1']}, index=['1567', '1967'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return (df_mols, df_frags)
#
# @pytest.fixture
# def df_overlap_frags():
#     """The set of fragments to use for overlap cases."""
#     df_frags = pd.DataFrame({'smiles': ['SC1CCCCC1',
#                                         'NC1CCCCC1',
#                                         'OC1CCCCC1',
#                                         'OC1CCCC1',
#                                         'NC1CCCC1',
#                                         'C1CCC1',
#                                         'C1CC1',
#                                         'NC1CCC1',
#                                         'OC1CCC1',
#                                         ]}, index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])
#     df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
#     return df_frags
#
# @pytest.fixture
# def df_overlap_case_1():
#     """A single combination of two overlapping fragments."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CCCCC1O')]}, index=['MOV1'])
#
# @pytest.fixture
# def df_overlap_case_2():
#     """Two overlapping fragments giving 2 alternative fragment maps with no common combination."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CC(CCC1O)C1CCC1')]}, index=['MOV2'])
#
# @pytest.fixture
# def df_overlap_case_3():
#     """Two overlapping fragments giving 2 alternative fragment maps with a common combination."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CC(CCC1O)C1CC(CC2CC2)C1')]}, index=['MOV3'])
#
# @pytest.fixture
# def df_overlap_case_4():
#     """Two independant sets of two overlapping fragments giving 2 alternative fragment maps with no common combination."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CC(CCCC2CCC(O)C(N)C2)CC1O')]}, index=['MOV4'])
#
# @pytest.fixture
# def df_overlap_case_5():
#     """Two sets of two overlapping fragments giving 2 alternative fragment maps with no common combination."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1C(O)CCC1C1CCC(O)C(N)C1')]}, index=['MOV5'])
#
# @pytest.fixture
# def df_overlap_case_6():
#     """Two sets of two overlapping fragments giving 2 alternative fragment maps with one common combination."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1C(O)CCC1C1CC(N)C(O)CC1CCC1CC(CCC2CC2)C1')]}, index=['MOV6'])
#
# @pytest.fixture
# def df_overlap_case_7():
#     """Three overlapping fragments combined with a same other fragment."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CCC1')]}, index=['MOV7'])
#
# @pytest.fixture
# def df_overlap_case_8():
#     """Three overlapping fragments combined with a same other fragment linked to another fragment."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1CC(CC(S)C1O)C1CC(CCC2CC2)C1')]}, index=['MOV8'])
#
# @pytest.fixture
# def df_overlap_case_9():
#     """Three sets of overlapping fragments with a common fragment at the end."""
#     return pd.DataFrame({'mol': [Chem.MolFromSmiles('NC1C(O)C(C2CC2CCC2CC2)C1C1CC(C(N)C1O)C1CCC(O)C(N)C1')]}, index=['MOV9'])

@pytest.fixture
def df_mol_fusion_spiro():
    """Example molecule with the fusion spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2(CN1)CCCOC2')]}, index=['fsp'])


@pytest.fixture
def df_mol_fusion_edge():
    """Example molecule with the fusion edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2CCOCC2N1')]}, index=['fed'])


@pytest.fixture
def df_mol_fusion_bridged():
    """Example molecule with the fusion bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CC1COC2')]}, index=['fbr'])


@pytest.fixture
def df_mol_fusion_other():
    """Example molecule with the fusion other fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2COCC1CCCCCC2')]}, index=['fot'])


@pytest.fixture
def df_mol_connection_monopodal():
    """Example molecule with the connection monopodal fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC(CN1)C1CCCOC1')]}, index=['cmo'])


@pytest.fixture
def df_mol_connection_bipodal_spiro():
    """Example molecule with the connection bipodal spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2(CCCOC2)C2CNCC12')]}, index=['cbs'])


@pytest.fixture
def df_mol_connection_bipodal_edge():
    """Example molecule with the connection bipodal edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NCC2C1CCC1OCCCC21')]}, index=['cbe'])


@pytest.fixture
def df_mol_connection_bipodal_bridged():
    """Example molecule with the connection bipodal bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CC1C1CCCOC21')]}, index=['cbb'])


@pytest.fixture
def df_mol_connection_bipodal_linker():
    """Example molecule with the connection bipodal linker fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1COC2C3CCCCCCC(CN3)C2C1')]}, index=['cbl'])



@pytest.fixture
def df_mol_connection_tripodal_spiro():
    """Example molecule with the connection tripodal spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1C2CNC3CCC4CCOCC14C23')]}, index=['cts'])


@pytest.fixture
def df_mol_connection_tripodal_edge():
    """Example molecule with the connection tripodal edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CCC3CCOC4CCC1C2C34')]}, index=['cte'])


@pytest.fixture
def df_mol_connection_tripodal_bridged():
    """Example molecule with the connection tripodal bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1C2NC3CCC4CCOC2C4C13')]}, index=['ctb'])


@pytest.fixture
def df_mol_connection_other_spiro():
    """Example molecule with the connection other spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC23CC4COC2CCC2NC(C4)C1C32')]}, index=['cos'])


@pytest.fixture
def df_mol_connection_other_edge():
    """Example molecule with the connection other edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2OCC3CC4NC1C1C4CCC3C21')]}, index=['coe'])


@pytest.fixture
def df_mol_connection_other_bridged():
    """Example molecule with the connection other brdiged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2C3COC4CCC5OC(CC1C5C24)C3')]}, index=['cob'])


@pytest.fixture
def df_mol_fusion_false_positive_substructure():
    """Example molecule with the connection false_positive substructure fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('CCC1CCNC1')]}, index=['ffs'])


@pytest.fixture
def df_frags_fusion_false_positive_substructure():
    """Two overlapping fragments used to identify false_positive substructure fragment combinations."""
    df_frags = pd.DataFrame({'smiles': ['C1NCCC1', 'CC1CCNC1']}, index=['QA', 'QF'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

@pytest.fixture
def df_frags_fusion_linker():
    """Two overlapping fragments used to identify fusion linkers."""
    df_frags = pd.DataFrame({'smiles': ['C1CCCCNCCCC1', 'C1CCCCCOCCCC1']}, index=['QA', 'QF'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

@pytest.fixture
def df_frags_connection_bipodal_linker():
    """Two overlapping fragments used to identify fusion linkers."""
    df_frags = pd.DataFrame({'smiles': ['C1CCCCNCCCC1', 'C1CCOCC1']}, index=['QA', 'QF'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

@pytest.fixture
def df_mol_connection_tripodal_linker():
    """Example molecule with the connection false_positive substructure fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2CCC3NC4CCCCCCC3C2C4O1')]}, index=['ctl'])

@pytest.fixture
def df_frags_connection_tripodal_linker():
    """Two overlapping fragments used to identify fusion linkers."""
    df_frags = pd.DataFrame({'smiles': ['C1CCCCNCCCC1', 'C1CCOCC1']}, index=['QA', 'QF'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

@pytest.fixture
def df_mol_connection_other_linker():
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2C3COC4CCC5OC(CC1C5C24)CCCCCC3')]}, index=['col'])

@pytest.fixture
def df_frags_connection_other_linker():
    df_frags = pd.DataFrame({'smiles': ['C1CCOCC1']}, index=['QB'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags

# def test_fcc_connection_other_linker(df_mol_connection_other_linker, df_frags_connection_other_linker):
#     """Check if connection bipodal bridged fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'
#     logging.debug("")
#     df_aidxf = fragment.find(df_mol_connection_other_linker, df_frags_connection_other_linker)
#     print(f"\n{df_aidxf}\n")
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

def test_fcc_fusion_false_positive_overlap():
    # df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
    # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    # result = df_fcc.iloc[0]
    # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'

    df_mol_fusion_false_positive_overlap = pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NCC2CC3COCCC3CC12')]}, index=['ffo'])
    df_frags_fusion_false_positive_overlap = pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NCC2CCCCC12'), Chem.MolFromSmiles('C1CCC2COCCC2C1')]}, index=['QA', 'QC'])
    logging.debug("")
    df_aidxf = fragment.find(df_mol_fusion_false_positive_overlap, df_frags_fusion_false_positive_overlap)
    print(f"\n{df_aidxf}\n")
    df_fcc = fragment_combination.classify_df(df_aidxf)
    df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
    df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
    # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
    print(f"\n{df_aidxf}\n")
    print(f"\n{df_fcc}\n")

@pytest.fixture
def df_mol_connection_false_positive_cutoff():
    """Example molecule with the connection false_positive cutoff fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(CCC1CCCOC1)CC1CCNC1')]}, index=['cfc'])






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# def test_fcc_run_substructure_match(fcc, fm, df_mol_fusion_spiro, df_frags):
#     """Check that the df_aidxf obtained by substructure match is adequate."""
#     df_aidxf = fm.run(df_mol_fusion_spiro, df_frags)
#     assert list(df_aidxf.index) == [0, 1]
#     assert list(df_aidxf['idm']) == ["fsp"] * 2
#     assert list(df_aidxf['idf']) == ['QA', 'QB']
#     assert list(df_aidxf['_aidxf']) == [{0, 1, 2, 3, 4}, {2, 5, 6, 7, 8, 9}]
#     assert list(df_aidxf['mol_perc']) == [50.0, 60.0]
#     assert list(df_aidxf['mol'].map(Chem.MolToSmiles)) == ['C1COCC2(C1)CCNC2'] * 2
#

# def test_fcc_fusion_spiro(df_mol_fusion_spiro, df_frags):
#     """Check if fusion spiro fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_fusion_spiro, df_frags)
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
    # df_fcc = fragment_combination.classify_df(df_aidxf)
    # df_fcc['mol'] = df_fcc['mol'].map(Chem.MolToSmiles)
    # print(f"\n{df_fcc}\n")
    # result = df_fcc.iloc[0]
    # assert result['category'] == 'fusion' and result['type'] == 'spiro' and result['subtype'] == '' and result['abbrev'] == 'fsp'


# def test_fcc_fusion_edge(df_mol_fusion_edge, df_frags):
#     """Check if fusion edge fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_fusion_edge, df_frags)
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
# #     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
# #     result = df_fcc.iloc[0]
# #     assert result['category'] == 'fusion' and result['type'] == 'edge' and result['subtype'] == '' and result['abbrev'] == 'fed'
# #
#
# def test_fcc_fusion_bridged(df_mol_fusion_bridged, df_frags):
#     """Check if fusion bridged fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_fusion_bridged, df_frags)
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
    # result = df_fcc.iloc[0]
    # assert result['category'] == 'fusion' and result['type'] == 'bridged' and result['subtype'] == '' and result['abbrev'] == 'fbr'

#
# @pytest.mark.skip  # no example for fot yet!
# def test_fcc_fusion_other(df_mol_fusion_other, df_frags_fusion_linker):
#     """Check if fusion other fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_fusion_other, df_frags_fusion_linker)
#
#
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'fusion' and result['type'] == 'other' and result['subtype'] == '' and result['abbrev'] == 'fot'
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")


#
#
# def test_fcc_connection_monopodal(df_mol_connection_monopodal, df_frags):
#     """Check if connection monopodal fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_monopodal, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cmo'
#     df_aidxf = fragment.find(df_mol_connection_monopodal, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_fcc_connection_bipodal_spiro(df_mol_connection_bipodal_spiro, df_frags):
#     """Check if connection bipodal spiro fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_connection_bipodal_spiro, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

# def test_fcc_connection_bipodal_edge(df_mol_connection_bipodal_edge, df_frags):
#     """Check if connection bipodal edge fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_bipodal_edge, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cbe'
#     df_aidxf = fragment.find(df_mol_connection_bipodal_edge, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_fcc_connection_bipodal_bridged(df_mol_connection_bipodal_bridged, df_frags):
#     """Check if connection bipodal bridged fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'
#     df_aidxf = fragment.find(df_mol_connection_bipodal_bridged, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

# def test_fcc_connection_bipodal_linker(df_mol_connection_bipodal_linker, df_frags_connection_bipodal_linker):
#     """Check if connection bipodal bridged fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'
#     df_aidxf = fragment.find(df_mol_connection_bipodal_linker, df_frags_connection_bipodal_linker)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_fcc_connection_tripodal_linker(df_mol_connection_tripodal_linker, df_frags_connection_tripodal_linker):
#     """Check if connection bipodal bridged fragment combinations are identified."""
#     logging.debug("")
#     # df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'
#     df_aidxf = fragment.find(df_mol_connection_tripodal_linker, df_frags_connection_tripodal_linker)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")



# def test_fcc_connection_tripodal_spiro(df_mol_connection_tripodal_spiro, df_frags):
#     """Check if connection tripodal spiro fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_tripodal_spiro, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cts'
#     df_aidxf = fragment.find(df_mol_connection_tripodal_spiro, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
# #
# def test_fcc_connection_tripodal_edge(df_mol_connection_tripodal_edge, df_frags):
#     """Check if connection tripodal edge fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_tripodal_edge, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cte'
#     # df_aidxf = fm.run(df_mol_connection_tripodal_spiro, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cts'
#     df_aidxf = fragment.find(df_mol_connection_tripodal_edge, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_fcc_connection_tripodal_bridged(df_mol_connection_tripodal_bridged, df_frags):
#     """Check if connection tripodal bridged fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_tripodal_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'ctb'
#     df_aidxf = fragment.find(df_mol_connection_tripodal_bridged, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_fcc_connection_other_spiro(df_mol_connection_other_spiro, df_frags):
#     """Check if connection other spiro fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_other_spiro, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'spiro' and result['abbrev'] == 'cos'
#     df_aidxf = fragment.find(df_mol_connection_other_spiro, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# #
# # @pytest.mark.skip  # classified as bridged now...
# def test_fcc_connection_other_edge(df_mol_connection_other_edge, df_frags):
#     """Check if connection other edge fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_other_edge, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'edge' and result['abbrev'] == 'coe'
#     df_aidxf = fragment.find(df_mol_connection_other_edge, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

# def test_fcc_connection_other_bridged(df_mol_connection_other_bridged, df_frags):
#     """Check if connection other bridged fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_other_bridged, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'bridged' and result['abbrev'] == 'cob'
#     df_aidxf = fragment.find(df_mol_connection_other_bridged, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf)
#     result = df_fcc.iloc[0]
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     # assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

# def test_fcc_fusion_false_positive_substructure(df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure):
#     """Check if fusion false_positive substructure fragment combinations are identified."""
#     df_aidxf = fragment.find(df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure)
#     df_fcc = fragment_combination.classify_df(df_aidxf, clean_combinations=False)
#     result = df_fcc.iloc[0]
#     # assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'substructure' and result['abbrev'] == 'ffs'
#     # assert len(fcc.clean(df_fcc).index) == 0
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")

# def test_fcc_connection_false_positive_cutoff(df_mol_connection_false_positive_cutoff, df_frags):
#     """Check if fusion false_positive cutoff fragment combinations are identified."""
#     # df_aidxf = fm.run(df_mol_connection_false_positive_cutoff, df_frags)
#     # df_fcc = fcc.classify_fragment_combinations(df_aidxf, clean=False)
#     # result = df_fcc.iloc[0]
#     # assert result['category'] == 'connection' and result['type'] == 'false_positive' and result['subtype'] == 'cutoff' and result['abbrev'] == 'cfc'
#     # assert len(fcc.clean(df_fcc).index) == 0
#     print(f"df_frags\n{df_frags}\n")
#     df_aidxf = fragment.find(df_mol_connection_false_positive_cutoff, df_frags)
#     df_fcc = fragment_combination.classify_df(df_aidxf, clean_combinations=False)
#     result = df_fcc.iloc[0]
#     # assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'substructure' and result['abbrev'] == 'ffs'
#     # assert len(fcc.clean(df_fcc).index) == 0
#     df_aidxf['mol'] = df_aidxf['mol'].map(Chem.MolToSmiles)
#     df_aidxf['mol_frag'] = df_aidxf['mol_frag'].map(Chem.MolToSmiles)
#     print(f"\n{df_aidxf}\n")
#     print(f"\n{df_fcc}\n")
#
# def test_case_repeated_frag(fcc, fm, df_case_repeated_frag):
#     """Test to see how repeated fragments are represented only once."""
#     df_aidxf = fm.run(df_case_repeated_frag[0], df_case_repeated_frag[1])
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for df_case_repeated_frag:\n{df_fcc.drop('mol', axis=1)}\n")
#     assert list(df_fcc['abbrev'].values) == ['fed', 'cmo', 'cmo']
#     df_map = fcc.map_frags(df_fcc)
#     assert len(df_map.index) == 1
#     logging.debug(f"\nFragment map for df_case_repeated_frag:\n{df_map.drop('mol', axis=1)}\n")
#     g = df_map.iloc[0]["_fmap"]
#     assert list(g.edges(data=True)) == [('A', 'A', {'abbrev': 'fed', 'n_abbrev': 1, 'idm': 'REPEATEDFRAGS'}), ('A', 'B', {'abbrev': 'cmo', 'n_abbrev': 2, 'idm': 'REPEATEDFRAGS'})]
#
#
# def test_case_chembl_1(fcc, fm, df_case_chembl_1):
#     """Test to see how ffo and ffs combinations are dealt with."""
#     df_aidxf = fm.run(df_case_chembl_1[0], df_case_chembl_1[1])
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf, clean=False)
#     logging.debug(f"\nRaw results for chembl_1:\n{df_fcc}\n")
#     assert list(df_fcc['abbrev'].values) == ['ffs', 'cmo', 'ffo', 'ffs', 'cmo',
#                                              'fed', 'fed', 'cmo', 'cmo', 'ffs']
#     df_fcc = fcc.clean(df_fcc)
#     logging.debug(f"\nClean results for chembl_1:\n{df_fcc}\n")
#     assert list(df_fcc['abbrev'].values) == ['cmo', 'ffo', 'cmo']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for chembl_1:\n{df_map}\n")
#     assert list(df_map['fmap_str'] == ["2:0[cmo]32:1", "32:1[cmo]320:3"])
#
#
# def test_case_chembl_2(fcc, fm, df_case_chembl_2):
#     """Test to see if aidxfs are processed correctly."""
#     df_aidxf = fm.run(df_case_chembl_2[0], df_case_chembl_2[1])
#     logging.debug(f"\nSubstructure hits for chembl_2:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for chembl_2:\n{df_fcc}\n")
#     assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for chembl_1:\n{df_map}\n")
#     assert list(df_map['fmap_str'] == ["678:1[fbr]1141:0"])


# def test_case_chembl_3(fcc, fm, df_case_chembl_3):
#     """Test to see if circular maps are identified correctly."""
#     df_aidxf = fm.run(df_case_chembl_3[0], df_case_chembl_3[1])
#     logging.debug(f"\nSubstructure hits for chembl_3:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for chembl_3:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for chembl_3:\n{df_map}\n")
#     # assert list(df_map['fmap_str'] == ["678:1[fbr]1141:0"])


# def test_case_overlap_1(fcc, fm, df_overlap_case_1, df_overlap_frags):
#     """Two overlapping fragments should not yield any fragment map."""
#     df_aidxf = fm.run(df_overlap_case_1, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_1:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_1:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_1:\n{df_map}\n")
#     assert len(df_map.index) == 0


# def test_case_overlap_2(fcc, fm, df_overlap_case_2, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_2, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_2:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_2:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_2:\n{df_map}\n")

# def test_case_overlap_3(fcc, fm, df_overlap_case_3, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_3, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_3:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_3:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_3:\n{df_map}\n")
#
# def test_case_overlap_4(fcc, fm, df_overlap_case_4, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_4, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_4:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_4:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_4:\n{df_map}\n")


# def test_case_overlap_5(fcc, fm, df_overlap_case_5, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_5, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_5:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_5:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_5:\n{df_map}\n")

# def test_case_overlap_6(fcc, fm, df_overlap_case_6, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_6, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_6:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_6:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_6:\n{df_map}\n")

#
# def test_case_overlap_7(fcc, fm, df_overlap_case_7, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_7, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_7:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_7:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_7:\n{df_map}\n")


# def test_case_overlap_8(fcc, fm, df_overlap_case_8, df_overlap_frags):
#     """Two overlapping fragments with a common other fragment should return two alternative fragment maps."""
#     df_aidxf = fm.run(df_overlap_case_8, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_8:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_8:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for test_case_overlap_8:\n{df_map}\n")


# def test_case_overlap_9(fcc, fm, df_overlap_case_9, df_overlap_frags):
#     """Three sets of overlapping fragments to confirm that all alternatives are considered at once."""
#     df_aidxf = fm.run(df_overlap_case_9, df_overlap_frags)
#     logging.debug(f"\nSubstructure hits for test_case_overlap_9:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for test_case_overlap_9:\n{df_fcc}\n")
    # assert list(df_fcc['abbrev'].values) == ['fbr']
    # df_map = fcc.map_frags(df_fcc)
    # logging.debug(f"\nFragment map for test_case_overlap_9:\n{df_map}\n")


#
# def test_case_dnp_1(fcc, fm, df_case_dnp_1):
#     """Test if some ffo are not mislabeled as cbr."""
#     df_aidxf = fm.run(df_case_dnp_1[0], df_case_dnp_1[1])
#     logging.debug(f"\nSubstructure hits for chembl_2:\n{df_aidxf}\n")
#     df_fcc = fcc.classify_fragment_combinations(df_aidxf)
#     logging.debug(f"\nClean results for chembl_2:\n{df_fcc}\n")
#     # assert list(df_fcc['abbrev'].values) == ['fbr']
#     df_map = fcc.map_frags(df_fcc)
#     logging.debug(f"\nFragment map for chembl_1:\n{df_map}\n")
#     # This case has no remaining fm after clearing ffo
#     assert len(df_map.index) == 0
