"""
tests.test_fcc
~~~~~~~~~~
Tests for the npfc.fcc module.
"""
# standard library
import logging
# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import fragment
import logging
# debug
# logging.basicConfig(level=logging.DEBUG)  # uncomment to debug


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def fcc():
    """Default instance of fcc."""
    return fragment.CombinationClassifier()


@pytest.fixture
def fm():
    """Default instance of fm."""
    return fragment.Matcher()


@pytest.fixture
def df_frags():
    """Two fragments used for most of substructue searches."""
    df_frags = pd.DataFrame({'smiles': ['C1NCCC1', 'C1CCCOC1']}, index=['QA', 'QB'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags


@pytest.fixture
def df_case_chembl_1():
    """A molecule from real-life study with false positive combinations such as ffo and ffs for testing clean and map functions."""
    df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('CS(=O)(=O)c1ccc2nc(O)c(c(O)c2c1)c3ccccc3')]}, index=['CHEMBL209576'])
    df_frags = pd.DataFrame({'smiles': ['c1ccc2ncccc2c1', 'c1ccccc1', 'Oc1ccccn1', 'c1ccncc1']}, index=['2', '32', '320', '328'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return (df_mols, df_frags)


@pytest.fixture
def df_case_chembl_2():
    """A molecule with wrong aidxfs in the final output."""
    df_mols = pd.DataFrame({'mol': [Chem.MolFromSmiles('COC1=CC=C2CC3C4C[C@]5(CCCCC6=CC=CC=C6)CO[C@@H]5C5OC1=C2C45CCN3C')]}, index=['CHEMBL10006'])
    df_frags = pd.DataFrame({'smiles': ['C1OC2CCCCC12', 'C1OC2=C3C(CC4CC13CCN4)=CC=C2']}, index=['678', '1141'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return (df_mols, df_frags)


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
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2NC1CO2')]}, index=['fot'])


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
def df_mol_connection_false_positive_cutoff():
    """Example molecule with the connection false_positive cutoff fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(CCC1CCCOC1)CC1CCNC1')]}, index=['cfc'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_fcc_run_substructure_match(fcc, fm, df_mol_fusion_spiro, df_frags):
    """Check that the df_aidxf obtained by substructure match is adequate."""
    df_aidxf = fm.run(df_mol_fusion_spiro, df_frags)
    assert list(df_aidxf.index) == [0, 1]
    assert list(df_aidxf['idm']) == ["fsp"] * 2
    assert list(df_aidxf['idf']) == ['QA', 'QB']
    assert list(df_aidxf['aidxf']) == [{0, 1, 2, 3, 4}, {2, 5, 6, 7, 8, 9}]
    assert list(df_aidxf['mol_perc']) == [50.0, 60.0]
    assert list(df_aidxf['mol'].map(Chem.MolToSmiles)) == ['C1COCC2(C1)CCNC2'] * 2


def test_fcc_fusion_spiro(fcc, fm, df_mol_fusion_spiro, df_frags):
    """Check if fusion spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'spiro' and result['subtype'] == '' and result['abbrev'] == 'fsp'


def test_fcc_fusion_edge(fcc, fm, df_mol_fusion_edge, df_frags):
    """Check if fusion edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'edge' and result['subtype'] == '' and result['abbrev'] == 'fed'


def test_fcc_fusion_bridged(fcc, fm, df_mol_fusion_bridged, df_frags):
    """Check if fusion bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'bridged' and result['subtype'] == '' and result['abbrev'] == 'fbr'


@pytest.mark.skip  # no example for fot yet!
def test_fcc_fusion_other(fcc, fm, df_mol_fusion_other, df_frags):
    """Check if fusion other fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_other, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'other' and result['subtype'] == '' and result['abbrev'] == 'fot'


def test_fcc_connection_monopodal(fcc, fm, df_mol_connection_monopodal, df_frags):
    """Check if connection monopodal fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_monopodal, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cmo'


def test_fcc_connection_bipodal_spiro(fcc, fm, df_mol_connection_bipodal_spiro, df_frags):
    """Check if connection bipodal spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'


def test_fcc_connection_bipodal_edge(fcc, fm, df_mol_connection_bipodal_edge, df_frags):
    """Check if connection bipodal edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cbe'


def test_fcc_connection_bipodal_bridged(fcc, fm, df_mol_connection_bipodal_bridged, df_frags):
    """Check if connection bipodal bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'


def test_fcc_connection_tripodal_spiro(fcc, fm, df_mol_connection_tripodal_spiro, df_frags):
    """Check if connection tripodal spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cts'


def test_fcc_connection_tripodal_edge(fcc, fm, df_mol_connection_tripodal_edge, df_frags):
    """Check if connection tripodal edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cte'


def test_fcc_connection_tripodal_bridged(fcc, fm, df_mol_connection_tripodal_bridged, df_frags):
    """Check if connection tripodal bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'ctb'


def test_fcc_connection_other_spiro(fcc, fm, df_mol_connection_other_spiro, df_frags):
    """Check if connection other spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'spiro' and result['abbrev'] == 'cos'


@pytest.mark.skip  # classified as bridged now...
def test_fcc_connection_other_edge(fcc, fm, df_mol_connection_other_edge, df_frags):
    """Check if connection other edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'edge' and result['abbrev'] == 'coe'


def test_fcc_connection_other_bridged(fcc, fm, df_mol_connection_other_bridged, df_frags):
    """Check if connection other bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'bridged' and result['abbrev'] == 'cob'


def test_fcc_fusion_false_positive_substructure(fcc, fm, df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure):
    """Check if fusion false_positive substructure fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf, clean=False)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'substructure' and result['abbrev'] == 'ffs'
    assert len(fcc.clean(df_fcc).index) == 0


def test_fcc_connection_false_positive_cutoff(fcc, fm, df_mol_connection_false_positive_cutoff, df_frags):
    """Check if fusion false_positive cutoff fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_false_positive_cutoff, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_aidxf, clean=False)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'false_positive' and result['subtype'] == 'cutoff' and result['abbrev'] == 'cfc'
    assert len(fcc.clean(df_fcc).index) == 0


def test_case_chembl_1(fcc, fm, df_case_chembl_1):
    """Test to see how ffo and ffs combinations are deal with."""
    df_aidxf = fm.run(df_case_chembl_1[0], df_case_chembl_1[1])
    df_fcc = fcc.classify_fragment_combinations(df_aidxf, clean=False)
    logging.debug(f"\nRaw results for chembl_1:\n{df_fcc}\n")
    assert list(df_fcc['abbrev'].values) == ['ffs', 'cmo', 'ffo', 'ffs', 'cmo',
                                             'fed', 'fed', 'cmo', 'cmo', 'ffs']
    df_fcc = fcc.clean(df_fcc)
    logging.debug(f"\nClean results for chembl_1:\n{df_fcc}\n")
    assert list(df_fcc['abbrev'].values) == ['cmo', 'ffo', 'cmo']
    df_map = fcc.map_frags(df_fcc)
    logging.debug(f"\nFragment map for chembl_1:\n{df_map}\n")
    assert list(df_map['map_str'] == ["2:0[cmo]32:1", "32:1[cmo]320:3"])


def test_case_chembl_2(fcc, fm, df_case_chembl_2):
    """Test to see if aidxfs are processed correctly."""
    df_aidxf = fm.run(df_case_chembl_2[0], df_case_chembl_2[1])
    logging.debug(f"\nSubstructure hits for chembl_2:\n{df_aidxf}\n")
    df_fcc = fcc.classify_fragment_combinations(df_aidxf)
    logging.debug(f"\nClean results for chembl_1:\n{df_fcc}\n")
    assert list(df_fcc['abbrev'].values) == ['fbr']
    df_map = fcc.map_frags(df_fcc)
    logging.debug(f"\nFragment map for chembl_1:\n{df_map}\n")
    assert list(df_map['map_str'] == ["678:1[fbr]1141:0"])
