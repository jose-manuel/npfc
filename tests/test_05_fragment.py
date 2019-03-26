"""
tests.test_fcc
~~~~~~~~~~
Tests for the npfc.fcc module.
"""
# standard library
# data handling
import pandas as pd
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc import fragment

# import logging
# logging.basicConfig(level=logging.DEBUG)


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
def df_mol_fusion_spiro():
    """Example molecule with the fusion spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2(CN1)CCCOC2')]}, index=['m_fs'])


@pytest.fixture
def df_mol_fusion_edge():
    """Example molecule with the fusion edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2CCOCC2N1')]}, index=['m_fe'])


@pytest.fixture
def df_mol_fusion_bridged():
    """Example molecule with the fusion bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CC1COC2')]}, index=['m_fb'])


@pytest.fixture
def df_mol_fusion_unknown():
    """Example molecule with the fusion unknown fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2NC1CO2')]}, index=['m_fu'])


@pytest.fixture
def df_mol_connection_monopodal():
    """Example molecule with the connection monopodal fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC(CN1)C1CCCOC1')]}, index=['m_cm'])


@pytest.fixture
def df_mol_connection_bipodal_spiro():
    """Example molecule with the connection bipodal spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2(CCCOC2)C2CNCC12')]}, index=['m_cbs'])


@pytest.fixture
def df_mol_connection_bipodal_edge():
    """Example molecule with the connection bipodal edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NCC2C1CCC1OCCCC21')]}, index=['m_cbe'])


@pytest.fixture
def df_mol_connection_bipodal_bridged():
    """Example molecule with the connection bipodal bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CC1C1CCCOC21')]}, index=['m_cbb'])


@pytest.fixture
def df_mol_connection_tripodal_spiro():
    """Example molecule with the connection tripodal spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1C2CNC3CCC4CCOCC14C23')]}, index=['m_cts'])


@pytest.fixture
def df_mol_connection_tripodal_edge():
    """Example molecule with the connection tripodal edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1NC2CCC3CCOC4CCC1C2C34')]}, index=['m_cte'])


@pytest.fixture
def df_mol_connection_tripodal_bridged():
    """Example molecule with the connection tripodal bridged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1C2NC3CCC4CCOC2C4C13')]}, index=['m_ctb'])


@pytest.fixture
def df_mol_connection_unknown_spiro():
    """Example molecule with the connection unknown spiro fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1C2C3C4CCC5OCC(CC135)CC2N4')]}, index=['m_cus'])


@pytest.fixture
def df_mol_connection_unknown_edge():
    """Example molecule with the connection unknown edge fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2OCC3CC4NC1C1C4CCC3C21')]}, index=['m_cue'])


@pytest.fixture
def df_mol_connection_unknown_bridged():
    """Example molecule with the connection unknown brdiged fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C1CC2C3COC4CCC5OC(CC1C5C24)C3')]}, index=['m_cub'])


@pytest.fixture
def df_mol_fusion_false_positive_substructure():
    """Example molecule with the connection false_positive substructure fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('CCC1CCNC1')]}, index=['m_cfps'])


@pytest.fixture
def df_frags_fusion_false_positive_substructure():
    """Two overlapping fragments used to identify false_positive substructure fragment combinations."""
    df_frags = pd.DataFrame({'smiles': ['C1NCCC1', 'CC1CCNC1']}, index=['QA', 'QF'])
    df_frags['mol'] = df_frags['smiles'].map(Chem.MolFromSmiles)
    return df_frags


@pytest.fixture
def df_mol_connection_false_positive_cutoff():
    """Example molecule with the connection false_positive cutoff fragment combination."""
    return pd.DataFrame({'mol': [Chem.MolFromSmiles('C(CCC1CCCOC1)CC1CCNC1')]}, index=['m_cfpc'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_fcc_run_substructure_match(fcc, fm, df_mol_fusion_spiro, df_frags):
    """Check that the df_aidxf obtained by substructure match is adequate."""
    df_aidxf = fm.run(df_mol_fusion_spiro, df_frags)
    df_aidxf_ref = pd.DataFrame({'idm': ['m_fs'] * 2,
                                 'idf': ['QA', 'QB'],
                                 'aidxf': [{0, 1, 2, 3, 4}, {2, 5, 6, 7, 8, 9}],
                                 'mol_perc': [50.0, 60.0]},
                                index=[0, 1])
    assert df_aidxf.equals(df_aidxf_ref)


def test_fcc_fusion_spiro(fcc, fm, df_mol_fusion_spiro, df_frags):
    """Check if fusion spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_spiro, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'spiro' and result['subtype'] == '' and result['abbrev'] == 'fsp'


def test_fcc_fusion_edge(fcc, fm, df_mol_fusion_edge, df_frags):
    """Check if fusion edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_edge, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'edge' and result['subtype'] == '' and result['abbrev'] == 'fed'


def test_fcc_fusion_bridged(fcc, fm, df_mol_fusion_bridged, df_frags):
    """Check if fusion bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_bridged, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'bridged' and result['subtype'] == '' and result['abbrev'] == 'fbr'


def test_fcc_fusion_unknown(fcc, fm, df_mol_fusion_unknown, df_frags):
    """Check if fusion unknown fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_unknown, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_unknown, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'unknown' and result['subtype'] == '' and result['abbrev'] == 'fun'


def test_fcc_connection_monopodal(fcc, fm, df_mol_connection_monopodal, df_frags):
    """Check if connection monopodal fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_monopodal, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_monopodal, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'monopodal' and result['subtype'] == '' and result['abbrev'] == 'cmo'


def test_fcc_connection_bipodal_spiro(fcc, fm, df_mol_connection_bipodal_spiro, df_frags):
    """Check if connection bipodal spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_bipodal_spiro, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cbs'


def test_fcc_connection_bipodal_edge(fcc, fm, df_mol_connection_bipodal_edge, df_frags):
    """Check if connection bipodal edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_bipodal_edge, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cbe'


def test_fcc_connection_bipodal_bridged(fcc, fm, df_mol_connection_bipodal_bridged, df_frags):
    """Check if connection bipodal bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_bipodal_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_bipodal_bridged, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'bipodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'cbb'


def test_fcc_connection_tripodal_spiro(fcc, fm, df_mol_connection_tripodal_spiro, df_frags):
    """Check if connection tripodal spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_tripodal_spiro, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'spiro' and result['abbrev'] == 'cts'


def test_fcc_connection_tripodal_edge(fcc, fm, df_mol_connection_tripodal_edge, df_frags):
    """Check if connection tripodal edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_tripodal_edge, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'edge' and result['abbrev'] == 'cte'


def test_fcc_connection_tripodal_bridged(fcc, fm, df_mol_connection_tripodal_bridged, df_frags):
    """Check if connection tripodal bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_tripodal_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_tripodal_bridged, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'tripodal' and result['subtype'] == 'bridged' and result['abbrev'] == 'ctb'


def test_fcc_connection_unknown_spiro(fcc, fm, df_mol_connection_unknown_spiro, df_frags):
    """Check if connection unknown spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_unknown_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_unknown_spiro, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'unknown' and result['subtype'] == 'spiro' and result['abbrev'] == 'cus'


@pytest.mark.skip   # could not come up with an example!
def test_fcc_connection_unknown_edge(fcc, fm, df_mol_connection_unknown_edge, df_frags):
    """Check if connection unknown edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_unknown_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_unknown_edge, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'unknown' and result['subtype'] == 'edge' and result['abbrev'] == 'cue'


def test_fcc_connection_unknown_bridged(fcc, fm, df_mol_connection_unknown_bridged, df_frags):
    """Check if connection unknown bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_unknown_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_unknown_bridged, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'unknown' and result['subtype'] == 'bridged' and result['abbrev'] == 'cub'


def test_fcc_fusion_false_positive_substructure(fcc, fm, df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure):
    """Check if fusion false_positive substructure fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_false_positive_substructure, df_frags_fusion_false_positive_substructure)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_false_positive_substructure, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'false_positive' and result['subtype'] == 'substructure' and result['abbrev'] == 'ffs'


def test_fcc_connection_false_positive_cutoff(fcc, fm, df_mol_connection_false_positive_cutoff, df_frags):
    df_aidxf = fm.run(df_mol_connection_false_positive_cutoff, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_false_positive_cutoff, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'false_positive' and result['subtype'] == 'cutoff' and result['abbrev'] == 'cfc'
