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
def df_fcc():
    """A simple fragment combination classification but complicated enough (cutoff) for testing fragment maps"""
    # raw data extracted from chembl_031_passed_synth_crm_fcc.csv.gz:
    # 3|CHEMBL1421|32|139|connection|false_positive|cutoff|{4, 5, 6, 22, 23, 24}|{8, 9, 10, 11, 27, 28}
    # 4|CHEMBL1421|32|547|connection|monopodal||{4, 5, 6, 22, 23, 24}|{0, 1, 18, 19, 20}
    # 5|CHEMBL1421|32|718|connection|false_positive|cutoff|{4, 5, 6, 22, 23, 24}|{3, 13, 14, 29, 30, 31}
    # 6|CHEMBL1421|32|818|connection|false_positive|cutoff|{4, 5, 6, 22, 23, 24}|{8, 9, 10, 11, 17, 27, 28}
    # 7|CHEMBL1421|139|547|connection|monopodal||{8, 9, 10, 11, 27, 28}|{0, 1, 18, 19, 20}
    # 8|CHEMBL1421|139|718|connection|monopodal||{8, 9, 10, 11, 27, 28}|{3, 13, 14, 29, 30, 31}
    # 9|CHEMBL1421|139|818|fusion|false_positive|substructure|{8, 9, 10, 11, 27, 28}|{8, 9, 10, 11, 17, 27, 28}
    # 10|CHEMBL1421|547|718|connection|false_positive|cutoff|{0, 1, 18, 19, 20}|{3, 13, 14, 29, 30, 31}
    # 11|CHEMBL1421|547|818|connection|monopodal||{0, 1, 18, 19, 20}|{8, 9, 10, 11, 17, 27, 28}
    # 12|CHEMBL1421|718|818|connection|monopodal||{3, 13, 14, 29, 30, 31}|{8, 9, 10, 11, 17, 27, 28}

    return pd.DataFrame({'idm': ['CHEMBL1421'] * 10,
                         'idf1': [32, 32, 32, 32, 139, 139, 139, 547, 547, 718],
                         'idxf1': [0, 0, 0, 0, 1, 1, 1, 2, 2, 3],
                         'idf2': [139, 547, 718, 818, 547, 718, 818, 718, 818, 818],
                         'idxf2': [1, 2, 3, 4, 2, 3, 4, 3, 4, 4],
                         'abbrev': ['cfc', 'cmo', 'cfc', 'cfc', 'cmo', 'cmo', 'ffs', 'cfc', 'cmo', 'cmo'],
                         'category': ['connection'] * 6 + ['fusion'] + ['connection'] * 3,
                         'type': ['false_positive', 'monopodal', 'false_positive', 'false_positive',
                                  'monopodal', 'monopodal', 'false_positive', 'false_positive',
                                  'monopodal', 'monopodal',
                                  ],
                         'subtype': ['cutoff', '', 'cutoff', 'cutoff', '', '', 'substructure',
                                     'cutoff', '', ''],
                         'aidxf1': [{4, 5, 6, 22, 23, 24}, {4, 5, 6, 22, 23, 24},
                                    {4, 5, 6, 22, 23, 24}, {4, 5, 6, 22, 23, 24},
                                    {8, 9, 10, 11, 27, 28}, {8, 9, 10, 11, 27, 28},
                                    {8, 9, 10, 11, 27, 28}, {0, 1, 18, 19, 20},
                                    {0, 1, 18, 19, 20}, {3, 13, 14, 29, 30, 31},
                                    ],
                         'aidxf2': [{8, 9, 10, 11, 27, 28}, {0, 1, 18, 19, 20},
                                    {3, 13, 14, 29, 30, 31}, {8, 9, 10, 11, 17, 27, 28},
                                    {0, 1, 18, 19, 20}, {3, 13, 14, 29, 30, 31},
                                    {8, 9, 10, 11, 17, 27, 28}, {3, 13, 14, 29, 30, 31},
                                    {8, 9, 10, 11, 17, 27, 28}, {8, 9, 10, 11, 17, 27, 28},
                                    ],
                         })


@pytest.fixture
def df_fcc_2():
    """Another simple fragment combination classification but complicated enough (overlap) for testing fragment maps"""

    return pd.DataFrame({'idm': ['CHEMBL209576'] * 9,
                         'idf1': [2, 2, 2, 2, 32, 32, 32, 32, 32],
                         'idxf1': [0, 0, 0, 0, 1, 1, 1, 2, 2],
                         'idf2': [32, 32, 320, 328, 32, 320, 328, 320, 328],
                         'idxf2': [1, 2, 3, 4, 2, 3, 4, 3, 4],
                         'abbrev': ['ffs', 'cmo', 'ffo', 'ffs', 'cmo', 'fed', 'fed', 'cmo', 'cmo'],
                         'category': ['fusion', 'connection', 'fusion', 'fusion', 'connection',
                                      'fusion', 'fusion', 'connection', 'connection',
                                      ],
                         'type': ['false_positive', 'monopodal', 'false_positive', 'false_positive',
                                  'monopodal', 'edge', 'edge', 'monopodal', 'monopodal',
                                  ],
                         'subtype': ['substructure', '', 'overlap', 'substructure', '', '', '',
                                     '', ''],
                         'aidxf1': [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                    {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5},
                                    {0, 1, 2, 3, 4, 5}, {12, 13, 14, 15, 16, 17},
                                    {12, 13, 14, 15, 16, 17},
                                    ],
                         'aidxf2': [{0, 1, 2, 3, 4, 5}, {12, 13, 14, 15, 16, 17},
                                    {4, 5, 6, 7, 8, 9, 10}, {4, 5, 6, 7, 8, 9},
                                    {12, 13, 14, 15, 16, 17}, {4, 5, 6, 7, 8, 9, 10},
                                    {4, 5, 6, 7, 8, 9}, {4, 5, 6, 7, 8, 9, 10},
                                    {4, 5, 6, 7, 8, 9},
                                    ],
                         })


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
    """Check if fusion false_positive cutoff fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_false_positive_cutoff, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_false_positive_cutoff, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'false_positive' and result['subtype'] == 'cutoff' and result['abbrev'] == 'cfc'


def test_fcc_fragmap(fcc, df_fcc, df_fcc_2):
    """Check the fragment map functionality"""

    print(fcc.map_frags(df_fcc))
    print(fcc.map_frags(df_fcc_2))
