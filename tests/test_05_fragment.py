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
def df_fcc_clean_5():
    """An example of fcc_clean with 5 combinations"""
    # raw data extracted from chembl_182_passed_synth_crm_fcc.csv.gz:
    # idm|idf1|idxf1|fid1|idf2|idxf2|fid2|abbrev|category|type|subtype|aidxf1|aidxf2
    # CHEMBL3991441|2|11|2:11|178|15|178:15|fed|fusion|edge||{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}|{0, 1, 10, 11, 12, 13, 14, 15, 16}
    # CHEMBL3991441|2|11|2:11|718|17|718:17|cmo|connection|monopodal||{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}|{19, 20, 21, 22, 23, 24}
    # CHEMBL3991441|2|11|2:11|718|18|718:18|cmo|connection|monopodal||{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}|{32, 27, 28, 29, 30, 31}
    # CHEMBL3991441|178|15|178:15|718|17|718:17|cmo|connection|monopodal||{0, 1, 10, 11, 12, 13, 14, 15, 16}|{19, 20, 21, 22, 23, 24}
    # CHEMBL3991441|178|15|178:15|718|18|718:18|cmo|connection|monopodal||{0, 1, 10, 11, 12, 13, 14, 15, 16}|{32, 27, 28, 29, 30, 31}

    return pd.DataFrame({'idm': ['CHEMBL3991441'] * 5,
                         'fid1': ['2:11'] * 3 + ['178:15'] * 2,
                         'fid2': ['178:15', '718:17', '718:18', '718:17', '718:18'],
                         'abbrev': ['fed'] + ['cmo'] * 4,

                         'aidxf1': ['{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}'] * 3 + ['{0, 1, 10, 11, 12, 13, 14, 15, 16}'] * 2,
                         'aidxf2': ['{0, 1, 10, 11, 12, 13, 14, 15, 16}', '{19, 20, 21, 22, 23, 24}',
                                    '{32, 27, 28, 29, 30, 31}', '{19, 20, 21, 22, 23, 24}',
                                    '{32, 27, 28, 29, 30, 31}',
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
                         'aidxf1': ['{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}', '{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}',
                                    '{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}', '{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}',
                                    '{0, 1, 2, 3, 4, 5}', '{0, 1, 2, 3, 4, 5}',
                                    '{0, 1, 2, 3, 4, 5}', '{12, 13, 14, 15, 16, 17}',
                                    '{12, 13, 14, 15, 16, 17}',
                                    ],
                         'aidxf2': ['{0, 1, 2, 3, 4, 5}', '{12, 13, 14, 15, 16, 17}',
                                    '{4, 5, 6, 7, 8, 9, 10}', '{4, 5, 6, 7, 8, 9}',
                                    '{12, 13, 14, 15, 16, 17}', '{4, 5, 6, 7, 8, 9, 10}',
                                    '{4, 5, 6, 7, 8, 9}', '{4, 5, 6, 7, 8, 9, 10}',
                                    '{4, 5, 6, 7, 8, 9}',
                                    ],
                         })


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
    df_aidxf_ref = pd.DataFrame({'idm': ['fsp'] * 2,
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


@pytest.mark.skip  # no example for fot yet!
def test_fcc_fusion_other(fcc, fm, df_mol_fusion_other, df_frags):
    """Check if fusion other fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_fusion_other, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_fusion_other, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'fusion' and result['type'] == 'other' and result['subtype'] == '' and result['abbrev'] == 'fot'


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


def test_fcc_connection_other_spiro(fcc, fm, df_mol_connection_other_spiro, df_frags):
    """Check if connection other spiro fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_spiro, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_other_spiro, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'spiro' and result['abbrev'] == 'cos'


@pytest.mark.skip  # classified as bridged now...
def test_fcc_connection_other_edge(fcc, fm, df_mol_connection_other_edge, df_frags):
    """Check if connection other edge fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_edge, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_other_edge, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'edge' and result['abbrev'] == 'coe'


def test_fcc_connection_other_bridged(fcc, fm, df_mol_connection_other_bridged, df_frags):
    """Check if connection other bridged fragment combinations are identified."""
    df_aidxf = fm.run(df_mol_connection_other_bridged, df_frags)
    df_fcc = fcc.classify_fragment_combinations(df_mol_connection_other_bridged, df_aidxf)
    result = df_fcc.iloc[0]
    assert result['category'] == 'connection' and result['type'] == 'other' and result['subtype'] == 'bridged' and result['abbrev'] == 'cob'


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


def test_fcc_fragmap(fcc, df_fcc_clean_5):
    """Check the fragment map functionality"""

    df_map = fcc.map_frags(df_fcc_clean_5)
    logging.debug("")
    logging.debug(df_map)
