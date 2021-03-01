"""
Module test_04_standardize
===========================
Tests for the standardize module.
"""
# standard
from pathlib import Path
import warnings
from copy import deepcopy
# data science
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Mol
from rdkit.Chem import rdinchi
from rdkit.Chem.MolStandardize.metal import MetalDisconnector
from rdkit.Chem.MolStandardize.normalize import Normalizer
from rdkit.Chem.MolStandardize.tautomer import TautomerCanonicalizer
# tests
import pytest
from npfc.standardize import Standardizer
from npfc.standardize import FullUncharger
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def input_files_dupl():
    p = Path('tests/tmp')
    return [str(f) for f in list(p.glob('test_save_dupl_00[1-4].csv.gz'))]


@pytest.fixture
def standardizer():
    return Standardizer()


@pytest.fixture
def full_uncharger():
    return FullUncharger()


@pytest.fixture
def mols():
    """An example of a DataFrame with molecules."""
    d = {'empty': '',
         'metal': 'CC(C)(C)[N+](=Cc1ccc(cc1S(=O)(=O)O[Na])S(=O)(=O)O[Na])[O-]',
         'mixture_1': 'Cl.Cl.Cl.NCCCCN(CCCN)Cc1ccc(cc1)B(O)O',
         'mixture_2': '[Rn].C1CCC1.C1CCCCC1',
         'mixture_3': 'C1CCCC1.CCCCCCCCCC',
         'isotope': 'F[14C](F)(F)C(Cl)C1CCCC1',
         'normalize': '[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1',
         'neutral_all': 'C1CCC1C[C@H]([NH3+])C([O-])=C',
         'neutral_neg': 'CC(C)(C)C1CCc2nnc[n+]([O-])c2C1',
         'tautomer_1': 'O=C1CC=CC=C1',
         'tautomer_2': 'OC1=CC=CC=C1',
         'stereo_chiral': 'CN1CCC[C@H]1CN2CCN(Cc3ccncc3)CC2',
         'stereo_doublebond': r'Br\C=C\1/CCC(C(=O)O1)c2cccc3ccccc23',
         'inorganic': 'OP1(N=P(N=P(N=P(N=1)(O)O)(O)O)(O)O)O',
         'mw_too_large': 'CCCN(NC(=O)[C@H]1CCCN1C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)C)[C@@H](C)O)C(C)C)C(C)C)C(=O)N[C@@H](CO)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](Cc3c[nH]c4ccccc34)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)O',
         'hac_too_small': 'C1CC1',
         'non_medchem': 'C1CCC(CC1)[Si](CCCN2CCCCC2)(C3=CC=CC=C3)O',
         'linear': 'CCCCCC',
         'dupl_1': 'C1CCCCC1',
         'dupl_2': 'C1CCCCC1',
         }
    for k in d.keys():
        d[k] = Chem.MolFromSmiles(d[k])
    return d


@pytest.fixture
def mols_timeout():
    """An example of a DataFrame with a molecule filtered because of timeout (>10s)."""
    d = {'timeout': 'Nc1nc(N)c2nc(CNc3ccc(cc3)C(=O)NC(CCCNC(=O)c4ccccc4)C(=O)O)cnc2n1',
         }
    for k in d.keys():
        d[k] = Chem.MolFromSmiles(d[k])
    return d


@pytest.fixture
def mols_bad():
    """An example of a DataFrame with molecules with errors."""
    d = {'incorrect_valence': 'C1CCCC1O=C',
         'bad_minor_cpd': 'C1CCCCC1.Cl=O(=Cl)(=Cl)(=Cl)',
         'incorrect_smiles': 'C1CCCC',
         'no_smiles': '',
         }
    for k in d.keys():
        d[k] = Chem.MolFromSmiles(d[k], sanitize=False)
    return d


@pytest.fixture
def df_fragments():
    """An example of a DataFrame with molecules from which to extract Murcko Scaffolds."""
    df = pd.DataFrame([['simple', 'Oc1ccccc1'],
                       ['minor_cpd', 'Oc1ccccc1.O'],
                       ['protB_crms_19', 'CCC12CCCCC1[N+](=O)CCC2']
                       ], columns=['idm', 'mol'])
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    df.index = df['idm']

    return df


@pytest.fixture
def df_deglyco_scaffold_hunter():
    """Examples of molecules to deglycosylate as described in the SCONP paper SI."""
    df = pd.DataFrame([['CC(C)=CC1CC(C)(O)C2C3CCC4C5(C)CCC(O)C(C)(C)C5CCC4(C)C33COC2(C3)O1'],
                       ['CC(C)=CC1CC(C)(O)C2C3CCC4C5(C)CCC(OC6OCC(O)C(O)C6O)C(C)(C)C5CCC4(C)C33COC2(C3)O1'],
                       ['CC1OC(OC2(C)CC(OC34CC5(CO3)C(CCC3C6(C)CCC(OC7OC(CO)C(O)C(O)C7O)C(C)(C)C6CCC53C)C24)C=C(C)C)C(O)C(O)C1O'],
                       ['CC1OC(OC2C(O)C(CO)OC(OC3C(O)COC(OC4CCC5(C)C(CCC6(C)C5CCC5C7C8(CC65CO8)OC(CC7(C)O)C=C(C)C)C4(C)C)C3O)C2O)C(O)C(OC(C)=O)C1O'],
                       ['CC(=O)OC1C(OC2CCC3(C)C(CCC4(C)C3CCC3C5C6(CC43CO6)OC(CC5(C)O)C=C(C)C)C2(C)C)OCC(O)C1OC1OC(CO)C(O)C(O)C1OC1OCC(O)C(O)C1O'],
                       ['CC1OC(OC2C(OC3CCC4(C)C(CCC5(C)C4CCC4C6C7(CC54CO7)OC(CC6(C)O)C=C(C)C)C3(C)C)OCC(O)C2OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'],
                       ['CC1OC(OC2C(OC3CCC4(C)C(CCC5(C)C4CCC4C6C7(CC54CO7)OC(CC6(C)OC4OC(C)C(O)C(O)C4O)C=C(C)C)C3(C)C)OCC(O)C2O)C(O)C(O)C1O'],
                       ['CC(C)=CC1CC(C)(OC2OCC(O)C(O)C2O)C2C3CCC4C5(C)CCC(OC6OCC(O)C(OC7OC(CO)C(O)C(O)C7O)C6OC6OC(CO)C(O)C6O)C(C)(C)C5CCC4(C)C33COC2(C3)O1'],
                       ['CC(=O)OC1C(O)C(CO)OC(OC2C(O)COC(OC3CCC4(C)C(CCC5(C)C4CCC4C6C7(CC54CO7)OC(CC6(C)O)C=C(C)C)C3(C)C)C2OC(C)=O)C1OC1OCC(O)C(O)C1O'],
                       ], columns=['mol'])
    df['smiles_ref'] = 'CC(C)=CC1CC(C)(O)C2C3CCC4C5(C)CCC(O)C(C)(C)C5CCC4(C)C34COC2(C4)O1'
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    return df


@pytest.fixture
def df_deglyco_extra_tests():
    """Examples of molecules to deglycosylate."""
    df = pd.DataFrame([['OC1COC(CC1O)OC1CC2C(CCCC2)CC1', 'OC1CCC2CCCCC2C1', '1 aglycan, 1 glycan'],
                       ['OC1CCC(OC2CC3C(CCCC3)CC2)OC1', 'OC1CCC(OC2CCC3CCCCC3C2)OC1', '2 aglycans'],
                       ['OCC1OC(OC2CC(O)C(O)CO2)C(O)C(O)C1O', 'OCC1OC(OC2CC(O)C(O)CO2)C(O)C(O)C1O', '2 glycans'],
                       ['OC1COC(CC2CC3C(CCCC3)CC2)CC1O', 'OC1COC(CC2CCC3CCCCC3C2)CC1O', '1 aglycan, 1 glycan but no glycosidic bond'],
                       ['OC1COC(OC2CCC3CCCCC3C2)C1O', 'OC1CCC2CCCCC2C1', '1 aglycan, 1 5-glycan'],
                       ['OC1CCOC1OC1CCC2CCCCC2C1', 'OC1CCC2CCCCC2C1', '1 aglycan, 1 5-glycan'],
                       ['OCC1CCOC1OC1CCC2CCCCC2C1', 'OCC1CCOC1OC1CCC2CCCCC2C1', '2 aglycans'],
                       ['OC1COCC(O)C1OC1CCC2CCCCC2C1', 'OC1COCC(O)C1OC1CCC2CCCCC2C1', '2 aglycans, but one could be matched for a glycan if one looks only for side chains'],
                       ['OC1CC(OC2CCC3CCCCC3C2)OCC1OC1CC2C(CCCC2)CC1', 'OC1CC(OC2CCC3CCCCC3C2)OCC1OC1CCC2CCCCC2C1', '2 aglycans, 1 glycan with glycosidic bonds but glycan is not terminal'],
                       ['OC1CC(COCC2CCC3CCCCC3C2)OC(O)C1O', 'OCC1CCC2CCCCC2C1', '1 glycan with glycosidic bond (COC)'],
                       ['OC1OC(COC2CCC3CCCCC3C2)C(O)C(O)C1O', 'OC1CCC2CCCCC2C1', '1 glycan with glycosidic bond (CO)'],
                       ['OC1COC(CC1O)OCC1CCC2CCCCC2C1', 'OCC1CCC2CCCCC2C1', '1 glycan with glycosidic bond (OC)'],
                       ['N[C@H]1C(O)O[C@H](COC2CCC3CCCCC3C2)[C@@H](O)[C@@H]1O', 'N[C@H]1C(O)O[C@H](COC2CCC3CCCCC3C2)[C@@H](O)[C@@H]1O', 'aglycan with sugar derivative (glycosamine)'],
                       ['CC(=O)C[C@H]1C(O)O[C@H](COC2CCC3CCCCC3C2)[C@@H](O)[C@@H]1O', 'CC(=O)C[C@H]1C(O)O[C@H](COC2CCC3CCCCC3C2)[C@@H](O)[C@@H]1O', 'aglycan with sugar another derivative'],
                       ['CCCC(CC)OC[C@H]1OC(O)C(O)C(O)C1O', 'CCCC(CC)OC[C@H]1OC(O)C(O)C(O)C1O', 'linear glycan with 1 glycan'],
                       ['CC(=O)OC1C(O)OC(COC2CCC3CCCCC3C2)C(OC(C)=O)C1OC(C)=O', 'OC1CCC2CCCCC2C1', 'aglycan with a less typical glycan (maybe false positive?)'],
                       ['CC(=O)OCC1OC(OC2CCC3(C)C(CCC4(C)C3CC(O)C3C(C(O)(CCC=C(C)C)COC5OC(CO)C(O)C(O)C5O)CCC34C)C2(C)C)C(OC2OC(C)C(O)C(O)C2O)C(OC2OCC(O)C(O)C2O)C1O', 'CC(C)=CCCC(O)(CO)C1CCC2(C)C1C(O)CC1C3(C)CCC(O)C(C)(C)C3CCC12C', 'aglycan with 4 glycans, one of those having a CCO linker'],
                       ], columns=['mol', 'smiles_deglyco', 'description'])
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    return df

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_fu_full_uncharge(full_uncharger, mols):
    """Test if all charges are removed when possible."""
    # neutral_all
    assert Chem.MolToSmiles(mols['neutral_all']) == "C=C([O-])[C@@H]([NH3+])CC1CCC1"
    mol_clean = full_uncharger.full_uncharge(mols['neutral_all'])
    assert Chem.MolToSmiles(mol_clean) == "C=C(O)[C@@H](N)CC1CCC1"
    # neutral_neg
    assert Chem.MolToSmiles(mols['neutral_neg']) == "CC(C)(C)C1CCc2nnc[n+]([O-])c2C1"
    mol_clean = full_uncharger.full_uncharge(mols['neutral_neg'])
    assert Chem.MolToSmiles(mol_clean) == "CC(C)(C)C1CCc2nnc[n+](O)c2C1"


def test_std_disconnect_metal(standardizer, mols):
    """Test if the MetalDisconnector works as expected to disconnect metals."""
    mol_clean = standardizer.metal_disconnector.disconnect(mols['metal'])
    assert len(Chem.GetMolFrags(mol_clean)) == 3


def test_std_init(standardizer):

    # default parameters
    standardizer = deepcopy(standardizer)  # if the standardizer is modified in place, then other tests will fail
    assert set(standardizer.protocol.keys()) == set(['tasks', 'filter_num_heavy_atoms', 'filter_molecular_weight', 'filter_num_rings', 'filter_elements'])
    assert standardizer.protocol['tasks'] == ['filter_empty',
                                              'disconnect_metal',
                                              'clear_mixtures',
                                              'deglycosylate',
                                              'filter_num_heavy_atoms',
                                              'filter_molecular_weight',
                                              'filter_num_rings',
                                              'filter_elements',
                                              'sanitize',
                                              'clear_isotopes',
                                              'normalize',
                                              'uncharge',
                                              'canonicalize',
                                              ]
    # workers
    assert isinstance(standardizer.metal_disconnector, MetalDisconnector)
    assert isinstance(standardizer.normalizer, Normalizer)
    assert isinstance(standardizer.canonicalizer, TautomerCanonicalizer)

    # from a json file
    json_config = 'tests/tmp/std_protocol.json'
    standardizer.protocol = json_config
    assert standardizer.protocol['tasks'] == ["sanitize", "filter_molecular_weight"]
    assert standardizer.protocol['filter_molecular_weight'] == "100.0 <= molecular_weight <= 1000.0"


def test_std_clear_mixtures(standardizer, mols):
    """Test if best fragments are extracted from mixtures."""
    # mol with smaller submols
    mol_clean = standardizer.clear_mixtures(mols['mixture_1'])
    assert Chem.MolToSmiles(mol_clean) == "NCCCCN(CCCN)Cc1ccc(B(O)O)cc1"
    # mol with a larger non-medchem submol
    mol_clean = standardizer.clear_mixtures(mols['mixture_2'])
    assert Chem.MolToSmiles(mol_clean) == "C1CCCCC1"
    # mol with a larger linear submol
    mol_clean = standardizer.clear_mixtures(mols['mixture_3'])
    assert Chem.MolToSmiles(mol_clean) == "C1CCCC1"
    # mol with only one submol
    mol_clean = standardizer.clear_mixtures(mols['tautomer_1'])
    assert Chem.MolToSmiles(mol_clean) == "O=C1C=CC=CC1"


def test_std_normalize(mols, standardizer):
    """Test if wrongly defined functional groups are normalized."""
    assert Chem.MolToSmiles(mols['normalize']) == "O=C(O[Na])c1ccc(C[S+2]([O-])[O-])cc1"
    mol_clean = standardizer.normalizer.normalize(mols['normalize'])
    assert Chem.MolToSmiles(mol_clean) == "O=C(O[Na])c1ccc(C[S](=O)=O)cc1"


def test_std_canonicalize(mols, standardizer):
    """Test if inchikeys obtained from different tautomers are the same after canonicalization."""
    # case 1
    assert rdinchi.MolToInchiKey(mols['tautomer_1']) == "WQPDQJCBHQPNCZ-UHFFFAOYSA-N"
    mol_clean = standardizer.canonicalizer.canonicalize(mols['tautomer_1'])
    assert rdinchi.MolToInchiKey(mol_clean) == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"

    # case 2
    assert rdinchi.MolToInchiKey(mols['tautomer_2']) == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"
    mol_clean = standardizer.canonicalizer.canonicalize(mols['tautomer_2'])
    assert rdinchi.MolToInchiKey(mol_clean) == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"


def test_std_remove_stereochemistry(mols):
    """Test if all stereochemistry centers (chiral and double bonds) are removed."""
    # chirality
    assert Chem.FindMolChiralCenters(mols['stereo_chiral'], includeUnassigned=True) == [(5, 'S')]
    Chem.RemoveStereochemistry(mols['stereo_chiral'])  # mol_ini is modified inplace
    assert Chem.FindMolChiralCenters(mols['stereo_chiral'], includeUnassigned=True) == [(5, '?')]
    # doublebond
    stereo_doublebond = [b.GetStereo() for b in mols['stereo_doublebond'].GetBonds() if b.GetStereo() != Chem.rdchem.BondStereo.STEREONONE]
    assert stereo_doublebond == [Chem.rdchem.BondStereo.STEREOE]
    Chem.RemoveStereochemistry(mols['stereo_doublebond'])  # mol_ini is modified inplace
    stereo_doublebond = [b.GetStereo() for b in mols['stereo_doublebond'].GetBonds() if b.GetStereo() != Chem.rdchem.BondStereo.STEREONONE]
    assert stereo_doublebond == []  # ideally it should be set to Chem.rdchem.BondStereo.STEREOANY, but whatever obscure reason it is set to STEREONONE...


def test_run_protocol(standardizer, mols, mols_bad):
    """Run default standardization protocol."""
    # initiate a global df_mols
    d = {}
    d['mol'] = []
    d['idm'] = []
    for k, m in mols.items():
        d['idm'].append(k)
        d['mol'].append(m)
    for k, m in mols_bad.items():
        d['idm'].append(k)
        d['mol'].append(m)
    df_mols = pd.DataFrame(d)
    assert len(df_mols.index) == 24
    # run default protocol
    df_passed, df_filtered, df_error = standardizer.run_df(df_mols)
    assert len(df_passed.index) == 15
    assert len(df_error.index) == 3
    assert len(df_filtered.index) == 6  # duplicates are not found


def test_standardizer_timeout(mols_timeout, standardizer):
    """Test if timeout is enforced. It is set to 10s and cannot be changed without reinstalling the library."""
    mol, status, task = standardizer.run(mols_timeout['timeout'], timeout=1)
    assert isinstance(mol, Mol) is True
    assert status == 'filtered'
    assert task == 'timeout'


def test_standardizer_murcko_scaffolds(df_fragments, standardizer):
    """Run Murcko Scaffold Extraction on fragment dataset. Protocols A and B are performed within the same function."""
    standardizer.protocol = {"tasks": ['clear_mixtures', 'extract_murcko']}
    df_passed, df_filtered, df_error = standardizer.run_df(df_fragments)
    assert len(df_filtered) == 0
    assert len(df_error) == 0

    # simple
    assert Chem.MolToSmiles(df_passed.loc['simple']['mol']) == 'c1ccccc1'
    # minor_cpd
    assert Chem.MolToSmiles(df_passed.loc['minor_cpd']['mol']) == 'c1ccccc1'
    # protB
    assert Chem.MolToSmiles(df_passed.loc['protB_crms_19']['mol']) == 'O=[N+]1CCCC2CCCCC21'


def test_standardizer_deglycosylate_scaffold_hunter(df_deglyco_scaffold_hunter, standardizer):
    """Deglycosylate examples molecules from the SCONP paper. All molecules should give the same result."""
    for i in range(len(df_deglyco_scaffold_hunter)):
        row = df_deglyco_scaffold_hunter.iloc[i]
        assert Chem.MolToSmiles(standardizer.deglycosylate(row['mol'])) == row['smiles_ref']


def test_standardizer_deglycosylate_extra_tests(df_deglyco_extra_tests, standardizer):
    """Deglycosylate extra examples molecules. Results should match expected smiles."""
    for i in range(len(df_deglyco_extra_tests)):
        row = df_deglyco_extra_tests.iloc[i]
        assert Chem.MolToSmiles(standardizer.deglycosylate(row['mol'])) == row['smiles_deglyco']
