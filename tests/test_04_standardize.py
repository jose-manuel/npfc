"""
Module test_04_standardize
====================
Tests for the standardize module.
"""
# standard
from pathlib import Path
import warnings
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
# configure logging
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
# debug
# import logging
# logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def input_files_dupl():
    p = Path('tests/tmp')
    return [str(f) for f in list(p.glob('test_save_dupl_00[1-4].csv.gz'))]


# @pytest.fixture
# def ref_file():
#     return 'tests/tmp/test_dupl_ref.hdf'


@pytest.fixture
def standardizer():
    return Standardizer()


@pytest.fixture
def full_uncharger():
    return FullUncharger()

#
# @pytest.fixture
# def duplicate_filter():
#     return DuplicateFilter()


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
    d = {'valence': 'C1CCCC1O=C',
         'bad_minor_cpd': 'C1CCCCC1.Cl=O(=Cl)(=Cl)(=Cl)',
         'incorrect_smiles': 'C1CCCC',
         'no_smiles': '',
         }
    for k in d.keys():
        d[k] = Chem.MolFromSmiles(d[k], sanitize=False)
    return d


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
    assert set(standardizer.protocol.keys()) == set(['tasks', 'filter_hac', 'filter_molweight', 'filter_nrings', 'filter_medchem'])
    assert standardizer.protocol['tasks'] == ['filter_empty',
                                              'disconnect_metal',
                                              'keep_best',
                                              'filter_hac',
                                              'filter_molweight',
                                              'filter_nrings',
                                              'filter_medchem',
                                              'remove_isotopes',
                                              'normalize',
                                              'uncharge',
                                              'canonicalize',
                                              'remove_stereo',
                                              ]
    # workers
    assert isinstance(standardizer.metal_disconnector, MetalDisconnector)
    assert isinstance(standardizer.normalizer, Normalizer)
    assert isinstance(standardizer.canonicalizer, TautomerCanonicalizer)

    # from a json file
    json_config = 'tests/tmp/std_protocol.json'

    standardizer.protocol = json_config
    assert standardizer.protocol['tasks'] == ["sanitize", "filter_molweight"]
    assert standardizer.protocol['filter_molweight'] == "100.0 <= molweight <= 1000.0"


def test_std_keep_best(standardizer, mols):
    """Test if best fragments are extracted from mixtures."""
    # mol with smaller submols
    mol_clean = standardizer.keep_best(mols['mixture_1'])
    assert Chem.MolToSmiles(mol_clean) == "NCCCCN(CCCN)Cc1ccc(B(O)O)cc1"
    # mol with a larger non-medchem submol
    mol_clean = standardizer.keep_best(mols['mixture_2'])
    assert Chem.MolToSmiles(mol_clean) == "C1CCCCC1"
    # mol with a larger linear submol
    mol_clean = standardizer.keep_best(mols['mixture_3'])
    assert Chem.MolToSmiles(mol_clean) == "C1CCCC1"
    # mol with only one submol
    mol_clean = standardizer.keep_best(mols['tautomer_1'])
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


@pytest.mark.skip  # skip this test to avoid 10s waiting time
def test_standardizer_timeout(mols_timeout, standardizer):
    """Test if timeout is enforced. It is set to 10s and cannot be changed without reinstalling the library."""
    mol, status, task = standardizer.run(mols_timeout['timeout'])
    assert isinstance(mol, Mol) is True
    assert status == 'filtered'
    assert task == 'timeout'


# def test_init_ref(ref_file):
#     """Make sure ref file is computed during this test."""
#     p = Path(ref_file)
#     if p.exists():
#         p.unlink()
#     assert p.exists() is False


# def test_remove_dupl(standardizer, input_files_dupl):
#     """Test the DuplicateFilter class from context of Standardizer"""
#     standardizer.protocol = {'tasks': []}  # remove duplicates only
#     ref_file = 'tests/tmp/test_dupl_ref.hdf'
#     # without ref file
#     passed = 0
#     filtered = 0
#     error = 0
#     for f in input_files_dupl:
#         df = load.file(f)
#         df_passed, df_filtered, df_error = standardizer.run_df(df)
#         passed += len(df_passed.index)
#         filtered += len(df_filtered.index)
#         error += len(df_error.index)
#     assert passed == 6 and filtered == 1 and error == 0
#
#     # with ref_file
#     passed = 0
#     filtered = 0
#     error = 0
#     standardizer.ref_file = ref_file
#     for f in input_files_dupl:
#         df = load.file(f)
#         df_passed, df_filtered, df_error = standardizer.run_df(df)
#         passed += len(df_passed.index)
#         filtered += len(df_filtered.index)
#         error += len(df_error.index)
#
#     assert passed == 4 and filtered == 3 and error == 0
