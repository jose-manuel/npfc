"""
Module test_04_standardize
====================
Tests for the standardize module.
"""
# data science
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdinchi
from rdkit.Chem.MolStandardize.metal import MetalDisconnector
from rdkit.Chem.MolStandardize.normalize import Normalizer
from rdkit.Chem.MolStandardize.tautomer import TautomerCanonicalizer
# tests
import pytest
from npfc.standardize import Standardizer
from npfc.standardize import FullUncharger
from npfc.standardize import DuplicateFilter

# configure logging
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
# import logging
# logging.basicConfig(level=logging.DEBUG)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def standardizer():
    return Standardizer()


@pytest.fixture
def full_uncharger():
    return FullUncharger()


@pytest.fixture
def duplicate_filter():
    return DuplicateFilter()


@pytest.fixture
def mols():
    """An example of a DataFrame with molecules."""
    d = {'metal': 'CC(C)(C)[N+](=Cc1ccc(cc1S(=O)(=O)O[Na])S(=O)(=O)O[Na])[O-]',
         'mixture_1': 'Cl.Cl.Cl.NCCCCN(CCCN)Cc1ccc(cc1)B(O)O',
         'mixture_2': '[Rn].C1CCC1.C1CCCCC1',
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
    assert standardizer.protocol['tasks'] == ['sanitize',
                                              'disconnect_metal',
                                              'keep_largest',
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
    with open(json_config, 'w') as JSON:
        JSON.write('''{
                    "tasks": ["sanitize", "filter_molweight"],
                    "filter_molweight": "100.0 <= molweight <= 1000.0"\n}''')
    standardizer.protocol = json_config
    assert standardizer.protocol['tasks'] == ["sanitize", "filter_molweight"]
    assert standardizer.protocol['filter_molweight'] == "100.0 <= molweight <= 1000.0"


def test_std_keep_largest(standardizer, mols):
    """Test if largest fragments are extracted from mixtures."""
    # mol with smaller fragments and store molweight property
    mol_clean = standardizer.keep_largest(mols['mixture_1'])
    assert Chem.MolToSmiles(mol_clean) == "NCCCCN(CCCN)Cc1ccc(B(O)O)cc1"
    # mol with larger non-medchem fragments, do not store molweight property
    mol_clean = standardizer.keep_largest(mols['mixture_2'])
    assert Chem.MolToSmiles(mol_clean) == "C1CCCCC1"
    # mol with only one fragment without sanitize
    mol_clean = standardizer.keep_largest(mols['tautomer_1'])
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


def test_remove_dupl(standardizer):
    """Remove duplicates accross chunks using a syn file"""
    pass


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
    assert len(df_mols.index) == 22
    # run default protocol
    df_passed, df_filtered, df_error = standardizer.run_df(df_mols)
    assert len(df_passed.index) == 11
    assert len(df_error.index) == 3
    assert len(df_filtered.index) == 8


@pytest.mark.skip   # decorators are set only once, so the TIMEOUT value cannot be changed inside of the decorators. no idea for a work-around now so just skip the test as it is long.
def test_standardizer_timeout(mols_timeout, standardizer_fast):
    """Test if timeout is enforced and can be configured."""
    mol, status, task = standardizer_fast.run(mols_timeout['timeout'])
    assert isinstance(mol, Chem.Mol) is True
    assert status == 'filtered'
    assert task == 'timeout'
