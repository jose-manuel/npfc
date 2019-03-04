"""
Module test_04_standardize
====================
Tests for the standardize module.
"""
# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc.standardize import Standardizer
from npfc.standardize import FullUncharger
from npfc.standardize import DuplicateFilter


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


def test_standardizer_full_uncharge(full_uncharger, mols):
    """Test if all charges are removed when possible."""
    # neutral_all
    assert Chem.MolToSmiles(mols['neutral_all']) == "C=C([O-])[C@@H]([NH3+])CC1CCC1"
    mol_clean = full_uncharger.full_uncharge(mols['neutral_all'])
    assert Chem.MolToSmiles(mol_clean) == "C=C(O)[C@@H](N)CC1CCC1"
    # neutral_neg
    assert Chem.MolToSmiles(mols['neutral_neg']) == "CC(C)(C)C1CCc2nnc[n+]([O-])c2C1"
    mol_clean = full_uncharger.full_uncharge(mols['neutral_neg'])
    assert Chem.MolToSmiles(mol_clean) == "CC(C)(C)C1CCc2nnc[n+](O)c2C1"
