"""
Module test_03_filter
====================
Tests for the filter module.
"""

# chemoinformatics
from rdkit import Chem
# tests
import pytest
from npfc.filter import Filter


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

@pytest.fixture
def filter():
    return Filter()


@pytest.fixture
def mol():
    return Chem.MolFromSmiles('C1CCCCC1')  # hac=6, molweight=84.1595, nrings=1, elements={'C'}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_filter(filter, mol):
    assert filter.filter_mol(mol, 'hac > 10') is True
