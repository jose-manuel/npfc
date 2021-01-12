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
# debug
import logging
logging.basicConfig(level=logging.ERROR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@pytest.fixture
def filter():
    return Filter()


@pytest.fixture
def mol():
    return Chem.MolFromSmiles('C1CCCCC1')  # hac=6, molweight=84.0939, nrings=1, elements={'C'}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_filter_parse_expression(filter, mol):
    """Test the parsing of the expression used to define the filter"""
    assert filter.filter_mol(mol, 'num_heavy_atoms > 10') is False
    assert filter.filter_mol(mol, 'num_heavy_atoms >= 10') is False
    assert filter.filter_mol(mol, 'num_heavy_atoms == 6') is True
    assert filter.filter_mol(mol, '100 <= molecular_weight <= 1000') is False
    assert filter.filter_mol(mol, 'elements in C') is True
    assert filter.filter_mol(mol, 'elements not in O') is True
    assert filter.filter_mol(mol, 'num_rings != 0') is True


def test_filter_compute_descriptors(filter, mol):
    """Test the computation of descriptors (all or just subset)"""
    # all descriptors
    result = filter.compute_descriptors(mol)
    assert sorted(list(result.keys())) == ['elements',
                                           'molecular_formula',
                                           'molecular_weight',
                                           'num_atoms_nitrogen',
                                           'num_atoms_oxygen',
                                           'num_hba',
                                           'num_hbd',
                                           'num_heavy_atoms',
                                           'num_rings',
                                           'num_rings_arom',
                                           'num_rotatable_bonds',
                                           'ring_size_max',
                                           'ring_size_min',
                                           'slogp',
                                           'tpsa',
                                           ]
    # only subset of descriptors
    result = filter.compute_descriptors(mol, descriptors=['num_heavy_atoms', 'molecular_weight'])
    assert sorted(list(result.keys())) == ['molecular_weight', 'num_heavy_atoms']

    # no descriptor
    try:
        result = filter.compute_descriptors(mol, descriptors=[])
        print('There should have been a ValueError!')
    except ValueError:
        pass

    # unknown descriptor
    try:
        result = filter.compute_descriptors(mol, descriptors=['ultimate_activity_prediction'])
        print('There should have been a KeyError!')
    except KeyError:
        pass
