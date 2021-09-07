"""
Module fragment_combination_point
==================================
This modules contains the functions for annotating and parsing fragment
combination point labels in fragments.
"""


# standard
import logging
import re
import string
from collections import defaultdict
# chemoinformatics
from rdkit.Chem import Mol
# docs
from typing import Union
from typing import Tuple
from typing import List


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def idx_to_label(n: int) -> str:
    """Convert a number to a corresponding letter in the alphabet.
    In case the number is higher than the number of letters in the english alphabet, then
    a second character is appended.

    >>> For instance:
    >>>idx_to_label(0)
    >>> 'a'
    >>> idx_to_label(25)
    >>> 'z'
    >>> idx_to_label(26)
    >>> 'aa'

    This function was inspired after:
    https://stackoverflow.com/questions/2267362/how-to-convert-an-integer-to-a-string-in-any-base

    :param n: the input number
    :return: the corresponding string
    """
    alphabet_size = 26
    digits = []
    n += 1
    while n:
        digits.append(int(n % alphabet_size - 1))
        n //= alphabet_size

    digits.reverse()
    return ''.join([string.ascii_lowercase[i] for i in digits])


def find_symmetry_groups(mol: Mol) -> Tuple[Tuple[int]]:
    """Search for symmetry groups within a molecule.
    Symmetry groups are defined as atoms that can have equivalent
    positions in substructure matches.

    This function has been copied from
    https://sourceforge.net/p/rdkit/mailman/message/27897393/

    :param mol: the input molecule
    :return: a tuple of tuples of atom indices. Each tuple represents a symmetry group and inside are the atom indices.
    """
    equivs = defaultdict(set)
    matches = mol.GetSubstructMatches(mol, uniquify=False)
    for match in matches:
        for idx1, idx2 in enumerate(match):
            equivs[idx1].add(idx2)
    classes = set()
    for s in equivs.values():
        classes.add(tuple(s))
    return tuple(classes)


def get_fcp_labels(mol: Mol) -> dict:
    """Search for symmetry groups within the input molecule and defines
    accordingly a label each atom. This label is composed of 2 parts:

        - a prefix: sequential number that gets incremented for each group
        - a suffix: sequential letter or string that gets incremented for each element in group


    This means that equivalent atoms (within the same group) are given the same prefix but different suffixes.

    >>> For instance:
    >>> '1a', '1b', '2'  # 2 equivalent atoms in position 1
    >>> '1', '2', '3', '4', '5', '6'  # a 6-atoms fragment with no symmetry group
    >>> '1a', '1b', '1c', '1d', '1e', '1f'  # a benzene or cyclohexane

    .. note:: labels prefixes begin at 1, not 0 like atom indices.

    :param mol: the input molecule
    :return: a dictionary with the correspondance between atom indices (keys) and atom labels (values)
    """
    symmetry_groups = find_symmetry_groups(mol)
    logging.debug('Symmetry groups: %s', symmetry_groups)
    d_fcp_labels = {}
    for i, aidx_eq in enumerate(symmetry_groups):
        for j, aidx in enumerate(aidx_eq):
            if len(aidx_eq) > 1:
                suffix = idx_to_label(j)
            else:
                suffix = ''
            fcp_label = f"{i+1}{suffix}"
            d_fcp_labels[aidx] = fcp_label

    logging.debug('FCP labels: %s', d_fcp_labels)
    return d_fcp_labels


# def simplify_fcp(fcp: List[str]) -> List[str]:
#     """Remove suffixes from each fcp within provided list.
#
#     >>> For instance:
#     >>> '1a', '1b', '2' becomes: '1', '1', '2'
#
#     :param fcp: a list of Fragment Connection Points
#     :return: the simplified list
#     """
#     return [re.sub("[^0-9]", "", x) for x in fcp]


def clear_fcp_suffixes(fcp_str: str) -> List[str]:
    """Remove suffixes from each fcp within provided string.

    >>> For instance:
    >>> '1a,1b,2' becomes: '1,1,2'

    :param fcp: a list of Fragment Connection Points
    :return: the simplified list
    """
    return re.sub("[^0-9,]", "", fcp_str)


def clear_fcp_suffixes_in_edges(edges: list) -> list:
    """Remove suffixes from each fcp within the attributes found in a list of edges

    :param edges: a list of edges ([[u1, v1, d1], [u2, v2, d2]]), with u and v the nodes and d a dictionary of attributes
    :return: the updated edges
    """
    for e in edges:
        e[2]['fcp_1'] = clear_fcp_suffixes(e[2]['fcp_1'])
        e[2]['fcp_2'] = clear_fcp_suffixes(e[2]['fcp_2'])

    return edges


def count_symmetry_groups(val: Union[Mol, List[str]]) -> int:
    """Return the number of symmetry groups found in a fragment.

    :param val: either a fragment as Mol object or the fragment connection points as dict
    :return: the count of symmetry groups
    """
    # init
    if isinstance(val, Mol):
        val = get_fcp_labels(val)
    elif not isinstance(val, dict):
        raise ValueError(f"Error! Input value is neither a Molecule nor a dictionary: '{val}' ({type(val)})")
    fcp = list(val.values())
    # regroup by symmetry group
    d = {}
    # populate each groups with the equivalent atom indices
    for s in fcp:
        m = re.search(r"[a-z]", s)
        if m:
            idx = m.start()
            prefix = s[:idx]
            if prefix in d.keys():
                d[prefix] += 1
            else:
                d[prefix] = 1

    # count only groups with more than 1 element
    return len([x for x in d.values() if x > 1])
