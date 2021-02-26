"""
Module filter
==============
This modules contains the class Filter, which is used to filter molecules using
molecular descriptors.
"""

# data handling
import logging
import re
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
# docs
from typing import List


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def count_violations_lipinski(molecular_weight, slogp, num_hbd, num_hba):
    """Lipinski, J Pharmacol Toxicol Methods. 2000 Jul-Aug;44(1):235-49.
    """
    n = 0
    if molecular_weight < 150 or molecular_weight > 500:
        n += 1
    if slogp > 5:
        n += 1
    if num_hbd > 5:
        n += 1
    if num_hba > 10:
        n += 1
    return n


def count_violations_veber(num_rotatable_bonds, tpsa):
    """Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD (June 2002).
    "Molecular properties that influence the oral bioavailability of drug candidates".
    J. Med. Chem. 45 (12): 2615–23.
    """
    n = 0
    if num_rotatable_bonds > 10:
        n += 1
    if tpsa > 140:
        n += 1
    return n


def count_violations_lead_like(molecular_weight, slogp, num_rotatable_bonds):
    """http://zinc.docking.org/browse/subsets/
    Teague, Davis, Leeson, Oprea, Angew Chem Int Ed Engl. 1999 Dec 16;38(24):3743-3748.

    """
    n = 0
    if molecular_weight < 250 or molecular_weight > 350:
        n += 1
    if slogp > 3.5:
        n += 1
    if num_rotatable_bonds > 7:
        n += 1
    return n


def count_violations_ppi_like(molecular_weight, slogp, num_hba, num_rings):
    """Hamon, V., Bourgeas, R., Ducrot, P., Theret, I., Xuereb, L., Basse, M.J., Brunel, J.M., Combes, S., Morelli, X., Roche, P., 2013.
    2P2IHUNTER: a tool for filtering orthosteric protein–protein interaction modulators via a dedicated support vector machine.
    Journal of The Royal Society Interface 11. doi:10.1098/rsif.2013.0860
    """
    n = 0
    if molecular_weight > 400:
        n += 1
    if slogp < 4:
        n += 1
    if num_hba < 4:
        n += 1
    if num_rings < 4:
        n += 1
    return n


def count_violations_fragment_like(molecular_weight, slogp, num_hba, num_hbd):
    """Congreve, M., Carr, R., Murray, C., Jhoti, H., 2003. A “Rule of Three” for fragment-based lead discovery?
    Drug Discovery Today 8, 876–877. doi:10.1016/S1359-6446(03)02831-9
    """
    n = 0
    if molecular_weight >= 300:
        n += 1
    if slogp > 3:
        n += 1
    if num_hba > 3:
        n += 1
    if num_hbd > 3:
        n += 1
    return n


def count_violations_fragment_like_ext(num_fragment_like_violations, tpsa, num_rotatable_bonds):
    """Congreve, M., Carr, R., Murray, C., Jhoti, H., 2003.
    A “Rule of Three” for fragment-based lead discovery? Drug Discovery Today 8, 876–877.
    doi:10.1016/S1359-6446(03)02831-9
    """
    n = num_fragment_like_violations
    if tpsa > 60:
        n += 1
    if num_rotatable_bonds:
        n += 1
    return n


def get_min_max_ring_sizes(mol):
    """Return a tuple wih (minimum, maximum) ring sizes of the input molecule.
    In case the molecule is linear, (0, 0) is returned.
    """
    ring_sizes = [len(x) for x in mol.GetRingInfo().AtomRings()]
    if len(ring_sizes) > 0:
        min_ring_size = min(ring_sizes)
        max_ring_size = max(ring_sizes)
    else:
        min_ring_size = 0
        max_ring_size = 0

    return (min_ring_size, max_ring_size)


DESCRIPTORS = {
              # classical molecular descriptors
              'num_heavy_atoms': lambda x: x.GetNumAtoms(),
              'molecular_weight': lambda x: round(Descriptors.ExactMolWt(x), 4),
              'num_rings': lambda x: rdMolDescriptors.CalcNumRings(x),
              'num_rings_arom': lambda x: rdMolDescriptors.CalcNumAromaticRings(x),
              'elements': lambda x: set([a.GetSymbol() for a in x.GetAtoms()]),
              'molecular_formula': lambda x: rdMolDescriptors.CalcMolFormula(x),
              'num_hbd': lambda x: rdMolDescriptors.CalcNumLipinskiHBD(x),
              'num_hba': lambda x: rdMolDescriptors.CalcNumLipinskiHBA(x),
              'slogp': lambda x: round(Crippen.MolLogP(x), 4),
              'tpsa': lambda x: round(rdMolDescriptors.CalcTPSA(x), 4),
              'num_rotatable_bonds': lambda x: rdMolDescriptors.CalcNumRotatableBonds(x),
              'num_atoms_oxygen': lambda x: len([a for a in x.GetAtoms() if a.GetAtomicNum() == 8]),
              'num_atoms_nitrogen': lambda x: len([a for a in x.GetAtoms() if a.GetAtomicNum() == 7]),
              # custom molecular descriptors
              # ring_sizes:
              # it would have been faster to access only once RingInfo for both min and max,
              # but this is tricky because I would have to start making exceptions in the way
              # the functions are accessed or more complicated downstream process.
              # Indeed, I do not think it is possible to set two dict keys at once
              # within a dict comprehension and it is not that bad for
              # performance to call it twice anyway.
              'ring_size_min': lambda x: min([len(y) for y in x.GetRingInfo().AtomRings()]),
              'ring_size_max': lambda x: max([len(y) for y in x.GetRingInfo().AtomRings()]),
              }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Filter:
    """A class for filtering molecules based on molecular descriptors."""

    def __init__(self, descriptors: list = DESCRIPTORS):
        """Create a Filter object."""

        self.descriptors = descriptors

    def compute_descriptors(self, mol: Mol, descriptors: List = None) -> dict:
        """Compute descriptors. A subset of descriptors can be computed if a list
        of descriptor names is provided. To get an idea of what descriptors can be computed,
        the method get_possible_descriptors can be used.

        :param mol: the input molecule
        :param descriptors: the list of descriptors to compute. If none is provided, all possible descriptors are computed.
        :return: a dictionary with all descriptors
        """
        # if no descriptor is specified, compute them all
        if descriptors is None:
            descriptors = list(self.descriptors.keys())

        if len(descriptors) == 0:
            raise ValueError('Error! No descriptor is specified for computation!')

        return {descriptors[i]: self.descriptors[descriptors[i]](mol) for i in range(len(descriptors))}

    def get_possible_descriptors(self) -> List:
        """Return a list of all descriptors that can be computed using this module.

        :return: the list of descriptors that can be computed
        """
        return sorted(list((self.descriptors.keys())))

    def filter_mol(self, mol: Mol, expr: str) -> bool:
        """Filter a molecule based on an expression.
        Two types of expressions are currently supported:

            - inclusion/exclusion
                - 'elements not in C, N, O'
                - 'elements in C, N, O'

            - numeric
                - 'num_heavy_atoms > 3'
                - '100.0 < molecular_weight <= 1000.0'
                - 'num_rings' != 0'
                - 'num_rings == 0'

        :param mol: the input molecule
        :param expr: the filter to apply
        :return: True if the molecule passes the filter, False otherwise
        """
        # init
        split_expr = [s.lower() for s in expr.lower().split()]
        # filters of type: 'elements in C, N, O'
        if 'in' in split_expr:  # 'in' or 'not in'
            return self._eval_set_expr(mol, expr)
        # filters of type: 'num_heavy_atoms > 3'
        return self._eval_numeric_expr(mol, expr.lower())

    def _eval_numeric_expr(self, mol, expr):
        """
        Evaluate if the statements stored in the expression are True or False.
        For now statement is composed of either 3 elements (['molweiht', '<=', '1000'])
        or 5 elements: (['0', '<=', 'molecular_weight', '<=', '1000']).
        ### No check has been added on this number because there might be an expanded functionality
        later on (combining statements with ';'?).
        Descriptors used for the comparisons need to be provided as a dictionary (name: value).

        Possible values for how: numeric, set or literal.
        """
        mol = Mol(mol)
        expr = expr.replace(" ", "")
        split_expr = self._split_expr(expr)  # something like 'molecular_weight', '<=', '1000'
        # replace descriptor names by their values
        split_expr = [self.descriptors[k](mol) if k in self.descriptors.keys() else k for k in split_expr]  # now it is '250.0', '<=', '1000'
        logging.debug("Applying numeric filter: %s", ' '.join(str(v) for v in split_expr))
        # convert all values extracted as string into their type
        split_expr = [float(x) if x not in split_expr[1::2] else x for x in split_expr]  # and now it is 250.0, '<=', 1000.0
        # operators are always at odd positions, whereas values are at even positions
        # and there is always a value on the left and on the right of an operator
        for i in range(1, len(split_expr), 2):
            operator = split_expr[i]
            left = split_expr[i-1]
            right = float(split_expr[i+1])
            if operator == "<=":
                if not left <= right:
                    return False
            elif operator == "<":
                if not left < right:
                    return False
            elif operator == "==":
                if not left == right:
                    return False
            elif operator == "!=":
                if not left != right:
                    return False
            elif operator == ">=":
                if not left >= right:
                    return False
            elif operator == ">":
                if not left > right:
                    return False
        return True

    def _eval_set_expr(self, mol, expr):
        """Helper function for _eval_expr.
        Look for keywords ' in ' and ' not in ' in expression and check the condition
        by defining left as the descriptor and right as the values, i.e.:
        descriptor in values ('elements in H, C, N, O')
        """
        for op in [' not in ', ' in ']:
            pattern = re.compile(op)  # raw string
            hits = [(m.start(0), m.end(0)) for m in re.finditer(pattern, expr)]
            if len(hits) > 0:
                break  # leave asap with op still set to the correct operator
        # in case we did not find anything, just stop
        if len(hits) == 0:
            raise ValueError(f"expected ' not in ' or ' in ' in expr ({expr})")
        expr_split = [e.replace(" ", "") for e in expr.split(op)]
        descriptor = self.descriptors[expr_split[0]](mol)  # left
        values = set(expr_split[1].split(","))  # right
        logging.debug("Applying inclusion/exclusion filter: %s", ''.join(str(v) for v in [descriptor, op, values]))
        if (op == ' in ' and descriptor.issubset(values)) or (op == ' not in ' and not descriptor.issubset(values)):
            return True
        else:
            return False

    def _split_expr(self, expr):
        """Helper function for _eval_expr.
        From a string containing an expression (i.e. 'molecular_weight < 1000'), return
        a list of values and operators (['molecular_weight', '<', '1000']).
        """
        opidx_eq = self._find_opidx("==", expr)
        opidx_diff = self._find_opidx("!=", expr)
        opidx_supeq = self._find_opidx(">=", expr)
        opidx_infeq = self._find_opidx("<=", expr)
        opidx_sup = self._find_opidx(">", expr)
        opidx_inf = self._find_opidx("<", expr)

        # filter sup and inf with supeq and infeq
        opidx_sup = self._filter_wrong_matches(opidx_supeq, opidx_sup)
        opidx_inf = self._filter_wrong_matches(opidx_infeq, opidx_inf)
        # split expr into values and operators
        # sorted operators so we can iterate over the expr from left to right
        opidx_all = sorted(opidx_eq + opidx_diff + opidx_supeq + opidx_infeq + opidx_sup + opidx_inf, key=lambda x: x[0])
        split_expr = []
        split_expr.append(expr[:opidx_all[0][0]])
        for i in range(len(opidx_all) - 1):
            # always take on the value on the right side of the op, so init the first part outside of the loop
            opidx_curr = opidx_all[i]
            opidx_next = opidx_all[i+1]
            operator = expr[opidx_curr[0]:opidx_curr[1]]
            split_expr.append(operator)
            value = expr[opidx_curr[1]:opidx_next[0]]
            split_expr.append(value)
        split_expr.append(expr[opidx_all[-1][0]:opidx_all[-1][1]])
        split_expr.append(expr[opidx_all[-1][1]:])
        return split_expr

    def _find_opidx(self, op, expr):
        """ Helper function for _split_expr.
        Return all occurrences indices of a comparison operator (op) within an expr.
        """
        # init possible operator symbols
        pattern = re.compile(op)  # raw string
        return [(m.start(0), m.end(0)) for m in re.finditer(pattern, expr)]

    def _filter_wrong_matches(self, opidx_larger, opidx_smaller):
        """Helper function for __split_expr.
        Filter out false positives of comparison operators. For instance,
        '<' beginning at the same position as '<=' should be discarded.
        """
        invalid = []
        for smaller in opidx_smaller:
            for larger in opidx_larger:
                if smaller[0] == larger[0]:
                    invalid.append(smaller)
        return [smaller for smaller in opidx_smaller if smaller not in invalid]
