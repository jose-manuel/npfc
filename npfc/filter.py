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
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


class Filter:
    """A class for filtering molecules based on molecular descriptors."""

    def __init__(self):
        """Create a Filter object with following descriptors:"""
        self.descriptors = {'hac': lambda x: x.GetNumAtoms(),
                            'molweight': lambda x: round(Descriptors.ExactMolWt(x), 4),
                            'nrings': lambda x: rdMolDescriptors.CalcNumRings(x),
                            'elements': lambda x: set([a.GetSymbol() for a in x.GetAtoms()]),
                            }

    def filter_mol(self, mol: Mol, expr: str) -> bool:
        f"""Filter a molecule based on an expression.
        Two types of expressions are currently supported:

            - inclusion/exclusions: 'elements not in C, N, O', 'elements in C, N, O'
            - numeric: 'hac > 3', '100.0 < molweight <= 1000.0', 'nrings' != 0, 'nrings == 0'

        List of currently supported descriptors:
        {', '.join([k for k in self.descriptors.keys()])}


        :param mol: the input molecule
        :param expr: the filter to apply
        :return: True if the molecule passes the filter, False otherwise
        """

        expr = expr.lower()
        split_expr = expr.lower().split()
        # filters of type: 'elements in C, N, O'
        if 'in' in split_expr:  # 'in' or 'not in'
            return self._eval_set_expr(mol, expr)
        # filters of type: 'hac > 3'
        return self._eval_numeric_expr(mol, expr)

    def _eval_numeric_expr(self, mol, expr):
        """
        Evaluate if the statements stored in the expression are True or False.
        For now statement is composed of either 3 elements (['molweiht', '<=', '1000'])
        or 5 elements: (['0', '<=', 'molweight', '<=', '1000']).
        ### No check has been added on this number because there might be an expanded functionality
        later on (combining statements with ';'?).
        Descriptors used for the comparisons need to be provided as a dictionary (name: value).

        Possible values for how: numeric, set or literal.
        """
        expr = expr.replace(" ", "")
        split_expr = self._split_expr(expr)  # something like 'molweight', '<=', '1000'
        # replace descriptor names by their values
        split_expr = [self.descriptors[k](mol) if k in self.descriptors.keys() else k for k in split_expr]  # now it is '250.0', '<=', '1000'
        logging.debug(f"applying numeric filter: {split_expr}")
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
        logging.debug(f"applying inclusion/exclusion filter: {' '.join([descriptor, op, values])}")
        if (op == ' in ' and descriptor.issubset(values)) or (op == ' not in ' and not descriptor.issubset(values)):
            return True
        else:
            return False

    def _split_expr(self, expr):
        """Helper function for _eval_expr.
        From a string containing an expression (i.e. 'molweight < 1000'), return
        a list of values and operators (['molweight', '<', '1000']).
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
