"""
Module fragment
====================
This modules contains the function to run substructure searches.
"""


# standard
import logging
# data handling
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import Mol
# docs
from typing import List


def find(df_mols: DataFrame,
        df_frags: DataFrame,
        col_mol_mols: str = 'mol',
        col_mol_frags: str = 'mol',
        ) -> DataFrame:
    """Create a DataFrame recording every substructure (fragment) match in the
    input molecule DataFrame.

    A record is composed of 6 fields:

    1) idm: the id of the molecule (rowid from df_mols)
    2) idf: the id of the fragment (rowid from df_frags)
    3) aidxf: the atom indices of the fragment found in the molecule
    4) mol_perc: the percentage of the molecule the fragment represents (based on hac)
    5) mol: the molecule
    6) mol_frag: the fragment

    :param df_mols: the input DataFrame with molecules
    :param df_frags: the input DataFrame with fragments to use for substructure search
    :param col_mol_frags: the input DataFrame column name with the fragments
    :return: the substructure matches as a DataFrame

    .. note:: Rowids are used for recording the ids of substructure hits: mols => idm, frags => idf
    """
    # init
    d = {}
    d['idm'] = []
    d['idf'] = []
    d['idxf'] = []
    d['_aidxf'] = []
    d['mol_perc'] = []  # proportion of the molecule the substructure represents
    d['mol'] = []  # encode the molecule here so we don't have to combine multiple files when trying to have a look at the results
    d['mol_frag'] = []  # strucutre of the fragment
    # begin
    for i in range(len(df_mols.index)):
        rowm = df_mols.iloc[i]
        mol = rowm[col_mol_mols]
        hac = mol.GetNumAtoms()
        for j in range(len(df_frags.index)):
            rowq = df_frags.iloc[j]
            # perform the substructure search on mol so the latest matching fragment does not get highlighted
            matches = mol.GetSubstructMatches(rowq[col_mol_frags])
            if len(matches) > 0:
                logging.debug(f"MOL {rowm.name} + FRAG {rowq.name} ==> {matches}")
                for i, m in enumerate(matches):
                    d['idm'].append(rowm.name)
                    d['idf'].append(rowq.name)
                    d['_aidxf'].append(frozenset(m))  # frozenset so we can use intersection, etc. and still remove dupl. easily
                    d['idxf'].append(str(i))
                    d['mol_perc'].append(round(len(m)/hac, 2) * 100)
                    d['mol'].append(rowm[col_mol_mols])
                    d['mol_frag'].append(rowq[col_mol_frags])

    return DataFrame(d)
