"""
Module fragment_search
==========================
This modules contains the function to run substructure searches.
"""


# standard
import logging
# data handling
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import rdTautomerQuery


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def get_fragment_hits(df_mols: DataFrame,
                      df_frags: DataFrame,
                      col_mol_mols: str = 'mol',
                      col_mol_frags: str = 'mol',
                      col_mol_inchikey: str = 'inchikey',
                      fcp_labels: str = None,
                      tautomer: bool = False,
                      col_to_index_mols: str = 'idm',
                      col_to_index_frags: str = 'idm',
                      ) -> DataFrame:
    """Create a DataFrame recording every Fragment Hit in the
    input molecule DataFrame.

    A Fragment Hit is composed of 6 fields:

    1) idm: the id of the molecule (rowid from df_mols)
    2) idf: the id of the fragment (rowid from df_frags)
    3) aidxf: the atom indices of the fragment found in the molecule
    4) mol_perc: the percentage of the molecule the fragment represents (based on hac)
    5) mol: the molecule as RDKit Mol object
    6) mol_frag: the fragment as RDKit Mol object

    :param df_mols: the input DataFrame with the molecules (df_mols)
    :param df_frags: the input DataFrame with fragments to use for substructure search (df_frags)
    :param col_mol_mols: the column name in df_mols with the molecules
    :param col_mol_frags: the column name in df_frags with the fragments
    :param col_mol_inchikey: the input DataFrame column name with the inchikey of the molecule. If it does not exist, then an empty column is created in the output.
    :param fcp_labels: the column name in the fragments dataframe with the fcp labels
    :param tautomer: if set to True, tautomers will be taken into account during fragment search (warning, tautomer-independant search is much slower!)
    :param col_to_index_mols: set the row indices of the DataFrame with the molecules to probe to the specified column. If empty (''), indices are left untouched.
    :param col_to_index_frags: set the row indices of the DataFrame with the fragments to seach for to the specified column. If empty (''), indices are left untouched.

    :return: the substructure matches as a DataFrame

    .. warning:: Row indices are used for recording the ids of substructure hits and are therefore required to be set to the molecule identifiers (i.e. idm).

    """
    # init
    d = {}
    d['idm'] = []
    d['idf'] = []
    d['_aidxf'] = []
    d['mol_perc'] = []  # proportion of the molecule the substructure represents
    d['mol'] = []  # encode the molecule here so we don't have to combine multiple files when trying to have a look at the results
    d['mol_frag'] = []  # strucutre of the fragment
    d['inchikey'] = []  # inchikey of the molecule

    # row ids are used for storing the correct molecule and fragment ids in the
    #  output table, so they need to be set before the computation occurs
    if col_to_index_mols != '':
        df_mols.index = list(df_mols[col_to_index_mols])
    if col_to_index_frags != '':
        df_frags.index = list(df_frags[col_to_index_frags])

    # if no inchikey was computed, generate an empty column
    if col_mol_inchikey not in df_mols.columns:
        df_mols['inchikey'] = ''

    # tautomers
    if tautomer:
        df_mols[col_mol_mols + '_taut'] = df_mols[col_mol_mols].map(lambda x: rdTautomerQuery.TautomerQuery(x).GetTemplateMolecule())

    # begin
    for i in range(len(df_mols.index)):
        rowm = df_mols.iloc[i]
        if tautomer:
            mol = rowm[col_mol_mols + '_taut']
        else:
            mol = rowm[col_mol_mols]
        hac = mol.GetNumAtoms()
        for j in range(len(df_frags.index)):
            rowq = df_frags.iloc[j]
            # perform the substructure search on mol so the latest matching fragment does not get highlighted
            matches = mol.GetSubstructMatches(rowq[col_mol_frags])
            if len(matches) > 0:
                logging.debug("MOL %s + FRAG %s ==> %s", rowm.name, rowq.name, matches)
                for m in matches:
                    d['idm'].append(rowm.name)
                    d['idf'].append(rowq.name)
                    d['_aidxf'].append(m)  # frozenset so we can use intersection, etc. and still remove dupl. easily
                    d['mol_perc'].append(round(len(m)/hac, 2) * 100)
                    d['mol'].append(rowm[col_mol_mols])
                    d['inchikey'].append(rowm[col_mol_inchikey])
                    d['mol_frag'].append(rowq[col_mol_frags])

    df_fs = DataFrame(d)
    df_fs['idf_idx'] = df_fs.groupby(['idm', 'idf']).cumcount()  # rank seems to be working with np.float types only...

    # add fcp labels
    if fcp_labels is not None and fcp_labels != '':
        df_frags = df_frags.rename({'idm': 'idf'}, axis=1)
        df_fs = df_fs.merge(df_frags[['idf', fcp_labels]], on='idf', how='left')

    return df_fs.rename({fcp_labels: '_fcp_labels'}, axis=1)
