"""
Module duplicate
===================
This modules is used to identify duplicate molecules within and accross multiple files.
"""

# standard
import logging
from pathlib import Path
import time
# data handling
import pandas as pd
from pandas import DataFrame
from tables.exceptions import HDF5ExtError
# chemoinformatics
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdinchi import MolToInchiKey
# docs
from typing import Union
from typing import Tuple
# dev
from npfc import save
# from npfc import save  # use of df.to_hdf with append and table modes for now


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def init_ref_file(ref_file: str, group_on: str, col_id: str) -> bool:
    """Initiate an empty reference hdf for identifying duplicates.

    :param ref_file: the reference file location
    :param group_on: the column name for grouping entries
    :param col_id: the column id for tracking kept entries
    :return: True if the reference file could be initialized, False otherwise.
    """
    try:
        # delete file if it already exists
        if Path(ref_file).is_file():
            Path.unlink(ref_file)
        # init
        key = Path(ref_file).stem
        # create a new ref file
        df_ref = pd.DataFrame({group_on: [], col_id: []})
        df_ref.index = df_ref[group_on]
        df_ref.to_hdf(ref_file, key=key, format="table")
        logging.debug("Created new ref_file at '%s'", ref_file)
        return True
    except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
        logging.critical("Could not create a new ref_file at '%s'", ref_file)
        return False


def filter_duplicates(df: DataFrame,
                      group_on: str = "inchikey",
                      col_id: str = "idm",
                      col_mol: str = "mol",
                      ref_file: str = None,
                      get_df_dupl: bool = False,
                      max_attempts: int = 5) -> Union[DataFrame, Tuple]:
    """Filter out duplicate molecules from a DataFrame.

    The comparison is performed by grouping entries on a column (inchikey or smiles) and the first entry of a group is kept.
    If no reference file is defined, then duplicates are removed independently for each input chunk.
    If a reference file is defined, then all passing entries of a chunk at a time are recorded in it so that further chunks molecules
    can be grouped with all previously recorded molecules.

    .. note:: the more chunks there are, the larger the reference file grows. For ChEMBL, I got something around 16Go (1.8M molecules) and for ZINC >100Go (14M molecules). To improve IO performances, I switched to append mode and table format, so the file does not to be rewritten entirely for every chunk.

    :param df: the DataFrame with the molecules to filter for duplicates
    :param col_mol: the column name for the molecules
    :param ref_file: the reference file location
    :param group_on: the column name for grouping entries
    :param col_id: the column id for tracking kept entries
    :param get_df_dupl: return the filered duplicate entries alongside the kept molecules, as another DataFrame
    :param max_attempts: in case two instances of this function runs very close in time to produce a new ref file, one will fail. This parameter gives the failing instance another attempt(s) every second until success or maximum is reached.
    :return: either the DataFrame with kept molecules only, or a tuple of DataFrames for kept and filtered molecules.
    """
    # avoid pandas warnings
    df = df.copy()  # ### this certainly is the worst way of doing this!

    # check on col_id
    if col_id not in df.columns:
        raise ValueError(f"Error! No column {col_id} found for identifying molecules.")

    # check on col_mol in case group_on is not present (needed for computing it!)
    if group_on not in df.columns and col_mol not in df.columns:
        raise ValueError(f"Error! No column {group_on} or {col_mol} found for grouping molecules.")

    # check on group_on
    if group_on not in ('inchikey', 'smiles'):
        raise ValueError(f"Error! Unauthorized value for on parameter ({group_on}).")
    # compute it if not present
    elif group_on not in df.columns:
        logging.warning("Column '%s' not found, so computing it.", group_on)
        if group_on == 'inchikey':
            df.loc[:, group_on] = df.loc[:, col_mol].map(MolToInchiKey)
        elif group_on == 'smiles':
            df.loc[:, group_on] = df.loc[:, col_mol].map(MolToSmiles)
        else:
            raise ValueError(f"Error! Expected values for group_on are: 'inchikey' or 'smiles' but got '{group_on}' instead!")

    # preprocess df into df_u (unique)
    df.set_index(group_on, inplace=True)
    df_u = df.loc[~df.index.duplicated(keep="first")]

    # filter duplicates found in the same input file
    dupl_inchikey = []
    dupl_idm_kept = []
    dupl_idm_filtered = []
    if len(df_u.index) == len(df.index):
        logging.debug("Number of duplicate molecule found in current chunk: 0")
    else:
        df_dupl = df[~df[col_id].isin(df_u[col_id])]
        logging.debug("Number of duplicate molecule found in current chunk: %s", len(df_dupl.index))
        for i in range(len(df_dupl)):
            row_dupl = df_dupl.iloc[i]
            dupl_inchikey.append(row_dupl.name)
            dupl_idm_kept.append(df_u.loc[row_dupl.name][col_id])  # this is why I have a loop
            dupl_idm_filtered.append(row_dupl[col_id])
    # df with entries that were removed because of another entry in the same chunk
    df_filtered = DataFrame({'group_on': dupl_inchikey, 'idm_kept': dupl_idm_kept, 'idm_filtered': dupl_idm_filtered})

    # load reference file
    if ref_file is not None:
        key = Path(ref_file).stem
        # init ref file if does not exist already
        num_attempts = 0
        ref_file_already_exists = False
        ref_file_could_be_created = False
        while num_attempts < max_attempts:
            try:
                if not Path(ref_file).is_file():
                    init_ref_file(ref_file, group_on=group_on, col_id=col_id)
                    ref_file_could_be_created = True
                    break
                else:
                    ref_file_already_exists = True
                    break
            except HDF5ExtError:
                logging.error("Error! Init File was not found to be existing but could not be created either. Retrying in 1s...")
                time.sleep(1)
            finally:
                if not ref_file_already_exists:
                    logging.debug("Reference file was already present")
                else:
                    if ref_file_could_be_created:
                        logging.debug("Reference file was created after %s attempt(s)!", num_attempts)
                    else:
                        logging.debug("Reference file could not be created, even after %s attempt(s)!", num_attempts)

            num_attempts += 1

        # open it with a lock as we'll need to update it at the end
        with save.SafeHDF5Store(ref_file) as store:
            try:
                df_ref = store[key]
            except KeyError:
                df_ref = pd.DataFrame({group_on: [], col_id: []})
            # use group_on as index for faster comparison
            df_ref.set_index(group_on, inplace=True)
            logging.debug("Number of entries in ref: %s", len(df_ref.index))

            # filter out already referenced compounds
            dupl_inchikey = []
            dupl_idm_kept = []
            dupl_idm_filtered = []
            df_u2 = df_u[~df_u.index.isin(df_ref.index)]
            # reset indices (feather does not support strings as rowids...)
            df_u2.reset_index(inplace=True)
            df_u2.drop([c for c in df_u.columns if c not in (group_on, col_id)], axis=1).to_hdf(ref_file, key=key, mode="a", format="table", append=True)
            if len(df_u2.index) == len(df_u.index):
                logging.debug("Number of duplicate molecules found by using ref_file: 0")
            else:
                df_dupl = df_u[~df_u[col_id].isin(df_u2[col_id])]
                logging.debug("Number of duplicate molecule found by using ref_file: %s", len(df_dupl.index))
                for i in range(len(df_dupl)):
                    row_dupl = df_dupl.iloc[i]
                    dupl_inchikey.append(row_dupl.name)
                    dupl_idm_kept.append(df_ref.loc[row_dupl.name][col_id])  # this is why I have a loop
                    dupl_idm_filtered.append(row_dupl[col_id])

            # return updated output in case of ref file
            df_u = df_u2
            df_ref.reset_index(inplace=True)
            df_filtered = pd.concat([df_filtered, DataFrame({'group_on': dupl_inchikey, 'idm_kept': dupl_idm_kept, 'idm_filtered': dupl_idm_filtered})])

    if get_df_dupl:
        return (df_u, df_filtered)

    return df_u
