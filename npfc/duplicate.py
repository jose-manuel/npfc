"""
Module duplicate
===================
This modules is used to identify duplicate molecules within and accross multiple files.
"""

# standard
import logging
from time import sleep
from pathlib import Path
from filelock import FileLock
# data handling
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdinchi import MolToInchiKey
# dev
from npfc import load
from npfc import save


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def init_ref_file(ref_file, group_on, col_id):
    """Initiate an empty reference file for identifying duplicates.

    :param ref_file: name of the reference file. Format will be deduced from extension.
    :param group_on: property to use for indexing data
    :param col_id: property to use for labelling data
    :return: True if the reference file could be initialized, False otherwise.
    """
    try:
        # delete file if it already exists
        if Path(ref_file).is_file():
            Path.unlink(ref_file)
        # init
        # create a new ref file
        df_ref = pd.DataFrame({group_on: [], col_id: []})
        # df_ref.set_index(group_on, inplace=True)
        save.file(df_ref, ref_file)
        logging.debug(f"Created new ref_file at '{ref_file}'")
        return True
    except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
        logging.critical(f"Could not create a new ref_file at '{ref_file}'")
        return False


def filter_duplicates(df: DataFrame, group_on: str = "inchikey", col_id: str = "idm", col_mol: str = "mol", ref_file: str = None):

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
        logging.warning(f"Column {group_on} not found, so computing it.")
        if group_on == 'inchikey':
            df.loc[:, group_on] = df.loc[:, col_mol].map(MolToInchiKey)
        elif group_on == 'smiles':
            df.loc[:, group_on] = df.loc[:, col_mol].map(MolToSmiles)
        else:
            raise ValueError(f"Error! Expected values for group_on are: 'inchikey' or 'smiles' but got '{group_on}' instead!")

    # preprocess df into df_u (unique)
    df.set_index(group_on, inplace=True)
    df_u = df.loc[~df.index.duplicated(keep="first")]

    # load reference file
    if ref_file:

        # define a lock file
        lock_file = ref_file + ".lock"

        # if the lock file exists, wait for it to be removed by its current job
        while Path(lock_file).exists():
            logging.debug(f"Waiting for lock file to be lifted at '{lock_file}'")
            sleep(1)  # wait 1s before trying to acces the file again

        # lock is not here, we can create one
        lock = FileLock(lock_file)
        logging.debug(f"Set lock file at '{lock_file}'. Reference file will not be accessible for read/write until this lock is lifted.")

        # work with lock on, it will be automatically lifted when work is done
        with lock:

            # load references
            if not Path(ref_file).exists():
                init_ref_file(ref_file, group_on, col_id)
            df_ref = load.file(ref_file, out_mol=None)

            # use group_on as index for faster comparison
            df_ref.set_index(group_on, inplace=True)

            # filter out already referenced compounds
            df_u = df_u[~df_u.index.isin(df_ref.index)]

            # reset indices (feather does not support strings as rowids...)
            df_ref.reset_index(inplace=True)
            df_u.reset_index(inplace=True)

            # record new references
            df_ref = pd.concat([df_ref, df_u.drop([c for c in df_u.columns if c not in (group_on, col_id)], axis=1)], axis=0, join='inner')  # no need to align both dataframes as they should always have the columns in the same order

            # overwrite reference file  ### might be worth to look into using pandas with table mode for appending new data, for now feathers works just fine
            save.file(df_ref, ref_file)

        # once lock is lifted, we can delete the lock file
        Path(lock_file).unlink()

    # This function has a lot of changes made to the index of the dataframes. This is because the feather format,
    # although super fast to read/write and memory efficient, cannot have strings as row indices.
    # Also, (un)setting the index is a fast operation, so it should not be a problem even for larger DataFrames.
    # The "only" real problem with this function is that it needs to rewrite the full reference file for each chunk.

    return df_u
