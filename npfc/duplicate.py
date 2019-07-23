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
# from npfc import save  # use of df.to_hdf with append and table modes for now


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
        p = Path(ref_file)
        if p.is_file():
            p.unlink()
        # init
        # create a new ref file
        df_ref = pd.DataFrame({group_on: [], col_id: []})
        # df_ref.set_index(group_on, inplace=True)
        df_ref.to_hdf(ref_file, key=p.stem, table=True)
        logging.debug(f"Created new ref_file at '{ref_file}'")
        return True
    except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
        logging.critical(f"Could not create a new ref_file at '{ref_file}'")
        return False


def filter_duplicates(df: DataFrame, group_on: str = "inchikey", col_id: str = "idm", col_mol: str = "mol", ref_file: str = None):
    """Filter out duplicate molecules from a DataFrame.

    The comparison is performed by grouping entries on a column (inchikey or smiles) and the first entry of a group is kept.
    If no reference file is defined, then duplicates are removed independently for each input chunk.
    If a reference file is defined, then all passing entries of a chunk at a time are recorded in it so that further chunks molecules
    can be grouped with all previously recorded molecules.

    .. note:: the more chunks there are, the larger the reference file grows. For ChEMBL, I got something around 16Go (1.8M molecules) and for ZINC >100Go (14M molecules). To improve IO performances, I switched to append mode and table format, so the file does not to be rewritten entirely for every chunk.

    .. warning:: the feather format looks very promising for very fast IO times, but it does not allow for anything below simple variable types or customized row indices. For now I don't use it because there would be much testing to do. This could be addressed in a further release of npfc.

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
            p = Path(ref_file)
            if not p.exists():
                init_ref_file(ref_file, group_on, col_id)
                logging.debug(f"Initialized new reference file at '{ref_file}'")
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

            # df_ref.index = range(len(df.index))  # error
            # update reference file
            df_ref.to_hdf(ref_file, key=p.stem, table=True, mode='a')

        # once lock is lifted, we can delete the lock file
        Path(lock_file).unlink()

    # This function has a lot of changes made to the index of the dataframes. This is because the feather format,
    # although super fast to read/write and memory efficient, cannot have strings as row indices.
    # Also, (un)setting the index is a fast operation, so it should not be a problem even for larger DataFrames.
    # The "only" real problem with this function is that it needs to rewrite the full reference file for each chunk.

    return df_u
