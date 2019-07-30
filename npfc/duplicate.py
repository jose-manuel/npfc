"""
Module duplicate
===================
This modules is used to identify duplicate molecules within and accross multiple files.
"""

# standard
import logging
from pathlib import Path
# data handling
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdinchi import MolToInchiKey
# dev
from npfc import save
# from npfc import save  # use of df.to_hdf with append and table modes for now


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def init_ref_file(ref_file, group_on, col_id) -> bool:
    """Initiate an empty reference hdf for identifying duplicates.

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

    # in case of debug only, record the list of all duplicates from the chunk
    if logging.getLogger().level == logging.DEBUG:
        if len(df_u.index) == len(df.index):
            logging.debug("number of duplicate molecule found in current chunk: 0")
        else:
            df_dupl = df[~df[col_id].isin(df_u[col_id])]
            logging.debug(f"number of duplicate molecule found in current chunk: {len(df_dupl.index)}")
            # log what molecule "lost" to what other in same DataFrame: InChiKey, kept, filtered
            logging.debug(f"HEADER:group_on|id_kept|id_filtered")
            for i in range(len(df_dupl)):
                row_dupl = df_dupl.iloc[i]
                logging.debug(f"RESULT: {row_dupl.name}|{df_u.loc[row_dupl.name][col_id]}|{row_dupl[col_id]}")

    # load reference file
    if ref_file is not None:
        key = Path(ref_file).stem
        # init ref file if does not exist already
        if not Path(ref_file).is_file():
            init_ref_file(ref_file, group_on=group_on, col_id=col_id)

        # open it with a lock as we'll need to update it at the end
        with save.SafeHDF5Store(ref_file) as store:
            try:
                df_ref = store[key]
            except KeyError:
                df_ref = pd.DataFrame({group_on: [], col_id: []})
            # use group_on as index for faster comparison
            df_ref.set_index(group_on, inplace=True)

            # filter out already referenced compounds
            df_u2 = df_u[~df_u.index.isin(df_ref.index)]

            # reset indices (feather does not support strings as rowids...)
            df_u2.reset_index(inplace=True)

            df_u2.drop([c for c in df_u.columns if c not in (group_on, col_id)], axis=1).to_hdf(ref_file, key=key, mode="a", format="table", append=True)

            # in case of debug only, record the list of all duplicates by using the ref file
            if logging.getLogger().level == logging.DEBUG:
                if len(df_u2.index) == len(df_u.index):
                    logging.debug("Number of duplicate molecule found by using ref_file: 0")
                else:
                    df_dupl = df_u[~df_u[col_id].isin(df_u2[col_id])]
                    logging.debug(f"Number of duplicate molecule found by using ref_file: {len(df_dupl.index)}")
                    logging.debug(f"HEADER:group_on|id_kept|id_filtered")
                    # log what molecule "lost" to what other in same DataFrame: InChiKey, kept, filtered
                    for i in range(len(df_dupl)):
                        row_dupl = df_dupl.iloc[i]
                        logging.debug(f"RESULT:{row_dupl.name}|{df_ref.loc[row_dupl.name][col_id]}|{row_dupl[col_id]}")

            # return updated output in case of ref file
            df_u = df_u2
            df_ref.reset_index(inplace=True)

    return df_u
