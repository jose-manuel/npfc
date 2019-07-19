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
from itertools import chain
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdinchi import MolToInchiKey
from rdkit.Chem import rdinchi
# dev
from npfc import utils
from npfc import load
from npfc import save


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# def init_ref_file(ref_file, group_on, col_id):
#     """Initiate an empty reference hdf for identifying duplicates.
#
#     :param ref_file: name of the reference file. Format will be deduced from extension.
#     :param group_on: property to use for indexing data
#     :param col_id: property to use for labelling data
#     :return: True if the reference file could be initialized, False otherwise.
#     """
#     try:
#         # delete file if it already exists
#         if Path(ref_file).is_file():
#             Path.unlink(ref_file)
#         # init
#         # create a new ref file
#         df_ref = pd.DataFrame({group_on: [], col_id: []})
#         df_ref.set_index(group_on, inplace=True)
#         save.file(df_ref)
#         logging.debug(f"Created new ref_file at '{ref_file}'")
#         return True
#     except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
#         logging.critical(f"Could not create a new ref_file at '{ref_file}'")
#         return False
#
#
# def filter_duplicates(df: DataFrame, group_on: str = "inchikey", col_id: str = "idm", col_mol: str = "mol", ref_file: str = None):
#
#     # avoid pandas warnings
#     df = df.copy()  # ### this certainly is the worst way of doing this!
#
#     # check on col_id
#     if col_id not in df.columns:
#         raise ValueError(f"Error! No column {col_id} found for identifying molecules.")
#
#     # check on col_mol in case group_on is not present (needed for computing it!)
#     if group_on not in df.columns and col_mol not in df.columns:
#         raise ValueError(f"Error! No column {group_on} or {col_mol} found for grouping molecules.")
#
#     # check on group_on
#     if group_on not in ('inchikey', 'smiles'):
#         raise ValueError(f"Error! Unauthorized value for on parameter ({group_on}).")
#     # compute it if not present
#     elif group_on not in df.columns:
#         logging.warning(f"Column {group_on} not found, so computing it.")
#         if group_on == 'inchikey':
#             df.loc[:, group_on] = df.loc[:, col_mol].map(MolToInchiKey)
#         else:
#             df.loc[:, group_on] = df.loc[:, col_mol].map(MolToSmiles)
#
#     # get df ready
#     df.set_index(group_on)
#
#     # load reference file
#     if ref_file:
#
#         # define a lock file
#         lock_file = ref_file + ".lock"
#
#         # if the lock file exists, wait for it to be removed by its current job
#         while Path(lock_file).exists():
#             sleep(1)
#
#         # lock is not here, we can create one
#         lock = FileLock(lock_file)
#
#         # work with lock on, it will be automatically lifted when work is done
#         with lock:
#
#             # load references
#             df_ref = load.file(ref_file)
#             df_ref.set_index(group_on)
#
#             # identify duplicate compounds





class DuplicateFilter:
    """A class used for filtering duplicate molecules in on or several DataFrame(s)."""

    def __init__(self, on: str = 'inchikey',
                 ref_file: str = None,
                 col_mol: str = 'mol',
                 col_id: str = 'idm'):
        """Create an instance of DuplicateFilter with following parameters:

        :param on: The property to use for identifying duplicate molecules.
        :param ref_file: The path to the reference file used to store synonyms.
        :param col_mol: The DataFrame column name where molecules are stored.
        :param col_id: The DataFrame column name where molecule identifiers are stored.

        .. note:: For now, only the 'inchikey' property works for the 'on' parameter. Moreover InchIKey are computed wether provided or not.

        """
        self._on = on
        self._col_mol = col_mol
        self._col_id = col_id
        self._col_id_synonyms = self.col_id + "_synonyms"
        self._ref_file = ref_file
        logging.debug(f"Initialized a new DuplicateFilter object")

    @property
    def col_id(self) -> str:
        return self._col_id

    @col_id.setter
    def col_id(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_id cannot be '{value}'.")
        self._col_id = value

    @property
    def col_mol(self) -> str:
        return self._col_mol

    @col_mol.setter
    def col_mol(self, value: str) -> None:
        if value is None:
            raise ValueError(f"Error! col_mol cannot be '{value}'.")
        self._col_mol = value

    @property
    def ref_file(self) -> str:
        return self._ref_file

    @ref_file.setter
    def ref_file(self, value: str) -> None:
        if value is not None and not isinstance(value, str):
            raise ValueError(f"Error! Either None or a str are expected for ref_file, not '{value}' ({type(value)}).")
        self._ref_file = value

    @property
    def on(self):
        return self._on

    @on.setter
    def on(self, value):
        if value != 'inchikey':
            raise ValueError(f"for now argument 'on' only accepts 'inchikey' as value, not '{value}'")
        self._on = value

    def find_synonyms(self, df: DataFrame) -> DataFrame:
        """Find synonyms within a DataFrame.

        The Synonyms DataFrame follows the follow syntax:

        - rowid: The InChIKey of the molecule ('inchikey').
        - idm_synonyms: The list of all molecule ids sharing the same InChIKey.
        - idm: The id of the molecule that was kept (first element of idm_synonyms).

        :param df: The input DataFrame.
        :return: The output DataFrame containg synonyms for identifying duplicate molecules.

        .. warning:: make sure there is no duplicate idm before using this function, as the final filtering of duplicates is performed with a whitelist based on idm.

        .. todo:: optimize code (table format for reference file, numpy for setting up lists when grouping?)
        """
        df = df.copy()  # this certainly removes the pandas warnings but is it the right way of doing this?
        # check on col_mol
        if self.col_mol not in df.columns:
            raise ValueError(f"Error! No column {self.col_mol} found for col_mol parameter.")
        # check on on
        if self.on != 'inchikey':
            raise ValueError(f"Error! Unauthorized value for on parameter ({self.on}).")
        elif self.on not in df.columns:
            logging.warning(f"Column {self.on} not found, so computing it.")
            df.loc[:, self.on] = df.loc[:, self.col_mol].map(rdinchi.MolToInchiKey)
        # init
        df.index = df[self.on]
        df.drop(self.on, axis=1, inplace=True)
        # define synonyms in current dataframe
        df_synonyms = pd.DataFrame(df.groupby(self.on)[self.col_id].apply(list))
        df_synonyms.rename({self.col_id: self._col_id_synonyms}, axis=1, inplace=True)
        df_synonyms[self.col_id] = df_synonyms[self._col_id_synonyms].map(lambda x: x[0])

        # df_synonyms[self.col_id] = df_synonyms[self._col_id_synonyms].map(lambda x: x[0])
        # use information stored in ref file as well, if provided
        if self.ref_file is not None:
            key = Path(self.ref_file).stem
            # init ref file if does not exist already
            if not Path(self.ref_file).is_file():
                self.init_ref_file()
            # open it with a lock as we'll need to update it at the end
            with utils.SafeHDF5Store(self.ref_file) as store:
                df_ref = store[key]
                df_ref = pd.concat([df_ref, df_synonyms])
                df_ref = pd.DataFrame(df_ref.groupby(self.on)[self._col_id_synonyms].apply(list))
                df_ref[self._col_id_synonyms] = df_ref[self._col_id_synonyms].map(lambda x: list(chain.from_iterable(x)))
                df_ref[self.col_id] = df_ref[self._col_id_synonyms].map(lambda x: x[0])
                df_ref.to_hdf(self.ref_file, key=key)
            return df_ref
        else:
            return df_synonyms

    def init_ref_file(self) -> bool:
        """Initiate an empty reference hdf for identifying duplicates.

        :return: True if the reference file could be initialized, False otherwise.
        """
        try:
            # delete file if it already exists
            if Path(self.ref_file).is_file():
                Path.unlink(self.ref_file)
            # init
            key = Path(self.ref_file).stem
            # create a new ref file
            df_ref = pd.DataFrame({self.on: [], self._col_id_synonyms: [], self.col_id: []})
            df_ref.index = df_ref[self.on]
            df_ref.drop(self.on, axis=1, inplace=True)
            df_ref.to_hdf(self.ref_file, key=key)
            logging.debug(f"Created new ref_file at '{self.ref_file}'")
            return True
        except ValueError:  # certainly not the only kind of error but PEP8 is against using just plain except.
            logging.critical(f"Could not create a new ref_file at '{self.ref_file}'")
            return False

    def mark_dupl(self, df: DataFrame) -> DataFrame:
        """Mark duplicate entries found in df or in ref. Designed to work from
        within Standardizer.run only as it updates columns 'status' and 'task' to
        respectively 'filtered' and 'filter_duplicates'. Might still be useful if exposed.

        Moreover, update the reference file (if specified) by using a lock (from utils module) so that
        no other processes can modify it at the same time. This enables a safe duplicate
        molecules filtering accross chunks without having to gather everything in a single file.

        :param df: the input DataFrame
        :return: the output DataFrame
        """
        # any molecule whose id is not in df_synonyms[col_id] is a duplicate of another
        df_synonyms = self.find_synonyms(df)
        if "status" not in df.columns:
            df["status"] = "passed"
        if "task" not in df.columns:
            df["task"] = "filter_dupl"
        # df = df.copy()  # supress warnings, but might significantly slow down the process..?
        df.loc[~df.loc[:, self.col_id].isin(df_synonyms.loc[:, self.col_id]), ['status', 'task']] = ('filtered', 'filter_dupl')
        return df
