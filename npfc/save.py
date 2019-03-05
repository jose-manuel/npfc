"""
Module save
================

A module containing the Saver class, used for storing DataFrames with molecules
on disk.
"""

# standard
import logging
import base64
from pathlib import Path
# data science
import pandas as pd
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import PandasTools
# docs
from typing import List
# dev
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def encode_mol(mol: Mol) -> str:
    """Convert a molecule to binary and then represent this binary as base64 string.

    :param mol: the input molecule
    :return: the molecule in base64
    """
    try:
        return base64.b64encode(mol.ToBinary()).decode()
    except AttributeError:
        return None


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def save(df: pd.DataFrame, output_file: str,
         shuffle: bool = False, random_seed: int = None,
         chunk_size: int = None, encode_mols: bool = True,
         col_mol: str = 'mol', col_id: str = 'idm')-> utils.Output_files:
        """A method for saving DataFrames with molecules to different file types.
        This is handy way of using the Saver class without having to keep a Saver object.

        .. note:: I did this because it was the most easiest way, I guess I'll have to investigate class methods for a better implementation.

        :param df: the input DataFrame
        :param output_file: the output file
        :param shuffle: randomize records
        :param random_seed: a number for reproducing the shuffling
        :param chunk_size: the maximum number of records per chunk. If this value is unset, no chunking is performed, otherwise each chunk filename gets appended with a suffix: file_XXX.ext.
        :param encode_mols: convert rdkit.Chem.Mol objects to base64 string representation. For HDF format, pandas stops complaining about PerformanceWarning, for csv molecules do not need to parsed again.
        :param col_mol: if molecules need to be encoded, then the encoding is perfomed on this column.
        :return: the list of output files with their number of records

        """
        s = Saver(shuffle=shuffle,
                  random_seed=random_seed,
                  chunk_size=chunk_size,
                  encode_mols=encode_mols,
                  col_mol=col_mol,
                  col_id=col_id,
                  )
        return s.save(df, output_file)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Saver:
    """A class for saving DataFrames with molecules to different file types."""

    def __init__(self,
                 shuffle: bool = False,
                 random_seed: int = None,
                 chunk_size: int = None,
                 encode_mols: bool = True,
                 col_mol: str = 'mol',
                 col_id: str = 'idm'):
        """Create a Saver object with below parameters.

        .. note:: The reason I have a class for this is because I wanted to avoid at lot of redundant code. But a class is not that perfect if one Saver object has to instanciated everytime. As a work-around, I also added a save function using this class.

        :param shuffle: randomize records
        :param random_seed: a number for reproducing the shuffling
        :param chunk_size: the maximum number of records per chunk. If this value is unset, no chunking is performed, otherwise each chunk filename gets appended with a suffix: file_XXX.ext.
        :param encode_mols: convert rdkit.Chem.Mol objects to base64 string representation. For HDF format, pandas stops complaining about PerformanceWarning, for csv molecules do not need to parsed again.
        :param col_mol: if molecules need to be encoded, then the encoding is perfomed on this column.
        """
        self._shuffle = shuffle
        self._random_seed = random_seed
        self._chunk_size = chunk_size
        self._encode_mols = encode_mols
        self._col_mol = col_mol
        self._col_id = col_id

    @property
    def shuffle(self):
        return self._shuffle

    @shuffle.setter
    def shuffle(self, value: bool):
        if utils.check_arg_bool(value):
            self._shuffle = value

    @property
    def random_seed(self):
        return self._random_seed

    @random_seed.setter
    def random_seed(self, value: int):
        if utils.check_arg_positive_number(value):
            self._random_seed = value

    @property
    def chunk_size(self):
        return self._chunk_size

    @chunk_size.setter
    def chunk_size(self, value: int):
        if utils.check_arg_positive_number(value):
            self._chunk_size = value

    @property
    def encode_mols(self):
        return self._encode_mols

    @encode_mols.setter
    def encode_mols(self, value: bool):
        if utils.check_arg_bool(value):
            self._encode_mols = value

    @property
    def col_mol(self):
        return self._col_mol

    @col_mol.setter
    def col_mol(self, value: str):
        self._col_mol = str(value)

    @property
    def col_id(self):
        return self._col_id

    @col_id.setter
    def col_id(self, value: str):
        self._col_id = str(value)

    def _save(self, df: pd.DataFrame, output_file: str, suffixes: List[str], key: str, sep: str):
        """Helper function for the save method.
        Does the actual export to the output file and picks a format based on provided infos.

        :param df: the input DataFrame
        :param suffixes: the suffixes of the output file
        :param key: the key for a HDF file
        :param sep: the separator for a CSV file
        """
        if suffixes[0] == '.csv':
            df.to_csv(output_file, sep=sep)
        elif suffixes[0] == '.hdf':
            df.to_hdf(output_file, key=key)
        elif suffixes[0] == '.sdf':
            PandasTools.WriteSDF(df, output_file, molColName=self.col_mol, idName=self.col_id, properties=list(df.columns))
        else:
            raise ValueError(f"Error! Cannot save DataFrame to unexpected format '{suffixes[0]}'.")
        logging.debug(f"Saved {len(df.index)} records at '{output_file}'.")

    def save(self, df: pd.DataFrame, output_file: str) -> utils.Output_files:
        """Save the input DataFrame on disk using Saver object parameters.

        :param df: the input DataFrame
        :param output_file: the output file
        :return: the list of output files with their number of records
        """
        # check output_file
        utils.check_arg_output_file(output_file)
        # init
        path_output_file = Path(output_file)
        ext_output_file = path_output_file.suffixes
        output_dir = path_output_file.resolve().parent
        output_files = []
        # for sdf, molecules cannot be encoded
        if ext_output_file[0] == '.sdf' and self.encode_mols:
            logging.warning(f"Format is SDF, so molecules are not encoded.")
        # avoid pandas warnings
        df = df.copy()
        # shuffle
        if self.shuffle:
            df = df.sample(frac=1, random_state=self.random_seed)
        # encode molecules
        if self.encode_mols and ext_output_file[0] != '.sdf':
            df[self.col_mol] = df[self.col_mol].map(encode_mol)
        # chunking
        if self.chunk_size is None:
            # single output
            self._save(df=df, output_file=output_file, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], sep='|')
            output_files.append([output_file, len(df.index)])
        else:
            # chunks
            start = 0
            j = 0
            for start in range(0, len(df.index), self.chunk_size):
                end = start + self.chunk_size
                output_chunk = str(output_dir) + "/" + path_output_file.stem.split('.')[0] + "_" + str(j).zfill(3) + ''.join(ext_output_file)  # stem returns file.csv for file.csv.gz
                self._save(df=df.iloc[start:end], output_file=output_chunk, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], sep='|')
                output_files.append([output_chunk, len(df.iloc[start:end].index)])
                j += 1
            logging.debug(f"{len(output_files)} chunks were created")

        return output_files
