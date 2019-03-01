"""
Module save
================
"""

# standard
import logging
import base64
from pathlib import Path
# data science
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import PandasTools
# docs
from typing import List
# dev
from npfc import utils

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


class Saver:
    """A class for saving DataFrames with molecules to different file types."""

    def __init__(self,
                 shuffle: bool = False,
                 random_seed: int = None,
                 chunk_size: int = None,
                 encode_mols: bool = True,
                 col_mol: str = 'mol',
                 col_id: str = 'idm'):
        """Create a Saver object with following parameters:

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
            self._random_seed = value

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
        return self._col_mol

    @col_id.setter
    def col_id(self, value: str):
        self._col_id = str(value)

    def encode_mol(self, mol: Mol) -> str:
        """Convert a molecule to binary and then represent this binary as base64 string.

        :param mol: the input molecule
        :return: the molecule in base64
        """
        try:
            return base64.b64encode(mol.ToBinary()).decode()
        except AttributeError:
            return None

    def _save(self, df: pd.DataFrame, output_file: str, suffixes: List[str], key: str, sep: str):
        """Helper function for the save method.
        Does the actual export to the output file and picks a format based on provided infos.

        :param df: the input DataFrame
        :suffixes: the suffixes of the output file
        :key: the key for a HDF file
        :sep: the separator for a CSV file
        """
        if suffixes[0] == '.csv':
            logging.debug(f"Saved {len(df.index)}")
            df.to_csv(output_file, sep=sep)
        elif suffixes[0] == '.hdf':
            df.to_hdf(output_file, key=key)
        elif suffixes[0] == '.sdf':
            PandasTools.WriteSDF(df, output_file, molColName=self.col_mol, idName=self.col_id, properties=list(df.columns))
        else:
            raise ValueError(f"Error! Cannot save DataFrame to unexpected format '{suffixes[0]}'.")

    def save(self, df: pd.DataFrame, output_file: str) -> List[List[str, int]]:
        """Save the input DataFrame on disk using Saver object parameters.

        :param df: the input DataFrame
        :output_file: the target file
        :return: the list of output files with their number of records
        """
        # check output_file
        utils.check_arg_output_file(output_file)
        # init
        path_output_file = Path(output_file)
        ext_output_file = path_output_file.suffixes
        output_dir = path_output_file.resolve().parent
        output_files = []
        # avoid pandas warnings
        df = df.copy()
        # shuffle
        if self.shuffle:
            df = df.sample(frac=1, random_state=self.random_seed)
        # encode molecules
        if self.encode_mols:
            df[self.col_mol] = df[self.col_mol].map(self.encode_mol_base64)
        # chunking
        if self.chunk_size is None:
            # single output
            self._save(df, output_file, suffixes=ext_output_file, key=path_output_file.stem, sep='|')
            output_files.append(output_file, len(df.index))
            logging.debug(f"Saved {len(df.index)} records to '{output_file}'.")
        else:
            # chunks
            start = 0
            j = 0
            for start in range(0, len(df.index), self.chunk_size):
                end = start + self.chunk_size
                output_chunk = str(output_dir) + "/" + output_file.stem + "_" + str(j).zfill(3) + ''.join(ext_output_file)
                output_files.append(self._save(df.iloc[start:end], ext_output_file[0], output_chunk))
                j += 1
            logging.debug(f"{len(output_files)} chunks were created")
