"""
Module save
================

A module containing the Saver class, used for storing DataFrames with molecules
on disk.
"""

# standard
import logging
import gzip
import shutil
import os
import numpy as np
import time
from pathlib import Path
from random import random
# data science
import pandas as pd
from pandas import DataFrame
from pandas import HDFStore
# chemoinformatics
from rdkit import Chem
from rdkit.Chem import SDWriter
# docs
from typing import List
# dev
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def file(df: pd.DataFrame,
         output_file: str,
         shuffle: bool = False,
         random_seed: int = None,
         chunk_size: int = None,
         encode: bool = True,
         col_mol: str = 'mol',
         col_id: str = 'idm',
         csv_sep: str = '|'):
    """A method for saving DataFrames with molecules to different file types.
    This is handy way of using the Saver class without having to keep a Saver object.

    :param df: the input DataFrame
    :param output_file: the output file
    :param shuffle: randomize records
    :param random_seed: a number for reproducing the shuffling
    :param chunk_size: the maximum number of records per chunk. If this value is unset, no chunking is performed, otherwise each chunk filename gets appended with a suffix: file_XXX.ext.
    :param encode: encode RDKit Mol objects and other objects in predefined columns as base64 strings.
    :param col_mol: if molecules need to be encoded, then the encoding is perfomed on this column.
    :param csv_sep: separator to use in case of csv output
    :return: the list of output files with their number of records
    """
    # check some arguments
    utils.check_arg_output_file(output_file)
    utils.check_arg_bool(shuffle)
    utils.check_arg_bool(encode)

    logging.debug("Excerpt of the data as provided to file function:\n\n%s\n", df.head(5))

    # init
    path_output_file = Path(output_file)
    ext_output_file = path_output_file.suffixes
    output_dir = path_output_file.resolve().parent
    output_files = []
    # for sdf, molecules cannot be encoded
    if ext_output_file[0] == '.sdf' and encode:
        logging.warning(f"Format is SDF, so column '{col_mol}' is not encoded.")
    # avoid pandas warnings
    df = df.copy()
    # shuffle
    if shuffle:
        logging.debug('Shuffling rows before saving file')
        df = df.sample(frac=1, random_state=random_seed)
    # encode predefined data
    # if nothing to encode, just don't
    if len(df.index) == 0:
        logging.warning("DataFrame is empty, skip encoding.")
    # in case there is stuff to encode, encode it:
    elif encode:
        # for SDF files, RDKit Mol objects to use for MolBlocks should not be encoded
        if ext_output_file[0] != '.sdf' and col_mol in df.columns:
            # df[col_mol] = df[col_mol].map(utils.encode_mol_smiles)
            df[col_mol] = df[col_mol].map(utils.encode_mol)
        # other RDKit Mol objects can be encoded though
        for col in ("mol", "mol_frag", "mol_frag_1", "mol_frag_2", "mol_rdkit"):
            if col in df.columns and col != col_mol:
                # df[col] = df[col].map(utils.encode_mol_smiles)
                df[col] = df[col].map(utils.encode_mol)
        # other objects are labelled with leading '_'
        for col in df.columns:
            if col.startswith('_') and col != '_Name':
                df[col] = df[col].map(utils.encode_object)

    logging.debug("Excerpt of the data to save before chuking:\n\n%s\n", df.head(3))

    # chunking
    if chunk_size is None:
        # single output
        _save(df=df, output_file=output_file, col_mol=col_mol, col_id=col_id, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], csv_sep=csv_sep)
        output_files.append([output_file, len(df.index)])
    else:
        # chunks
        start = 0
        j = 1
        for start in range(0, len(df.index), chunk_size):
            end = start + chunk_size
            output_chunk = str(output_dir) + "/" + path_output_file.stem.split('.')[0] + "_" + str(j).zfill(3) + ''.join(ext_output_file)  # stem returns file.csv for file.csv.gz
            _save(df=df.iloc[start:end], output_file=output_chunk, col_mol=col_mol, col_id=col_id, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], csv_sep=csv_sep)
            output_files.append([output_chunk, len(df.iloc[start:end].index)])
            j += 1
        logging.debug("%s chunks were created", len(output_files))

    return output_files


def _save(df: DataFrame,
          output_file: str,
          col_mol: str,
          col_id: str,
          suffixes: List[str],
          key: str,
          csv_sep: str):
    """Helper function for the save method.
    Does the actual export to the output file and picks a format based on provided infos.

    :param df: the input DataFrame
    :param suffixes: the suffixes of the output file
    :param key: the key for a HDF file
    :param csv_sep: the separator for a CSV file
    """
    # infer from pandas.to_csv does not work as expected (no compression!)
    # so I need to specify the compression type manually.
    utils.check_arg_output_file(output_file)
    out_format, out_compression = utils.get_file_format(output_file)
    if out_format == 'CSV':
        if out_compression == 'gzip':
            df.to_csv(output_file, sep=csv_sep, compression=out_compression, index=False)
        else:
            df.to_csv(output_file, sep=csv_sep, index=False)
    elif out_format == 'HDF':
        df.to_hdf(output_file, key=key)
    elif out_format == 'SDF':
        # write the uncompressed file
        if out_compression == 'gzip':
            # init
            output_file_base = '.'.join(output_file.split('.')[:-1])
            logging.debug("Output_file_base: %s", output_file_base)
            # write the file uncompressed
            write_sdf(df, output_file_base, molColName=col_mol, idName=col_id, properties=list(df.columns))
            # compress the file
            with open(output_file_base, 'rb') as OUTPUT:
                with gzip.open(output_file, 'wb') as ARCHIVE:
                    shutil.copyfileobj(OUTPUT, ARCHIVE)
                # delete the uncompressed file as it is only a byproduct
                Path(output_file_base).unlink()
        else:
            write_sdf(df, output_file, molColName=col_mol, idName=col_id, properties=list(df.columns))
    elif out_format == 'FEATHER':
        df.to_feather(output_file)
    else:
        raise ValueError(f"Error! Cannot save DataFrame to unexpected format '{suffixes[0]}'.")
    logging.debug("Saved %s records at '%s'.", len(df.index), output_file)


def write_sdf(df, out, molColName='ROMol', idName=None, properties=None, allNumeric=False):
    """
    Redefinition of PandasTools.WriteSDF because RDKit 2019.03.1 is incompatible with Pandas 25.1.

    Write an SD file for the molecules in the dataframe. Dataframe columns can be exported as
    SDF tags if specified in the "properties" list. "properties=list(df.columns)" would export
    all columns.
    The "allNumeric" flag allows to automatically include all numeric columns in the output.
    User has to make sure that correct data type is assigned to column.
    "idName" can be used to select a column to serve as molecule title. It can be set to
    "RowID" to use the dataframe row key as title.
    """
    close = None
    if isinstance(out, str):
        if out.lower()[-3:] == ".gz":
            out = gzip.open(out, "wt")
            close = out.close

    writer = SDWriter(out)
    if properties is None:
        properties = []
    else:
        properties = list(properties)
    if allNumeric:
        properties.extend([
          dt for dt in df.dtypes.keys()
          if (np.issubdtype(df.dtypes[dt], np.floating) or np.issubdtype(df.dtypes[dt], np.integer))
        ])

    if molColName in properties:
        properties.remove(molColName)
    if idName in properties:
        properties.remove(idName)
    writer.SetProps(properties)
    for row in df.iterrows():
        # make a local copy I can modify
        mol = Chem.Mol(row[1][molColName])

        if idName is not None:
            if idName == 'RowID':
                mol.SetProp('_Name', str(row[0]))
            else:
                mol.SetProp('_Name', str(row[1][idName]))
        for p in properties:
            cell_value = row[1][p]
            # Make sure float does not get formatted in E notation
            if np.issubdtype(type(cell_value), np.floating):
                s = '{:f}'.format(cell_value).rstrip("0")  # "f" will show 7.0 as 7.00000
                if s[-1] == ".":
                    s += "0"  # put the "0" back on if it's something like "7."
                mol.SetProp(p, s)
            else:
                mol.SetProp(p, str(cell_value))
        writer.write(mol)
    writer.close()
    if close is not None:
        close()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class SafeHDF5Store(HDFStore):
    """Implement safe HDFStore by obtaining file lock. Multiple writes will queue if lock is not obtained.

    Edited after:
    https://stackoverflow.com/questions/41231678/obtaining-a-exclusive-lock-when-writing-to-an-hdf5-file
    """

    def __init__(self, *args, **kwargs):
        """Initialize and obtain file lock."""

        interval = kwargs.pop('probe_interval', random())
        self._lock = f"{args[0]}.lock"
        while True:
            try:
                self._flock = os.open(self._lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                break
            except (IOError, OSError):
                time.sleep(interval)

        HDFStore.__init__(self, *args, **kwargs)

    def __exit__(self, *args, **kwargs):
        """Exit and remove file lock."""

        HDFStore.__exit__(self, *args, **kwargs)
        os.close(self._flock)
        os.remove(self._lock)
