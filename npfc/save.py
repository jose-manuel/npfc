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
import base64
from pathlib import Path
# data science
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit.Chem import Mol
from rdkit.Chem import PandasTools
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
         sep: str = '|'):
    """A method for saving DataFrames with molecules to different file types.
    This is handy way of using the Saver class without having to keep a Saver object.

    :param df: the input DataFrame
    :param output_file: the output file
    :param shuffle: randomize records
    :param random_seed: a number for reproducing the shuffling
    :param chunk_size: the maximum number of records per chunk. If this value is unset, no chunking is performed, otherwise each chunk filename gets appended with a suffix: file_XXX.ext.
    :param encode: encode RDKit Mol objects and other objects in predefined columns as base64 strings.
    :param col_mol: if molecules need to be encoded, then the encoding is perfomed on this column.
    :return: the list of output files with their number of records
    """
    # check some arguments
    utils.check_arg_output_file(output_file)
    utils.check_arg_bool(shuffle)
    utils.check_arg_bool(encode)

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
        df = df.sample(frac=1, random_state=random_seed)
    # encode predefined data
    # if nothing to encode, just don't
    if len(df.index) == 0:
        logging.warning("DataFrame is empty, skip encoding.")
    # in case there is stuff to encode, encode it:
    elif encode:
        # for SDF files, RDKit Mol objects to use for MolBlocks should not be encoded
        if ext_output_file[0] != '.sdf':
            df[col_mol] = df[col_mol].map(utils.encode_mol)
        # other RDKit Mol objects can be though
        for col in ("mol", "mol_frag"):
            if col in df.columns and col != col_mol:
                df[col] = df[col].map(utils.encode_mol)
        # encode other predefined objects
        for col in ('graph', 'colormap', 'aidxf', 'aidxf1', 'aidxf2', 'd_aidxs'):
            if col in df.columns:
                df[col] = df[col].map(utils.encode_object)

    # chunking
    if chunk_size is None:
        # single output
        _save(df=df, output_file=output_file, col_mol=col_mol, col_id=col_id, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], sep=sep)
        output_files.append([output_file, len(df.index)])
    else:
        # chunks
        start = 0
        j = 1
        for start in range(0, len(df.index), chunk_size):
            end = start + chunk_size
            output_chunk = str(output_dir) + "/" + path_output_file.stem.split('.')[0] + "_" + str(j).zfill(3) + ''.join(ext_output_file)  # stem returns file.csv for file.csv.gz
            _save(df=df.iloc[start:end], output_file=output_chunk, col_mol=col_mol, col_id=col_id, suffixes=ext_output_file, key=path_output_file.stem.split('.')[0], sep=sep)
            output_files.append([output_chunk, len(df.iloc[start:end].index)])
            j += 1
        logging.debug(f"{len(output_files)} chunks were created")

    return output_files


def _save(df: DataFrame,
          output_file: str,
          col_mol: str,
          col_id: str,
          suffixes: List[str],
          key: str,
          sep: str):
    """Helper function for the save method.
    Does the actual export to the output file and picks a format based on provided infos.

    :param df: the input DataFrame
    :param suffixes: the suffixes of the output file
    :param key: the key for a HDF file
    :param sep: the separator for a CSV file
    """
    # infer from pandas.to_csv does not work as expected (no compression!)
    # so I need to specify the compression type manually.
    utils.check_arg_output_file(output_file)
    format, compression = utils.get_file_format(output_file)
    if format == 'CSV':
        if compression == 'gzip':
            df.to_csv(output_file, sep=sep, compression=compression)
        else:
            df.to_csv(output_file, sep=sep)
    elif format == 'HDF':
        df.to_hdf(output_file, key=key)
    elif format == 'SDF':
        # write the uncompressed file
        if compression == 'gzip':
            # init
            output_file_base = '.'.join(output_file.split('.')[:-1])
            logging.debug(f"Output_file_base: {output_file_base}")
            # write the file uncompressed
            PandasTools.WriteSDF(df, output_file_base, molColName=col_mol, idName=col_id, properties=list(df.columns))
            # compress the file
            with open(output_file_base, 'rb') as OUTPUT:
                with gzip.open(output_file, 'wb') as ARCHIVE:
                    shutil.copyfileobj(OUTPUT, ARCHIVE)
                # delete the uncompressed file as it is only a byproduct
                Path(output_file_base).unlink()
        else:
            PandasTools.WriteSDF(df, output_file, molColName=col_mol, idName=col_id, properties=list(df.columns))
    else:
        raise ValueError(f"Error! Cannot save DataFrame to unexpected format '{suffixes[0]}'.")
    logging.debug(f"Saved {len(df.index)} records at '{output_file}'.")
