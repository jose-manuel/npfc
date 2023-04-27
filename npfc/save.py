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


"""
Module save
===========

A module for saving DataFrames into files in different formats.
"""

# standard
import logging
import gzip
import shutil
import numpy as np
from pathlib import Path
# data science
from pandas import DataFrame
# chemoinformatics
from rdkit import Chem
from rdkit.Chem import SDWriter
# docs
from typing import List
from typing import Tuple
from typing import Union

# dev
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


from npfc.utils import FORMATS_IO
from npfc.utils import COLUMNS_MOL
from npfc.utils import COLUMNS_ENCODED


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _save_sdf(df: DataFrame,
              output_file: Union[str, Path],
              col_mol: str,
              col_idm: str,
              ) -> Tuple[str, int]:
    """
    Save molecules in input DF as a SDF, using RDKit.

    In case of the desired output file is an archive (i.e. file.sdf.gz), an uncompressed file will first
    be written (i.e. file.sdf), which may overwrite any pre-existing file with the same path.

    :param df: input DataFrame
    :param output_file: output SDF path
    :param col_mol: column with the RDKit Mol objects
    :param col_idm: column with the molecule ids, if not None, info is saved as property and as molecule title
    :return: a tuple containing the output file path and the number of saved molecules
    """
    def _format_value(value):
        """
        Helper function to format float values for export to avoid
        awkward scientific notation or flaots without trailing 0, i.e. '7.'.
        
        :param value: the value to format
        :return: the formatted value
        """
        # floats
        if np.issubdtype(type(value), np.floating):
            s = '{:f}'.format(value).rstrip("0")  # remove long trailing 0s
            if s[-1] == ".":
                s += "0"  # put the "0" back on if the result is something like "7."
            return s
        
        # others
        else:
            return str(value)

    # init
    output_file = str(output_file)
    path_output_file = Path(output_file)
    path_output_file.parent.mkdir(parents=True, exist_ok=True)
    format, compression = utils.get_file_format(output_file)

    # check input
    if format != 'SDF':
        raise ValueError("Error! Input file '{input_sdf}' has a format different from 'SDF' ('{format}')")

    # in case output file has to be compressed, remove the .gz from the filename for the output file for now
    if compression == 'gzip':
        output_archive = output_file
        output_file = str(path_output_file.parent / path_output_file.stem)
        path_output_file = Path(output_file)

    # format properties (floats)
    properties = [x for x in df.columns if x != col_mol]
    for p in properties:
        df[p] = df[p].map(_format_value)

    # init writer
    writer = SDWriter(output_file)
    writer.SetProps(properties)

    # export molecule to SDF
    for rowid, row in df.iterrows():
        # copy mol so it can be edited without consequence
        mol = Chem.Mol(row[col_mol])
        # set molecule title
        if col_idm is not None:
            mol.SetProp('_Name', str(row[col_idm]))
        else:
            mol.SetProp('_Name', '')
        # set properties
        for p in properties:
            mol.SetProp(p, row[p])
        # save mol
        writer.write(mol)
    
    # end the writing of mols
    writer.close()

    # compress the now-closed output file
    if compression == 'gzip':
        with open(output_file, 'rb') as OUTPUT:
            with gzip.open(output_archive, 'wb') as ARCHIVE:
                shutil.copyfileobj(OUTPUT, ARCHIVE)
            # delete the uncompressed file as it is only a byproduct
            Path(output_file).unlink()
            # consider output file to be the archive for return value 
            output_file = output_archive

    return output_file, len(df)


def file(df: DataFrame,
         output_file: Union[str, Path],
         col_mol: str = 'mol',
         col_idm: str = 'idm',
         csv_sep: str = '|',
         encode: bool = True,
         ):
    """
    Save an input DataFrame into a single file.

    :param df: input DataFrame
    :param output_file: output file path
    :param col_mol: for SDF format only, column with the RDKit Mol objects
    :param col_idm: for SDF format only, column with the molecule ids, if not None, info is saved as property and as molecule title
    :param csv_sep: for CSV format only, delimiter to use
    :param encode: encode mols and objects into base64 strings based on predefined column names.
    :return: a tuple containing the output file name and its number of records
    """
    # check some arguments
    utils.check_arg_output_file(output_file)
    format, compression = utils.get_file_format(output_file)
    if format not in FORMATS_IO:
        raise ValueError(f"Error! Output file format is unknown ('{format}'). Authorized formats are: {', '.join(FORMATS_IO)}.")


    # init
    logging.debug("Excerpt of the data as provided to save.file function:\n\n%s\n", df.head(5))
    path_output_file = Path(output_file)

    # avoid pandas warnings
    df = df.copy()

    # encode predefined data
    if encode:
        cols_to_encode = [x for x in COLUMNS_ENCODED if x in df.columns]
        cols_mol = [x for x in COLUMNS_MOL if x in df.columns]
    else:
        cols_to_encode = []
        cols_mol = []

    # other objects
    for c in cols_to_encode:
        df[c] = df[c].map(utils.encode_object)
    
    if format == 'SDF':
        # encode all mols except the 'mol' column
        if len(cols_mol) > 1:
            cols_mol_to_encode = [x for x in cols_mol if x != col_mol]
            logging.warning("More than 1 RDKit Mol objects columns detected, using CTAB from column '%s' and encoding the rest: %s", col_mol, ','.join(cols_mol_to_encode))
            for c in cols_mol_to_encode:
                df[c] = df[c].map(utils.encode_mol)
        return _save_sdf(df, output_file, col_mol, col_idm)
    else: # CSV/HDF
        # RDKit Mol objects
        for c in cols_mol:
            print(f"save.file: encoding column {c}")
            df[c] = df[c].map(utils.encode_mol)
        # export to other formats
        if format == 'CSV':
            df.to_csv(output_file, index=False, sep=csv_sep)
        elif format == 'HDF':
            key = path_output_file.stem.split('.')[0]
            df.to_hdf(output_file, key=key)
        else:
            raise ValueError("Error! Output file format '{format}' is not supported.")
            
    return (output_file, len(df))


def chunk(df: DataFrame,
          chunk_name_template: str,
          chunk_size: int,
          shuffle: bool = False,
          random_seed: int = None,
          col_mol: str = 'mol',
          col_idm: str = 'idm',
          csv_sep: str = '|',
          encode : bool = True,
          ) -> List[Tuple[str, int]]:
    """
    Save an input DataFrame into several chunks, written on disk.
    The input data has to be converted as a whole to a DF first.

    :param df: input DataFrame
    :param chunk_name_template: path of the output file if there was only one (i.e. dir/file.csv). Is modified to add chunk IDs (i.e. dir/file_001.csv, dir/file_002.csv, etc.).
    :param chunk_size: number of record for each chunk (last chunk might contain less)
    :param shuffle: shuffle the records before splitting into chunks
    :param random_seed: random seed to use for shuffling records
    :param col_mol: for SDF format only, column with the RDKit Mol objects
    :param col_idm: for SDF format only, column with the molecule ids, if not None, info is saved as property and as molecule title
    :param csv_sep: for CSV format only, delimiter to use
    :param encode: encode RDKit Mol and other predefined objects into base64 strings.
    :return: a list of tuples containing each chunk name and its number of records
    """
    # check arguments
    if not isinstance(chunk_size, int) and chunk_size < 1:
        raise ValueError("Error! Argument chunk_size needs to be defined as a non-zero positive integer. (value was: '{chunk_size}')")
    
    # shuffle records to remove bias from input order 
    if shuffle:
        df = df.sample(frac=1, random_state=random_seed)

    # init iteration
    num_records = len(df)
    output_chunks = []
    start = 0
    output_chunk_template = Path(chunk_name_template)
    output_dir_path = output_chunk_template.parent
    output_chunk_basename = Path(output_chunk_template.stem).stem  # 2-layer stem for cases like 'file.sdf.gz'
    output_chunk_suffixes = output_chunk_template.suffixes
    
    # iteration
    for i, start in enumerate(range(0, num_records, chunk_size)):
        end = start + chunk_size
        output_chunk = output_dir_path / Path(output_chunk_basename + "_" + str(i+1).zfill(3) + ''.join(output_chunk_suffixes))
        results = file(df=df.iloc[start:end],
                       output_file=output_chunk,
                       col_mol=col_mol,
                       col_idm=col_idm,
                       csv_sep=csv_sep,
                       encode=encode,
                       )
        output_chunks.append(results)

    return output_chunks


def chunk_sdf(input_sdf: str,
              output_dir: str,
              chunk_size: int = None,
              prefix: str = None,
              keep_uncompressed: bool = False,
              ) -> List[Tuple[str, int]]:
    """
    Split an input SDF file into SDF chunks using memory-efficient line by line text parsing, suitable for large files.
    Molecules are not parsed, no change is made to the molblocks.

    :param input_sdf: input SDF
    :param output_dir: 
    :param chunk_name_template: path of the output file as if there was only one (i.e. dir/file.csv). It is modified to add chunk IDs (i.e. dir/file_001.csv, dir/file_002.csv, etc.).
    :param chunk_size: number of record for each chunk (last chunk might contain less)
    :param prefix: prefix to use for chunks. If left to None, the input SDF filename will be used.
    :param keep_uncompressed: in case of gzip input, keep the uncompressed file instead of deleting it as a temp file
    :return: a list of tuples containing each chunk name and its number of records

    TODO: support for gzip outputs
    TODO: support for gzip input
    """

    def _save_chunk(current_lines, output_dir_path, output_chunk_basename, current_chunk_idx, format, compression):
        """Helper function for saving currently gathered text into the current uncompressed chunk."""
        # save to uncompressed  output file
        output_chunk = output_dir_path / Path(f"{output_chunk_basename}_{str(current_chunk_idx).zfill(3)}.{format.lower()}")
        with open(output_chunk, 'w+') as OFH:
            OFH.write(''.join(current_lines))

        # compress the now-closed output file
        if compression == 'gzip':
            output_archive = f"{output_chunk}.gz"
            with open(output_chunk, 'rb') as OUTPUT:
                with gzip.open(output_archive, 'wb') as ARCHIVE:
                    shutil.copyfileobj(OUTPUT, ARCHIVE)
                # delete the uncompressed file as it is only a byproduct
                Path(output_chunk).unlink()
                # consider output file to be the archive for return value 
                output_chunk = output_archive
    
        return output_chunk

    # define I/O
    input_dir_path = Path(input_sdf).parent
    output_dir_path = Path(output_dir)
    # output_chunk_basename = Path(Path(input_sdf).stem).stem

    # path_chunk_template = Path(chunk_name_template)
    # output_dir_path = path_chunk_template.parent
    format, compression = utils.get_file_format(input_sdf)
    path_chunk_template = output_dir_path / Path(input_sdf).name
    # format_out, compression_out = utils.get_file_format(chunk_name_template)

    # in case input file is compressed, uncompress it
    if compression == 'gzip':
        input_sdf_uncompressed = str(input_dir_path / Path(input_sdf).stem)
        with open(input_sdf_uncompressed, 'wb') as FH_UNCOMPRESSED, gzip.open(input_sdf, 'rb') as FH_COMPRESSED:
            bindata = FH_COMPRESSED.read()
            FH_UNCOMPRESSED.write(bindata)
        input_sdf = input_sdf_uncompressed

        # in case output file has to be compressed, remove the .gz from the filename for the output file for now
        output_archive = path_chunk_template
        path_chunk_template = output_dir_path / path_chunk_template.stem

    if prefix is None:
        output_chunk_basename = path_chunk_template.stem
    else:
        output_chunk_basename = prefix

    # check format and compression
    if format != 'SDF' or compression not in (None, 'gzip'):
        raise ValueError("Error! Input SDF should be of format SDF and, if compressed, be a gzip archive.")
    
    # begin        
    with open(input_sdf,'r') as FH:

        # init iteration
        current_num_record = 0
        current_chunk_idx = 1
        output_chunks = []
        current_lines = []
        file_read = False

        # iteration line by line
        while True:
            line = FH.readline()

            # exit when reached the end of the file
            if not line:
                FH.close()
                file_read = True
                break
            
            # identify end of record
            if line.strip() == '$$$$':
                current_num_record += 1 

            # append the lines to the current chunk
            current_lines.append(line)

            # save chunk when enough records are gathered
            if current_num_record >= chunk_size:
                output_chunk = _save_chunk(current_lines, output_dir_path, output_chunk_basename, current_chunk_idx, format, compression)
                output_chunks.append((output_chunk, current_num_record))
                
                # clean-up
                current_lines = []
                current_num_record = 0
                current_chunk_idx += 1
        
        # saving last chunk with current records
        if file_read and len(current_lines) > 0:
            output_chunk = _save_chunk(current_lines, output_dir_path, output_chunk_basename, current_chunk_idx, format, compression)
            output_chunks.append((output_chunk, current_num_record))
            
    # clean-up
    if compression == 'gzip' and not keep_uncompressed:
        Path(input_sdf).unlink()  # remove the temporary uncompressed input SDF

    return output_chunks


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
