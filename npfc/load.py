"""
Module load
============

A module for loading files in different formats into DataFrames.
"""

# standard
import logging
import gzip
from pathlib import Path
# data science
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit import Chem
# dev
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

from npfc.utils import COLUMN_IDM
from npfc.utils import COLUMNS_MOL
from npfc.utils import COLUMNS_ENCODED
from npfc.utils import FORMATS_IO

CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'rdkit': lambda x: x,  # nothing to do here, but removes the needs for more if/elif
              }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _from_sdf(input_sdf: str, col_mol: str = 'mol'):
    """
    Parses a SDF and load molecules into a DataFrame.
    Contrarly to PandasTools.LoadSDF function, empty molecules are not silently
    filtered out, so their properties can be accessed for better error tracking.

    :param input_sdf: input SDF
    :param col_mol: column name with the RDKit Mol objects
    :param compression: 
    """
    # determine compression
    compression = None
    if str(input_sdf).endswith('.gz'):
        compression = 'gzip'
    
    # create two file handlers, one for mols and one for properties (raw)
    if compression is None:
        logging.debug("Read sdf from uncompressed file")
        FH_mols = open(input_sdf, 'rb')
        FH_props = open(input_sdf, 'rb')
    elif compression == 'gzip':
        logging.debug("Read sdf from compressed file (%s)", compression)
        FH_mols = gzip.open(input_sdf)
        FH_props = gzip.open(input_sdf)
    else:
        raise ValueError(f"Error! Unknown compression type for SDF: '{compression}'")
    # init
    i = 0
    rows = []
    row_idx = []
    # double iteration over sdf but generators so it should be ok memory-wise
    for mol, mol_raw in zip(Chem.ForwardSDMolSupplier(FH_mols), Chem.ForwardSDMolSupplier(FH_props, sanitize=False)):
        try:
            # properties
            row = dict((k, mol_raw.GetProp(k)) for k in mol_raw.GetPropNames())
            # molecule title
            row['_Name'] = mol_raw.GetProp('_Name')
            # mol
            row[col_mol] = mol
        except AttributeError:
            logging.warning(f"Molecule #{i} could not be parsed and was skipped!")
            row = None
        # record entry
        if row is not None:
            rows.append(row)
            row_idx.append(i)
        i += 1

    # clean-up
    FH_mols.close()
    FH_props.close()

    return DataFrame(rows, index=row_idx)


def file(input_file: str,
         col_mol: str = None,
         col_idm: str = None,
         mol_format: str = 'rdkit',
         keep_props: bool = True,
         decode: bool = True,
         csv_sep = '|',
         ):
    """Load a file into a DataFrame.

    :param input_file: the input file to load
    :param col_idm: the column/property to use for molecule ids. If left by default and no idm col is found, then _Name is used instead. If this property is not set, then a sequential idm will be generated (MOL_0000001, etc.).
    :param col_mol: the column to use for molecules (irrerlevant for SDF)
    :param csv_sep: the column separator to use for parsing the input file (CSV)
    :param mol_format: the input format for molecules
    :param out_id: the column name used for storing molecule ids
    :param out_mol: the column name used for storing molecules
    :param keep_props: keep all properties found in the input file. If False, then only out_id and out_mol are kept.
    :param decode: decode base64 strings into objects. Columns with encoded objects are labelled with a leading '_'. For molecules, reserved names are 'mol' and 'mol_frag'.
    :return: a DataFrame
    
    ..warning:: if a 'idm' property exists in the input file but the user picks another property for in_id, the pre-existing 'idm' will be renamed into 'idm.1' (and overwritten if already present). 
    """
    # init
    out_id = 'idm'
    out_mol = 'mol'

    # check arguments
    utils.check_arg_input_file(input_file)
    path_input_file = Path(input_file)
    format, compression = utils.get_file_format(input_file)
    if format not in FORMATS_IO:
        raise ValueError("Error! Unsupported format for input file '{input_file}' ('{format}').")
    logging.debug("Loading file '%s' (format=%s, compression=%s)", input_file, format, compression)
    
    # read mols
    if format == 'SDF':
        df = _from_sdf(input_file)
    elif format == 'HDF':
        df = pd.read_hdf(input_file)  #, key=path_input_file.stem)  # does not work anymore if pytables is not used?
    elif format == 'CSV':
        df = pd.read_csv(input_file, sep=csv_sep)
    else:
        raise ValueError("Error! Unsupported format for file '{input_file}' ('{format}')")
    
    # process data
    # logging.debug("EXCERPT OF PARSED DATA BEFORE MODIFICATION:\n%s", df.head(3))

    # automatically use idm and mol columns, when possible, until stated otherwise
    if col_idm is None:   
        if 'idm' in df.columns:
            col_idm = 'idm'
        elif format == 'SDF':
            col_idm = '_Name'

    if col_mol is None and 'mol' in df.columns:
        col_mol = 'mol'
    
    # back up idm col if already existing in props but not used as idm
    if keep_props and col_idm != out_id and out_id in df.columns:
        df = df.rename({out_id: f"{out_id}.1"}, axis=1)
        logging.warning("Warning! Reserved name col_idm '%s' is already present in DF, renaming it to '%s' to avoid overwriting column.", out_id, f"{out_id}.1")

    # set idm/mol
    if col_idm is not None:
        df['idm'] = df[col_idm]
    if 'idm' in df.columns:
        df['idm'] = df['idm'].astype(str)  # ids should never be numbers, everything is simpler as str to compare between datasets, etc.
    if col_mol is not None and col_mol != 'mol':
        df['mol'] = df[col_mol]

    # keep_props
    if not keep_props:
        df = df.drop([c for c in df.columns if c not in (out_id, out_mol)], axis=1)

    # mol_format: create RDKit Mol objects for SMILEs and MolBlocks, but leave (encoded) RDKit as found 
    if out_mol in df.columns:
        df[out_mol] = df[out_mol].map(CONVERTERS[mol_format])  

    # decode
    if decode:

        # handle mols differently becasue it is faster this way
        # alread in Mol format for SDF files
        cols_to_ignore = []
        # in case of SDF, mols are already RDKit Mol objects
        if format == 'SDF':
            cols_to_ignore.append('mol')

        # decode mols 
        cols_to_decode = [c for c in COLUMNS_MOL if c not in cols_to_ignore and c in df.columns]
        logging.debug("COLUMNS WITH MOLECULES TO DECODE: %s", ','.join(cols_to_decode))
        for col in cols_to_decode:
            df[col] = df[col].map(utils.decode_mol)

        # decode other objects
        cols_to_decode = [c for c in COLUMNS_ENCODED if c in df.columns]
        logging.debug("COLUMNS WITH OBJECTS TO DECODE: %s", ','.join(cols_to_decode))
        for col in cols_to_decode:
            df[col] = df[col].map(utils.decode_object)
        
    logging.debug("First 3 rows loaded:\n\n%s\n", df.head(3))

    return df


def count_mols(input_file: str, buffer_size: int = 10240, keep_uncompressed: bool = False):
    """
    Count the number of molecules in an input file.
    The method varies depending on the format:

        - SDF: count the $$$$ pattern
        - CSV: count the number of lines, minus 1 for column headers
        - HDF: load file into memory using Pandas and then count number of rows
        - PARQUET: not implemented yet
    
    This function is optmized for memory, so it should handle very large files, apart from HDF files.
    
    :param input_file: input file
    :param buffer_size: buffer size in bytes to use for scanning the input text file (SDF and CSV). Default is 10Mb.
    :param keep_uncompressed: in case of gzip file, leave the uncompressed file after execution
    :return: counf of molecules
    """

    def _count_generator(reader, buffer_size):
        """Use a generator it iterate over the raw data for larger files:
        https://pynative.com/python-count-number-of-lines-in-file/
        """
        b = reader(buffer_size * buffer_size)
        while b:
            yield b
            b = reader(buffer_size * buffer_size)

    # check arguments
    utils.check_arg_input_file(input_file)
    path_input_file = Path(input_file)
    format, compression = utils.get_file_format(input_file)
    if format not in FORMATS_IO:
        raise ValueError("Error! Unsupported format for input file '{input_file}' ('{format}').")
    
    # in case of compressed file, uncompress it
    if compression == 'gzip':
        uncompressed_file = input_file.split('.gz')[0]
        with open(uncompressed_file, 'wb') as uncompressed, gzip.open(input_file, 'rb') as compressed:
            bindata = compressed.read()
            uncompressed.write(bindata)
        input_file = uncompressed_file

    # count mols
    with open(input_file, 'rb') as FH:

        # create an iterator 
        c_generator = _count_generator(FH.raw.read, buffer_size)

        # apply pattern to count molecules
        if format == 'SDF':
            count = sum(buffer.count(b'$$$$') for buffer in c_generator)
        elif format == 'CSV':
            count = sum(buffer.count(b'\n') for buffer in c_generator) - 1  # ignore headers
        elif format == 'HDF':
            count = len(file(input_file))
        elif format == 'PARQUET':
            count = -1
            logging.warning("PARQUET FORMAT NOT IMPLEMENTED YET")
        else:
            raise ValueError("Error! Unsupported format for input file '{input_file}' ('{format}').")

    # remove uncompressed file
    if compression == 'gzip' and not keep_uncompressed:
        Path(uncompressed_file).unlink()

    return count
