"""
Module load
================

A module containing the Loader class, used for storing DataFrames with molecules
on disk.
"""

# standard
import logging
import gzip
from pathlib import Path
import base64
import json
# data science
import psycopg2
import pandas as pd
# chemoinformatics
from rdkit import Chem
from rdkit.Chem import Mol
# docs
from typing import List
# dev
from npfc import utils

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'rdkit': lambda x: x,  # nothing to do here, but removes the needs for more if/elif
              }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def decode_mol_base64(string: str) -> Mol:
    """Convert a string to a RDKit Mol object.

    :param string: a string with a Mol object in bytes with a base64 string representation
    :return: a Mol object upon success, None otherwise

    """
    try:
        return Chem.Mol(base64.b64decode(string))
    except TypeError:
        return None


def from_pgsql(dbname: str, user: str,
               psql: str, src_id: str,
               src_mol: str, mol_format: str = None,
               col_mol: str = 'mol', col_id: str = 'idm',
               keep_db_cols: bool = False) -> pd.DataFrame:
    """Load molecules from a PGSQL query.
    The col_mol will is parsed by RDKit depending on the mol_format argument.
    If no mol_format is specified, then the column is returned untouched.

    .. note:: For this function to work, you need to be able to log into the target database without prompt. Tested only with ChEMBL24 with default ports.

    :param dbname: the input postgres database name
    :param user: the user name for logging into database
    :param pgsql: the posgresql statement to execute
    :param mol_format: the molecular format to use to parse the molecules. If none is specified then no parsing is performed. Currently only molblock and smiles are allowed
    :param col_mol: the column with the molecules to parse
    :param keep_db_cols: keep src_id and src_mol columns in output DataFrame. This does not impact any other column extracted from the psql query
    :return: a DataFrame with Mol objects
    """
    logging.debug(f"Retrieving data from dbname={dbname} wiht user={user} with psql:\n{psql}\nwith src_id={src_id}, src_mol={src_mol}, col_mol={col_mol}, col_id={col_id} and keep_db_cols={keep_db_cols}\n")
    if mol_format not in ('molblock', 'smiles', None):
        raise ValueError(f"Format '{mol_format}' is currently not supported. Please use either 'molblock' or 'sdf'.")
    if src_id is None:
        raise ValueError(f"src_id needs to be defined but values None.")
    # establish connection
    conn = psycopg2.connect(dbname=dbname, user=user)
    cur = conn.cursor()
    # execute query
    cur.execute(psql)
    # generate DataFrame
    df = pd.DataFrame(cur.fetchall())
    # apply table column names to DataFrame
    df.columns = [c[0] for c in cur.description]
    # check for columns
    if mol_format is not None:
        if src_mol not in df.columns:
            raise ValueError(f"src_mol '{src_mol}' could not be found in pgsql query results.\nAvailable columns are: {list(df.columns.values)}")
        # convert mol_col to RDKit
        df[col_mol] = df[src_mol].map(CONVERTERS[mol_format])

    # cleaning up
    if col_id != src_id:
        df[col_id] = df[src_id]
    if not keep_db_cols:
        df.drop([src_id, src_mol], axis=1, inplace=True)

    return df


def from_file(input_file: str,
              in_id: str = 'idm',
              in_mol: str = 'mol',
              in_sep: str = '|',
              mol_format: str = 'rdkit',
              out_id: str = 'idm',
              out_mol: str = 'mol',
              keep_props: bool = True,
              decode: bool = True,
              cols_list: List[str] = [],
              ):
    # check arguments
    utils.check_arg_input_file(input_file)
    format, compression = utils.get_file_format(input_file)
    logging.debug(f"Format: {format}, compression: {compression}")
    # read mols
    if format == 'SDF':
        df = _from_sdf(input_file, compression=compression)
        logging.debug("Decode is set to False")
        decode = False
    elif format == 'HDF':
        df = _from_hdf(input_file)
    else:  # has to be CSV
        df = _from_csv(input_file, in_sep=in_sep, compression=compression)

    # process data

    # out_id
    if out_id is None:
        out_id = in_id
    elif out_id != in_id:
        df.rename({in_id: out_id}, axis=1, inplace=True)

    # out_mol
    if out_mol is None:
        out_mol = in_mol
    elif out_mol != in_mol:
        df.rename({in_mol: out_mol}, axis=1, inplace=True)

    # keep_props
    if not keep_props:
        logging.debug(f"Found properties: {[c for c in df.columns]}")
        logging.debug(f"Keeping properties: {[out_id, out_mol]}")
        logging.debug(f"Dropping properties: {[c for c in df.columns if c not in (out_id, out_mol)]}")
        df.drop([c for c in df.columns if c not in (out_id, out_mol)], axis=1, inplace=True)

    # mol_format
    df[out_mol] = df[out_mol].map(CONVERTERS[mol_format])

    # decode
    if decode:
        logging.debug(f"Decoding structures")
        df[in_mol] = df[in_mol].map(decode_mol_base64)

    return df


def _from_sdf(input_sdf: str, col_mol: str = 'mol', compression: str = None):
    if compression is None:
        logging.debug("Read sdf from uncompressed file")
        FH = open(input_sdf, 'rb')
        FH_raw = open(input_sdf, 'rb')
    elif compression == 'gzip':
        logging.debug(f"Read sdf from compressed file ({compression})")
        FH = gzip.open(input_sdf)
        FH_raw = gzip.open(input_sdf)
    else:
        raise ValueError(f"Error! Unknown compression type for SDF: '{compression}'")
    # init
    i = 0
    rows = []
    row_idx = []
    # double iteration over sdf but generators so it should be ok memory-wise
    for mol, mol_raw in zip(Chem.ForwardSDMolSupplier(FH), Chem.ForwardSDMolSupplier(FH_raw, sanitize=False)):
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

    return pd.DataFrame(rows, index=row_idx)


def _from_hdf(input_hdf: str):
    return pd.read_hdf(input_hdf, key=Path(input_hdf).stem)


def _from_csv(input_csv: str, in_sep: str = '|', compression: str = None):
    if compression is None:
        try:
            return pd.read_csv(input_csv, sep=in_sep, index_col='Unnamed: 0')  # define rowidx with rowids, if any
        except KeyError:
            return pd.read_csv(input_csv, sep=in_sep)
    elif compression == 'gzip':
        try:
            return pd.read_csv(input_csv, sep=in_sep, index_col='Unnamed: 0', compression=compression)  # define rowidx with rowids, if any
        except KeyError:
            return pd.read_csv(input_csv, sep=in_sep, compression=compression)

# def from_csv(input_csv: str,
#              mol_format: str = 'rdkit',
#              decode_mols: bool = True,
#              col_mol: str = 'mol',
#              col_id: str = 'idm',
#              sep: str = '|',
#              keep_props: bool = True,
#              cols_list: List[str] = []) -> pd.DataFrame:
#     """Load molecules from a CSV file.
#     In case molecules were stored with encoding (base64), they need to be decoded
#     before using them as RDKit molecules.
#
#     :param input_csv: the input hdf filename
#     :param decode_mols: use base64 to decode molecules encoded as strings
#     :param col_mol: the DataFrame column name where molecules are stored
#     :return: a DataFrame with Mol objects
#     """
#     # check for filename
#     logging.debug(f"Checking if input_hdf exists at {input_csv}")
#     utils.check_arg_input_file(input_csv)
#     # read file according to its compression
#     suffixes = Path(input_csv).suffixes
#     if len(suffixes) == 1:
#         # read data
#         df = pd.read_csv(input_csv, sep=sep, index_col='Unnamed: 0')
#     else:
#         compression = utils.get_conversion(suffixes[1])
#         try:
#             df = pd.read_csv(input_csv, sep=sep, index_col='Unnamed: 0', compression=compression)  # in case csv comes with rowidx, use rowidx
#         except KeyError:
#             df = pd.read_csv(input_csv, sep=sep, compression=compression)
#     logging.debug(f"Found {len(df.index)} records in input_hdf")
#     # process data
#     df[col_mol] = df[col_mol].map(CONVERTERS[format])
#     num_failed = df[col_mol].isna().sum()
#     if num_failed > 0:
#         logging.warning(f"{num_failed} structures could not be initialized from format '{format}'")
#     if decode_mols:
#         logging.debug(f"Decoding structures")
#         if col_mol not in df.columns:
#             raise ValueError(f"col_mol {col_mol} could not be found in DataFrame, available columns are: {list(df.columns.values)}")
#         else:
#             df[col_mol] = df[col_mol].map(decode_mol_base64)
#     # restaure list format
#     if len(cols_list) > 0:
#         for c in cols_list:
#             df[c] = df[c].map(json.loads)
#     # keep_props
#     if not keep_props:
#         df.drop([c for c in df.columns if c not in (col_id, col_mol)], axis=1, inplace=True)
#
#     return df
#
#
# def from_file(input_file: str,
#               decode_mols: bool = True,
#               src_id: str = '_Name',
#               col_mol: str = 'mol',
#               col_id: str = 'idm',
#               sep: str = '|',
#               keep_props: bool = True,
#               cols_list: List[str] = []):
#     """Load molecules from a file with molecules (i.e. SDF, HDF or CSV).
#     In case molecules were stored with encoding (base64), they need to be decoded
#     before using them as RDKit molecules.
#
#     :param input_file: the input filename (sdf, hdf or csv)
#     :param decode_mols: use base64 to decode molecules encoded as strings
#     :param col_mol: the DataFrame column name where molecules are stored
#     :return: a DataFrame with Mol objects
#     """
#     # check input file
#     utils.check_arg_input_file(input_file)
#     # check format and compression
#     input_format = Path(input_file).suffixes[0]
#
#     if input_format == '.sdf':
#         df_mols = from_sdf(input_file,
#                            src_id=src_id,
#                            col_mol=col_mol,
#                            col_id=col_id,
#                            keep_props=keep_props,
#                            )
#     elif input_format == '.csv':
#         df_mols = from_csv(input_file,
#                            col_id=col_id,
#                            col_mol=col_mol,
#                            keep_props=keep_props,
#                            decode_mols=decode_mols,
#                            sep=sep,
#                            )
#     else:   # check on argument for input does not leave any other option than sdf, csv or hdf
#         df_mols = from_hdf(input_file,
#                            decode_mols=decode_mols,
#                            col_mol=col_mol,
#                            )
#
#         return df_mols
