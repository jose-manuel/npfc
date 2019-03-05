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
              'base64': lambda x: decode_mol_base64,
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


def from_sdf(input_sdf: str, src_id: str = '_Name',
             col_id: str = 'idm',  col_mol: str = 'mol',
             keep_props: bool = True, cols_list: List[str] = []) -> pd.DataFrame:
    """Load molecules from a SDF file.
    It is heavily inspired from PandasTools but differs as molecules with errors are
    not silently discarded but stored as None instead. Also, each molecule is read first
    without parsing to retrieve all properties and then with parsing to retrieve
    the molecular structure, if possible.

    .. warning:: If keep_props is set to True, then only the last molecule is used for determing what properties should be retrieved. This could be an issue if the input SDF has unconsistent properties.

    :param input_sdf: the input sdf filename
    :param src_id: the SDF property to use for defining the col_id
    :param col_id: the column to use to generate the 'idm' column, used as single identifier
    :param keep_props: keep all properties found in the SDF. If set to False, then only 'idm' and 'mol' columns are present in the DataFrame
    :return: a DataFrame with Mol objects
    """
    logging.debug(f"Checking if input_sdf exists at {input_sdf}")
    # check for filename
    utils.check_arg_input_file(input_sdf)
    # get suffixes
    suffixes = Path(input_sdf).suffixes
    # get the correct file handlers
    if len(suffixes) == 1:
        FH = open(input_sdf, 'rb')
        FH_raw = open(input_sdf, 'rb')
    else:
        if utils.get_conversion(suffixes[1]) == 'gzip':
            FH = gzip.open(input_sdf)
            FH_raw = gzip.open(input_sdf)
        else:
            raise ValueError(f"Error! Unknown compression type: '{suffixes}'.")
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
            if mol_raw.GetProp('_Name') != '':
                row['_Name'] = mol_raw.GetProp('_Name')
            else:
                row['_Name'] = f'NONAME_{i}'
            # mol
            row[col_mol] = mol
            # src_id
            try:
                row[col_id] = row[src_id]
            except KeyError:
                logging.warning(f"No property {src_id} could be found for record #{i}, setting {col_id} to 'NOID_{i}'")
                row[col_id] = f"NOID_{i}"
            # delete properties
            if not keep_props:
                tmp = row.copy()
                for k in tmp.keys():
                    if k != col_mol and k != col_id:
                        del row[k]
        except AttributeError:
            logging.warning("Molecule #{i} could not be parsed and was skipped!")
            row = None
        # record entry
        if row is not None:
            rows.append(row)
            row_idx.append(i)
        i += 1
    df = pd.DataFrame(rows, index=row_idx, columns=[col_mol, col_id] + [k for k in row.keys() if k not in (col_mol, col_id)])
    # restaure list format
    if len(cols_list) > 0:
        for c in cols_list:
            df[c] = df[c].map(json.loads)

    return df


def from_hdf(input_hdf: str, decode_mols: bool = True, col_mol: str = 'mol')-> pd.DataFrame:
    """Load molecules from a HDF file.
    In case molecules were stored with encoding (base64), they need to be decoded
    before using them as RDKit molecules.

    :param input_hdf: the input hdf filename
    :param decode_mols: use base64 to decode molecules encoded as strings
    :param col_mol: the DataFrame column name where molecules are stored
    :return: a DataFrame with Mol objects
    """
    # check for filename
    logging.debug(f"Checking if input_hdf exists at {input_hdf}")
    if not Path(input_hdf).is_file():
        raise ValueError(f"Error! File input_sdf could not be found at '{input_hdf}'.")
    # init
    key = Path(input_hdf).stem.split('.')[0]
    logging.debug(f"Key for hdf is {key}")
    # read data
    df = pd.read_hdf(input_hdf, key=key)
    logging.debug(f"Found {len(df.index)} records in input_hdf")
    if decode_mols:
        logging.debug(f"Decoding structures")
        if col_mol not in df.columns:
            raise ValueError(f"Error! Parameter col_mol {col_mol} could not be found in DataFrame, available columns are: '{list(df.columns.values)}'.")
        else:
            df[col_mol] = df[col_mol].map(decode_mol_base64)

    return df


def from_csv(input_csv: str, decode_mols: bool = True,
             col_mol: str = 'mol', col_id: str = 'idm',
             sep: str = '|',, keep_props: bool = True,
             cols_list: List[str] = []) -> pd.DataFrame:
    """Load molecules from a CSV file.
    In case molecules were stored with encoding (base64), they need to be decoded
    before using them as RDKit molecules.

    :param input_csv: the input hdf filename
    :param decode_mols: use base64 to decode molecules encoded as strings
    :param col_mol: the DataFrame column name where molecules are stored
    :return: a DataFrame with Mol objects
    """
    # check for filename
    logging.debug(f"Checking if input_hdf exists at {input_csv}")
    utils.check_arg_input_file(input_csv)
    # read file according to its compression
    suffixes = Path(input_csv).suffixes
    if len(suffixes) == 1:
        # read data
        df = pd.read_csv(input_csv, sep=sep, index_col='Unnamed: 0')
    else:
        compression = utils.get_conversion(suffixes[1])
        df = pd.read_csv(input_csv, sep=sep, index_col='Unnamed: 0', compression=compression)
    logging.debug(f"Found {len(df.index)} records in input_hdf")
    # process data
    if decode_mols:
        logging.debug(f"Decoding structures")
        if col_mol not in df.columns:
            raise ValueError(f"col_mol {col_mol} could not be found in DataFrame, available columns are: {list(df.columns.values)}")
        else:
            df[col_mol] = df[col_mol].map(decode_mol_base64)
    # restaure list format
    if len(cols_list) > 0:
        for c in cols_list:
            df[c] = df[c].map(json.loads)
    # keep_props
    if not keep_props:
        df.drop([c for c in df.columns if c not in (col_id, col_mol)], axis=1, inplace=True)

    return df
