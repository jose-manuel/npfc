"""
Module load
================

A module containing the Loader class, used for storing DataFrames with molecules
on disk.
"""

# standard
import logging
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'base64': lambda x: decode_mol_base64,
              }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def decode_mol_base64(string: str) -> Mol:
    """Convert a string to a RDKit Mol object.

    :param string: a string with a Mol object in bytes with a base64 string representation
    :return: a RDKit Mol object upon success, None otherwise

    """
    try:
        return Chem.Mol(base64.b64decode(string))
    except TypeError:
        return None


def from_pgsql(dbname: str,
               user: str,
               psql: str,
               db_id: str,
               db_mol: str,
               mol_format: str = None,
               col_mol: str = 'mol',
               col_id: str = 'idm',
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
    :param keep_db_cols: keep db_id and db_mol columns in output DataFrame. This does not impact any other column extracted from the psql query
    :return: the output DataFrame containing molecules as RDKit Mol objects
    """
    logging.debug(f"Retrieving data from dbname={dbname} wiht user={user} with psql:\n{psql}\nwith db_id={db_id}, db_mol={db_mol}, col_mol={col_mol}, col_id={col_id} and keep_db_cols={keep_db_cols}\n")
    if mol_format not in ('molblock', 'smiles', None):
        raise ValueError(f"Format '{mol_format}' is currently not supported. Please use either 'molblock' or 'sdf'.")
    if db_id is None:
        raise ValueError(f"db_id needs to be defined but values None.")
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
        if db_mol not in df.columns:
            raise ValueError(f"db_mol '{db_mol}' could not be found in pgsql query results.\nAvailable columns are: {list(df.columns.values)}")
        # convert mol_col to RDKit
        df[col_mol] = df[db_mol].map(CONVERTERS[mol_format])

    # cleaning up
    if col_id != db_id:
        df[col_id] = df[db_id]
    if not keep_db_cols:
        df.drop([db_id, db_mol], axis=1, inplace=True)

    return df


def from_sdf(input_sdf: str, id_sdf: str = '_Name',
             col_id: str = 'idm',  col_mol: str = 'mol',
             keep_props: bool = True, cols_list: List[str] = []) -> pd.DataFrame:
    """Load molecules from a SDF file.
    It is heavily inspired from PandasTools but differs as molecules with errors are
    not silently discarded but stored as None instead. Also, each molecule is read first
    without parsing to retrieve all properties and then with parsing to retrieve
    the molecular structure, if possible.

    :param input_sdf: the input sdf filename
    :param id_sdf: the SDF property to use for defining the col_id
    :param col_id: the column to use to generate the 'idm' column, used as single identifier
    :param keep_props: keep all properties found in the SDF. If set to False, then only 'idm' and 'mol' columns are present in the DataFrame
    :return: the output DataFrame containing molecules as RDKit Mol objects
    """
    logging.debug(f"Checking if input_sdf exists at {input_sdf}")
    # check for filename
    if not Path(input_sdf).is_file():
        raise ValueError(f"Error! File input_sdf could not be found at {input_sdf}.")
    # begin
    i = 0
    rows = []
    row_idx = []
    # double iteration over sdf but generators so it should be ok memory-wise
    for mol, mol_raw in zip(Chem.ForwardSDMolSupplier(input_sdf), Chem.ForwardSDMolSupplier(input_sdf, sanitize=False)):
        # properties
        row = dict((k, mol_raw.GetProp(k)) for k in mol_raw.GetPropNames())
        # molecule title
        if mol_raw.GetProp('_Name') != '':
            row['_Name'] = mol_raw.GetProp('_Name')
        else:
            row['_Name'] = f'NONAME_{i}'

        row[col_mol] = mol
        try:
            row[col_id] = row[id_sdf]
        except KeyError:
            logging.warning(f"No property {id_sdf} could be found for record #{i}, setting {col_id} to 'NOID_{i}'")
            row[col_id] = f"NOID_{i}"
        # delete properties
        if not keep_props:
            tmp = row.copy()
            for k in tmp.keys():
                if k != col_mol and k != col_id:
                    del row[k]
        # record entry
        rows.append(row)
        row_idx.append(i)
        i += 1
    df = pd.DataFrame(rows, index=row_idx)
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
    :return: the output DataFrame
    """
    # check for filename
    logging.debug(f"Checking if input_hdf exists at {input_hdf}")
    if not Path(input_hdf).is_file():
        raise ValueError(f"Error! File input_sdf could not be found at '{input_hdf}'.")
    # init
    key = Path(input_hdf).stem
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
             col_mol: str = 'mol', sep: str = '|',
             cols_list: List[str] = []) -> pd.DataFrame:
    """Load molecules from a CSV file.
    In case molecules were stored with encoding (base64), they need to be decoded
    before using them as RDKit molecules.

    :param input_csv: the input hdf filename
    :param decode_mols: use base64 to decode molecules encoded as strings
    :param col_mol: the DataFrame column name where molecules are stored
    :return: the output DataFrame
    """
    # check for filename
    logging.debug(f"Checking if input_hdf exists at {input_csv}")
    if not Path(input_csv).is_file():
        raise ValueError(f"File input_csv could not be found at {input_csv}")
    # read data
    df = pd.read_csv(input_csv, sep=sep, index_col='Unnamed: 0')
    logging.debug(f"Found {len(df.index)} records in input_hdf")
    if decode_mols:
        logging.debug(f"Decoding structures")
        if col_mol not in df.columns:
            raise ValueError(f"col_mol {col_mol} could not be found in DataFrame, available columns are: {list(df.columns.values)}")
        else:
            df[col_mol] = df[col_mol].map(_decode_mol_base64)
    # restaure list format
    if len(cols_list) > 0:
        for c in cols_list:
            df[c] = df[c].map(json.loads)

    return df
