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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'base64': lambda x: decode_mol_base64,
              }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def decode_mol_base64(string: str) -> Mol:
    """Convert a string to a RDKit Mol object."""
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

    :param dbname: The input postgres database name.
    :param user: The user name for logging into database.
    :param pgsql: The posgresql statement to execute.
    :param mol_format: The molecular format to use to parse the molecules. If none is specified then no parsing is performed. Currently only molblock and smiles are allowed.
    :param col_mol: The column with the molecules to parse.
    :param keep_db_cols: Keep db_id and db_mol columns in output DataFrame. This does not impact any other column extracted from the psql query.
    :return df: The output DataFrame containing molecules as RDKit Mol objects.
    :rtype: pandas.DataFrame
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
