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
# data science
import psycopg2
import pandas as pd
from pandas import DataFrame
# chemoinformatics
from rdkit import Chem
# docs
from typing import List
from typing import Union
# dev
from npfc import utils

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'rdkit': lambda x: x,  # nothing to do here, but removes the needs for more if/elif
              }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def pgsql(dbname: str,
          user: str,
          psql: str,
          src_id: str,
          src_mol: str,
          mol_format: str = None,
          col_mol: str = 'mol',
          col_id: str = 'idm',
          keep_db_cols: bool = False) -> DataFrame:
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
    df = DataFrame(cur.fetchall())
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


def file(input_file: str,
         in_id: str = 'idm',
         in_mol: str = 'mol',
         csv_sep: str = '|',
         mol_format: str = 'rdkit',
         out_id: str = 'idm',
         out_mol: str = 'mol',
         keep_props: bool = True,
         decode: Union[bool, List[str]] = True,
         ) -> DataFrame:
    """Load a file into a DataFrame.

    :param input_file: the input file to load
    :param in_id: the column/property to use for molecule ids
    :param in_mol: the column to use for molecules (irrerlevant for SDF)
    :param csv_sep: the column separator to use for parsing the input file (CSV)
    :param mol_format: the input format for molecules
    :param out_id: the column name used for storing molecule ids
    :param out_mol: the column name used for storing molecules
    :param keep_props: keep all properties found in the input file. If False, then only out_id and out_mol are kept.
    :param decode: decode base64 strings into objects. It is not implemented yet, but later it would be possible to define just column names with objects to parse. The expected behavior would be: list -> specific column names; True: all predefined columns; False: no decoding at all. Predefined columns are: mol, mol_frag, graph, colormap, aidxf, aidxf1, aidxf2, d_aidxs.
    :return: a DataFrame
    """
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
    elif format == "FEATHER":
        df = pd.read_feather(input_file)
    else:  # has to be CSV
        df = _from_csv(input_file, csv_sep=csv_sep, compression=compression)

    # process data

    # out_id
    if out_id is None and in_id is not None:
        out_id = in_id
    elif out_id != in_id:
        df.rename({in_id: out_id}, axis=1, inplace=True)

    # out_mol
    if out_mol is None and in_mol is not None:
        out_mol = in_mol
    elif out_mol != in_mol and in_mol in df.columns:
        df.rename({in_mol: out_mol}, axis=1, inplace=True)

    # keep_props
    if not keep_props:
        logging.debug(f"Found properties: {[c for c in df.columns]}")
        logging.debug(f"Keeping properties: {[out_id, out_mol]}")
        logging.debug(f"Dropping properties: {[c for c in df.columns if c not in (out_id, out_mol)]}")
        df.drop([c for c in df.columns if c not in (out_id, out_mol)], axis=1, inplace=True)

    # mol_format
    if out_mol in df.columns:
        df[out_mol] = df[out_mol].map(CONVERTERS[mol_format])

    # decode
    if decode:
        logging.debug(f"Decoding structures")
        # mols
        for col in ('mol', 'mol_frag'):
            if col in df.columns:
                df[col] = df[col].map(utils.decode_mol)
        # other objects
        for col in ('graph', 'colormap', 'aidxf', 'aidxf1', 'aidxf2', 'd_aidxs'):
            if col in df.columns:
                df[col] = df[col].map(utils.decode_object)

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

    return DataFrame(rows, index=row_idx)


def _from_hdf(input_hdf: str):
    return pd.read_hdf(input_hdf, key=Path(input_hdf).stem)


def _from_csv(input_csv: str, csv_sep: str = '|', compression: str = None):
    if compression is None:
        try:
            return pd.read_csv(input_csv, sep=csv_sep, index_col='Unnamed: 0')  # define rowidx with rowids, if any
        except KeyError:
            return pd.read_csv(input_csv, sep=csv_sep)
    elif compression == 'gzip':
        try:
            return pd.read_csv(input_csv, sep=csv_sep, index_col='Unnamed: 0', compression=compression)  # define rowidx with rowids, if any
        except (ValueError, KeyError):
            return pd.read_csv(input_csv, sep=csv_sep, compression=compression)
