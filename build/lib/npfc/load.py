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
import subprocess
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
import warnings


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


CONVERTERS = {'molblock': lambda x: Chem.MolFromMolBlock(x),
              'smiles': lambda x: Chem.MolFromSmiles(x),
              'rdkit': lambda x: x,  # nothing to do here, but removes the needs for more if/elif
              }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


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
    :param decode: decode base64 strings into objects. Columns with encoded objects are labelled with a leading '_'. For molecules, reserved names are 'mol' and 'mol_frag'.
    :return: a DataFrame
    
    ..warning:: if a 'idm' property exists in the input file but the user picks another property for in_id, the pre-existing 'idm' will be renamed into 'idm_2'. 
    """
    # check arguments
    utils.check_arg_input_file(input_file)
    format, compression = utils.get_file_format(input_file)
    logging.debug("Loading file '%s' (format=%s, compression=%s)", input_file, format, compression)

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
    
    # in case in_id is different than default but some column with the same name as out_id already exists
    if in_id != out_id and out_id in df.columns:
        df.rename({out_id: f"{out_id}.2"}, axis=1, inplace=True)
        logging.debug('in_id=%s, out_id=%s, but out_id already present. Renaming existing out_id into %s', in_id, out_id, f"{out_id}.2.")

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
        logging.debug("Found properties: %s", ', '.join([c for c in df.columns]))
        logging.debug("Keeping properties: %s", ', '.join([out_id, out_mol]))
        logging.debug("Dropping properties: %s", ', '.join([c for c in df.columns if c not in (out_id, out_mol)]))
        df.drop([c for c in df.columns if c not in (out_id, out_mol)], axis=1, inplace=True)

    # mol_format
    if out_mol in df.columns:
        df[out_mol] = df[out_mol].map(CONVERTERS[mol_format])

    # decode
    if decode:
        # faster way for mols
        for col in ('mol', 'mol_frag', 'mol_frag_1', 'mol_frag_2', 'mol_rdkit'):
            if col in df.columns:
                logging.debug("Decoding column='%s'", col)
                df[col] = df[col].map(utils.decode_mol)

        # other objects are labelled with leading '_'
        for col in df.columns:
            if col.startswith('_') and col != '_Name':
                logging.debug("Decoding column='%s'", col)
                df[col] = df[col].map(utils.decode_object)

    # consistent id comparison accross datasets
    # columns to convert to str
    cols_str = [out_id, 'idf1', 'idf2']
    for col in cols_str:
        if col in df.columns:
            df[col] = df[col].astype(str)

    logging.debug("First 3 rows loaded:\n\n%s\n", df.head(3))

    return df


def _from_sdf(input_sdf: str, col_mol: str = 'mol', compression: str = None):
    if compression is None:
        logging.debug("Read sdf from uncompressed file")
        FH = open(input_sdf, 'rb')
        FH_raw = open(input_sdf, 'rb')
    elif compression == 'gzip':
        logging.debug("Read sdf from compressed file (%s)", compression)
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
        return pd.read_csv(input_csv, sep=csv_sep)
    elif compression == 'gzip':
        return pd.read_csv(input_csv, sep=csv_sep, compression=compression)


def count_mols(input_file: str, keep_uncompressed: bool = False):
    """Count the number of molecules found in a text file.
    In case the file is compressed (gzip), it is uncompressed first. The resulting
    uncompressed file can be kept for further use.

    #### this function failed to the ZINC (9,902,598)

    :param input_file: the input file
    :param keep_uncompressed: if the input file is compressed (gzip), do not remove the uncompressed file when finished
    """

    utils.check_arg_input_file(input_file)
    format, compression = utils.get_file_format(input_file)
    logging.debug("input_file='%s' (format: %s, compression: %s)", input_file, format, compression)

    # in case of compressed file, uncompress it
    if compression == 'gzip':
        uncompressed_file = input_file.split('.gz')[0]
        with open(uncompressed_file, 'wb') as uncompressed, gzip.open(input_file, 'rb') as compressed:
            bindata = compressed.read()
            uncompressed.write(bindata)
        input_file = uncompressed_file

    # count mols
    if format == 'SDF':
        # look for the '$$$$' pattern in the SDF file
        result = subprocess.check_output(f"grep -c '$$$$' {input_file}", shell=True).decode('utf-8').strip()
    elif format == 'HDF':
        # get the size of the dataframe
        result = len(file(input_file, decode=False).index)
    else:
        # any text-based file with 1 molecule per row (csv, smiles, inchi, etc.)
        result = subprocess.check_output(f"wc -l {input_file}", shell=True).decode('utf-8').split()[0].strip()

    # read count
    count = int(result)

    # remove uncompressed file
    if compression and not keep_uncompressed:
        Path(uncompressed_file).unlink()

    return count
