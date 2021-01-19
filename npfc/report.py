"""
Module report
================

A module with helper functions for computing pre-defined plots for the analysis
of fragment combinations.
"""

import warnings
import logging
import argparse
import sys
from shutil import rmtree
from datetime import datetime
import re
from pathlib import Path
from collections import OrderedDict
# data handling
import numpy as np
import json
import pandas as pd
from pandas import DataFrame
import networkx as nx
# data visualization
from matplotlib import pyplot as plt
import seaborn as sns
# from pylab import savefig
from adjustText import adjust_text
# chemoinformatics
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Mol
from rdkit.Chem import rdChemReactions
# docs
from typing import List
from typing import Tuple
from rdkit import RDLogger
# custom libraries
import npfc
from npfc import utils
from npfc import load
from npfc import save
from npfc import fragment_combination


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# test


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class ReporterProcess:

    def __init__(self,
                 chunk_load: str,
                 chunk_std_passed: str,
                 chunk_std_filtered: str,
                 chunk_std_error: str,
                 chunk_dedupl: str,
                 WD_out: str,
                 max_examples: int = 1,
                 ):
        pass


class ReporterFragmentSearch:

    def __init__(self):
        pass

class ReporterFragmentCombination:

    def __init__(self):
        pass

class ReporterFragmentCombinationGraph:

    def __init__(self):
        pass


class ReporterPNP:

    def __init__(self):
        pass




def _parse_std_chunks(chunks: List[str]) -> DataFrame:
    """Parse all output files of a category (passed, filtered or error) for the std step and return a corresponding a results summary.

    :param chunks: output files for a category of std results
    :return: summary DF with counts
    """
    # parse all files
    dfs = []
    for c in chunks:
        df = pd.read_csv(c, sep="|", compression="gzip").groupby("task").count()[['status']].rename({'status': 'Count'}, axis=1)
        if len(df.index) > 0:
            dfs.append(df)

    # if no case was found, return an empty dataframe
    if len(dfs) == 0:
        df = pd.DataFrame([], columns=["Count", "Category"])
        return df

    # concatenate all dataframes and compute the sum of all counts
    df = pd.concat(dfs)
    df["Category"] = df.index  # I don't know how to group by index!
    df = df.groupby("Category").sum()
    df["Category"] = df.index.map(lambda x: x.replace('filter_', ''))
    # df['Perc_status'] = df['Count'].map(lambda x: f"{x/tot_mols:.2%}")

    return df


def preprocess(input_load: str,
               output_file: str,
               input_std_passed: str = None,
               input_std_filtered: str = None,
               input_std_error: str = None,
               input_dedupl: str = None,
               input_depict: str = None,
               num_examples: int = 0,
               ):
    """The information is not looked everywhere using the same logic:

        - input_load: the log file from the chunk being loaded
    """
    # load
    df = pd.read_csv(input_load, sep="@", header=None)  # char not found in the log file
    records = df[df[0].str.contains("FAILURE")].iloc[0][0].split()
    num_total = int(records[6])
    num_errors = int(records[9])
    num_passed = int(df[df[0].str.contains("SAVED")].iloc[0][0].split()[6])
    if num_total != num_errors + num_passed:
        raise ValueError(f"Error during parsing of log file: '{input_load}': {num_passed} + {num_errors} != {num_total}")
    df_load = DataFrame({'Category': ['loaded', 'cannot_load'], 'Count': [num_passed, num_errors]})

    # standardize
    df_std_passed = load.file(input_std_passed, decode=False)[['task', 'status']].groupby("task").count()[['status']].reset_index().rename({'task': 'Category', 'status': 'Count'}, axis=1)
    df_std_filtered = load.file(input_std_filtered, decode=False)[['task', 'status']].groupby("task").count()[['status']].rename({'status': 'Count'}, axis=1)
    df_std_error = load.file(input_std_error, decode=False)[['task', 'status']].groupby("task").count()[['status']].rename({'status': 'Count'}, axis=1)

    # dedupl
    df = pd.read_csv(input_dedupl, sep="@", header=None)  # char not found in the log file so we can extract all lines as one column
    num_passed, num_total = [int(x) for x in df[df[0].str.contains("REMAINING MOLECULES")].iloc[0][0].split("MOLECULES:")[1].split("/")]
    num_filtered = num_total - num_passed
    df_dedupl = pd.DataFrame({'Category': ['unique', 'duplicate'], 'Count': [num_passed, num_filtered]})

    # depict
    pass  # nothing to do here at the moment since I never saw any problem at this step

    # merge data


def get_df_dedupl(WD: str) -> DataFrame:
    """Get a DF summarizing the results of the deduplication step.

    :param WD: the directory of the std step
    :return: a DF summarizing results of the deduplication step
    """
    logger.info("PREP -- COMPUTING DEDUPL RESULTS")
    # iterate over the log files to count status
    pattern = ".*([0-9]{3})?_dedupl.log"
    chunks = _get_chunks(f"{WD}/log", pattern)
    chunks = [c for c in chunks if c.split('.')[-1] == 'log']  # ###### quick and dirty
    # print(f"{chunks=}")
    logger.info(f"PREP -- FOUND {len(chunks):,d} CHUNKS")
    # initiate counts
    num_tot = 0
    num_passed = 0
    num_filtered = 0
    for c in chunks:
        df = pd.read_csv(c, sep="@", header=None)  # char not found in the log file so we can extract all lines as one column
        passed, total = [int(x) for x in df[df[0].str.contains("REMAINING MOLECULES")].iloc[0][0].split("MOLECULES:")[1].split("/")]
        num_passed += passed
        num_filtered += total - passed
        num_tot += total
    # create a dataframe with counts
    df = pd.DataFrame({'Category': ['unique', 'duplicate'], 'Count': [num_passed, num_filtered]})
    logger.info(f"PREP -- RESULTS FOR DEDUPL:\n\n{df}\n")

    return df


def get_dfs_prep_frags(WD: str) -> Tuple[DataFrame]:
    """
    Same as dfs_prep but adapted to the fragments protocol.
    The latter differs from other protocols in several ways:

        - no deglycoslyation step
        - 2 competiting murcko scaffold extraction schemes:
            - A: load > murcko > std
            - B: load > std > murcko > stdms
        - concatenation of A and B to prioritize B over A

        (with std without filters)
    """

    p = Path(WD)
    WD_LOAD = [str(x) for x in list(p.glob("*_load"))][0]
    WD_STD = [str(x) for x in list(p.glob("*_std*"))]
    print(WD_STD)
    WD_DEDUPL = [str(x) for x in list(p.glob("*_dedupl"))][0]

    sys.exit(0)


def get_dfs_prep(WD: str) -> Tuple[DataFrame]:
    """Get a list of DFs summarizing the whole preprocess superstep: load, deglyco, std and dedupl.

    - DF_deglyco is the detailed summary of deglycosylation appended with the number of mols that did not get processed because they could not be loaded (NA).
    - DF_prep_filtered is the detailed summary of std and dedupl
    - DF_prep_error is the detailed summary of std and load
    - DF_prep_all is the general summary with the final number of passed, filtered and error molecules.

    :param WD: the main directory of the dataset data (i.e. 'natural/dnp/data')
    :return: a list of DFs of interest: [DF_deglyco, DF_prep_filtered, DF_prep_error, DF_prep_all]
    """

    logger.info("PREP -- COMPUTE RESULTS FOR PREPROCESS")
    logger.info("PREP -- PROPRESS CONTAINS LOAD, DEGLYCO, STD AND DEDUPL STEPS")

    # define subfolders
    p = Path(WD)
    WD_LOAD = [str(x) for x in list(p.glob("*_load"))][0]
    # there is no deglyco step anymore for fragments
    WD_DEGLYCO = [str(x) for x in list(Path(WD).glob("*_deglyco"))]
    if len(WD_DEGLYCO) < 1:
        logging.debug(f"NO DEGLYCO STEP COULD BE FOUND AT '{WD}'")
        WD_DEGLYCO = ''
    else:
        WD_DEGLYCO = WD_DEGLYCO[0]
    WD_STD = [str(x) for x in list(p.glob("*_std"))][0]
    WD_DEDUPL = [str(x) for x in list(p.glob("*_dedupl"))][0]

    # get dfs
    df_load = get_df_load(WD_LOAD)
    df_deglyco = get_df_deglyco(WD_DEGLYCO)
    df_std_passed = get_df_std_passed(WD_STD)
    df_std_filtered = get_df_std_filtered(WD_STD)
    df_std_error = get_df_std_error(WD_STD)
    df_dedupl = get_df_dedupl(WD_DEDUPL)

    # get total of molecules in input
    num_mols_tot = df_load['Count'].sum()
    if WD_DEGLYCO != '':
        # count not loaded molecules as well in deglyco
        num_mols_deglyco_na = df_load[df_load['Category'] == 'cannot_load']['Count']
        df_deglyco = pd.concat([df_deglyco, pd.DataFrame({'Category': ['NA'], 'Count': [num_mols_deglyco_na]})])
        df_deglyco.reset_index(inplace=True, drop=True)
        df_deglyco['Count'] = df_deglyco['Count'].astype(int)
        df_deglyco['Perc_Status'] = df_deglyco['Count'].map(lambda x: f"{x/num_mols_tot:.2%}")
        logger.info(f"PREP -- RESULTS FOR DEGLYCOSYLATION:\n\n{df_deglyco}\n")
    else:
        df_deglyco = None

    # gather all filtered molecules
    df_dedupl_dupl = df_dedupl[df_dedupl['Category'] == 'duplicate']
    num_dedupl_dupl = df_dedupl_dupl['Count'].sum()
    df_std_filtered = pd.concat([df_std_filtered, df_dedupl_dupl], sort=True)
    # count even unoccurred cases in df_std_filtered
    filters = ['empty', 'hac', 'molweight', 'nrings', 'medchem', 'timeout', 'duplicate']
    df_std_filtered.set_index('Category', inplace=True)
    df_std_filtered = df_std_filtered.reindex(filters)
    df_std_filtered.reset_index(inplace=True)
    df_std_filtered.fillna(0, inplace=True)
    df_std_filtered['Count'] = df_std_filtered['Count'].astype(int)
    df_std_filtered['Perc_Status'] = df_std_filtered['Count'].map(lambda x: f"{x/num_mols_tot:.2%}")
    logger.info(f"PREP -- RESULTS FOR STD_FILTERED:\n\n{df_std_filtered}\n")

    # gather all molecules that raised an error
    df_std_error = pd.concat([df_std_error, df_load[df_load['Category'] == 'cannot_load']], sort=True)
    # count even unoccurred cases in df_std_error
    errors = ['cannot_load', 'initiate_mol', 'disconnect_metal', 'sanitize', 'remove_isotopes', 'normalize', 'uncharge', 'canonicalize', 'remove_stereo']
    df_std_error.set_index('Category', inplace=True)
    df_std_error = df_std_error.reindex(errors)
    df_std_error.reset_index(inplace=True)
    df_std_error.fillna(0, inplace=True)
    df_std_error['Count'] = df_std_error['Count'].astype(int)
    df_std_error['Perc_Status'] = df_std_error['Count'].map(lambda x: f"{x/num_mols_tot:.2%}")
    logger.info(f"PREP -- RESULTS FOR STD_ERRORS:\n\n{df_std_error}\n")

    # general count for passed/filtered/errors
    num_tot_filtered = df_std_filtered['Count'].sum()
    num_tot_passed = df_std_passed['Count'].sum() - num_dedupl_dupl  # dedupl happens after std, so std contains passsed mols that get filtered
    num_tot_errors = df_std_error['Count'].sum()
    df_std_all = pd.DataFrame({'Category': ['passed', 'filtered', 'errors'], 'Count': [num_tot_passed, num_tot_filtered, num_tot_errors]})
    df_std_all['Perc_Status'] = df_std_all['Count'].map(lambda x: f"{x/num_mols_tot:.2%}")
    logger.info(f"PREP -- RESULTS FOR STD_ALL:\n\n{df_std_all}\n")

    return {'df_prep_overview': df_std_all, 'df_prep_filtered': df_std_filtered, 'df_prep_error': df_std_error, 'df_prep_deglyco': df_deglyco}
