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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
