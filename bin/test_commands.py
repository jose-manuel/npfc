"""
Module test_commands
=====================
Tests for the npfc executables:

    - load_mols
    - standardize_mols
    - substruct_mols
    - fcc

This command uses the installed library, so output are written in /tmp folder by default.
"""
# standard
from pathlib import Path
import subprocess as sp
# data science
import pandas as pd
# chemoinformatics
from rdkit import Chem
# pytest
import pytest
# dev library
from npfc import load
from npfc import save
from npfc import standardize
from npfc import fragment
# logging
import logging
logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
WD = '/tmp/'


@pytest.fixture
def df_mols():
    """Example of a DataFrame with some molecules with a phenyl."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'C1CCCCC1', 'NC1CCCC1', 'FC1CCCCC1', '[Cl]C1CCCCC1'],
                       })
    df['idm'] = [f"mol_{i}" for i in range(len(df.index))]
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    # export the mols
    print()
    save.save(df, WD + 'test_commands_mols_in.sdf')
    return df


@pytest.fixture
def df_fragsl():
    """Example of a DataFrame with some fragments, including a phenyl."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'C1CCCC1'],
                       })
    df['idm'] = [f"mol_{i}" for i in range(len(df.index))]
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    print()
    save.save(df, WD + 'test_commands_frags_in.sdf')
    return df

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init_files(df_mols, df_fragsl):
    pass

def test_load_mols():
    print()
    command = f"""load_mols {WD + 'test_commands_mols_in.sdf'} {WD + 'test_commands_mols_out.sdf'} -n 2 --in_id idm"""
    returncode = sp.call(command, shell=True)
    print(returncode)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


if __name__ == '__main__':
    test_load_mols()
