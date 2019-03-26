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
from npfc import save
# logging
import logging
logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIXTURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
WD = '/tmp/npfc_test/'  # would be cool if WD could be specified so it could work on cluster too


@pytest.fixture
def df_mols():
    """Example of a DataFrame with some molecules with a phenyl."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'C1CCCCC1', 'NC1CCCC1', 'FC1CCC1', 'C1CCC(CC1)C1CCCCC1'],
                       })
    df['idm'] = [f"mol_{i}" for i in range(len(df.index))]
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    # export the mols
    print()
    print(df.columns)
    save.save(df, WD + 'test_commands_mols_in.sdf', col_id='idm')
    return df


@pytest.fixture
def df_frags():
    """Example of a DataFrame with some fragments, including a phenyl."""
    df = pd.DataFrame({'mol': ['C1CCCCC1', 'C1CCCC1'],
                       })
    df['idm'] = [f"mol_{i}" for i in range(len(df.index))]
    df['mol'] = df['mol'].map(Chem.MolFromSmiles)
    print()
    save.save(df, WD + 'test_commands_frags_in.sdf')
    return df

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init_files(df_mols, df_frags):
    pass


def test_load_mols():
    """Load molecules in SDF format and export them as base64 in a csv.gz file."""
    print()
    input_file = WD + 'test_commands_mols_in.sdf'
    output_file = WD + 'test_commands_mols_out.csv.gz'
    command = f"""load_mols {input_file} {output_file} --in_id _Name --out_id idm"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(output_file).is_file() is True


def test_load_frags():
    """Load frags in SDF format and export them as base64 in a csv.gz file."""
    print()
    input_file = WD + 'test_commands_frags_in.sdf'
    output_file = WD + 'test_commands_frags_out.csv.gz'
    command = f"""load_mols {input_file} {output_file} --in_id _Name --out_id idm"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(output_file).is_file() is True


def test_standardize_mols():
    """Standardize structures for postprocessing passed molecules."""
    print()
    input_file = WD + 'test_commands_mols_out.csv.gz'
    ref_file = WD + 'test_commands_mols_ref.hdf'
    output_passed = WD + 'test_commands_mols_out_passed.csv.gz'
    output_filtered = WD + 'test_commands_mols_out_filtered.csv.gz'
    output_error = WD + 'test_commands_mols_out_error.csv.gz'
    command = f"""standardize_mols {input_file} -r {ref_file}"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(ref_file).is_file() is True
    assert Path(output_passed).is_file() is True
    assert Path(output_filtered).is_file() is True
    assert Path(output_error).is_file() is True


def test_standardize_frags():
    """Standardize fragment structures."""
    print()
    input_file = WD + 'test_commands_frags_out.csv.gz'
    ref_file = WD + 'test_commands_frags_ref.hdf'
    output_passed = WD + 'test_commands_frags_out_passed.csv.gz'
    output_filtered = WD + 'test_commands_frags_out_filtered.csv.gz'
    output_error = WD + 'test_commands_frags_out_error.csv.gz'
    command = f"""standardize_mols {input_file} -r {ref_file}"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(ref_file).is_file() is True
    assert Path(output_passed).is_file() is True
    assert Path(output_filtered).is_file() is True
    assert Path(output_error).is_file() is True


def test_substruct_mols():
    """Run a substructure search on standardized mols and frags"""
    print()
    input_mols = WD + 'test_commands_mols_out_passed.csv.gz'
    input_frags = WD + 'test_commands_frags_out_passed.csv.gz'
    output_sub = WD + 'test_commands_frags_out_sub.csv.gz'
    command = f"""substruct_mols {input_mols} {input_frags} {output_sub}"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(output_sub).is_file() is True


def test_classify_frags_comb():
    """Classify identified fragment combinations"""
    print()
    input_mols = WD + 'test_commands_mols_out_passed.csv.gz'
    input_sub = WD + 'test_commands_frags_out_sub.csv.gz'
    output_fcc = WD + 'test_commands_frags_out_fcc.csv.gz'
    command = f"""classify_frags_comb {input_mols} {input_sub} {output_fcc}"""
    return_code = sp.call(command, shell=True)
    assert return_code == 0
    assert Path(output_fcc).is_file() is True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


if __name__ == '__main__':
    test_load_mols()
