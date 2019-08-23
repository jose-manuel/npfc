"""
tests.test_06_local_workflow
~~~~~~~~~~
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.
"""
# standard library
import subprocess
from math import ceil
from pathlib import Path
import shutil
import pkg_resources
# dev
from npfc import load
# debug
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_00_init_folder():
    """Reset the folder where outputs are computed to its initial state."""

    data_ori = "tests/data/"
    data_tgt = "tests/tmp/"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_01_fragments():
    """Run the 'fragments' protocol applied to the cr dataset."""

    output_files = ["tests/tmp/fragments/crms/data/08_gen2D/data/crms_gen2D.csv.gz"]
    output_svg = 'tests/tmp/fragments/crms/crms_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol fragments --chunksize 100'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])


def test_02_natural():
    """Run the 'natural' protocol applied to a subset of the DNP."""

    output_files = [f"tests/tmp/natural/dnp/data/09_fmap/data/dnp_{str(cid+1).zfill(3)}_fmap.csv.gz" for cid in range(3)]
    output_svg = 'tests/tmp/natural/dnp/dnp_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol natural --chunksize 100'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])


def test_03_synthetic():
    """Run the 'synthetic' protocol applied to a subset of the ChEMBL."""

    output_files = [f"tests/tmp/synthetic/chembl/data/12_pnp/data/chembl_{str(cid+1).zfill(3)}_pnp.csv.gz" for cid in range(3)]
    output_svg = 'tests/tmp/synthetic/chembl/chembl_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol synthetic --chunksize 100'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
