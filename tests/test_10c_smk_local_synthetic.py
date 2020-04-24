"""
tests.test_08c_smk_local_synthetic
~~~~~~~~~~
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_08a_smk_local_fragments and test_08b_smk_local_natural
completed successfully so their outputs can be used.
"""
# standard library
import subprocess
from pathlib import Path
import shutil
# debug
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init():
    """Reset the folder where outputs are computed to its initial state."""

    data_ori = "tests/data/synthetic"
    data_tgt = "tests/tmp/synthetic"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_run():
    """Run the 'synthetic' protocol applied to a subset of the ChEMBL."""

    output_files = [f"tests/tmp/synthetic/chembl/data/natref_dnp/frags_crms/10_pnp/data/chembl_{str(cid+1).zfill(3)}_pnp.csv.gz" for cid in range(3)]
    output_svg = 'tests/tmp/synthetic/chembl/synthetic_chembl_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol synthetic --chunksize 20'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
