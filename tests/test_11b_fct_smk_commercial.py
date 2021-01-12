"""
tests.test_11b_smk_commercial
=====================================
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_10b_fc_smk_commercial completed successfully so its outputs
can be used.
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

    data_ori = "tests/data/fct/crms_dnp_chembl_zinc/data/02_commercial"
    data_tgt = "tests/tmp/fct/crms_dnp_chembl_zinc/data/02_commercial"

    # delete any previous jobs from workflows
    if Path(data_tgt).is_dir():
        shutil.rmtree(data_tgt)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    if Path(data_ori).is_dir():
        shutil.copytree(str(data_ori), str(data_tgt))


def test_run():
    """Run the 'commercial' protocol applied to a subset of the ZINC ('AAAAML.xaa.sdf')."""
    # define outputs
    WD = 'tests/tmp/fct/crms_dnp_chembl_zinc/data/02_commercial'
    cids = [f"{str(cid+1).zfill(3)}" for cid in range(2)]
    output_files = [f"{WD}/molecule/data/molecule_{cid}.csv.gz" for cid in cids]
    output_svg = f"{WD}/test_commercial_zinc_tasktree.svg"
    # run protocol
    command_smk = 'run_protocol_fct commercial -c fct/crms_dnp_chembl_zinc/data/02_commercial/test_commercial_zinc.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
