"""
tests.test_11c_fct_smk_natural
======================================
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_10c_fc_smk_natural completed successfully so its outputs
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

    data_ori = "tests/data/fct/crms_dnp_chembl_zinc/data/03_natural"
    data_tgt = "tests/tmp/fct/crms_dnp_chembl_zinc/data/03_natural"

    # delete any previous jobs from workflows
    if Path(data_tgt).is_dir():
        shutil.rmtree(data_tgt)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    if Path(data_ori).is_dir():
        shutil.copytree(str(data_ori), str(data_tgt))


def test_run():
    """Run the 'natural' protocol applied to a subset of the DNP."""
    # define outputs
    WD = 'tests/tmp/fct/crms_dnp_chembl_zinc/data/03_natural'
    cids = [f"{str(cid+1).zfill(3)}" for cid in range(2)]
    output_files = [f"{WD}/dataset/data/dataset.csv.gz"] + \
                   [f"{WD}/fragment_fragment/data/fragment_fragment_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/molecule/data/molecule_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/molecule_molecule/data/molecule_molecule_{cid}.csv.gz" for cid in cids]

    output_svg = f"{WD}/test_natural_dnp_fragscrms_tasktree.svg"
    # run protocol
    command_smk = 'run_protocol_fct natural -c fct/crms_dnp_chembl_zinc/data/03_natural/test_natural_dnp_fragscrms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
