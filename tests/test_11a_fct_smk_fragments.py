"""
tests.test_11a_fct_smk_fragments
=====================================
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_10a_fc_smk_fragments completed successfully so its outputs
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

    data_ori = "tests/data/fct/crms_dnp_chembl_zinc/data/01_fragments"
    data_tgt = "tests/tmp/fct/crms_dnp_chembl_zinc/data/01_fragments"
    # delete any previous jobs from workflows
    if Path(data_tgt).is_dir():
        shutil.rmtree(data_tgt)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    if Path(data_ori).is_dir():
        shutil.copytree(str(data_ori), str(data_tgt))


def test_run():
    """Run the 'fragments' protocol applied to the cr dataset."""

    WD = 'tests/tmp/fct/crms_dnp_chembl_zinc/data/01_fragments'
    output_files = [f"{WD}/dataset/data/dataset.csv.gz",
                    f"{WD}/fragment/data/fragment.csv.gz",
                    f"{WD}/molecule/data/molecule.csv.gz",
                    f"{WD}/molecule_dataset/data/molecule_dataset.csv.gz",
                    f"{WD}/molecule_molecule/data/molecule_molecule.csv.gz",
                    ]
    output_svg = f"{WD}/test_fragments_crms_tasktree.svg"
    # run protocol
    command_smk = 'run_protocol_fct fragments -c fct/crms_dnp_chembl_zinc/data/01_fragments/test_fragments_crms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
