"""
tests.test_08c_smk_local_commercial
=====================================
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.
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

    data_ori = "tests/data/commercial"
    data_tgt = "tests/tmp/commercial"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_run():
    """Run the 'commercial' protocol applied to a subset of the ZINC ('AAAAML.xaa.sdf')."""

    output_files = [f"tests/tmp/commercial/zinc/data/prep/05_depict/data/zinc_{str(cid+1).zfill(3)}_depict.csv.gz" for cid in range(2)]
    output_svg = 'tests/tmp/commercial/zinc/commercial_zinc_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol commercial -c commercial/zinc/test_zinc.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
