"""
tests.test_08a_smk_local_fragments
~~~~~~~~~~
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_08a_smk_local_fragments completed successfully so its outputs
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

    data_ori = "tests/data/fc/01_fragments"
    data_tgt = "tests/tmp/fc/01_fragments"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_run():
    """Run the 'fragments' protocol applied to the cr dataset."""

    output_files = ["tests/tmp/fc/01_fragments/crms/data/prep/05_fcp/data/crms_fcp.csv.gz"]
    output_svg = 'tests/tmp/fc/01_fragments/crms/fragments_crms_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol_fc fragments -c fc/01_fragments/crms/test_fragments_crms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
