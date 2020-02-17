"""
tests.test_08b_smk_local_natural
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

    data_ori = "tests/data/natural"
    data_tgt = "tests/tmp/natural"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_run():
    """Run the 'natural' protocol applied to a subset of the DNP."""

    output_files = [f"tests/tmp/natural/dnp/data/frags_crms/09_fgraph/data/dnp_{str(cid+1).zfill(3)}_fgraph.csv.gz" for cid in range(3)]
    output_svg = 'tests/tmp/natural/dnp/natural_dnp_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol natural --chunksize 50'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
