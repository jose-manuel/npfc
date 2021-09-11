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

    data_ori = "tests/data/fc/03_natural"
    data_tgt = "tests/tmp/fc/03_natural"
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

    output_files = [f"tests/tmp/fc/03_natural/coconut/data/prep/frags_crms/08_fcg/data/coconut_{str(cid+1).zfill(3)}_fcg.csv.gz" for cid in range(1)]
    output_svg = 'tests/tmp/fc/03_natural/coconut/natural_coconut_tasktree.svg'
    # run protocol
    command_smk = 'run_protocol_fc natural -c fc/03_natural/coconut/test_natural_coconut_fragscrms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
