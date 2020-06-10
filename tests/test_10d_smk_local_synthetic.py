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

    # delete any previous jobs from workflows
    data_tgt = "tests/tmp/synthetic/chembl/prep/natref_dnp/frags_crms_nobn"
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)    

    data_tgt = "tests/tmp/natural/dnp/prep/frags_crms_nobn"
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir) 


def test_run():
    """Run the 'synthetic' protocol applied to a subset of the ChEMBL."""

    output_files = [f"tests/tmp/synthetic/chembl/data/prep/natref_dnp/frags_crms_nobn/10_pnp/data/chembl_{str(cid+1).zfill(3)}_pnp.csv.gz" for cid in range(3)]

    # run protocol
    fs_filter_fhits tests/tmp/synthetic/chembl/data/prep/natref_dnp/frags_crms_nobn
    command_smk = 'run_protocol synthetic -c synthetic/chembl/test_synthetic_chembl_natrefdnp_fragscrms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')

    assert all([Path(f).exists() for f in output_files])
