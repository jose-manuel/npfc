"""
tests.test_11d_fct_smk_synthetic
==================================
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.

This requires that the test_10d_fc_smk_synthetic completed successfully so its outputs
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

    data_ori = "tests/data/fct/crms_dnp_chembl_zinc/data/04_synthetic"
    data_tgt = "tests/tmp/fct/crms_dnp_chembl_zinc/data/04_synthetic"

    # delete any previous jobs from workflows
    if Path(data_tgt).is_dir():
        shutil.rmtree(data_tgt)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    if Path(data_ori).is_dir():
        shutil.copytree(str(data_ori), str(data_tgt))


def test_run():
    """Run the 'synthetic' protocol applied to a subset of the ChEMBL."""

    WD = 'tests/tmp/fct/crms_dnp_chembl_zinc/data/04_synthetic'
    cids = [f"{str(cid+1).zfill(3)}" for cid in range(2)]
    output_files = [f"{WD}/data/assay/data/assay_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/assay_document/data/assay_document_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/assay_species/data/assay_species_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/dataset/data/dataset.csv.gz"] + \
                   [f"{WD}/data/document/data/document_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/fragment_fragment/data/fragment_fragment_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/molecule/data/molecule_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/molecule_assay/data/molecule_assay_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/molecule_dataset/data/molecule_dataset_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/molecule_document/data/molecule_document_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/molecule_molecule/data/molecule_molecule_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/species/data/species_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/target/data/target_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/target_assay/data/target_assay_{cid}.csv.gz" for cid in cids] + \
                   [f"{WD}/data/target_species/data/target_species_{cid}.csv.gz" for cid in cids]

    output_svg = f"{WD}/test_synthetic_chembl_natrefdnp_fragscrms_tasktree.svg"
    # run protocol
    command_smk = 'run_protocol_fct synthetic -c fct/crms_dnp_chembl_zinc/data/04_synthetic/test_synthetic_chembl_natrefdnp_fragscrms.json'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
