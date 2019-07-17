"""
tests.test_06_local_workflow
~~~~~~~~~~
Functional test to make sure all step input/outputs are compatible
for a local run of the npfc workflow.
"""
# standard library
import subprocess
from pathlib import Path
import pkg_resources
# debug
# import logging
# logging.basicConfig(level=logging.DEBUG)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def test_init_folder():
    """Removes all files in tmp worklow output folder"""
    p = Path('tests/tmp/local_workflow')
    p.mkdir(parents=True, exist_ok=True)
    for f in [x for x in p.iterdir() if x.is_file()]:
        f.unlink()


def test_01_prep_crms_small():
    """Run the pipeline for preparing fragments (crms) on a small dataset"""
    smk_file = pkg_resources.resource_filename('npfc', 'data/prep_crms_small.smk')
    subprocess.run("snakemake -s " + smk_file, shell=True, check=True)
    # assert Path("tests/tmp/local_workflow/local_workflow.png").exists()
    # assert Path("tests/tmp/local_workflow/data/chembl_small_001_map_crms.csv.gz").exists()


def test_02_fcc_dnp_small():
    """Run the pipeline for executing the FCC protocol on a small dataset from the DNP"""
    smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_dnp_small.smk')
    subprocess.run("snakemake -j 4 -s " + smk_file, shell=True, check=True)
    # assert Path("tests/tmp/local_workflow/local_workflow.png").exists()
    # assert Path("tests/tmp/local_workflow/data/chembl_small_001_map_crms.csv.gz").exists()


def test_02_fcc_chembl_small():
    """Run the pipeline for executing the FCC protocol on a small dataset from the ChEMBL, using synth subsetting"""
    smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_chembl_small.smk')
    subprocess.run("snakemake -j 4 -s " + smk_file, shell=True, check=True)
    # assert Path("tests/tmp/local_workflow/local_workflow.png").exists()
    # assert Path("tests/tmp/local_workflow/data/chembl_small_001_map_crms.csv.gz").exists()


# def test_03_fcc_zinc_small():
#     """Run the pipeline for executing the FCC protocol on a small dataset from the ChEMBL, using synth subsetting"""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_zinc_small.smk')
#     subprocess.run("snakemake -j 4 -s " + smk_file, shell=True, check=True)
#     # assert Path("tests/tmp/local_workflow/local_workflow.png").exists()
#     # assert Path("tests/tmp/local_workflow/data/chembl_small_001_map_crms.csv.gz").exists()
