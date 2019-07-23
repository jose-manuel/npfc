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

#
# def test_00_init_folder():
#     """Removes all files in tmp worklow output folder."""
#     p = Path('tests/tmp/local_workflow')
#     p.mkdir(parents=True, exist_ok=True)
#     for f in [x for x in p.iterdir() if x.is_file()]:
#         f.unlink()


# def test_01a_prep_crms_small_tasktree():
#     """Compute the task tree for the prep_crms_small pipeline."""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/prep_crms_small.smk')
#     smk_svg = "tests/tmp/prep_crms_small.svg"
#     subprocess.run(f"snakemake -s {smk_file} --dag | dot -Tsvg > {smk_svg}", shell=True, check=True)
#     assert Path(smk_svg).exists()
#
#
# def test_01b_prep_crms_small():
#     """Run the pipeline for preparing fragments (crms) on a small dataset."""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/prep_crms_small.smk')
#     output_files = ["tests/tmp/scaffolds/crms/data/8_gen2d/data/crms_2d.csv.gz"]
#     subprocess.run("snakemake -s " + smk_file, shell=True, check=True)
#     assert all([Path(f).exists() for f in output_files])
#
#
# def test_02a_fcc_dnp_small_tasktree():
#     """Compute the task tree for the fcc_dnp_small pipeline."""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_dnp_small.smk')
#     smk_svg = "tests/tmp/fcc_dnp_small.svg"
#     subprocess.run(f"snakemake -s {smk_file} --dag | dot -Tsvg > {smk_svg}", shell=True, check=True)
#     assert Path(smk_svg).exists()
#
#
# def test_02b_fcc_dnp_small():
#     """Run the pipeline for executing the FCC protocol on a small dataset from the DNP."""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_dnp_small.smk')
#     output_files = [f"tests/tmp/dnp/data/8_map/data/dnp_00{i+1}_map.csv.gz" for i in range(3)]
#     subprocess.run("snakemake -j 4 -s " + smk_file, shell=True, check=True)
#     assert all([Path(f).exists() for f in output_files])


# def test_03a_fcc_chembl_small_tasktree():
#     """Compute the task tree for the fcc_chembl_small pipeline."""
#     smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_chembl_small.smk')
#     smk_svg = "tests/tmp/fcc_chembl_small.svg"
#     subprocess.run(f"snakemake -s {smk_file} --dag | dot -Tsvg > {smk_svg}", shell=True, check=True)
#     assert Path(smk_svg).exists()


def test_03b_fcc_chembl_small():
    """Run the pipeline for executing the FCC protocol on a small dataset from the ChEMBL, using synth subsetting."""
    smk_file = pkg_resources.resource_filename('npfc', 'data/fcc_chembl_small.smk')
    output_files = [f"tests/tmp/chembl/data/11_pnp/data/chembl_00{i+1}_pnp.csv.gz" for i in range(3)]
    subprocess.run("snakemake -j 4 -s " + smk_file, shell=True, check=True)
    assert all([Path(f).exists() for f in output_files])
