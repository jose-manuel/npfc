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
import shutil
# debug
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_00_init_folder():
    """Reset the folder where outputs are computed to its initial state."""

    data_ori = "tests/data/"
    data_tgt = "tests/tmp/"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_01_fragments():
    """Create the tasktree and run the pipeline for preparing fragments on a small dataset."""
    smk_file = pkg_resources.resource_filename('npfc', 'data/fragments.smk')
    smk_svg = "tests/tmp/fragments/fragments.svg"
    output_files = ["tests/tmp/fragments/crms/data/08_gen2D/data/crms_gen2D.csv.gz"]
    prefix = 'crms'
    molid = 'Cluster'
    WD = 'tests/tmp/fragments/crms/data/'
    input_file = 'tests/tmp/fragments/crms/data/00_raw/data/cr.sdf.gz'

    command_base = f"snakemake -j 4 -s {smk_file}  \
                   --config  prefix='{prefix}' molid='{molid}' WD='{WD}' input_file='{input_file}'"
    # compute task tree svg
    command_svg = command_base + f" --dag | dot -Tsvg > {smk_svg}"
    subprocess.run(command_svg, shell=True, check=True)
    assert Path(smk_svg).exists()
    # run protocol
    command_run = command_base  # no need for filtering stdout with fragments
    subprocess.run(command_run, shell=True, check=True)
    assert all([Path(f).exists() for f in output_files])


def test_02_natural():
    """Create the tasktree and run the pipeline for executing the FCC protocol on a small dataset from the DNP."""
    smk_file = pkg_resources.resource_filename('npfc', 'data/natural.smk')
    smk_svg = "tests/tmp/natural/natural.svg"
    output_files = [f"tests/tmp/natural/dnp/data/09_fmap/data/dnp_{str(cid+1).zfill(3)}_fmap.csv.gz" for cid in range(3)]
    prefix = 'dnp'
    molid = 'UKEY'
    WD = 'tests/tmp/natural/dnp/data/'
    input_file = 'tests/tmp/natural/dnp/data/00_raw/data/dnp.sdf.gz'
    frags_file = 'tests/tmp/fragments/crms/data/08_gen2D/data/crms_gen2D.csv.gz'
    chunksize = 100
    command_base = f"snakemake -j 4 -s {smk_file}  \
                   --config  prefix='{prefix}' molid='{molid}' chunksize={chunksize} \
                   WD='{WD}' input_file='{input_file}' frags_file={frags_file}"
    # compute task tree svg
    command_svg = command_base + f" --dag | dot -Tsvg > {smk_svg}"
    subprocess.run(command_svg, shell=True, check=True)
    assert Path(smk_svg).exists()
    # run protocol
    command_run = command_base + f" 2>&1 | grep -v INFO:"
    # I have a bug with snakemake when one rule runs a subprocess.
    # All log outputs are doubled, one in the usual green and one in white that starts with 'INFO:'.
    # The best comporomise I found is to filter out these extra INFO lines from the stderr stream,
    # but it cancels the snakemake color scheme.
    subprocess.run(command_run, shell=True, check=True)
    assert all([Path(f).exists() for f in output_files])


def test_03_synthetic():
    """Create the tasktree and run the pipeline for executing the FCC protocol on a small dataset from the DNP."""
    # init
    smk_file = pkg_resources.resource_filename('npfc', 'data/synthetic.smk')
    smk_svg = "tests/tmp/synthetic/synthetic.svg"
    output_files = [f"tests/tmp/synthetic/chembl/data/12_pnp/data/chembl_{str(cid+1).zfill(3)}_pnp.csv.gz" for cid in range(3)]
    prefix = 'chembl'
    molid = 'chembl_id'
    chunksize = 100
    WD = 'tests/tmp/synthetic/chembl/data/'
    input_file = 'tests/tmp/synthetic/chembl/data/00_raw/data/chembl.sdf.gz'
    frags_file = 'tests/tmp/fragments/crms/data/08_gen2D/data/crms_gen2D.csv.gz'
    natref_uni_reffile = 'tests/tmp/natural/dnp/data/05_uni/dnp_ref.hdf'
    natref_fmap_dir = 'tests/tmp/natural/dnp/data/09_fmap/data/'
    act_file = 'tests/data/synthetic/chembl/data/00_raw/data/chembl_act_raw.csv.gz'
    command_base = f"snakemake -k -j 4 -s {smk_file}  \
                   --config  prefix='{prefix}' molid='{molid}' chunksize={chunksize} \
                   WD='{WD}' input_file='{input_file}' frags_file={frags_file} \
                   natref_uni_reffile='{natref_uni_reffile}' natref_fmap_dir='{natref_fmap_dir}' act_file='{act_file}'"
    # compute task tree svg
    command_svg = command_base + f" --dag | dot -Tsvg > {smk_svg}"
    subprocess.run(command_svg, shell=True, check=True)
    assert Path(smk_svg).exists()
    # run protocol
    command_run = command_base + f" 2>&1 | grep -v INFO:"
    subprocess.run(command_run, shell=True, check=True)
    assert all([Path(f).exists() for f in output_files])
