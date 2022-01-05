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
# data
import pandas as pd
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
    report_count = 'tests/tmp/fc/01_fragments/crms/data/prep/report/data/crms_count_mols.csv'
    # run protocol
    command_smk = 'run_protocol_fc fragments -c fc/01_fragments/crms/test_fragments_crms.json > fc/01_fragments/crms/test_fragments_crms.log 2>&1'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    # check output files
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])
    # check number of molecules (no diff with entries for fragments)
    df = pd.read_csv(report_count, sep='|')
    assert df.iloc[0]['01_load_num_mols'] == 117 # first step
    assert df.iloc[0]['05_fcp_num_mols'] == 109 # last step
    # display counts on the terminal
    df = df.T.reset_index()
    df['step'] = df['index'].map(lambda x: str(x.split('_num_')[0]))
    df_mols = df[df['index'].str.contains('_num_mols')].rename({0: 'count_mols'}, axis=1).drop('index', axis=1)
    df_entries = df[df['index'].str.contains('_num_entries')].rename({0: 'count_entries'}, axis=1).drop('index', axis=1)
    df_mols['count_mols'] = df_mols['count_mols'].astype(int)
    df_entries['count_entries'] = df_entries['count_entries'].astype(int)
    df = df_mols.merge(df_entries, on='step')[['step', 'count_mols', 'count_entries']]
    print(f"\nResults:\n{df}\n")
