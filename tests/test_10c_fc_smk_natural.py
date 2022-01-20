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
# data
import pandas as pd
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
    report_count = 'tests/tmp/fc/03_natural/coconut/data/prep/frags_crms/report/data/coconut_count_mols.csv'

    # run protocol
    command_smk = 'run_protocol_fc natural -c fc/03_natural/coconut/test_natural_coconut_fragscrms.json > fc/03_natural/coconut/test_natural_coconut_fragscrms.log 2>&1'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])

    # check number of entries/molecules
    df = pd.read_csv(report_count, sep='|')
    df = df[df['subset'] == 'total']
    assert df.iloc[0]['01_chunk_num_mols'] == 20 # init
    assert df.iloc[0]['05_depict_num_mols'] == 17 # prep
    assert df.iloc[0]['06_fs_num_mols'] == 17 and df.iloc[0]['06_fs_num_entries'] == 55 # fragment search
    assert df.iloc[0]['07_fcc_num_mols'] == 14 and df.iloc[0]['07_fcc_num_entries'] == 70 # fragment combination
    assert df.iloc[0]['08_fcg_num_mols'] == 12 and df.iloc[0]['08_fcg_num_entries'] == 12 # final

    # display counts on the terminal
    df = df.T.reset_index()
    name_last_col = [c for c in df.columns][-1]
    df['step'] = df['index'].map(lambda x: str(x.split('_num_')[0]))
    df_mols = df[df['index'].str.contains('_num_mols')].rename({name_last_col: 'count_mols'}, axis=1).drop('index', axis=1)
    df_entries = df[df['index'].str.contains('_num_entries')].rename({name_last_col: 'count_entries'}, axis=1).drop('index', axis=1)
    df_mols['count_mols'] = df_mols['count_mols'].astype(int)
    df_entries['count_entries'] = df_entries['count_entries'].astype(int)
    df = df_mols.merge(df_entries, on='step')[['step', 'count_mols', 'count_entries']]
    print(f"\nResults:\n{df}\n")
