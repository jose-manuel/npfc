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
# data
import pandas as pd
# debug
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_init():
    """Reset the folder where outputs are computed to its initial state."""

    data_ori = "tests/data/fc/04_synthetic"
    data_tgt = "tests/tmp/fc/04_synthetic"
    # delete any previous jobs from workflows
    for subdir in Path(data_tgt).glob("*"):
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # copy reference subdirectories for reseting the workflows folder to their initial state
    for subdir in Path(data_ori).glob("*"):
        if subdir.is_dir():
            shutil.copytree(str(subdir), f"{data_tgt}/{subdir.stem}")


def test_run():
    """Run the 'synthetic' protocol applied to a subset of the ChEMBL."""

    output_files = [f"tests/tmp/fc/04_synthetic/chembl/data/prep/natref_coconut/frags_crms/10_pnp/data/chembl_{str(cid+1).zfill(3)}_pnp.csv.gz" for cid in range(2)]
    output_svg = 'tests/tmp/fc/04_synthetic/chembl/synthetic_chembl_tasktree.svg'
    report_count = 'tests/tmp/fc/04_synthetic/chembl/data/prep/natref_coconut/frags_crms/report/data/chembl_count_mols.csv'

    # run protocol
    command_smk = 'run_protocol_fc synthetic -c fc/04_synthetic/chembl/test_synthetic_chembl_natrefdnp_fragscrms.json > fc/04_synthetic/chembl/test_synthetic_chembl_natrefdnp_fragscrms.log 2>&1'
    subprocess.run(command_smk, shell=True, check=True, cwd='tests/tmp')
    assert Path(output_svg).exists()
    assert all([Path(f).exists() for f in output_files])

    # check number of entries/molecules
    df = pd.read_csv(report_count, sep='|')
    df = df[df['subset'] == 'total']
    assert df.iloc[0]['01_chunk_num_mols'] == 200 # init
    assert df.iloc[0]['05_depict_num_mols'] == 139 # prep
    assert df.iloc[0]['06_subset_num_mols'] == 139 # subset
    assert df.iloc[0]['07_fs_num_mols'] == 126 and df.iloc[0]['07_fs_num_entries'] == 381 # fragment search
    assert df.iloc[0]['08_fcc_num_mols'] == 77 and df.iloc[0]['08_fcc_num_entries'] == 384 # fragment combination
    assert df.iloc[0]['09_fcg_num_mols'] == 72 and df.iloc[0]['09_fcg_num_entries'] == 84 # fragment combination graphs
    assert df.iloc[0]['10_pnp_num_mols'] == 54 and df.iloc[0]['10_pnp_num_entries'] == 64 # final

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
