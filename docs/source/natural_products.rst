================
Natural Products
================

Natural products are used a reference dataset to identify PNPs in synthetic compounds.

************
The workflow
************

The natural workflow consists of 5 main steps:

    1. Preparation (Chunk, Load, Standardize, Deduplicate, Depict)
    2. Fragment Search (FS)
    3. Fragment Combination Classification (FCC)
    4. Fragment Combination Graph (FCG) 
    5. Reporting

The task tree below illustrates the workflow's logic:

.. image:: _images/workflow_natural_tasktree.svg
    :align: center
    :width: 100%

********************
Running the workflow
********************

As for fragments, two configuration files are required:

    - :download:`workflow configuration <_data/test_natural_coconut_fragscrms.json>`
    - :download:`standardization protocol <_data/test_natural_coconut_std.json>`

To run the natural workflow, run the following command:

>>> run_protocol_fc natural -c fc/03_natural/coconut/test_natural_coconut_fragscrms.json > fc/03_natural/coconut/test_natural_coconut_fragscrms.log 2>&1

***********
Folder tree
***********

The folder tree is the same as for fragments, but has an extra frags subfolder, where the results involving fragments are stored:

::

    fc
    ├── 01_fragments
    ├── 03_natural
    │   └── coconut
    │       ├── data
    │       │   ├── 00_raw
    │       │   │   └── data
    │       │   │       ├── coconut_num_mols.json
    │       │   │       └── coconut_test.sdf.gz
    │       │   └── prep
    │       │       ├── 01_chunk
    │       │       │   ├── data
    │       │       │   │   └── coconut_001.sdf.gz
    │       │       │   └── log
    │       │       │       └── coconut_chunk.log
    │       │       ├── 02_load
    │       │       │   ├── data
    │       │       │   │   └── coconut_001.csv.gz
    │       │       │   └── log
    │       │       │       └── coconut_001_load.log
    │       │       ├── 03_std
    │       │       │   ├── data
    │       │       │   │   └── coconut_001_std.csv.gz
    │       │       │   └── log
    │       │       │       ├── coconut_001_error.csv.gz
    │       │       │       ├── coconut_001_filtered.csv.gz
    │       │       │       └── coconut_001_std.log
    │       │       ├── 04_dedupl
    │       │       │   ├── coconut_ref.hdf
    │       │       │   ├── data
    │       │       │   │   └── coconut_001_dedupl.csv.gz
    │       │       │   └── log
    │       │       │       ├── coconut_001_dedupl.log
    │       │       │       ├── coconut_001_filtered.csv.gz
    │       │       │       └── coconut_001_synonyms.csv.gz
    │       │       ├── 05_depict
    │       │       │   ├── data
    │       │       │   │   └── coconut_001_depict.csv.gz
    │       │       │   └── log
    │       │       │       └── coconut_001_depict.log
    │       │       ├── frags_crms
    │       │       │   ├── 06_fs
    │       │       │   │   ├── data
    │       │       │   │   │   └── coconut_001_fs.csv.gz
    │       │       │   │   └── log
    │       │       │   │       └── coconut_001_fs.log
    │       │       │   ├── 07_fcc
    │       │       │   │   ├── data
    │       │       │   │   │   └── coconut_001_fcc.csv.gz
    │       │       │   │   └── log
    │       │       │   │       └── coconut_001_fcc.log
    │       │       │   ├── 08_fcg
    │       │       │   │   ├── data
    │       │       │   │   │   └── coconut_001_fcg.csv.gz
    │       │       │   │   └── log
    │       │       │   │       └── coconut_001_fcg.log
    │       │       │   └── report
    │       │       │       ├── data
    │       │       │       │   ├── 08_fcg
    │       │       │       │   │   ├── coconut_001_fcg_counts.csv
    │       │       │       │   │   ├── coconut_001_fcg_fcc.csv
    │       │       │       │   │   ├── coconut_001_fcg_fc.csv
    │       │       │       │   │   ├── coconut_001_fcg_fragratio.csv
    │       │       │       │   │   ├── coconut_001_fcg_nfcgpermol.csv
    │       │       │       │   │   ├── coconut_001_fcg_nhits.csv
    │       │       │       │   │   ├── coconut_001_fcg_nhits_u.csv
    │       │       │       │   │   ├── coconut_001_fcg_topfrags.csv
    │       │       │       │   │   ├── coconut_001_fcg_topfrags_u.csv
    │       │       │       │   │   └── coconut_001.log
    │       │       │       │   ├── coconut_count_mols.csv
    │       │       │       │   ├── coconut_fcg_fcc.csv
    │       │       │       │   ├── coconut_fcg_fc.csv
    │       │       │       │   ├── coconut_fcg_fragmolcov.csv
    │       │       │       │   ├── coconut_fcg_nfcgpermol.csv
    │       │       │       │   ├── coconut_fcg_nhits.csv
    │       │       │       │   ├── coconut_fcg_nhits_u.csv
    │       │       │       │   ├── coconut_fcg_top10frags.csv
    │       │       │       │   ├── coconut_fcg_top10frags_u.csv
    │       │       │       │   └── coconut_time.csv
    │       │       │       ├── log
    │       │       │       │   ├── coconut_001_count_mols.log
    │       │       │       │   ├── coconut_001_time.log
    │       │       │       │   └── report_fcg_coconut.log
    │       │       │       ├── plot
    │       │       │       │   ├── coconut_fcg_fcc.svg
    │       │       │       │   ├── coconut_fcg_fc.svg
    │       │       │       │   ├── coconut_fcg_fragmolcov.svg
    │       │       │       │   ├── coconut_fcg_fragmolcov_wo_side_chain.svg
    │       │       │       │   ├── coconut_fcg_fragmolcov._zoom.svg
    │       │       │       │   ├── coconut_fcg_nfcgpermol.svg
    │       │       │       │   ├── coconut_fcg_nhits.svg
    │       │       │       │   ├── coconut_fcg_nhits._zoom.svg
    │       │       │       │   ├── coconut_fcg_top10frags.svg
    │       │       │       │   └── coconut_fcg_top10frags_u.svg
    │       │       │       └── report_fcg_coconut.log
    │       │       └── report
    │       │           ├── data
    │       │           │   ├── coconut_prep_error.csv
    │       │           │   ├── coconut_prep_filtered.csv
    │       │           │   └── coconut_prep_overview.csv
    │       │           ├── log
    │       │           │   └── report_prep_coconut.log
    │       │           ├── plot
    │       │           │   ├── coconut_prep_error.svg
    │       │           │   ├── coconut_prep_filtered.svg
    │       │           │   └── coconut_prep_overview.svg
    │       │           └── report_prep_coconut.log
    │       ├── natural_coconut_tasktree.svg
    │       ├── test_natural_coconut_fragscrms.json
    │       ├── test_natural_coconut_fragscrms.log
    │       └── test_natural_coconut_std.json
