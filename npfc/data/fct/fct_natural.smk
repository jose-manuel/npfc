#!/usr/bin/env python

"""
Script fragment_combination_tree.smk
=====================================
The pipeline to apply for generating the fragment combination tree (FCT).
It uses 4 types of inputs:

    - pre-defined raw tables extracted from the postgres ChEMBL database
    - results from the fragments pipeline
    - results from the natural pipeline
    - results from the synthetic pipeline

At this point in the development, there is no intent to merge runs from different conditions,
but the graph database structure is likeley flexible enough to do so.
"""

from pathlib import Path
import shutil
import pkg_resources
from npfc import load
from npfc import utils


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# common
WD = config['WD']
# natural
root_dir = config['root_dir']
prep_subdir = config['prep_subdir']
prefix = config['prefix']
config_file = config['config_file']
# fragments
frags_subdir = config['frags_subdir']
commercial_ref = config['commercial_ref']
color = config['color']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from rowids for consistency
for wd in [WD,
           root_dir, prep_subdir,
           # natural_root, natural_prep_subdir, natural_frags_subdir,
           # synthetic_root, synthetic_prep_subdir, synthetic_frags_subdir, synthetic_frags_subdir,
           ]:
    if wd.endswith('/'):
        wd = wd[:-1]

chunk_ids = [str(i+1).zfill(3) for i in range(config['num_chunks'])]

# # count the number of chunks for naturral and synthetic runs
# natural_num_chunks = len(list(Path(natural_prep_subdir).glob('01_chunk/data/*')))
# synthetic_num_chunks = len(list(Path(synthetic_prep_subdir).glob('01_chunk/data/*')))
#
# # define chunk_ids for wildcard expansion
# natural_chunk_ids = [str(i+1).zfill(3) for i in range(natural_num_chunks)]
# synthetic_chunk_ids = [str(i+1).zfill(3) for i in range(synthetic_num_chunks)]

DATA = f"{WD}/data"
REPORT = f"{WD}/report"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule END:
    input:
        dataset = DATA + "/dataset/data/dataset.csv.gz",
        fragment_fragment = expand("{DATA}/fragment_fragment/data/fragment_fragment_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule = expand("{DATA}/molecule/data/molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_dataset = expand("{DATA}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_molecule = expand("{DATA}/molecule_molecule/data/molecule_molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        report_molecule = REPORT + "/molecule/molecular_features.svg"


rule REPORT_MOL:
    input: expand("{DATA}/molecule/data/molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids)
    output: REPORT + "/molecule/molecular_features.svg"
    log: REPORT + "/molecule/molecular_features.log"
    shell: "report_fct_molecule " + DATA + "/molecule/data" + " {output} " + color + " " + prefix + " >{log} 2>&1"


rule FRAG_FRAG:
    input: root_dir + "/data/" + prep_subdir + "/" + frags_subdir + "/08_fcg/data/" + prefix + "_{cid}_fcg.csv.gz"
    output: DATA + "/fragment_fragment/data/fragment_fragment_{cid}.csv.gz"
    log: DATA + "/fragment_fragment/log/fragment_fragment_{cid}.log"
    shell: "fct_fragment_fragment {input} {output}  >{log} 2>&1"


rule MOL_DATASET:
    input:
        molecule = "{DATA}/molecule/data/molecule_{cid}.csv.gz",
        dataset = "{DATA}/dataset/data/dataset.csv.gz"
    output: "{DATA}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz"
    log: "{DATA}/molecule_dataset/log/molecule_dataset_{cid}.log"
    shell: "fct_molecule_dataset natural {input.molecule} {input.dataset} {output}  >{log} 2>&1"


rule DATASET:
    input: config_file
    output: DATA + "/dataset/data/dataset.csv.gz"
    log: DATA + "/dataset/log/dataset.log"
    shell: "fct_dataset natural {input} {output}  >{log} 2>&1"


rule MOL_MOL:
    input:
        synonyms = root_dir + "/data/" + prep_subdir + "/04_dedupl/log/" + prefix + "_{cid}_synonyms.csv.gz",
        molecule = DATA + "/molecule/data/molecule_{cid}.csv.gz"
    output: DATA + "/molecule_molecule/data/molecule_molecule_{cid}.csv.gz"
    log: DATA + "/molecule_molecule/log/molecule_molecule_{cid}.csv.gz"
    shell: "fct_molecule_molecule {input.synonyms} {input.molecule} {output}  >{log} 2>&1"


rule MOL:
    input:
        load_step = root_dir + "/data/" + prep_subdir + "/02_load/data/" + prefix + "_{cid}.csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/" + frags_subdir + "/08_fcg/data/" + prefix + "_{cid}_fcg.csv.gz",
        commercial_ref = commercial_ref
    output: DATA + "/molecule/data/molecule_{cid}.csv.gz"
    log: DATA + "/molecule/log/molecule_{cid}.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {input.commercial_ref} {output}  >{log} 2>&1"
