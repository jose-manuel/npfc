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
# fragments
root_dir = config['root_dir']
prep_subdir = config['prep_subdir']
prefix = config['prefix']
config_file = config['config_file']
commercial_ref = config['commercial_ref']
color = config['color']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from directory names for consistency
for wd in [WD,
           root_dir, prep_subdir,
           ]:
    if wd.endswith('/'):
        wd = wd[:-1]

DATA = f"{WD}/data"
REPORT = f"{WD}/report"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule END:
    input:
        dataset = DATA + "/dataset/data/dataset.csv.gz",
        fragment = DATA + "/fragment/data/fragment.csv.gz",
        molecule = DATA + "/molecule/data/molecule.csv.gz",
        molecule_dataset = DATA + "/molecule_dataset/data/molecule_dataset.csv.gz",
        molecule_molecule = DATA + "/molecule_molecule/data/molecule_molecule.csv.gz",
        report_molecule = REPORT + "/molecule/molecular_features.svg"


rule REPORT_MOL:
    input: DATA + "/molecule/data/molecule.csv.gz"
    output: REPORT + "/molecule/molecular_features.svg"
    log: REPORT + "/molecule/molecular_features.log"
    shell: "report_fct_molecule {input} {output} " + color + " " + prefix + " >{log} 2>&1"


rule FRAG:
    input: root_dir + "/data/" + prep_subdir + "/05_fcp/data/" + prefix + "_fcp.csv.gz"
    output: DATA + "/fragment/data/fragment.csv.gz"
    log: DATA + "/fragment/log/fragment.log"
    shell: "fct_fragment {input} {output} >{log} 2>&1"


rule MOL_DATASET:
    input:
        molecule = DATA + "/molecule/data/molecule.csv.gz",
        dataset = DATA + "/dataset/data/dataset.csv.gz"
    output: DATA + "/molecule_dataset/data/molecule_dataset.csv.gz"
    log: DATA + "/molecule_dataset/log/molecule_dataset.log"
    shell: "fct_molecule_dataset fragments {input.molecule} {input.dataset} {output} >{log} 2>&1"


rule DATASET:
    input: config_file
    output: DATA + "/dataset/data/dataset.csv.gz"
    log: DATA + "/dataset/log/dataset.log"
    shell: "fct_dataset fragments {input} {output} >{log} 2>&1"


rule MOL_MOL:
    input:
        synonyms = root_dir + "/data/" + prep_subdir + "/03_dedupl/log/" + prefix + "_synonyms.csv.gz",
        molecule = DATA + "/molecule/data/molecule.csv.gz"
    output: DATA + "/molecule_molecule/data/molecule_molecule.csv.gz"
    log: DATA + "/molecule_molecule/log/molecule_molecule.log"
    shell: "fct_molecule_molecule {input.synonyms} {input.molecule} {output} >{log} 2>&1"


rule MOL:
    input:
        load_step = root_dir + "/data/" + prep_subdir + "/01_load/data/" + prefix + ".csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/05_fcp/data/" + prefix + "_fcp.csv.gz",
        commercial_ref = commercial_ref
    output: DATA + "/molecule/data/molecule.csv.gz"
    log: DATA + "/molecule/log/molecule.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {input.commercial_ref} {output} >{log} 2>&1"
