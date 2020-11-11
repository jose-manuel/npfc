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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from directory names for consistency
for wd in [WD,
           root_dir, prep_subdir,
           ]:
    if wd.endswith('/'):
        wd = wd[:-1]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule END:
    input:
        dataset = WD + "/dataset/data/dataset.csv.gz",
        fragment = WD + "/fragment/data/fragment.csv.gz",
        molecule = WD + "/molecule/data/molecule.csv.gz",
        molecule_dataset = WD + "/molecule_dataset/data/molecule_dataset.csv.gz",
        molecule_molecule = WD + "/molecule_molecule/data/molecule_molecule.csv.gz"

rule FRAG:
    input: root_dir + "/data/" + prep_subdir + "/05_depict/data/" + prefix + "_depict.csv.gz"
    output: WD + "/fragment/data/fragment.csv.gz"
    log: WD + "/fragment/log/fragment.log"
    shell: "fct_fragment {input} {output}  >{log} 2>&1"

rule MOL_DATASET:
    input:
        molecule = WD + "/molecule/data/molecule.csv.gz",
        dataset = WD + "/dataset/data/dataset.csv.gz"
    output: WD + "/molecule_dataset/data/molecule_dataset.csv.gz"
    log: WD + "/molecule_dataset/log/molecule_dataset.log"
    shell: "fct_molecule_dataset fragments {input.molecule} {input.dataset} {output}  >{log} 2>&1"

rule DATASET:
    input: config_file
    output: WD + "/dataset/data/dataset.csv.gz"
    log: WD + "/dataset/log/dataset.log"
    shell: "fct_dataset fragments {input} {output}  >{log} 2>&1"

rule MOL_MOL:
    input:
        synonyms = root_dir + "/data/" + prep_subdir + "/03_dedupl/log/" + prefix + "_synonyms.csv.gz",
        molecule = WD + "/molecule/data/molecule.csv.gz"
    output: WD + "/molecule_molecule/data/molecule_molecule.csv.gz"
    log: WD + "/molecule_molecule/log/molecule_molecule.log"
    shell: "fct_molecule_molecule {input.synonyms} {input.molecule} {output}  >{log} 2>&1"

rule MOL:
    input:
        load_step = root_dir + "/data/" + prep_subdir + "/01_load/data/" + prefix + ".csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/05_depict/data/" + prefix + "_depict.csv.gz"
    output: WD + "/molecule/data/molecule.csv.gz"
    log: WD + "/molecule/log/molecule.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {output}  >{log} 2>&1"
