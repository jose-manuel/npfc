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

# # natural
# natural_root = config['root_dir']
# natural_prep_subdir = config['natural_prep_subdir']
# natural_frags_subdir = config['natural_frags_subdir']
# natural_prefix = config['natural_prefix']
# # synthetic
# synthetic_root = config['synthetic_root']
# synthetic_prep_subdir = config['synthetic_prep_subdir']
# synthetic_frags_subdir = config['synthetic_frags_subdir']
# synthetic_frags_subdir = config['synthetic_natref_subdir']
# synthetic_prefix = config['synthetic_prefix']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from rowids for consistency
for wd in [WD,
           root_dir, prep_subdir,
           # natural_root, natural_prep_subdir, natural_frags_subdir,
           # synthetic_root, synthetic_prep_subdir, synthetic_frags_subdir, synthetic_frags_subdir,
           ]:
    if wd.endswith('/'):
        wd = wd[:-1]

# # count the number of chunks for naturral and synthetic runs
# natural_num_chunks = len(list(Path(natural_prep_subdir).glob('01_chunk/data/*')))
# synthetic_num_chunks = len(list(Path(synthetic_prep_subdir).glob('01_chunk/data/*')))
#
# # define chunk_ids for wildcard expansion
# natural_chunk_ids = [str(i+1).zfill(3) for i in range(natural_num_chunks)]
# synthetic_chunk_ids = [str(i+1).zfill(3) for i in range(synthetic_num_chunks)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule END:
    input:
        dataset = WD + "/data/dataset.csv.gz",
        molecule = WD + "/data/molecule.csv.gz",
        molecule_dataset = WD + "/data/molecule_dataset.csv.gz",
        molecule_molecule = WD + "/data/molecule_molecule.csv.gz"

rule MOLECULE_DATASET:
    input:
        molecule = WD + "/data/molecule.csv.gz",
        dataset = WD + "/data/dataset.csv.gz"
    output: WD + "/data/molecule_dataset.csv.gz"
    log: WD + "/log/molecule_dataset.log"
    shell: "fct_molecule_dataset fragments {input.molecule} {input.dataset} {output}  >{log} 2>&1"

rule DATASET:
    input: config_file
    output: WD + "/data/dataset.csv.gz"
    log: WD + "/log/dataset.log"
    shell: "fct_dataset fragments {input} {output}  >{log} 2>&1"

rule MOL_MOL:
    input:
        synonyms = root_dir + "/data/" + prep_subdir + "/03_dedupl/log/" + prefix + "_synonyms.csv.gz",
        molecule = WD + "/data/molecule.csv.gz"
    output: WD + "/data/molecule_molecule.csv.gz"
    log: WD + "/log/molecule_molecule.log"
    shell: "fct_molecule_molecule {input.synonyms} {input.molecule} {output}  >{log} 2>&1"

rule MOL:
    input:
        # load_step = expand("{a}/data/{b}/01_load/data/{c}.csv.gz", a=[root_dir], b=[prep_subdir], c=[prefix]),
        # latest_step = expand("{a}/data/{b}/05_depict/data/{c}_depict.csv.gz", a=[root_dir], b=[prep_subdir], c=[prefix])
        load_step = root_dir + "/data/" + prep_subdir + "/01_load/data/" + prefix + ".csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/05_depict/data/" + prefix + "_depict.csv.gz"
    output: WD + "/data/molecule.csv.gz"
    log: WD + "/log/molecule.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {output}  >{log} 2>&1"

rule START:
    # This step has for unique goal to define the outputs that will be used later on.
    # But this is not without risk... Never run this pipeline while forcing it from the beginning
    # as it would certainly erase the outputs from the other pipelines (load and latest step)!
    output:
        load_step = root_dir + "/data/" + prep_subdir + "/01_load/data/" + prefix + ".csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/05_depict/data/" + prefix + "_depict.csv.gz",
        config_file = config_file
    shell:
        """
        mv {output.load_step} {output.load_step}.tmp; mv {output.load_step}.tmp {output.load_step}
        mv {output.latest_step} {output.latest_step}.tmp; mv {output.latest_step}.tmp {output.latest_step}
        mv {output.config_file} {output.config_file}.tmp; mv {output.config_file}.tmp {output.config_file}
        """
