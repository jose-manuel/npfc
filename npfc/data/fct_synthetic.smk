#!/usr/bin/env python

"""
Script fct_synthetic.smk
=====================================
The pipeline to apply for generating the fragment combination tree (FCT).
It uses 2 types of inputs:

    - pre-defined raw tables extracted from the postgres ChEMBL database
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
# natref
natref_subdir = config['natref_subdir']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from rowids for consistency
for wd in [WD,
           root_dir, prep_subdir,
           natref_subdir,
           frags_subdir,
           ]:
    if wd.endswith('/'):
        wd = wd[:-1]

chunk_ids = [str(i+1).zfill(3) for i in range(config['num_chunks'])]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule END:
    input:
        # nodes
        assay = expand("{WD}/assay/data/assay_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        dataset = WD + "/dataset/data/dataset.csv.gz",
        document = expand("{WD}/document/data/document_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        molecule = expand("{WD}/molecule/data/molecule_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        species = expand("{WD}/species/data/species_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        target = expand("{WD}/target/data/target_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        # relationships
        assay_document = expand("{WD}/assay_document/data/assay_document_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        assay_species = expand("{WD}/assay_species/data/assay_species_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        fragment_fragment = expand("{WD}/fragment_fragment/data/fragment_fragment_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        molecule_assay = expand("{WD}/molecule_assay/data/molecule_assay_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        molecule_dataset = expand("{WD}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        molecule_document = expand("{WD}/molecule_document/data/molecule_document_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        molecule_molecule = expand("{WD}/molecule_molecule/data/molecule_molecule_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        target_assay = expand("{WD}/target_assay/data/target_assay_{cid}.csv.gz", WD=WD, cid=chunk_ids),
        target_species = expand("{WD}/target_species/data/target_species_{cid}.csv.gz", WD=WD, cid=chunk_ids)

rule SPECIES:
    input:
        species_raw = "{WD}/raw/data/species.csv.gz",
        target_species = "{WD}/target_species/data/target_species_{cid}.csv.gz",
        assay_species = "{WD}/assay_species/data/assay_species_{cid}.csv.gz",
    output: "{WD}/species/data/species_{cid}.csv.gz"
    log: "{WD}/species/log/species_{cid}.log"
    shell: "fct_species {input.species_raw} {input.target_species} {input.assay_species} {output} >{log} 2>&1"


rule TARGET_SPECIES:
    input:
        target_species_raw = "{WD}/raw/data/target_species.csv.gz",
        target = "{WD}/target/data/target_{cid}.csv.gz"
    output: "{WD}/target_species/data/target_species_{cid}.csv.gz"
    log: "{WD}/target_species/log/target_species_{cid}.log"
    shell: "fct_target_species {input.target_species_raw} {input.target} {output} >{log} 2>&1"


rule TARGET:
    input:
        target_raw = "{WD}/raw/data/target.csv.gz",
        target_assay = "{WD}/target_assay/data/target_assay_{cid}.csv.gz"
    output: "{WD}/target/data/target_{cid}.csv.gz"
    log: "{WD}/target/log/target_{cid}.log"
    shell: "fct_target {input.target_raw} {input.target_assay} {output} >{log} 2>&1"

rule TARGET_ASSAY:
    input:
        target_assay_raw = "{WD}/raw/data/target_assay.csv.gz",
        assay = "{WD}/assay/data/assay_{cid}.csv.gz"
    output: "{WD}/target_assay/data/target_assay_{cid}.csv.gz"
    log: "{WD}/target_assay/log/target_assay_{cid}.log"
    shell: "fct_target_assay {input.target_assay_raw} {input.assay} {output} >{log} 2>&1"

rule ASSAY_SPECIES:
    input:
        assay_species_raw = "{WD}/raw/data/assay_species.csv.gz",
        assay = "{WD}/assay/data/assay_{cid}.csv.gz"
    output: "{WD}/assay_species/data/assay_species_{cid}.csv.gz"
    log: "{WD}/assay_species/log/assay_species_{cid}.log"
    shell: "fct_assay_species {input.assay_species_raw} {input.assay} {output} >{log} 2>&1"

rule DOC:
    input:
        document_raw = "{WD}/raw/data/document.csv.gz",
        molecule_document = "{WD}/molecule_document/data/molecule_document_{cid}.csv.gz",
        assay_document = "{WD}/assay_document/data/assay_document_{cid}.csv.gz"
    output: "{WD}/document/data/document_{cid}.csv.gz"
    log: "{WD}/document/log/document_{cid}.log"
    shell: "fct_document {input.document_raw} {input.molecule_document} {input.assay_document} {output} >{log} 2>&1"

rule ASSAY_DOC:
    input:
        assay_document_raw = "{WD}/raw/data/assay_document.csv.gz",
        assay = "{WD}/assay/data/assay_{cid}.csv.gz"
    output: "{WD}/assay_document/data/assay_document_{cid}.csv.gz"
    log: "{WD}/assay_document/log/assay_document_{cid}.log"
    shell: "fct_assay_document {input.assay_document_raw} {input.assay} {output} >{log} 2>&1"

rule MOL_DOC:
    input:
        molecule_document_raw = "{WD}/raw/data/molecule_document.csv.gz",
        molecule = "{WD}/molecule/data/molecule_{cid}.csv.gz"
    output: "{WD}/molecule_document/data/molecule_document_{cid}.csv.gz"
    log: "{WD}/molecule_document/log/molecule_document_{cid}.log"
    shell: "fct_molecule_document {input.molecule_document_raw} {input.molecule} {output} >{log} 2>&1"

rule ASSAY:
    input:
        assay_raw = "{WD}/raw/data/assay.csv.gz",
        molecule_assay = "{WD}/molecule_assay/data/molecule_assay_{cid}.csv.gz"
    output: "{WD}/assay/data/assay_{cid}.csv.gz"
    log: "{WD}/assay/log/assay_{cid}.log"
    shell: "fct_assay {input.assay_raw} {input.molecule_assay} {output} >{log} 2>&1"

rule MOL_ASSAY:
    input:
        mol_assay_raw = "{WD}/raw/data/molecule_assay.csv.gz",
        molecule = "{WD}/molecule/data/molecule_{cid}.csv.gz"
    output: "{WD}/molecule_assay/data/molecule_assay_{cid}.csv.gz"
    log: "{WD}/molecule_assay/log/molecule_assay_{cid}.log"
    shell: "fct_molecule_assay {input.mol_assay_raw} {input.molecule} {output} >{log} 2>&1"

rule FRAG_FRAG:
    input: root_dir + "/data/" + prep_subdir + "/" + natref_subdir + "/" + frags_subdir + "/10_pnp/data/" + prefix + "_{cid}_pnp.csv.gz"
    output: WD + "/fragment_fragment/data/fragment_fragment_{cid}.csv.gz"
    log: WD + "/fragment_fragment/log/fragment_fragment_{cid}.log"
    shell: "fct_fragment_fragment {input} {output}  >{log} 2>&1"

rule MOL_DATASET:
    input:
        molecule = "{WD}/molecule/data/molecule_{cid}.csv.gz",
        dataset = "{WD}/dataset/data/dataset.csv.gz"
    output: "{WD}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz"
    log: "{WD}/molecule_dataset/log/molecule_dataset_{cid}.log"
    shell: "fct_molecule_dataset synthetic {input.molecule} {input.dataset} {output}  >{log} 2>&1"

rule DATASET:
    input: config_file
    output: WD + "/dataset/data/dataset.csv.gz"
    log: WD + "/dataset/log/dataset.log"
    shell: "fct_dataset synthetic {input} {output}  >{log} 2>&1"

rule MOL_MOL:
    input:
        synonyms = root_dir + "/data/" + prep_subdir + "/04_dedupl/log/" + prefix + "_{cid}_synonyms.csv.gz",
        molecule = WD + "/molecule/data/molecule_{cid}.csv.gz"
    output: WD + "/molecule_molecule/data/molecule_molecule_{cid}.csv.gz"
    log: WD + "/molecule_molecule/log/molecule_molecule_{cid}.csv.gz"
    shell: "fct_molecule_molecule {input.synonyms} {input.molecule} {output}  >{log} 2>&1"

rule MOL:
    input:
        load_step = root_dir + "/data/" + prep_subdir + "/02_load/data/" + prefix + "_{cid}.csv.gz",
        latest_step = root_dir + "/data/" + prep_subdir + "/" + natref_subdir + "/" + frags_subdir + "/10_pnp/data/" + prefix + "_{cid}_pnp.csv.gz"
    output: WD + "/molecule/data/molecule_{cid}.csv.gz"
    log: WD + "/molecule/log/molecule_{cid}.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {output}  >{log} 2>&1"
