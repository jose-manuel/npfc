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
commercial_ref = config['commercial_ref']
color = config['color']


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
DATA = f"{WD}/data"
REPORT = f"{WD}/report"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rule END:
    input:
        # nodes
        assay = expand("{DATA}/assay/data/assay_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        dataset = DATA + "/dataset/data/dataset.csv.gz",
        document = expand("{DATA}/document/data/document_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule = expand("{DATA}/molecule/data/molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        species = expand("{DATA}/species/data/species_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        target = expand("{DATA}/target/data/target_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        # relationships
        assay_document = expand("{DATA}/assay_document/data/assay_document_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        assay_species = expand("{DATA}/assay_species/data/assay_species_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        fragment_fragment = expand("{DATA}/fragment_fragment/data/fragment_fragment_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_assay = expand("{DATA}/molecule_assay/data/molecule_assay_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_dataset = expand("{DATA}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_document = expand("{DATA}/molecule_document/data/molecule_document_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        molecule_molecule = expand("{DATA}/molecule_molecule/data/molecule_molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        target_assay = expand("{DATA}/target_assay/data/target_assay_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        target_species = expand("{DATA}/target_species/data/target_species_{cid}.csv.gz", DATA=DATA, cid=chunk_ids),
        report_molecule = REPORT + "/molecule/molecular_features.svg"


rule REPORT_MOL:
    input: expand("{DATA}/molecule/data/molecule_{cid}.csv.gz", DATA=DATA, cid=chunk_ids)
    output: REPORT + "/molecule/molecular_features.svg"
    log: REPORT + "/molecule/molecular_features.log"
    shell: "report_fct_molecule " + DATA + "/molecule/data" + " {output} " + color + " " + prefix + " >{log} 2>&1"


rule SPECIES:
    input:
        species_raw = WD + "/raw/data/species.csv.gz",
        target_species = "{DATA}/target_species/data/target_species_{cid}.csv.gz",
        assay_species = "{DATA}/assay_species/data/assay_species_{cid}.csv.gz",
    output: "{DATA}/species/data/species_{cid}.csv.gz"
    log: "{DATA}/species/log/species_{cid}.log"
    shell: "fct_species {input.species_raw} {input.target_species} {input.assay_species} {output} >{log} 2>&1"


rule TARGET_SPECIES:
    input:
        target_species_raw = WD + "/raw/data/target_species.csv.gz",
        target = "{DATA}/target/data/target_{cid}.csv.gz"
    output: "{DATA}/target_species/data/target_species_{cid}.csv.gz"
    log: "{DATA}/target_species/log/target_species_{cid}.log"
    shell: "fct_target_species {input.target_species_raw} {input.target} {output} >{log} 2>&1"


rule TARGET:
    input:
        target_raw = WD + "/raw/data/target.csv.gz",
        target_assay = "{DATA}/target_assay/data/target_assay_{cid}.csv.gz"
    output: "{DATA}/target/data/target_{cid}.csv.gz"
    log: "{DATA}/target/log/target_{cid}.log"
    shell: "fct_target {input.target_raw} {input.target_assay} {output} >{log} 2>&1"


rule TARGET_ASSAY:
    input:
        target_assay_raw = WD + "/raw/data/target_assay.csv.gz",
        assay = "{DATA}/assay/data/assay_{cid}.csv.gz"
    output: "{DATA}/target_assay/data/target_assay_{cid}.csv.gz"
    log: "{DATA}/target_assay/log/target_assay_{cid}.log"
    shell: "fct_target_assay {input.target_assay_raw} {input.assay} {output} >{log} 2>&1"


rule ASSAY_SPECIES:
    input:
        assay_species_raw = WD + "/raw/data/assay_species.csv.gz",
        assay = "{DATA}/assay/data/assay_{cid}.csv.gz"
    output: "{DATA}/assay_species/data/assay_species_{cid}.csv.gz"
    log: "{DATA}/assay_species/log/assay_species_{cid}.log"
    shell: "fct_assay_species {input.assay_species_raw} {input.assay} {output} >{log} 2>&1"


rule DOC:
    input:
        document_raw = WD + "/raw/data/document.csv.gz",
        molecule_document = "{DATA}/molecule_document/data/molecule_document_{cid}.csv.gz",
        assay_document = "{DATA}/assay_document/data/assay_document_{cid}.csv.gz"
    output: "{DATA}/document/data/document_{cid}.csv.gz"
    log: "{DATA}/document/log/document_{cid}.log"
    shell: "fct_document {input.document_raw} {input.molecule_document} {input.assay_document} {output} >{log} 2>&1"


rule ASSAY_DOC:
    input:
        assay_document_raw = WD + "/raw/data/assay_document.csv.gz",
        assay = "{DATA}/assay/data/assay_{cid}.csv.gz"
    output: "{DATA}/assay_document/data/assay_document_{cid}.csv.gz"
    log: "{DATA}/assay_document/log/assay_document_{cid}.log"
    shell: "fct_assay_document {input.assay_document_raw} {input.assay} {output} >{log} 2>&1"


rule MOL_DOC:
    input:
        molecule_document_raw = WD + "/raw/data/molecule_document.csv.gz",
        molecule = "{DATA}/molecule/data/molecule_{cid}.csv.gz"
    output: "{DATA}/molecule_document/data/molecule_document_{cid}.csv.gz"
    log: "{DATA}/molecule_document/log/molecule_document_{cid}.log"
    shell: "fct_molecule_document {input.molecule_document_raw} {input.molecule} {output} >{log} 2>&1"


rule ASSAY:
    input:
        assay_raw = WD + "/raw/data/assay.csv.gz",
        molecule_assay = "{DATA}/molecule_assay/data/molecule_assay_{cid}.csv.gz"
    output: "{DATA}/assay/data/assay_{cid}.csv.gz"
    log: "{DATA}/assay/log/assay_{cid}.log"
    shell: "fct_assay {input.assay_raw} {input.molecule_assay} {output} >{log} 2>&1"


rule MOL_ASSAY:
    input:
        mol_assay_raw = WD + "/raw/data/molecule_assay.csv.gz",
        molecule = "{DATA}/molecule/data/molecule_{cid}.csv.gz"
    output: "{DATA}/molecule_assay/data/molecule_assay_{cid}.csv.gz"
    log: "{DATA}/molecule_assay/log/molecule_assay_{cid}.log"
    shell: "fct_molecule_assay {input.mol_assay_raw} {input.molecule} {output} >{log} 2>&1"


rule FRAG_FRAG:
    input: root_dir + "/data/" + prep_subdir + "/" + natref_subdir + "/" + frags_subdir + "/10_pnp/data/" + prefix + "_{cid}_pnp.csv.gz"
    output: DATA + "/fragment_fragment/data/fragment_fragment_{cid}.csv.gz"
    log: DATA + "/fragment_fragment/log/fragment_fragment_{cid}.log"
    shell: "fct_fragment_fragment {input} {output}  >{log} 2>&1"


rule MOL_DATASET:
    input:
        molecule = "{DATA}/molecule/data/molecule_{cid}.csv.gz",
        dataset = "{DATA}/dataset/data/dataset.csv.gz"
    output: "{DATA}/molecule_dataset/data/molecule_dataset_{cid}.csv.gz"
    log: "{DATA}/molecule_dataset/log/molecule_dataset_{cid}.log"
    shell: "fct_molecule_dataset synthetic {input.molecule} {input.dataset} {output}  >{log} 2>&1"


rule DATASET:
    input: config_file
    output: DATA + "/dataset/data/dataset.csv.gz"
    log: DATA + "/dataset/log/dataset.log"
    shell: "fct_dataset synthetic {input} {output}  >{log} 2>&1"


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
        latest_step = root_dir + "/data/" + prep_subdir + "/" + natref_subdir + "/" + frags_subdir + "/10_pnp/data/" + prefix + "_{cid}_pnp.csv.gz",
        commercial_ref = commercial_ref
    output: DATA + "/molecule/data/molecule_{cid}.csv.gz"
    log: DATA + "/molecule/log/molecule_{cid}.log"
    shell: "fct_molecule {input.load_step} {input.latest_step} {input.commercial_ref} {output}  >{log} 2>&1"
