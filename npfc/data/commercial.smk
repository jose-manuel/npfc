#!/usr/bin/env python

"""
Script purchasable.smk
===========================
The pipeline to apply for purchasable scenario.
"""

from pathlib import Path
import shutil
import pkg_resources
from math import ceil
from npfc import load


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# get the arguments
WD = config['WD']
molid = config['molid']
try:
    prefix = config['prefix']
except KeyError:
    prefix = ''
input_file = config['input_file']
chunksize = config['chunksize']
# preprocess subdir
try:
    prep_subdir = config['prep_subdir']
except KeyError:
    prep_subdir = 'prep'

# protocol for standardizing molecules
fallback_default_std_mols = False
try:
    config_std_mols = config['std_protocol']
    if config_std_mols == '' or config_std_mols == 'DEFAULT':
        fallback_default_std_mols = True
except KeyError:
    fallback_default_std_mols = True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'
# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(config['num_chunks'])]

# fall back to default std configuration in case either missing from JSON file or empty string
if fallback_default_std_mols:
    config_std_mols = pkg_resources.resource_filename('npfc', 'data/std_mols.json')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(f"{WD}/{prep_subdir}/05_depict/data/{prefix}" + '_{cid}_depict.csv.gz', cid=chunk_ids)

rule DEPICT:
    priority: 14
    input: ancient("{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz")
    output: "{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz"
    log: "{WD}/{prep_subdir}/05_depict/log/{prefix}_{cid}_depict.log"
    shell: "mols_depict {input} {output} 2>{log}"

rule DEDUPL:
    priority: 15
    input: ancient("{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz")
    output:
        passed = "{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz",
        filtered = "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_filtered.csv.gz",
        synonyms = "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_synonyms.csv.gz"
    log: "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_dedupl.log"
    shell: "mols_dedupl {input} {output.passed} -d {output.filtered} -s {output.synonyms} -r {WD}/{prep_subdir}/04_dedupl/{prefix}_ref.hdf 2>{log}"

rule STD:
    priority: 16
    input: ancient(WD + "/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz")
    output:
        std = "{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz",
        filtered = "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} -p " + config_std_mols + " 2>{log}"  # mols_standardize takes a dir as output

rule LOAD:
    priority: 18
    input: ancient("{WD}/{prep_subdir}/01_chunk/data/{prefix}_{cid}.sdf.gz")
    output: "{WD}/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz"
    log: "{WD}/{prep_subdir}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 19
    input: ancient(input_file)
    output: expand(WD + '/' + prep_subdir + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + '/' + prep_subdir + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -p '{prefix}' -o {WD}/{prep_subdir}/01_chunk/data/ >{log} 2>&1"
