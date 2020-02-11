#!/usr/bin/env python

"""
Script prep_crms.smk
===========================
This script is used preparing the cluster representative fragments for the FCC
of the DNP and ChEMBL datasets.
"""

from pathlib import Path
import shutil
import pkg_resources


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


WD = config['WD']
molid = config['molid']
prefix = config['prefix']
input_file = config['input_file']

fallback_default_std_frags = False
try:
    config_std_frags = config['config_std_frags']
    if config_std_frags == '':
        fallback_default_std_frags = True
except KeyError:
    fallback_default_std_frags = True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')
# remove trailing / from WD if any
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'

# fall back to default std configuration in case either missing from JSON file or empty string
if fallback_default_std_frags:
    config_std_frags = pkg_resources.resource_filename('npfc', 'data/std_fragments.json')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: WD + "/08_depict/data/" + prefix + "_depict.csv.gz"  # rule all does not accept wildcards

rule DEPICT:
    priority: 100
    input: "{WD}/07_dedupl/data/{prefix}_dedupl.csv.gz"
    output: "{WD}/08_depict/data/{prefix}_depict.csv.gz"
    log: "{WD}/08_depict/log/{prefix}_depict.log"
    shell: "mols_depict {input} {output} 2>{log}"

rule DEDUPLMS:
    priority: 101
    input: "{WD}/06_concat/data/{prefix}_concat.csv.gz"
    output: "{WD}/07_dedupl/data/{prefix}_dedupl.csv.gz"
    log: "{WD}/07_dedupl/log/{prefix}_dedupl.log"
    shell: "mols_dedupl {input} {output} -r " + "{WD}/07_dedupl/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule CONCAT:
    priority: 102
    input:
        "{WD}/05b_std/data/{prefix}_passed.csv.gz",  # first input takes priority for scaffolds with same hac
        "{WD}/04a_std/data/{prefix}_passed.csv.gz"
    output: "{WD}/06_concat/data/{prefix}_concat.csv.gz"
    log: "{WD}/06_concat/log/{prefix}_concat.log"
    shell: "mols_concat {input[0]} {input[1]} {output} -s 'idm:a, hac:a, dataset:a' 2>{log}"

rule STDMS_A:
    priority: 103
    input: "{WD}/03a_murcko/data/{prefix}_murcko.csv.gz"
    output:
        "{WD}/04a_std/data/{prefix}_passed.csv.gz",
        "{WD}/04a_std/data/{prefix}_filtered.csv.gz",
        "{WD}/04a_std/data/{prefix}_error.csv.gz"
    log: "{WD}/04a_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {WD}/04a_std/data/{prefix}.csv.gz 2>{log}"  # mols_standardize takes a dir as output

rule MURCKO_A:
    priority: 104
    input: "{WD}/02_deglyco/data/{prefix}_deglyco.sdf.gz"
    output: "{WD}/03a_murcko/data/{prefix}_murcko.csv.gz"
    log: "{WD}/03a_murcko/log/{prefix}_murcko.log"
    shell: "mols_extract_murcko {input} {output} 2>{log}"

rule STDMS_B:
    priority: 105
    input: "{WD}/04b_murcko/data/{prefix}_murcko.csv.gz"
    output:
        "{WD}/05b_std/data/{prefix}_passed.csv.gz",
        "{WD}/05b_std/data/{prefix}_filtered.csv.gz",
        "{WD}/05b_std/data/{prefix}_error.csv.gz"
    log: "{WD}/05b_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {WD}/05b_std/data/{prefix}.csv.gz 2>{log}"  # mols_standardize takes a dir as output

rule MURCKO_B:
    priority: 106
    input: "{WD}/03b_std/data/{prefix}_passed.csv.gz"
    output: "{WD}/04b_murcko/data/{prefix}_murcko.csv.gz"
    log: "{WD}/04b_murcko/log/{prefix}_murcko.log"
    shell: "mols_extract_murcko {input} {output} 2>{log}"

rule STD_B:
    priority: 107
    input: "{WD}/02_deglyco/data/{prefix}_deglyco.sdf.gz"
    output:
        "{WD}/03b_std/data/{prefix}_passed.csv.gz",
        "{WD}/03b_std/data/{prefix}_filtered.csv.gz",
        "{WD}/03b_std/data/{prefix}_error.csv.gz"
    log: "{WD}/03b_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {WD}/03b_std/data/{prefix}.csv.gz -p " + config_std_frags + " 2>{log}"  # mols_standardize takes a dir as output

rule DGC:
    priority: 108
    input: "{WD}/01_load/data/{prefix}.sdf.gz"
    output: "{WD}/02_deglyco/data/{prefix}_deglyco.sdf.gz"
    log: "{WD}/02_deglyco/log/{prefix}_deglyco_knwf.log"  # log from the KNIME execution, another log is computed by the workflow
    shell: "mols_deglyco -s {input} -i {molid} -o {WD}/02_deglyco -w {file_knwf} >{log} 2>&1"

rule LOAD:
    priority: 109
    input: input_file
    output: "{WD}/01_load/data/{prefix}.sdf.gz"
    log: "{WD}/01_load/log/{prefix}.log"
    shell: "mols_load {input} {output} --in_id {molid} 2>{log}"
