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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD if any
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: WD + "/07_gen2D/data/" + prefix + "_gen2D.csv.gz"  # rule all does not accept wildcards

rule GEN2D:
    priority: 101
    input: "{WD}/06_unims/data/{prefix}_uni.csv.gz"
    output: "{WD}/07_gen2D/data/{prefix}_gen2D.csv.gz"
    log: "{WD}/07_gen2D/log/{prefix}_gen2D.log"
    shell: "mols_gen2D {input} {output} 2>{log}"

rule UNIMS:
    priority: 102
    input: "{WD}/05_std/data/{prefix}_passed.csv.gz"
    output: "{WD}/06_unims/data/{prefix}_uni.csv.gz"
    log: "{WD}/06_unims/log/{prefix}_uni.log"
    shell: "mols_filter_dupl {input} {output} -r " + "{WD}/06_unims/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 103
    input: "{WD}/04_murcko/data/{prefix}_murcko.csv.gz"
    output:
        "{WD}/05_std/data/{prefix}_passed.csv.gz",
        "{WD}/05_std/data/{prefix}_filtered.csv.gz",
        "{WD}/05_std/data/{prefix}_error.csv.gz"
    log: "{WD}/05_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {WD}/05_std/data/{prefix}.csv.gz 2>{log}"  # mols_standardize takes a dir as output

rule MURCKO:
    priority: 104
    input: "{WD}/03_uni/data/{prefix}_uni.csv.gz"
    output: "{WD}/04_murcko/data/{prefix}_murcko.csv.gz"
    log: "{WD}/04_murcko/log/{prefix}_murcko.log"
    shell: "mols_extract_murcko {input} {output} 2>{log}"

rule UNI:
    priority: 105
    input: "{WD}/02_deglyco/data/{prefix}_deglyco.sdf.gz"
    output: "{WD}/03_uni/data/{prefix}_uni.csv.gz"
    log: "{WD}/03_uni/log/{prefix}_uni.log"
    shell: "mols_filter_dupl {input} {output} -r {WD}/03_uni/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule DGC:
    priority: 107
    input: "{WD}/01_load/data/{prefix}.sdf.gz"
    output: "{WD}/02_deglyco/data/{prefix}_deglyco.sdf.gz"
    log: "{WD}/02_deglyco/log/{prefix}_deglyco_knwf.log"  # log from the KNIME execution, another log is computed by the workflow
    shell: "mols_deglyco -s {input} -i {molid} -o {WD}/02_deglyco -w {file_knwf} >{log} 2>&1"

rule LOAD:
    priority: 108
    input: input_file
    output: "{WD}/01_load/data/{prefix}.sdf.gz"
    log: "{WD}/01_load/log/{prefix}.log"
    shell: "mols_load {input} {output} --in_id {molid} 2>{log}"
