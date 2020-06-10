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
# preprocess subdir
try:
    prep_subdir = config['prep_subdir']
except KeyError:
    prep_subdir = 'prep'

# protocol for standardizing fragments
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
    input: WD + '/' + prep_subdir + "/04_depict/data/" + prefix + "_depict.csv.gz"  # rule all does not accept wildcards

rule DEPICT:
    priority: 100
    input: "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz"
    output: "{WD}/{prep_subdir}/04_depict/data/{prefix}_depict.csv.gz"
    log: "{WD}/{prep_subdir}/04_depict/log/{prefix}_depict.log"
    shell: "mols_depict {input} {output} 2>{log}"

rule DEDUPL:
    priority: 101
    input: "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz"
    output:
        passed = "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz",
        filtered = "{WD}/{prep_subdir}/03_dedupl/log/{prefix}_filtered.csv.gz"
    log: "{WD}/{prep_subdir}/03_dedupl/log/{prefix}_dedupl.log"
    shell: "mols_dedupl {input} {output.passed} -d {output.filtered} -r {WD}/{prep_subdir}/03_dedupl/{prefix}_ref.hdf 2>{log}"

rule STD_MURCKO:
    priority: 102
    input: "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz"
    output:
        std = "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz",
        filtered = "{WD}/{prep_subdir}/02_std/log/{prefix}_filtered.csv.gz",
        error = "{WD}/{prep_subdir}/02_std/log/{prefix}_error.csv.gz"
    log: "{WD}/{prep_subdir}/02_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} -p " + config_std_frags + " 2>{log}"

rule LOAD:
    priority: 103
    input: input_file
    output: "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz"
    log: "{WD}/{prep_subdir}/01_load/log/{prefix}.log"
    shell: "mols_load {input} {output} --in_id {molid} 2>{log}"
