#!/usr/bin/env python

"""
Script fc_fragments.smk
===========================
This script is used preparing the cluster representative fragments for the FCC
of the DNP and ChEMBL datasets.
"""

from pathlib import Path
import shutil
import pkg_resources


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

configfile: "/home/gally/Projects/NPFC/src/tests/tmp/fc/01_fragments/crms/test_fragments_crms.json"

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
    config_std_frags = config['protocol_std']
    if config_std_frags == '' or config_std_frags == 'DEFAULT':
        fallback_default_std_frags = True
except KeyError:
    fallback_default_std_frags = True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


configfile: "/home/gally/Projects/NPFC/src/tests/tmp/fc/01_fragments/crms/test_fragments_crms.json"

if config['WD'].endswith('/'):
    config['WD'] = config['WD'][:-1]
if config['WD'].split('/')[-1] != 'data':
    config['WD'] += '/data'

print(config)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: f"{config['WD']}/{config['prep_subdir']}/01_load/data/{config['prefix']}.csv.gz"  # rule all does not accept wildcards


rule STD_MURCKO:
    priority: 104
    input: "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz"
    output:
        std = "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz",
        filtered = "{WD}/{prep_subdir}/02_std/log/{prefix}_filtered.csv.gz",
        error = "{WD}/{prep_subdir}/02_std/log/{prefix}_error.csv.gz"
    log: "{WD}/{prep_subdir}/02_std/log/{prefix}_std.log"
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} -p " + config_std_frags + " 2>{log}"

rule LOAD:
    priority: 105
    input: config['input_file']
    output: f"{config['WD']}/{config['prep_subdir']}/01_load/data/{config['prefix']}.csv.gz"
    log: f"{config['WD']}/{config['prep_subdir']}/01_load/log/{config['prefix']}.log"
    shell: "echo '' > {output}"
