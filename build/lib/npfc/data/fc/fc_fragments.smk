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


# remove trailing / from WD if any
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'

# fall back to default std configuration in case either missing from JSON file or empty string
if fallback_default_std_frags:
    config_std_frags = pkg_resources.resource_filename('npfc', 'data/std_fragments.json')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input:
        mols = WD + '/' + prep_subdir + "/05_fcp/data/" + prefix + "_fcp.csv.gz",  # rule all does not accept wildcards
        count_mols = WD + '/' + prep_subdir + '/report/data/' + prefix + '_count_mols.csv',
        time = WD + '/' + prep_subdir + '/report/data/' + prefix + '_time.csv'


rule REPORT_TIME:
    priority: 100
    input:
        load = "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz",
        std_passed = "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz",
        dedupl = "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz",
        depict = "{WD}/{prep_subdir}/04_depict/data/{prefix}_depict.csv.gz",
        fcp = "{WD}/{prep_subdir}/05_fcp/data/{prefix}_fcp.csv.gz"
    output: "{WD}/{prep_subdir}/report/data/{prefix}_time.csv"
    log: "{WD}/{prep_subdir}/report/log/{prefix}_time.log"
    shell: "report_time {WD}/{prep_subdir} '{prefix}*' {output} -p {prep_subdir} 2>{log}"


rule COUNT_MOLS:
    priority: 100
    input:
        load = "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz",
        std_passed = "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz",
        dedupl = "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz",
        depict = "{WD}/{prep_subdir}/04_depict/data/{prefix}_depict.csv.gz",
        fcp = "{WD}/{prep_subdir}/05_fcp/data/{prefix}_fcp.csv.gz"
    output: "{WD}/{prep_subdir}/report/data/{prefix}_count_mols.csv"
    log: "{WD}/{prep_subdir}/report/log/{prefix}_count_mols.log"
    shell: "mols_count {WD}/{prep_subdir} {prefix}* {output} 2>{log}"


rule FCP:
    priority: 101
    input: frags = "{WD}/{prep_subdir}/04_depict/data/{prefix}_depict.csv.gz"
    output:
        frags = "{WD}/{prep_subdir}/05_fcp/data/{prefix}_fcp.csv.gz",
        counts = "{WD}/{prep_subdir}/05_fcp/log/{prefix}_fcp_symcounts.csv"
    log: "{WD}/{prep_subdir}/05_fcp/log/{prefix}_fcp.log"
    shell: "frags_annotate_fcp {input} {output.frags} -c {output.counts} 2>{log}"


rule DEPICT:
    priority: 102
    input: "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz"
    output: "{WD}/{prep_subdir}/04_depict/data/{prefix}_depict.csv.gz"
    log: "{WD}/{prep_subdir}/04_depict/log/{prefix}_depict.log"
    shell: "mols_depict {input} {output} -m rdDepictor 2>{log}"


rule DEDUPL:
    priority: 103
    input: "{WD}/{prep_subdir}/02_std/data/{prefix}_std.csv.gz"
    output:
        passed = "{WD}/{prep_subdir}/03_dedupl/data/{prefix}_dedupl.csv.gz",
        filtered = "{WD}/{prep_subdir}/03_dedupl/log/{prefix}_filtered.csv.gz",
        synonyms = "{WD}/{prep_subdir}/03_dedupl/log/{prefix}_synonyms.csv.gz"
    log: "{WD}/{prep_subdir}/03_dedupl/log/{prefix}_dedupl.log"
    shell: "mols_dedupl {input} {output.passed} -d {output.filtered} -s {output.synonyms} -r {WD}/{prep_subdir}/03_dedupl/{prefix}_ref.hdf 2>{log}"


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
    input: input_file
    output: "{WD}/{prep_subdir}/01_load/data/{prefix}.csv.gz"
    log: "{WD}/{prep_subdir}/01_load/log/{prefix}.log"
    shell: "mols_load {input} {output} --in_id {molid} 2>{log}"
