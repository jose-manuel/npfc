#!/usr/bin/env python

"""
Script natural.smk
===========================
The pipeline to apply for natural scenarii. This is similar to the Synthetic scenario:
molecules are standardized and FCs are categorized, but Natural Products are not filtered out and PNPs are not computed.
"""

from pathlib import Path
import shutil
import pkg_resources
from math import ceil
from npfc import load
import pandas as pd
import time


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# get the arguments
WD = config['WD']
molid = config['molid']
try:
    prefix = config['prefix']
except KeyError:
    prefix = ''
input_file = config['input_file']
frags_file = config['frags_file']
frags_subdir = config['frags_subdir']
chunksize = config['chunksize']
tautomer = config['tautomer']
# preprocess subdir
try:
    prep_subdir = config['prep_subdir']
except KeyError:
    prep_subdir = 'prep'

# protocol for standardizing molecules
fallback_default_std_mols = False
try:
    config_std_mols = config['protocol_std']
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
    input:
        expand(f"{WD}/{prep_subdir}/{frags_subdir}/08_fcg/data/{prefix}" + '_{cid}_fcg.csv.gz', cid=chunk_ids),
        count_mols = '/'.join([WD, prep_subdir, frags_subdir]) + '/report/data/' + prefix + '_count_mols.csv',
        time = WD + '/' + prep_subdir + '/' + frags_subdir + '/report/data/' + prefix + '_time.csv'


rule REPORT_TIME_SUM:
    priority: 20
    input: expand(f"{WD}/{prep_subdir}/{frags_subdir}/report/data/{prefix}" + '_{cid}_time.csv', cid=chunk_ids)
    output: WD + '/' + prep_subdir + '/' + frags_subdir + '/report/data/' + prefix + '_time.csv'
    run:
        # refined input/output files so that are just plain strings and not hidden smk magic
        input_files = list(input)
        output_file = str(output)
        # gather result
        df = pd.concat([pd.read_csv(x, sep='|') for x in input_files]).rename({'pattern': 'subset'}, axis=1)
        df['subset'] = df['subset'].map(lambda x: x.replace('*', ''))
        # total
        df_tot = pd.DataFrame(df.sum()).T
        df_tot['subset'] = 'total'
        # merge data
        df = pd.concat([df, df_tot])
        # time.sleep(1)
        print(df)
        # save data
        df.to_csv(output_file, sep='|', index=False)
        # in case the operation above succeeded, delere temporary files
        if Path(output_file).exists():
            [Path(f).unlink() for f in input_files]


rule REPORT_TIME:
    priority: 100
    input:
        load = "{WD}/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz",
        std_passed = "{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz",
        dedupl = "{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz",
        depict = "{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz",
        fsearch = "{WD}/{prep_subdir}/{frags_subdir}/06_fs/data/{prefix}_{cid}_fs.csv.gz",
        fcc = "{WD}/{prep_subdir}/{frags_subdir}/07_fcc/data/{prefix}_{cid}_fcc.csv.gz",
        fcg = "{WD}/{prep_subdir}/{frags_subdir}/08_fcg/data/{prefix}_{cid}_fcg.csv.gz"
    output: "{WD}/{prep_subdir}/{frags_subdir}/report/data/{prefix}_{cid}_time.csv"
    log: "{WD}/{prep_subdir}/{frags_subdir}/report/log/{prefix}_{cid}_time.log"
    shell: "report_time {WD}/{prep_subdir} '{prefix}_{wildcards.cid}*' {output} -p {prep_subdir} -f {frags_subdir} 2>{log}"


rule COUNT_MOLS_SUM:
    priority: 0
    input: expand('/'.join([WD, prep_subdir, frags_subdir, "report/data", prefix]) + '_{cid}_count_mols.csv', cid=chunk_ids)
    output: '/'.join([WD, prep_subdir, frags_subdir]) + '/report/data/' + prefix + '_count_mols.csv'
    run:
        # refined input/output files so that are just plain strings and not hidden smk magic
        input_files = list(input)
        output_file = str(output)
        # gather result
        df = pd.concat([pd.read_csv(x, sep='|') for x in input_files]).rename({'pattern': 'subset'}, axis=1)
        df['subset'] = df['subset'].map(lambda x: x.replace('*', ''))
        # total
        df_tot = pd.DataFrame(df.sum()).T
        df_tot['subset'] = 'total'
        # merge data
        df = pd.concat([df, df_tot])
        time.sleep(0.2)  # make sure that this rule result (print below) is displayed last
        print(df)
        # save data
        df.to_csv(output_file, sep='|', index=False)
        # in case the operation above succeeded, delere temporary files
        if Path(output_file).exists():
            [Path(f).unlink() for f in input_files]


rule COUNT_MOLS:
    priority: 0
    input:
        chunk = "{WD}/{prep_subdir}/01_chunk/data/{prefix}_{cid}.sdf.gz",
        load = "{WD}/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz",
        std_passed = "{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz",
        dedupl = "{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz",
        depict = "{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz",
        fsearch = "{WD}/{prep_subdir}/" + frags_subdir + "/06_fs/data/{prefix}_{cid}_fs.csv.gz",
        fcc = "{WD}/{prep_subdir}/" + frags_subdir + "/07_fcc/data/{prefix}_{cid}_fcc.csv.gz",
        fcg = "{WD}/{prep_subdir}/" + frags_subdir + "/08_fcg/data/{prefix}_{cid}_fcg.csv.gz",
    output: "{WD}/{prep_subdir}/" + frags_subdir + "/report/data/{prefix}_{cid}_count_mols.csv"
    log: "{WD}/{prep_subdir}/" + frags_subdir + "/report/log/{prefix}_{cid}_count_mols.log"
    shell: "mols_count {WD}/{prep_subdir} {prefix}_{wildcards.cid}* {output} 2>{log}"


rule FCG:
    priority: 11
    input: ancient("{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/07_fcc/data/{prefix}_{cid}_fcc.csv.gz")
    output: "{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/08_fcg/data/{prefix}_{cid}_fcg.csv.gz"
    log: "{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/08_fcg/log/{prefix}_{cid}_fcg.log"
    shell: "fcg_generate {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"


rule FCC:
    priority: 12
    input: ancient("{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/06_fs/data/{prefix}_{cid}_fs.csv.gz")
    output: "{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/07_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}" + f"/{prep_subdir}/{frags_subdir}" + "/07_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"


rule FS:
    priority: 13
    input:
        mols = ancient("{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz"),
        frags = ancient(frags_file)
    output: "{WD}/{prep_subdir}" + f"/{frags_subdir}" + "/06_fs/data/{prefix}_{cid}_fs.csv.gz"
    log: "{WD}/{prep_subdir}" + f"/{frags_subdir}" + "/06_fs/log/{prefix}_{cid}_fs.log"
    shell: "frags_search {input.mols} {input.frags} {output} -t " + f"{tautomer}" + " >{log} 2>&1"


rule DEPICT:
    priority: 14
    input: ancient("{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz")
    output: "{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz"
    log: "{WD}/{prep_subdir}/05_depict/log/{prefix}_{cid}_depict.log"
    shell: "mols_depict {input} {output} -m rdDepictor 2>{log}"


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
