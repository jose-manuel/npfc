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
# preprocess subdir
try:
    prep_subdir = config['prep_subdir']
except KeyError:
    prep_subdir = 'prep'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')
# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'
# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(config['num_chunks'])]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(f"{WD}/{prep_subdir}/{frags_subdir}/08_fcg/data/{prefix}" + '_{cid}_fcg.csv.gz', cid=chunk_ids)

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
    shell: "frags_search {input.mols} {input.frags} {output} >{log} 2>&1"

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
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} 2>{log}"  # mols_standardize takes a dir as output

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
