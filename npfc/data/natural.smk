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
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# get the arguments
WD = config['WD']
molid = config['molid']
prefix = config['prefix']
input_file = config['input_file']
frags_file = config['frags_file']
chunksize = config['chunksize']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]

# # for counting mols, need to process the unzipped file
# input_file_uncompressed = input_file.split('.gz')[0]
#
# # count mols + uncompress input file
# num_mols = load.count_mols(input_file, keep_uncompressed=True)
#
# # determine the number of chunks to generate
# num_chunks = ceil(num_mols / chunksize)
input_file_uncompressed = config['input_file_uncompressed']
num_chunks = config['num_chunks']
# define chunk_ids for wildcard expansionWcx7g5!Qu
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]

# expected outputs at key points of the pipeline
INPUT_OUT = [f"{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz" for cid in chunk_ids]  # chunk_sdf
MAP_OUT = [f"{WD}/09_fmap/data/{prefix}_{cid}_fmap.csv.gz" for cid in chunk_ids]  # fc_map


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: MAP_OUT

rule FMAP:
    priority: 11
    input: "{WD}/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    output: "{WD}/09_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    log: "{WD}/09_fmap/log/{prefix}_{cid}_fmap.log"
    shell: "fc_map {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 12
    input: "{WD}/07_sub/data/{prefix}_{cid}_sub.csv.gz"
    output: "{WD}/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}/08_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule SUB:
    priority: 13
    input:
        mols = "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz",
        frags = frags_file
    output: "{WD}/07_sub/data/{prefix}_{cid}_sub.csv.gz"
    log: "{WD}/07_sub/log/{prefix}_{cid}_sub.log"
    shell: "mols_substruct {input.mols} {input.frags} {output} >{log} 2>&1"

rule GEN2D:
    priority: 14
    input: "{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz"
    output: "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz"
    log: "{WD}/06_gen2D/log/{prefix}_{cid}_gen2D.log"
    shell: "mols_gen2D {input} {output} 2>{log}"

rule UNI:
    priority: 15
    input: "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz"
    output: "{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz"
    log: "{WD}/05_uni/log/{prefix}_{cid}_uni.log"
    shell: "mols_filter_dupl {input} {output} -r {WD}/05_uni/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 16
    input: WD + "/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    output:
        passed = "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz",
        filtered = "{WD}/04_std/data/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/04_std/data/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/04_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} $(echo {output.passed} | rev | cut -d_ -f2- | rev).csv.gz >{log}  2>&1"  # {cid} unkown on this context

rule DGC:
    priority: 17
    input: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    log: "{WD}/03_deglyco/log/{prefix}_{cid}_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "mols_deglyco -s {input} -i {molid} -o {WD}/03_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 18
    input: "{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    log: "{WD}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 19
    input: input_file_uncompressed
    output: INPUT_OUT
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1; rm {input_file_uncompressed}"
