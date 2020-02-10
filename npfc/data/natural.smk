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
prefix = config['prefix']
input_file = config['input_file']
frags_file = config['frags_file']
frags_subdir = "frags_" + config['frags_subdir']
chunksize = config['chunksize']


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
    input: expand(f"{WD}/{frags_subdir}/09_fmap/data/{prefix}" + '_{cid}_fmap.csv.gz', cid=chunk_ids)

rule FMAP:
    priority: 11
    input: "{WD}" + f"/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    output: "{WD}" + f"/{frags_subdir}" + "/09_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    log: "{WD}" + f"/{frags_subdir}" + "/09_fmap/log/{prefix}_{cid}_fmap.log"
    shell: "fc_map {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 --log DEBUG >{log} 2>&1"

rule FCC:
    priority: 12
    input: "{WD}" + f"/{frags_subdir}" + "/07_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    output: "{WD}" + f"/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}" + f"/{frags_subdir}" + "/08_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule FSEARCH:
    priority: 13
    input:
        mols = "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz",
        frags = frags_file
    output: "{WD}" + f"/{frags_subdir}" + "/07_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    log: "{WD}" + f"/{frags_subdir}" + "/07_fsearch/log/{prefix}_{cid}_fsearch.log"
    shell: "mols_substruct {input.mols} {input.frags} {output} >{log} 2>&1"

rule GEN2D:
    priority: 14
    input: "{WD}/05_dedupl/data/{prefix}_{cid}_dedupl.csv.gz"
    output: "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz"
    log: "{WD}/06_gen2D/log/{prefix}_{cid}_gen2D.log"
    shell: "mols_gen2D {input} {output} 2>{log}"

rule DEDUPL:
    priority: 15
    input: "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz"
    output: "{WD}/05_dedupl/data/{prefix}_{cid}_dedupl.csv.gz"
    log: "{WD}/05_dedupl/log/{prefix}_{cid}_dedupl.log"
    shell: "mols_dedupl {input} {output} -r {WD}/05_dedupl/{prefix}_ref.hdf --log DEBUG 2>{log}"

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
    input: input_file
    output: expand(WD + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1"
