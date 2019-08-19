#!/usr/bin/env python

"""
Script synthetic.smk
===========================
The pipeline to apply for synthetic scenarii.
Molecules are processed exactly like for the Natural scenario but Natural Products
are filtered before Substrure Search and PNPs are annotated as final step.
"""

from pathlib import Path
import shutil
import pkg_resources
from math import ceil
from npfc import load
from npfc import utils
import logging
logging.basicConfig(level=logging.WARNING)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARGUMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# common
WD = config['WD']
molid = config['molid']
prefix = config['prefix']
input_file = config['input_file']
# additional
frags_file = config['frags_file']  # fragment file to use for substructure search
chunksize = config['chunksize']  # maximum number of molecules per chunk
# specific to synthetic
natref_uni_reffile = config['natref_uni_reffile']  # ref file for duplicate removal from natural dataset so we can define a synthetic dataset
act_file_raw = config['act_file']  # raw file with activity for annotating fmaps. For now only works with the ChEMBL
natref_fmap_dir = config['natref_fmap_dir']  # fmaps from natural dataset so that we can identify PNPs


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]

# for counting mols, need to process the unzipped file
input_file_uncompressed = input_file.split('.gz')[0]

# count mols + uncompress input file
num_mols = load.count_mols(input_file, keep_uncompressed=True)

# determine the number of chunks to generate
num_chunks = ceil(num_mols / chunksize)

# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]

# expected outputs at key points of the pipeline
INPUT_OUT = [f"{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz" for cid in chunk_ids]  # chunk_sdf
UNI_OUT = [f"{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz" for cid in chunk_ids]  # chunk_sdf
PNP_OUT = [f"{WD}/12_pnp/data/{prefix}_{cid}_pnp.csv.gz" for cid in chunk_ids]  # annotate_pnp


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input:
        PNP_OUT  # final output of the pipeline

rule PNP:
    priority: 0
    input: "{WD}/11_map/data/{prefix}_{cid}_map.csv.gz"
    output: "{WD}/12_pnp/data/{prefix}_{cid}_pnp.csv.gz"
    log: "{WD}/12_pnp/log/{prefix}_{cid}_pnp.log"
    shell: "fmaps_annotate_pnp {input} {natref_fmap_dir} {output} >{log} 2>&1"

rule MAP:
    priority: 1
    input: "{WD}/10_act/data/{prefix}_{cid}_act.csv.gz"
    output: "{WD}/11_map/data/{prefix}_{cid}_map.csv.gz"
    log: "{WD}/11_map/log/{prefix}_{cid}_map.log"
    shell: "fc_map {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule ACT:
    priority: 2
    input:
        fcs = "{WD}/09_fcc/data/{prefix}_{cid}_fcc.csv.gz",
        act = "{WD}/10_act/{prefix}_act_prep.csv.gz"
    output: "{WD}/10_act/data/{prefix}_{cid}_act.csv.gz"
    log: "{WD}/10_act/log/{prefix}_{cid}_act.log"
    shell: "act_annotate_fc {input.fcs} {input.act} {output} >{log} 2>&1"

rule ACT_PREP:
    priority: 3
    input:
        act_file_raw,  # input file
        UNI_OUT  # wait until all chunks are processed
    output: "{WD}/10_act/{prefix}_act_prep.csv.gz"
    log: "{WD}/10_act/{prefix}_act_prep.log"
    shell: "act_preprocess {input[0]} {output} {WD}/05_uni/{prefix}_ref.hdf {WD}/05_uni/log >{log} 2>&1"

rule FCC:
    priority: 4
    input: "{WD}/08_sub/data/{prefix}_{cid}_sub.csv.gz"
    output: "{WD}/09_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}/09_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule SUB:
    priority: 5
    input:
        mols = "{WD}/07_gen2D/data/{prefix}_{cid}_gen2D.csv.gz",
        frags = frags_file
    output: "{WD}/08_sub/data/{prefix}_{cid}_sub.csv.gz"
    log: "{WD}/08_sub/log/{prefix}_{cid}_sub.log"
    shell: "mols_substruct {input.mols} {input.frags} {output} >{log} 2>&1"

rule GEN2D:
    priority: 6
    input: "{WD}/06_synth/data/{prefix}_{cid}_synth.csv.gz"
    output: "{WD}/07_gen2D/data/{prefix}_{cid}_gen2D.csv.gz"
    log: "{WD}/07_gen2D/log/{prefix}_{cid}_gen2D.log"
    shell: "mols_gen2D {input} {output} 2>{log}"

rule SYNTH:
    priority: 7
    input:  "{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz"
    output: "{WD}/06_synth/data/{prefix}_{cid}_synth.csv.gz"
    log: "{WD}/06_synth/log/{prefix}_{cid}_sub.log"
    shell: "mols_subset {input} {natref_uni_reffile} {output} >{log} 2>&1"

rule UNI:
    priority: 8
    input: "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz"
    output: "{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz"
    log: "{WD}/05_uni/log/{prefix}_{cid}_uni.log"
    shell: "mols_filter_dupl {input} {output} -r {WD}/05_uni/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 9
    input: WD + "/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    output:
        passed = "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz",
        filtered = "{WD}/04_std/data/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/04_std/data/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/04_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} $(echo {output.passed} | rev | cut -d_ -f2- | rev).csv.gz >{log}  2>&1"  # {cid} unkown on this context

rule DGC:
    priority: 10
    input: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    log: "{WD}/03_deglyco/log/{prefix}_{cid}_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "mols_deglyco -s {input} -i {molid} -o {WD}/03_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 11
    input: "{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    log: "{WD}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 12
    input: input_file_uncompressed
    output: INPUT_OUT
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1; rm {input_file_uncompressed}"
