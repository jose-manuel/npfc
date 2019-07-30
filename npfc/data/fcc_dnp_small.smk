#!/usr/bin/env python

"""
Script fcc_dnp.smk
===========================
This script is used running the FCC protocol on the DNP dataset.
"""

from pathlib import Path
import shutil
import pkg_resources


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


data_ori = "tests/data/dnp/"
data_tgt = "tests/tmp/dnp/"
if Path(data_tgt).exists():
    shutil.rmtree(data_tgt)
shutil.copytree(data_ori, data_tgt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USER INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# Pick either of ROOT folders:
# ROOT = "/home/gally/Projects/NPFC/data/"  # local (testing)
# ROOT = "/scratch/josemanuel.gally/"  # remote (production)
ROOT = "tests/tmp/"

# folder where the magic happens
WD = ROOT + "dnp/data/"

# id to use for tracking molecules
molid = "UKEY"


# define how many chunks and thus jobs will have to be produced
# WARNING! you have to calculate yourself how many chunks are expected!
num_records_per_chunk = 100
num_chunks = 3


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/deglyco_mols.knwf')

# define chunk_ids for wildcard expansionWcx7g5!Qu
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]

# expected outputs at key points of the pipeline
INPUT_OUT = [WD + "1_input/data/dnp_" + cid + ".sdf.gz" for cid in chunk_ids]  # chunk_sdf
MAP_OUT = [WD + "9_map/data/dnp_" + cid + "_map.csv.gz" for cid in chunk_ids]  # map_frags


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: MAP_OUT

rule MAP:
    priority: 11
    input: WD + "8_fcc/data/dnp_{cid}_fcc.csv.gz"
    output: WD + "9_map/data/dnp_{cid}_map.csv.gz"
    log: WD + "9_map/log/dnp_{cid}_map.log"
    shell: "map_frags {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 12
    input: WD + "7_sub/data/dnp_{cid}_sub.csv.gz"
    output: WD + "8_fcc/data/dnp_{cid}_fcc.csv.gz"
    log: WD + "8_fcc/log/dnp_{cid}_fcc.log"
    shell: "classify_frags {input} {output} -c 3 >{log} 2>&1"

rule SUB:
    priority: 13
    input:
        mols = WD + "6_gen2d/data/dnp_{cid}_gen2d.csv.gz",
        frags = ROOT + "scaffolds/crms/data/8_gen2d/data/crms_2d.csv.gz"
    output: WD + "7_sub/data/dnp_{cid}_sub.csv.gz"
    log: WD + "7_sub/log/dnp_{cid}_sub.log"
    shell: "substruct_mols {input.mols} {input.frags} {output} >{log} 2>&1"

rule GEN2D:
    priority: 14
    input: WD + "5_uni/data/dnp_{cid}_uni.csv.gz"
    output: WD + "6_gen2d/data/dnp_{cid}_gen2d.csv.gz"
    log: WD + "6_gen2d/log/dnp_{cid}_gen2d.log"
    shell: "compute2D_mols {input} {output} 2>{log}"

rule UNI:
    priority: 15
    input: WD + "4_std/data/dnp_{cid}_passed.csv.gz"
    output: WD + "5_uni/data/dnp_{cid}_uni.csv.gz"
    log: WD + "5_uni/log/dnp_{cid}_uni.log"
    shell: "filter_dupl_mols {input} {output} -r " + WD + "5_uni/dnp_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 16
    input: WD + "3_deglyco/data/dnp_{cid}_deglyco.sdf.gz"
    output:
        WD + "4_std/data/dnp_{cid}_passed.csv.gz",
        WD + "4_std/data/dnp_{cid}_filtered.csv.gz",
        WD + "4_std/data/dnp_{cid}_error.csv.gz"
    log: WD + "4_std/log/dnp_{cid}_std.log"
    shell: "standardize_mols {input} $(echo {output[0]} | rev | cut -d_ -f2- | rev).csv.gz >{log}  2>&1"

rule DGC:
    priority: 17
    input: WD + "2_load/data/dnp_{cid}.sdf.gz"
    output: WD + "3_deglyco/data/dnp_{cid}_deglyco.sdf.gz"
    log: WD + "3_deglyco/log/dnp_{cid}_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "deglyco_mols -s {input} -i " + molid + " -o " + WD + "3_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 18
    input: WD + "1_input/data/dnp_{cid}.sdf.gz"
    output: WD + "2_load/data/dnp_{cid}.sdf.gz"
    log: WD + "2_load/log/dnp_{cid}_load.log"
    shell: "load_mols {input} {output} --in_id " + molid + " >{log} 2>&1"

rule INPUT:
    priority: 19
    input: WD + "0_raw/data/dnp.sdf.gz"
    output: INPUT_OUT
    log: WD + "1_input/log/dnp_chunk.log"
    shell: "chunk_sdf -i {input} -n " + str(num_records_per_chunk) + " -o " + WD + "1_input/data/ >{log} 2>&1"
