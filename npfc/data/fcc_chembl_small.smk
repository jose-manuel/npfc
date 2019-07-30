#!/usr/bin/env python

"""
Script fcc_chembl.smk
===========================
This script is used running the FCC protocol on the chembl dataset.
"""

from pathlib import Path
import shutil
import pkg_resources


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


data_ori = "tests/data/chembl/"
data_tgt = "tests/tmp/chembl/"
if Path(data_tgt).exists():
    shutil.rmtree(data_tgt)
shutil.copytree(data_ori, data_tgt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USER INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Pick either of ROOT folders:
# ROOT = "/home/gally/Projects/NPFC/data/"  # local (testing)
# ROOT = "/scratch/josemanuel.gally/"  # remote (production)
ROOT = "tests/tmp/"

# folder where the magic happens
WD = ROOT + "chembl/data/"

# id to use for tracking molecules
molid = "chembl_id"

# define how many chunks and thus jobs will have to be produced
# WARNING! you have to calculate yourself how many chunks are expected!
num_records_per_chunk = 100
num_chunks = 3


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/deglyco_mols.knwf')

# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]

# because we have to filter ChEMBL using DNP, we also need DNP information:
num_chunks_dnp = 61
chunks_ids_dnp = [str(i+1).zfill(3) for i in range(num_chunks_dnp)]

# expected outputs at key points of the pipeline
ref_syn = ROOT + "dnp/data/5_uni/dnp_ref.hdf"
ref_dir = ROOT + "dnp/data/9_map/data/"

INPUT_OUT = [WD + "1_input/data/chembl_" + cid + ".sdf.gz" for cid in chunk_ids]  # chunk_sdf
PNP_OUT = [WD + "11_pnp/data/chembl_" + cid + "_pnp.csv.gz" for cid in chunk_ids]  # annotate_pnp


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input:
        PNP_OUT  # final output of the pipeline

rule PNP:
    priority: 0
    input: WD + "10_map/data/chembl_{cid}_map.csv.gz"
    output: WD + "11_pnp/data/chembl_{cid}_pnp.csv.gz"
    log: WD + "11_pnp/log/chembl_{cid}_pnp.log"
    shell: "annotate_pnp {input} " + ref_dir + " {output} >{log} 2>&1"

rule MAP:
    priority: 1
    input: WD + "9_fcc/data/chembl_{cid}_fcc.csv.gz"
    output: WD + "10_map/data/chembl_{cid}_map.csv.gz"
    log: WD + "10_map/log/chembl_{cid}_map.log"
    shell: "map_frags {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 2
    input: WD + "8_sub/data/chembl_{cid}_sub.csv.gz"
    output: WD + "9_fcc/data/chembl_{cid}_fcc.csv.gz"
    log: WD + "9_fcc/log/chembl_{cid}_fcc.log"
    shell: "classify_frags {input} {output} -c 3 >{log} 2>&1"

rule SUB:
    priority: 3
    input:
        mols = WD + "7_gen2d/data/chembl_{cid}_gen2d.csv.gz",
        frags = ROOT + "scaffolds/crms/data/8_gen2d/data/crms_2d.csv.gz"
    output: WD + "8_sub/data/chembl_{cid}_sub.csv.gz"
    log: WD + "8_sub/log/chembl_{cid}_sub.log"
    shell: "substruct_mols {input.mols} {input.frags} {output} >{log} 2>&1"

rule GEN2D:
    priority: 4
    input: WD + "6_synth/data/chembl_{cid}_synth.csv.gz"
    output: WD + "7_gen2d/data/chembl_{cid}_gen2d.csv.gz"
    log: WD + "7_gen2d/log/chembl_{cid}_gen2d.log"
    shell: "compute2D_mols {input} {output} 2>{log}"

rule SYNTH:
    priority: 5
    input:  WD + "5_uni/data/chembl_{cid}_uni.csv.gz"
    output: WD + "6_synth/data/chembl_{cid}_synth.csv.gz"
    log: WD + "6_synth/log/chembl_{cid}_sub.log"
    shell: "subset_mols {input} " + ref_syn + " {output} >{log} 2>&1"

rule UNI:
    priority: 6
    input: WD + "4_std/data/chembl_{cid}_passed.csv.gz"
    output: WD + "5_uni/data/chembl_{cid}_uni.csv.gz"
    log: WD + "5_uni/log/chembl_{cid}_uni.log"
    shell: "filter_dupl_mols {input} {output} -r " + WD + "5_uni/chembl_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 7
    input: WD + "3_deglyco/data/chembl_{cid}_deglyco.sdf.gz"
    output:
        WD + "4_std/data/chembl_{cid}_passed.csv.gz",
        WD + "4_std/data/chembl_{cid}_filtered.csv.gz",
        WD + "4_std/data/chembl_{cid}_error.csv.gz",

    log: WD + "4_std/log/chembl_{cid}_std.log"
    shell: "standardize_mols {input} $(echo {output[0]} | rev | cut -d_ -f2- | rev).csv.gz >{log}  2>&1"

rule DGC:
    priority: 8
    input: WD + "2_load/data/chembl_{cid}.sdf.gz"
    output: WD + "3_deglyco/data/chembl_{cid}_deglyco.sdf.gz"
    log: WD + "3_deglyco/log/chembl_{cid}_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "deglyco_mols -s {input} -i " + molid + " -o " + WD + "3_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 9
    input: WD + "1_input/data/chembl_{cid}.sdf.gz"
    output: WD + "2_load/data/chembl_{cid}.sdf.gz"
    log: WD + "2_load/log/chembl_{cid}_load.log"
    shell: "load_mols {input} {output} --in_id " + molid + " >{log} 2>&1"

rule INPUT:
    priority: 10
    input: WD + "0_raw/data/chembl.sdf.gz"
    output: INPUT_OUT
    log: WD + "1_input/log/chembl_chunk.log"
    shell: "chunk_sdf -i {input} -n " + str(num_records_per_chunk) + " -o " + WD + "1_input/data/ >{log} 2>&1"
