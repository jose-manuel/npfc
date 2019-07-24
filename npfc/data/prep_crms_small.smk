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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


ROOT = "tests/tmp/"
data_ori = "tests/data/scaffolds/"
data_tgt = ROOT + "scaffolds/"
# reset folder in case it already exists
if Path(data_tgt).exists():
    shutil.rmtree(data_tgt)
shutil.copytree(data_ori, data_tgt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/deglyco_mols.knwf')

# WD = "/home/gally/Projects/NPFC/data/scaffolds/crms/"  # local (testing)
# WD = "/scratch/josemanuel.gally/scaffolds/crms/data/"  # production on cluster
WD = data_tgt + "crms/data/"  # tests

# id to use for tracking molecules
molid = "Cluster"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: WD + "8_gen2d/data/crms_2d.csv.gz"

rule GEN2D:
    priority: 107
    input: WD + "7_unims/data/crms_uni.csv.gz"
    output: WD + "8_gen2d/data/crms_2d.csv.gz"
    log: WD + "8_gen2d/log/crms_2d.log"
    shell: "compute2D_mols {input} {output} 2>{log}"

rule UNIMS:
    priority: 106
    input: WD + "6_stdms/data/crms_passed.csv.gz"
    output: WD + "7_unims/data/crms_uni.csv.gz"
    log: WD + "7_unims/log/crms_uni.log"
    shell: "filter_dupl_mols {input} {output} -r " + WD + "7_unims/crms_ref.hdf 2>{log}"

rule STDMS:
    priority: 105
    input: WD + "5_murcko/data/cr_murcko.csv.gz"
    output:
        WD + "6_stdms/data/crms_passed.csv.gz",
        WD + "6_stdms/data/crms_filtered.csv.gz",
        WD + "6_stdms/data/crms_error.csv.gz"
    log: WD + "6_stdms/log/crms_std.log"
    shell: "standardize_mols {input} " + WD + "6_stdms/data/crms.csv.gz 2>{log}"

rule EMS:
    priority: 104
    input: WD + "4_uni/data/cr_uni.csv.gz"
    output: WD + "5_murcko/data/cr_murcko.csv.gz"
    log: WD + "5_murcko/log/cr_murcko.log"
    shell: "murcko_mols {input} {output} 2>{log}"

rule UNI:
    priority: 103
    input: WD + "3_std/data/cr_passed.csv.gz"
    output: WD + "4_uni/data/cr_uni.csv.gz"
    log: WD + "4_uni/log/cr_uni.log"
    shell: "filter_dupl_mols {input} {output} -r " + WD + "4_uni/cr_ref.hdf 2>{log}"

rule STD:
    priority: 102
    input: WD + "2_deglyco/data/cr_deglyco.sdf.gz"
    output:
        WD + "3_std/data/cr_passed.csv.gz",
        WD + "3_std/data/cr_filtered.csv.gz",
        WD + "3_std/data/cr_error.csv.gz"
    log: WD + "3_std/log/cr_std.log"
    shell: "standardize_mols {input} " + WD + "3_std/data/cr.csv.gz 2>{log}"

rule DGC:
    priority: 101
    input: WD + "1_load/data/cr.sdf.gz"
    output: WD + "2_deglyco/data/cr_deglyco.sdf.gz"
    log: WD + "2_deglyco/log/cr_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "deglyco_mols -s {input} -i Cluster -o " + WD + "2_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 100
    input: WD + "0_raw/data/cr.sdf.gz"
    output: WD + "1_load/data/cr.sdf.gz"
    log: WD + "1_load/log/cr_load.log"
    shell: "load_mols {input} {output} --in_id " + molid + " 2>{log}"
