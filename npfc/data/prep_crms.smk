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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: WD + "3_murcko/data/crms_passed.csv.gz"

rule MURCKO:
    priority: 102
    input: WD + "2_deglyco/data/cr_deglyco.sdf.gz"
    output: WD + "3_murcko/data/crms_passed.csv.gz"
    log: WD + "3_murcko/log/crms_murcko.log"
    shell: "standardize_mols {input} -o $(echo {output} | rev | cut -d_ -f2- | rev).csv.gz -r " + WD + "3_murcko/data/crms_ref.hdf -m True --compute_2D True 2>{log}"

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
    shell: "load_mols {input} {output} --in_id Cluster 2>{log}"
