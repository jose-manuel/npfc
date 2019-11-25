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
natref = config['natref']
# for activity
act_file_raw = config['act_file']  # raw file with activity for annotating fmaps. For now only works with the ChEMBL


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'
if natref.endswith('/'):
    natref = natref[:-1]

natref_uni_reffile = f"{config['natref']}/data/05_uni/dnp_ref.hdf"
natref_fmap_dir = f"{config['natref']}/data/09_fmap/data/"
natref_fmap_files = [str(x) for x in Path(natref_fmap_dir).glob('*')]


num_chunks = config['num_chunks']
# define chunk_ids for wildcard expansionWcx7g5!Qu
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(WD + '/12_pnp/data/' + prefix + '_{cid}_pnp.csv.gz', cid=chunk_ids)

rule PNP:
    priority: 0
    input: "{WD}/11_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    output: "{WD}/12_pnp/data/{prefix}_{cid}_pnp.csv.gz"
    log: "{WD}/12_pnp/log/{prefix}_{cid}_pnp.log"
    shell: "fmaps_annotate_pnp {input} {natref_fmap_dir} {output} >{log} 2>&1"

rule MAP:
    priority: 1
    input: "{WD}/10_act/data/{prefix}_{cid}_act.csv.gz"
    output: "{WD}/11_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    log: "{WD}/11_fmap/log/{prefix}_{cid}_fmap.log"
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
        expand(WD + '/05_uni/data/' + prefix + '_{cid}_uni.csv.gz', cid=chunk_ids)
    output: "{WD}/10_act/{prefix}_act_prep.csv.gz"
    log: "{WD}/10_act/{prefix}_act_prep.log"
    shell: "act_preprocess_cpa {input[0]} {output} {WD}/05_uni/{prefix}_ref.hdf {WD}/05_uni/log >{log} 2>&1"

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
    input:
        mols = "{WD}/05_uni/data/{prefix}_{cid}_uni.csv.gz",
        ref = natref + "/data/05_uni/dnp_ref.hdf"
    output: "{WD}/06_synth/data/{prefix}_{cid}_synth.csv.gz"
    log: "{WD}/06_synth/log/{prefix}_{cid}_sub.log"
    shell: "mols_subset {input.mols} {input.ref} {output} >{log} 2>&1"

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
    input: input_file
    output: expand(WD + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1"
