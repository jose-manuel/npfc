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
frags_subdir = config['frags_subdir']
chunksize = config['chunksize']  # maximum number of molecules per chunk
# specific to synthetic
natref = config['natref']
act_file_raw = config['act_file']  # raw file with activity for annotating fmaps. For now only works with the ChEMBL
# from master script
num_chunks = config['num_chunks']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'
if natref.endswith('/'):
    natref = natref[:-1]

# frags
frags_filename = Path(frags_file).stem.split('.')[0]   # basename from frags file
# natref
natref_filename = Path(natref).stem   # basename from natref dir
natref_subdir = Path(natref).stem.split('.')[0]  # subfolder where outputs using nat_ref are stored
# natref dedupl - for subset
natref_dedupl_dir = list(Path(f"{natref}/data").glob("[0-9][0-9]_dedupl"))[0]
natref_dedupl_reffile = list(natref_dedupl_dir.glob('*.hdf'))[0]
# natref fmap - for pnp
natref_fmap_dir = [str(x) for x in Path(f"{natref}/data/{frags_subdir}").glob("[0-9][0-9]_fmap")][0] + '/data'
# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(f"{WD}/{natref_subdir}/{frags_subdir}/11_pnp/data/{prefix}" + '_{cid}_pnp.csv.gz', cid=chunk_ids)

rule PNP:
    priority: 0
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/11_pnp/data/{prefix}_{cid}_pnp.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/11_pnp/log/{prefix}_{cid}_pnp.log"
    shell: "fmaps_annotate_pnp {input} {natref_fmap_dir} {output} >{log} 2>&1"

rule FMAP:
    priority: 1
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_fmap/data/{prefix}_{cid}_fmap.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_fmap/log/{prefix}_{cid}_fmap.log"
    shell: "fc_map {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 2
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule FSEARCH:
    priority: 3
    input:
        mols = "{WD}" + f"/{natref_subdir}" + "/07_subset/data/{prefix}_{cid}_subset.csv.gz",
        frags = frags_file
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fsearch/log/{prefix}_{cid}_fsearch.log"
    shell: "mols_substruct {input.mols} {input.frags} {output} >{log} 2>&1"

rule SUBSET:
    priority: 4
    input:
        mols = "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz",
        ref = natref + "/data/05_dedupl/dnp_ref.hdf"
    output: "{WD}" + f"/{natref_subdir}" + "/07_subset/data/{prefix}_{cid}_subset.csv.gz"
    log: "{WD}" + f"/{natref_subdir}" + "/07_subset/log/{prefix}_{cid}_subset.log"
    shell: "mols_subset {input.mols} {input.ref} {output} >{log} 2>&1"

rule GEN2D:
    priority: 5
    input: "{WD}/05_dedupl/data/{prefix}_{cid}_dedupl.csv.gz"
    output: "{WD}/06_gen2D/data/{prefix}_{cid}_gen2D.csv.gz"
    log: "{WD}/06_gen2D/log/{prefix}_{cid}_gen2D.log"
    shell: "mols_gen2D {input} {output} 2>{log}"

rule DEDUPL:
    priority: 6
    input: "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz"
    output: "{WD}/05_dedupl/data/{prefix}_{cid}_dedupl.csv.gz"
    log: "{WD}/05_dedupl/log/{prefix}_{cid}_dedupl.log"
    shell: "mols_dedupl {input} {output} -r {WD}/05_dedupl/{prefix}_ref.hdf --log DEBUG 2>{log}"

rule STD:
    priority: 7
    input: WD + "/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    output:
        passed = "{WD}/04_std/data/{prefix}_{cid}_passed.csv.gz",
        filtered = "{WD}/04_std/data/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/04_std/data/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/04_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} $(echo {output.passed} | rev | cut -d_ -f2- | rev).csv.gz >{log}  2>&1"  # {cid} unkown on this context

rule DGC:
    priority: 8
    input: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/03_deglyco/data/{prefix}_{cid}_deglyco.sdf.gz"
    log: "{WD}/03_deglyco/log/{prefix}_{cid}_deglyco_knwf.log"  # log not from the job but from the workflow execution, useful for checking if everything went fine
    shell: "mols_deglyco -s {input} -i {molid} -o {WD}/03_deglyco -w {file_knwf} >{log} 2>&1"   # 1>{log} 2>&1" #/dev/null"

rule LOAD:
    priority: 9
    input: "{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/02_load/data/{prefix}_{cid}.sdf.gz"
    log: "{WD}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 10
    input: input_file
    output: expand(WD + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1"
