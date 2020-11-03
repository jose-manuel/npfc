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
try:
    prefix = config['prefix']
except KeyError:
    prefix = ''
input_file = config['input_file']
# additional
frags_file = config['frags_file']  # fragment file to use for substructure search
frags_subdir = config['frags_subdir']
natref_dedupl_reffile = config['natref_dedupl_reffile']
natref_fcg_dir = config['natref_fcg_dir']
chunksize = config['chunksize']  # maximum number of molecules per chunk
# specific to synthetic
natref_subdir = config['natref_subdir']  # WD for defining natural compounds, subdir with same frags is also searched for pnp annotation
# by default consider fcc for pnp attributes, if one wants to only compare pairs of fragments, provide empty string in config file instead
try:
    pnp_attributes = config['pnp_attributes']
except KeyError:
    pnp_attributes = 'fcc'
# preprocess subdir
try:
    prep_subdir = config['prep_subdir']
except KeyError:
    prep_subdir = 'prep'
# from master script
num_chunks = config['num_chunks']

# protocol for standardizing molecules
fallback_default_std_mols = False
try:
    config_std_mols = config['std_protocol']
    if config_std_mols == '' or config_std_mols == 'DEFAULT':
        fallback_default_std_mols = True
except KeyError:
    fallback_default_std_mols = True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# find where the KNIME workflow is (installed as data of the npfc package)
file_knwf = pkg_resources.resource_filename('npfc', 'data/mols_deglyco.knwf')

# remove trailing / from WD for consistency
if WD.endswith('/'):
    WD = WD[:-1]
WD += '/data'
# remove trailing / from natref_fcg_dir for consistency
if natref_fcg_dir.endswith('/'):
    natref_fcg_dir = natref_fcg_dir[:-1]


# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]

# fall back to default std configuration in case either missing from JSON file or empty string
if fallback_default_std_mols:
    config_std_mols = pkg_resources.resource_filename('npfc', 'data/std_mols.json')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(f"{WD}/{prep_subdir}/{natref_subdir}/{frags_subdir}/10_pnp/data/{prefix}" + '_{cid}_pnp.csv.gz', cid=chunk_ids)

rule PNP:
    priority: 0
    input: ancient("{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/09_fcg/data/{prefix}_{cid}_fcg.csv.gz")
    output:
        fgraphs = "{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/10_pnp/data/{prefix}_{cid}_pnp.csv.gz",
        list_pnps = "{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/10_pnp/log/{prefix}_{cid}_list_pnp.csv.gz"
    log: "{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/10_pnp/log/{prefix}_{cid}_pnp.log"
    shell: "fcg_annotate_pnp {input} {natref_fcg_dir} {output.fgraphs} -l {output.list_pnps} -d '" + pnp_attributes + "' >{log} 2>&1"

rule FCG:
    priority: 1
    input: ancient("{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz")
    output: "{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/09_fcg/data/{prefix}_{cid}_fcg.csv.gz"
    log: "{WD}" + f"/{prep_subdir}/{natref_subdir}/{frags_subdir}" + "/09_fcg/log/{prefix}_{cid}_fcg.log"
    shell: "fcg_generate {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 2
    input: ancient("{WD}/{prep_subdir}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fs/data/{prefix}_{cid}_fs.csv.gz")
    output: "{WD}/{prep_subdir}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}/{prep_subdir}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule FS:
    priority: 3
    input:
        mols = ancient("{WD}/{prep_subdir}" + f"/{natref_subdir}" + "/06_subset/data/{prefix}_{cid}_subset.csv.gz"),
        frags = ancient(frags_file)
    output: "{WD}/{prep_subdir}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fs/data/{prefix}_{cid}_fs.csv.gz"
    log: "{WD}/{prep_subdir}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fs/log/{prefix}_{cid}_fs.log"
    shell: "frags_search {input.mols} {input.frags} {output} >{log} 2>&1"

rule SUBSET:
    priority: 4
    input:
        mols = ancient("{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz"),
        ref = ancient(natref_dedupl_reffile)
    output: "{WD}/{prep_subdir}" + f"/{natref_subdir}" + "/06_subset/data/{prefix}_{cid}_subset.csv.gz"
    log: "{WD}/{prep_subdir}" + f"/{natref_subdir}" + "/06_subset/log/{prefix}_{cid}_subset.log"
    shell: "mols_subset {input.mols} {input.ref} {output} >{log} 2>&1"

rule DEPICT:
    priority: 5
    input: ancient("{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz")
    output: "{WD}/{prep_subdir}/05_depict/data/{prefix}_{cid}_depict.csv.gz"
    log: "{WD}/{prep_subdir}/05_depict/log/{prefix}_{cid}_depict.log"
    shell: "mols_depict {input} {output} 2>{log}"

rule DEDUPL:
    priority: 6
    input: ancient("{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz")
    output:
        passed = "{WD}/{prep_subdir}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz",
        filtered = "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_filtered.csv.gz",
        synonyms = "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_synonyms.csv.gz"
    log: "{WD}/{prep_subdir}/04_dedupl/log/{prefix}_{cid}_dedupl.log"
    shell: "mols_dedupl {input} {output.passed} -d {output.filtered} -s {output.synonyms} -r {WD}/{prep_subdir}/04_dedupl/{prefix}_ref.hdf 2>{log}"

rule STD:
    priority: 7
    input: ancient(WD + "/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz")
    output:
        std = "{WD}/{prep_subdir}/03_std/data/{prefix}_{cid}_std.csv.gz",
        filtered = "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/{prep_subdir}/03_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} -p " + config_std_mols + " 2>{log}"  # mols_standardize takes a dir as output

rule LOAD:
    priority: 9
    input: ancient("{WD}/{prep_subdir}/01_chunk/data/{prefix}_{cid}.sdf.gz")
    output: "{WD}/{prep_subdir}/02_load/data/{prefix}_{cid}.csv.gz"
    log: "{WD}/{prep_subdir}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 10
    input: ancient(input_file)
    output: expand(WD + '/' + prep_subdir + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + '/' + prep_subdir + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -p '{prefix}' -o {WD}/{prep_subdir}/01_chunk/data/ >{log} 2>&1"
