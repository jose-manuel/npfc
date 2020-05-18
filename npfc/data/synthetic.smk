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
frags_subdir = "frags_" + config['frags_subdir']
chunksize = config['chunksize']  # maximum number of molecules per chunk
# specific to synthetic
natref = config['natref']  # WD for defining natural compounds, subdir with same frags is also searched for pnp annotation
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
natref_subdir = "natref_" + Path(natref).stem.split('.')[0]  # subfolder where outputs using nat_ref are stored
# natref dedupl - for subset
natref_dedupl_dir = list(Path(f"{natref}/data").glob("[0-9][0-9]_dedupl"))[0]
natref_dedupl_reffile = list(natref_dedupl_dir.glob('*.hdf'))[0]
# natref fg - for pnp
natref_fgraph_dir = [str(x) for x in Path(f"{natref}/data/{frags_subdir}").glob("[0-9][0-9]_fgraph")][0] + '/data'
# define chunk_ids for wildcard expansion
chunk_ids = [str(i+1).zfill(3) for i in range(num_chunks)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rule all:
    input: expand(f"{WD}/{natref_subdir}/{frags_subdir}/10_pnp/data/{prefix}" + '_{cid}_pnp.csv.gz', cid=chunk_ids)

rule PNP:
    priority: 0
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fgraph/data/{prefix}_{cid}_fgraph.csv.gz"
    output:
        fgraphs = "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_pnp/data/{prefix}_{cid}_pnp.csv.gz",
        list_pnps = "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_pnp/log/{prefix}_{cid}_list_pnp.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/10_pnp/log/{prefix}_{cid}_pnp.log"
    shell: "fgraph_annotate_pnp {input} {natref_fgraph_dir} {output.fgraphs} -l {output.list_pnps} >{log} 2>&1"

rule FGRAPH:
    priority: 1
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fgraph/data/{prefix}_{cid}_fgraph.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/09_fgraph/log/{prefix}_{cid}_fgraph.log"
    shell: "fgraph_generate {input} {output} --min-frags 2 --max-frags 9999 --max-overlaps 5 >{log} 2>&1"

rule FCC:
    priority: 2
    input: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fcc/data/{prefix}_{cid}_fcc.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/08_fcc/log/{prefix}_{cid}_fcc.log"
    shell: "fc_classify {input} {output} -c 3 >{log} 2>&1"

rule FSEARCH:
    priority: 3
    input:
        mols = "{WD}" + f"/{natref_subdir}" + "/06_subset/data/{prefix}_{cid}_subset.csv.gz",
        frags = frags_file
    output: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fsearch/data/{prefix}_{cid}_fsearch.csv.gz"
    log: "{WD}" + f"/{natref_subdir}/{frags_subdir}" + "/07_fsearch/log/{prefix}_{cid}_fsearch.log"
    shell: "mols_fsearch {input.mols} {input.frags} {output} >{log} 2>&1"

rule SUBSET:
    priority: 4
    input:
        mols = "{WD}/05_depict/data/{prefix}_{cid}_depict.csv.gz",
        ref = natref + "/data/04_dedupl/dnp_ref.hdf"
    output: "{WD}" + f"/{natref_subdir}" + "/06_subset/data/{prefix}_{cid}_subset.csv.gz"
    log: "{WD}" + f"/{natref_subdir}" + "/06_subset/log/{prefix}_{cid}_subset.log"
    shell: "mols_subset {input.mols} {input.ref} {output} >{log} 2>&1"

rule DEPICT:
    priority: 5
    input: "{WD}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz"
    output: "{WD}/05_depict/data/{prefix}_{cid}_depict.csv.gz"
    log: "{WD}/05_depict/log/{prefix}_{cid}_depict.log"
    shell: "mols_depict {input} {output} 2>{log}"

rule DEDUPL:
    priority: 6
    input: "{WD}/03_std/data/{prefix}_{cid}_std.csv.gz"
    output:
        passed = "{WD}/04_dedupl/data/{prefix}_{cid}_dedupl.csv.gz",
        filtered = "{WD}/04_dedupl/log/{prefix}_{cid}_filtered.csv.gz"
    log: "{WD}/04_dedupl/log/{prefix}_{cid}_dedupl.log"
    shell: "mols_dedupl {input} {output.passed} -d {output.filtered} -r {WD}/04_dedupl/{prefix}_ref.hdf 2>{log}"

rule STD:
    priority: 7
    input: WD + "/02_load/data/{prefix}_{cid}.csv.gz"
    output:
        std = "{WD}/03_std/data/{prefix}_{cid}_std.csv.gz",
        filtered = "{WD}/03_std/log/{prefix}_{cid}_filtered.csv.gz",
        error = "{WD}/03_std/log/{prefix}_{cid}_error.csv.gz"
    log: "{WD}/03_std/log/{prefix}_{cid}_std.log"
    shell: "mols_standardize {input} {output.std} -f {output.filtered} -e {output.error} 2>{log}"  # mols_standardize takes a dir as output

rule LOAD:
    priority: 9
    input: "{WD}/01_chunk/data/{prefix}_{cid}.sdf.gz"
    output: "{WD}/02_load/data/{prefix}_{cid}.csv.gz"
    log: "{WD}/02_load/log/{prefix}_{cid}_load.log"
    shell: "mols_load {input} {output} --in_id {molid} >{log} 2>&1"

rule CHUNK:
    priority: 10
    input: input_file
    output: expand(WD + '/01_chunk/data/' + prefix + '_{cid}.sdf.gz', cid=chunk_ids)
    log: WD + "/01_chunk/log/" + prefix + "_chunk.log"
    shell: "chunk_sdf -i {input} -n {chunksize} -o {WD}/01_chunk/data/ >{log} 2>&1"
