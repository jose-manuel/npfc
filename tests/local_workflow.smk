LOCAL = "tests/tmp/local_workflow/"

rule all:
    input:
        LOCAL + "data/chembl_small_001_map_crms.csv.gz",
        LOCAL + "local_workflow.png"


rule map:
    input: LOCAL + "data/chembl_small_001_fcc_crms.csv.gz"
    output: LOCAL + "data/chembl_small_001_map_crms.csv.gz"
    log: LOCAL + "log/chembl_small_001_map_crms.log"
    shell: "map_frags {input} {output} 2>{log}"


rule classify:
    input: LOCAL + "data/chembl_small_001_sub_crms.csv.gz"
    output: LOCAL + "data/chembl_small_001_fcc_crms.csv.gz"
    log: LOCAL + "log/chembl_small_001_fcc_crms.log"
    shell: "classify_frags {input} {output} 2>{log}"


rule substructure:
    input:
        LOCAL + "data/chembl_small_001_passed.csv.gz",  # molecules
        "tests/data/crms_passed.csv.gz"  # fragments
    output: LOCAL + "data/chembl_small_001_sub_crms.csv.gz"
    log: LOCAL + "log/chembl_small_001_sub_crms.log"
    shell: "substruct_mols {input[0]} {input[1]} {output} 2>{log}"


rule standardize:
    input: LOCAL + "data/chembl_small_001_deglyco.sdf.gz"
    output:
        LOCAL + "data/chembl_small_001_passed.csv.gz",
        LOCAL + "data/chembl_small_ref.hdf"
    log: LOCAL + "log/chembl_small_001_std.log"
    shell: "standardize_mols {input} -o $(echo {input} | rev | cut -d_ -f2- | rev).csv.gz -r {output[1]} 2>{log}"


rule deglyco:
    input: LOCAL + "data/chembl_small_001.sdf.gz"
    output: LOCAL + "data/chembl_small_001_deglyco.sdf.gz"
    log: LOCAL + "log/chembl_small_deglyco.log"
    shell: "deglyco_mols -s {input} -i chembl_id -o {LOCAL} -w /home/gally/Projects/NPFC/src/bin/deglyco_mols.knwf >{log} 2>/dev/null"


rule chunk:
    input:
        "tests/data/chembl_small.sdf.gz",
        LOCAL + "data",
        LOCAL + "log"
    output: LOCAL + "data/chembl_small_001.sdf.gz"
    log: LOCAL + "log/chembl_small.log"
    shell: "load_mols {input[0]} $(echo {output} | rev | cut -d_ -f2- | rev).sdf.gz -n 500 --in_id chembl_id 2>{log}"


rule tasktree:
    input: LOCAL
    output: LOCAL + "local_workflow.png"
    shell: "snakemake -s tests/local_workflow.smk --dag | dot -Tpng > {output}"


rule mkdir:
    output:
        LOCAL + "data",
        LOCAL + "log"
    shell: "rm -rf {output[0]} {output[1]}; mkdir {output[0]} {output[1]}"
