LOCAL = "tests/tmp/local_workflow/"

rule all:
    input: LOCAL + "data/chembl_small_001_sub.csv.gz"


rule substructure:
    input:
        LOCAL + "data/chembl_small_001_passed.csv.gz",  # molecules
        "tests/data/crms_passed.sdf.gz"  # fragments
    output: LOCAL + "data/chembl_small_001_sub.csv.gz"
    log: LOCAL + "log/chembl_small_001_sub.log"
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


rule mkdir:
    output:
        LOCAL + "data",
        LOCAL + "log"
    shell: "rm -rf {output[0]} {output[1]}; mkdir {output[0]} {output[1]}"
