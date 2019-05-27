LOCAL = "tests/tmp/local_workflow/"

rule all:
    input:
        LOCAL + "data/chembl_small_001_deglyco.sdf.gz"




# rule standardize:
#     input: LOCAL + "data/chembl_small_001_deglyco.sdf.gz"
#     output:
#         LOCAL + "data/chembl_small_001_passed.csv.gz",
#         LOCAL + "data/chembl_small_ref.hdf"
#     log: LOCAL + "log/chembl_small_001_std.log"
#     shell: "standardize_mols {input} -o $(echo {input} | cut -d. -f1).csv.gz -r {output[1]} 2>{log}"


rule deglyco:
    input: LOCAL + "data/chembl_small_001.sdf.gz"
    output: LOCAL + "data/chembl_small_001_deglyco.sdf.gz"
    log: LOCAL + "log/chembl_small_deglyco.log"
    shell: "run_deglyco -s ${input} -i chembl_id -o $(echo ${output} | rev | cut -d/ -f3- | rev)"


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
    shell:
        "mkdir {output[0]} {output[1]}"
