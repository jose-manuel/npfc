rule all:
    input:
        "tests/data/chembl_small_001_passed.csv.gz"

rule chunk:
    input: "tests/data/chembl_small.sdf.gz"
    output: "tests/data/chembl_small_001.sdf.gz"
    log: "tests/data/chembl_small.log"
    shell: "load_mols {input} {input} -n 500 --in_id chembl_id 2>{log}"


rule deglyco:
    input: "tests/data/chembl_small_001.sdf.gz"
    output: "tests/data/chembl_small_001_deglyco.sdf.gz"
    shell: "cp {input} {output}"  # not motivated right now to make the KNIME worklow executable in CLI


rule standardize:
    input: "tests/data/chembl_small_001_deglyco.sdf.gz"
    output:
        "tests/data/chembl_small_001_passed.csv.gz",
        "tests/data/chembl_small_ref.hdf"
    log: "tests/data/chembl_small_001_std.log"
    shell: "standardize_mols {input} -o $(echo {input} | cut -d. -f1).csv.gz -r {output[1]} 2>{log}"
