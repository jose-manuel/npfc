rule all:
    input:
        "tests/data/chembl_small_001.sdf.gz",
        "tests/data/chembl_small_002.sdf.gz"

rule chunk:
    input:
        "tests/data/chembl_small.sdf.gz"
    output:
        "tests/data/chembl_small_001.sdf.gz",
        "tests/data/chembl_small_002.sdf.gz"
    log:
        "tests/data/chembl_small.log"
    shell:
        "load_mols {input} {input} -n 250 --in_id chembl_id 2>{log}"
