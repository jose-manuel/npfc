rule all:
    input:
        "log/inputs_to_cluster"

rule chunk:
    input:
        mols = '0_raw/DNP272_ct_no8.sdf.gz'
    output:
        mols = dynamic("data/{chunk}.csv.gz")
    log:
        "log/chunk_dnp.log"
    shell:
        "load_mols {input} data/dnp.csv.gz -n 10000 --in_id UKEY 2>{log}"

rule standardize:
    input:
        dynamic("data/{chunk}.csv.gz")
    output:
        "log/inputs_to_cluster"
    shell:
        "ssh gwdg1 mkdir -p /scratch/josemanuel.gally/dnp/data; scp -r 1_input gwdg1:/scratch/josemanuel.gally/dnp/data; echo 'True' > {output}"
