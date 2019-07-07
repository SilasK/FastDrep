




rule pyani:
    input:
        filter_genome_folder,
        genome_folder
    output:
        "pyani/{method}/{method}_percentage_identity.tab"
    params:
        outdir= "pyani/{method}"
    threads:
        config['threads']
    conda:
        "../envs/pyani.yaml"
    log:
        "logs/pyani/{method}.log"
    shell:
        "average_nucleotide_identity.py "
        "-i {input[0]} -o {params.outdir} "
        "-m {wildcards.method} -g -l {log}"
