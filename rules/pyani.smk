



rule pyani:
    input:
        genome_folder
    output:
        directory("pyani/{method}") #{method}_percentage_identity.tab"
    threads:
        config['threads']
    conda:
        "../envs/pyani.yaml"
    log:
        "logs/pyani/{method}.log"
    threads:
        12
    shell:
        "average_nucleotide_identity.py --noclobber "
        "--workers {threads} "
        "-i {input[0]} -o {output[0]} "
        "-m {wildcards.method} -g -l {log}"
