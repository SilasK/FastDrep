


average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIm_output -m ANIm -g


rule pyani:
    input:
        filter_genome_folder,
        genome_folder
    output:
        directory("pyani/{method}")
    threads:
        config['threads']
    conda:
        "../envs/pyani.yaml"
    log:
        "logs/pyani/{method}.log"
    shell:
        "average_nucleotide_identity.py "
        "-i {input[0]} -o {output[0]} "
        "-m {wildcards.method} -g -l {log}"
