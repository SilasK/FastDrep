

import os

genome_folder= config['genome_folder']
genomes= glob_wildcards(os.path.join(genome_folder,'{genome}.fasta')).genome

sketch_folder= config['sketch_folder']


rule all:
    input:
        expand(os.path.join(sketch_folder,"{genome}.sketch.gz"),genome=genomes)



rule bbsketch:
    input:
        input=os.path.join(genome_folder,"{genome}.fasta")
    output:
        out=os.path.join(sketch_folder,"{genome}.sketch.gz")
    params:
        k= ",".join([str(k) for k in config['k'] ]),
        translate=config['amino'],
        overwrite=True,
        command="bbsketch.sh",
        name0="{genome}"
    resources:
        mem= 1,
    group:
        "bbsketch"
    log:
        f"logs/{sketch_folder}/{{genome}}.log"
    # conda:
    #     "../envs/bbmap.yaml"
    threads:
        1
    script:
        "../scripts/runBB.py"
