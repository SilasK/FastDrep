





rule bbsketch_mags:
    input:
        genome_folder
    wildcard_constraints:
        NTorAA="(aa|nt)",
    output:
        out="bbsketch/mags_{NTorAA}.sketch.gz"
    params:
        k= lambda wildcards: config['bbsketch'][wildcards.NTorAA]['k'],
        translate=lambda wildcards: wildcards.NTorAA=='aa',
        overwrite=True,
        command=f"bbsketch.sh perfile {genome_folder}/*.fasta",
    resources:
        time= 10
    log:
        "logs/bbsketch/sketch_mags_{NTorAA}.log"
    benchmark:
        "logs/benchmark/bbsketch/sketch_mags_{NTorAA}.txt"
    conda:
        "../envs/bbmap.yaml"
    threads:
        16
    script:
        "../scripts/runBB.py"




rule allvall:
    input:
        ref="bbsketch/mags_{NTorAA}.sketch.gz"
    output:
        out="tables/bbsketch_{NTorAA}.tsv"
    wildcard_constraints:
        NTorAA="(nt|aa)"
    params:
        amino=lambda wildcards: wildcards.NTorAA=='aa',
        overwrite=True,
        command="comparesketch.sh alltoall",
        prealloc=0.75,
        format=3,
        k=lambda wildcards: config['bbsketch'][wildcards.NTorAA]['k'],
    shadow:
        "minimal"
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= config['mem']['large']
    benchmark:
        "logs/benchmark/bbsketch/alltoall_{NTorAA}.txt"
    log:
        "logs/bbsketch/alltoall_{NTorAA}.log"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"


rule sendsketch:
    input:
        "bbsketch/mags_aa.sketch.gz"
    output:
        "tables/mapping2refseq_aa.sketch.gz"
    params:
        minid=0.9
    log:
        "logs/bbsketch/sendsketch.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        1
    resources:
        mem= config['mem']['large'],
    benchmark:
        "logs/benchmark/bbsketch/sendsketch.txt"
    shell:
        "sendsketch.sh in={input} out={output} protein format=3 minid={params.minid} 2> {log}"
