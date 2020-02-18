


# rule mergesketch:
#     input:
#         lambda wildcards: expand("bbsketch/sketches_{{NTorAA}}/{genome}.sketch.gz",
#                genome=get_representatives(wildcards))
#     wildcard_constraints:
#         NTorAA="(aa|nt)",
#         resolution_level="(species|strains)"
#     group:
#         "bbsketch"
#     output:
#         out="bbsketch/{resolution_level}_{NTorAA}.sketch.gz"
#     threads:
#         1
#     run:
#         io.cat_files(input,output[0],gzip=False)
#



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
        mem= 10,
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
        format=3,
        k=lambda wildcards: config['bbsketch'][wildcards.NTorAA]['k'],
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= 50
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
        mem= 50,
    benchmark:
        "logs/benchmark/bbsketch/sendsketch.txt"
    shell:
        "sendsketch.sh in={input} out={output} protein format=3 minid={params.minid} usetaxidname=t 2> {log}"
