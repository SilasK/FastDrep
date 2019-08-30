
rule bbsketch:
    input:
        input=os.path.join(genome_folder,"{genome}.fasta")
    output:
        out="bbsketch/sketches_{NTorAA}/{genome}.sketch.gz"
    params:
        k= lambda wildcards: config['bbsketch'][wildcards.NTorAA]['k'],
        translate=lambda wildcards: wildcards.NTorAA=='aa',
        overwrite=True,
        command="bbsketch.sh",
        name0="{genome}"
    resources:
        mem= 1,
        time= 5
    log:
        "logs/bbsketch/sketch_{NTorAA}/{genome}.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        1
    script:
        "../scripts/runBB.py"


localrules: mergesketch
rule mergesketch:
    input:
        lambda wildcards: expand("bbsketch/sketches_{{NTorAA}}/{genome}.sketch.gz",
               genome=get_representatives(wildcards))
    wildcard_constraints:
        NTorAA="(aa|nt)",
        resolution_level="(species|strains)"
    output:
        out="bbsketch/{resolution_level}_{NTorAA}.sketch.gz"
    threads:
        1
    run:
        io.cat_files(input,output[0],gzip=False)


def get_mags(wildcards):

    folder=checkpoints.rename_genomes.get().output.genome_folder

    genomes= glob_wildcards(os.path.join(folder,'{genome}.fasta')).genome
    return lsit(genomes)


rule mergesketch_mags:
    input:
        lambda wildcards: expand("bbsketch/sketches_{{NTorAA}}/{genome}.sketch.gz",
               genome=get_mags(wildcards))
    wildcard_constraints:
        NTorAA="(aa|nt)",
    output:
        out="bbsketch/mags_{NTorAA}.sketch.gz"
    threads:
        1
    run:
        io.cat_files(input,output[0],gzip=True)



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
        format=3
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= 50
    log:
        "logs/bbsketch/alltoall_{NTorAA}.log"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"


rule sendsketch:
    input:
        "bbsketch/{resolution_level}_aa.sketch.gz"
    output:
        "tables/refseq_mapping_{resolution_level}.tsv"
    params:
        minid=0.9
    log:
        "logs/bbsketch/sendsketch_{resolution_level}.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        1
    resources:
        mem= 1,
    shell:
        "sendsketch.sh in={input} out={output} protein format=3 minid={params.minid} 2> {log}"
