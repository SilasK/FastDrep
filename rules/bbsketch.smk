

checkpoint build_sketch:
    input:
        genome_folder= genome_folder,
    output:
        sketch_folder= directory("bbsketch/sketches_{NTorAA}"),
    wildcard_constraints:
        NTorAA="(nt|aa)"
    threads:
        config['threads']
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= 50
    log:
        "logs/bbsketch/workflow_{NTorAA}.log"
    params:
        path= os.path.dirname(workflow.snakefile),
        amino=lambda wildcards: wildcards.NTorAA=='aa',
        k=lambda wildcards: config['bbsketch'][wildcards.NTorAA]['k']
    shell:
        "snakemake -s {params.path}/rules/buildsketch.smk "
        "--config  "
        " k='{params.k}' "
        " genome_folder='{input.genome_folder}' "
        " sketch_folder='{output.sketch_folder}' "
        " amino={params.amino} "
        " --rerun-incomplete "
        "-j {threads} --nolock 2> {log}"


def get_mags(wildcards):

    folder=checkpoints.rename_genomes.get().output.genome_folder

    genomes= glob_wildcards(os.path.join(folder,'{genome}.fasta')).genome
    return list(genomes)


def mergesketch_input(wildcards):

    checkpoints.build_sketch.get(NTorAA=wildcards.NTorAA)

    if wildcards.resolution_level=='mags':
        return  expand("bbsketch/sketches_{NTorAA}/{genome}.sketch.gz",
                       genome=get_mags(wildcards),**wildcards)
    else:
        return expand("bbsketch/sketches_{NTorAA}/{genome}.sketch.gz",
               genome=get_representatives(wildcards),**wildcards)





localrules: mergesketch
rule mergesketch:
    input:
        lambda wildcards: mergesketch_input(wildcards)
    wildcard_constraints:
        NTorAA="(aa|nt)",
        resolution_level="(species|strains|mags)"
    output:
        out="bbsketch/{resolution_level}_{NTorAA}.sketch.gz"
    threads:
        1
    run:
        io.cat_files(input,output[0],gzip=False)







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
