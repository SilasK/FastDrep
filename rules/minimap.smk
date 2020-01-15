

# rule minimap_index:
#     input:
#         fasta= f"{filter_genome_folder}/{genome}.fasta",
#     output:
#         "minimap/index/{genome}.mmi"
#     conda:
#         "../envs/minimap2.yaml"
#     group:
#         "minimap2"
#     conda:
#         "../envs/minimap2.yaml"
#     threads:
#         config['threads']
#     params:
#         extra= "-x asm10" #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
#     shell:
#         "minimap2 {params.preset}  -t {threads} -d {output} {input.fasta}   2> {log}"
#


rule minimap:
    input:
        ref= f"{genome_folder}/{{genome2}}.fasta",
        querry=f"{genome_folder}/{{genome1}}.fasta",
    output:
        "minimap/paf/{genome1}-{genome2}.paf"
    conda:
        "../envs/minimap2.yaml"
    group:
        "minimap2"
    conda:
        "../envs/minimap2.yaml"
    threads:
        config['threads']
    params:
        preset= "-c --secondary=no -x asm10" #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    shell:
        "minimap2 {params.preset}  -t {threads} {input.querry} {input.ref}   > {output}"



rule many_minimap:
    input:
        genome_folder=genome_folder,
        alignment_list="minimap/subsets/{subset}.txt",
        genome_stats= "tables/genome_stats.tsv"
    output:
        alignments_stats="minimap/alignments_stats/{subset}.tsv",
    log:
        "logs/minimap2/{subset}.txt"
    threads:
        config['threads']
    conda:
        "../envs/minimap2.yaml"
    params:
        minimap_extra= "-c --secondary=no", #"-x asm10", #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
        paf_folder="minimap/paf",
        extension='.fasta'
    script:
        "../scripts/many_minimap.py"




def combine_paf_input(wildcards):

    subsets=get_mummer_subsets(wildcards)


    return expand("minimap/alignments_stats/{subset}.tsv",subset=subsets)


localrules: combine_paf
rule combine_paf:
    input:
        combine_paf_input
    output:
        "tables/minimap_dists.tsv"
    run:
        shell("head -n1 {input[0]} > {output}")
        for fin in input:
            shell("tail -n+2 {fin} >> {output}")

        in_dir= os.path.dirname(input[0])
        import shutil
        shutil.rmtree(in_dir)
