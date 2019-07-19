




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
        3
    params:
        preset= "-c --secondary=no" #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    shell:
        "minimap2 {params.preset}  -t {threads} {input.querry} {input.ref}   > {output}"



rule many_minimap:
    input:
        genome_folder=filter_genome_folder,
        alignment_list='minimap/alignment_lists/subset_{n}.txt',
        genome_stats= "genome_stats.tsv"
    output:
        alignments_stats="minimap/alignments_stats/subset_{n}.tsv",
    log:
        "logs/minimap2/subset_{n}.txt"
    threads:
        3
    conda:
        "../envs/minimap2.yaml"
    params:
        minimap_extra= "-c --secondary=no", #"-x asm10", #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
        paf_folder="minimap/paf",
        extension='.fasta'
    script:
        "../scripts/many_minimap.py"




def combine_paf_input(wildcards):

    alignment_list_folder = checkpoints.filter_mash.get().output[0]
    N= glob_wildcards(os.path.join(alignment_list_folder,"subset_{n}.txt")).n

    return expand("minimap/alignments_stats/subset_{n}.tsv",n=N)


localrules: combine_paf
rule combine_paf:
    input:
        combine_paf_input
    output:
        "alignments_stats.tsv"
    run:
        shell("head -n1 {input[0]} > {output}")
        for fin in input:
            shell("tail -n+2 {fin} >> {output}")

        in_dir= os.path.dirname(input[0])
        import shutil
        shutil.rmtree(in_dir)
