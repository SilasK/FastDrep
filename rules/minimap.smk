

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

localrules: minimap_subsets
checkpoint minimap_subsets:
    input:
        genome_sketch
    output:
        temp(directory('minimap/subsets'))
    params:
        treshold=config['mummer']['max_dist'],
        N=config['mummer']['subset_size']
    run:

        F= gd.load_mash(input[0])
        G= gd.to_graph(F.query(f"Distance<={params.treshold}"))
        if hasattr(G,'selfloop_edges'):
            G.remove_edges_from(G.selfloop_edges())

        os.makedirs(output[0])

        fout=None
        for i,e in enumerate(G.edges()):
            if (i % params.N) ==0:
                n= int(i // params.N )+1
                if fout is not None: fout.close()
                fout= open(f"{output[0]}/subset_{n}.txt",'w')

            fout.write("\t".join(sorted(e))+'\n')

def get_minimap_subsets(wildcards):
    subset_dir= checkpoints.minimap_subsets.get().output[0]

    subsets= glob_wildcards(os.path.join(subset_dir,'{subset}.txt')).subset

    return subsets


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
        minimap_extra= "-c --secondary=no -x asm20", #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
        paf_folder="minimap/paf",
        extension='.fasta'
    script:
        "../scripts/many_minimap.py"




def combine_paf_input(wildcards):

    subsets=get_minimap_subsets(wildcards)


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
