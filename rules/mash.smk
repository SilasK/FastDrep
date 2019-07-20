

rule mash_sketch_genome:
    input:
        filter_genome_folder,
        genome_folder
    output:
        "genomes.msh"
    params:
        k= config['mash']['k'],
        s= config['mash']['sketch_size'],
        out_name = lambda wc, output: os.path.splitext(output[0])[0]
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/sketch.log"
    shell:
        "mash sketch -o {params.out_name} -p {threads} -s {params.s} -k {params.k} {input[0]}/* 2> {log}"


rule mash_calculate_dist:
    input:
        genomes=rules.mash_sketch_genome.output
    output:
        "mash_dists.txt"
    params:
        d= config['mash']['max_dist']
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/dist.log"
    shell:
        "mash dist -p {threads} -d {params.d} {input.genomes} {input.genomes} > {output[0]} 2> {log}"




# rule group_species_mash:
#     input:
#         mash_dists=
#         quality=
#     output:
#
#     params:
#         threshold = 0.95,
#         fillna=0.8,
#         linkage_method='average',
#         square=False
#     run:
#         M= gd.load_mash(input.mash_dists)
#         labels= gd.group_species_linkage(M,**params)
#
#         Mapping_species= best_genome_from_table(labels,Q.QualityScore)

localrules: filter_mash
checkpoint filter_mash:
    input:
        rules.mash_calculate_dist.output[0]
    output:
        temp(directory('minimap/alignment_lists'))
    params:
        treshold=config['mash']['dist_treshold'],
        N= 1000
    run:

        F= gd.load_mash(input[0])
        G= gd.to_graph(F.query(f"Distance<={params.treshold}"))
        G.remove_edges_from(G.selfloop_edges())

        os.makedirs(output[0])


        fout=None
        for i,e in enumerate(G.edges()):
            if (i % params.N) ==0:
                n= int(i // params.N )
                if fout is not None: fout.close()
                fout= open(f"{output[0]}/subset_{n}.txt",'w')
            else:
                fout.write("\t".join(sorted(e))+'\n')
