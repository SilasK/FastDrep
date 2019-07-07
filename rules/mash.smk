

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
        genomes=rules.mash_sketch_genome.output,
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



localrules: filter_mash
checkpoint filter_mash:
    input:
        rules.mash_calculate_dist.output[0]
    output:
        temp("alignment_list.txt")
    params:
        treshold=config['mash']['dist_treshold']
    run:

        F= gd.load_mash(input[0])
        G= gd.to_graph(F.query(f"Distance<={params.treshold}"))
        G.remove_edges_from(G.selfloop_edges())

        with open(output[0],'w') as fout:
            for e in G.edges():
                fout.write("\t".join(sorted(e))+'\n')
