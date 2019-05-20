

rule mash_sketch_genome:
    input:
        filter_genome_folder,
        genome_folder
    output:
        "filtered_genomes.msh"
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

localrules: mash_cluster, combine_fastANI
checkpoint mash_cluster:
    input:
        rules.mash_calculate_dist.output
    output:
        clusters=directory("clusters/mash"),
        singeltons= "singletons.txt"
    params:
        treshold=config['mash']['cluster_dist']
    run:
        from common import genome_pdist as gd
        import shutil, os

        M=gd.load_mash(input[0],False)
        G= gd.to_graph(M.query(f"Distance<={params.treshold}"))
        connected_components= gd.get_connected_components(G)

        if os.path.exists(output[0]): shutil.rmtree(output[0])
        os.makedirs(output[0])

        Sinlge_clusters=set()

        i =0
        for cc in connected_components:
            if len(cc)>1:
                i+=1
                with open(os.path.join(output[0],f"Cluster{i}.txt"),'w') as f:
                    f.write("\n".join(cc)+'\n')
            else:
                Sinlge_clusters.add(cc.pop())

        with open(os.path.join(output.singeltons),'w') as f:
            f.write("\n".join(Sinlge_clusters)+'\n')
