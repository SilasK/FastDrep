if config.get('skip_filter',False):
    filter_genome_folder= genome_folder
else:
    filter_genome_folder='filtered_bins'

    rule filter_genomes:
        input:
            dir=os.path.abspath(genome_folder),
            quality=config['genome_qualities']
        output:
            directory(filter_genome_folder)
        params:
            filter=config['filter_criteria']
        run:
            import pandas as pd

            Q= pd.read_csv(input.quality,sep='\t',index_col=0)
            assert not Q.index.duplicated().any()
            filtered= Q.query(params.filter).index



            os.makedirs(output[0])
            for f in filtered:
                os.symlink(os.path.join(input.dir,f),
                           os.path.join(output[0],f))

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



rule fastANI:
    input:
        genome_list= "clusters/mash/Cluster{i}.txt",
        genomes= filter_genome_folder,
        all_Genomes= genome_folder
    output:
        "clusters/fastani/ANIcluster{i}.txt",
    resources:
        mem=30
    threads:
        config['threads']
    conda:
        "envs/fastANI.yaml"
    log:
        "logs/fastANI/cluster{i}.log"
    benchmark:
        "logs/benchmarks/fastANI/cluster{i}.log"
    shell:
        " fastANI "
        " --threads {threads} "
        " --queryList {input.genome_list} --refList {input.genome_list} "
        " -k {config[fastani][k]} "
        " --fragLen {config[fastani][fragLen]} "
        " --minFrag {config[fastani][minFrag]} "
        " -o {output} "
        " 2> {log}"


def get_allFastANI(wildcards):
    ALL_I=  glob_wildcards(os.path.join(checkpoints.mash_cluster.get(**wildcards).output.clusters,
                                        "Cluster{i}.txt")).i
    return expand("clusters/fastani/ANIcluster{i}.txt", i=ALL_I )


rule combine_fastANI:
    input:
        get_allFastANI
    output:
        "ANI.tsv"
    shell:
        "cat {input} > {output};"
        "rm -r clusters "
