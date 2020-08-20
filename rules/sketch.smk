localrules: list_genomes
rule list_genomes:
    input:
        genome_folder
    output:
        temp("sketches/genome_list_{sketcher}.txt")
    run:
        from glob import glob

        fasta_files=glob(f"{input[0]}/*{config['fasta_extension']}")
        if len(fasta_files)<=2:
            Exception(f"You have less than 2 fasta file with extension ({config['fasta_extension']})"
                      f" in the folder {genome_folder}"
                      )

        with open(output[0],'w') as f:
            f.write('\n'.join(fasta_files)+'\n' )


rule mash_sketch_genome:
    input:
        "sketches/genome_list_mash.txt"
    output:
        "sketches/genomes.msh"
    params:
        k= config['sketch_k'],
        s= config['sketch_size'],
        out_name = lambda wc, output: os.path.splitext(output[0])[0]
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/sketch.log"
    benchmark:
        "logs/benchmark/mash_sketch.txt"
    shell:
        "mash sketch -o {params.out_name} -p {threads} -s {params.s} -k {params.k} -l {input[0]} 2> {log}"


rule mash_calculate_dist:
    input:
        genomes=rules.mash_sketch_genome.output
    output:
        "tables/mash_dists.tsv"
    params:
        d= config['sketch_max_dist'],
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/dist.log"
    benchmark:
        "logs/benchmark/mash_dists.tsv"
    shell:
        "mash dist -p {threads} -d {params.d} "
        "{input.genomes} {input.genomes} > {output[0]} 2> {log}"


rule bindash_sketch_genome:
    input:
        "sketches/genome_list_bindash.txt"
    output:
        "sketches/genomes.bdsh",
        "sketches/genomes.bdsh.dat",
        "sketches/genomes.bdsh.txt"
    params:
        k= config['sketch_k'],
        sketchsize64= int(config['sketch_size'])//64,
        extra=config.get('bindash_extra',"")
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/sketch.log"
    benchmark:
        "logs/benchmark/bindash_sketch.txt"
    shell:
        "bindash sketch "
        "--outfname={output[0]} "
        "--nthreads={threads} "
        "--sketchsize64={params.sketchsize64} "
        "--kmerlen={params.k} "
        "{params.extra} "
        "--listfname={input[0]} 2> {log}"

rule bindash_dist:
    input:
        rules.bindash_sketch_genome.output[0]
    output:
        "tables/bindash_dists.tsv"
    params:
        d= config['sketch_max_dist']
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/dist.log"
    benchmark:
        "logs/benchmark/bindash_dist.txt"
    shell:
        "bindash dist "
        "--nthreads={threads} "
        "--mthres={params.d} "
        "--outfname={output} {input[0]} 2> {log}"



rule precluster:
    input:
        dists=f"tables/{config['sketcher']}_dists.tsv",
        quality="tables/Genome_quality.tsv",
    output:
        preclustering="tables/preclustering.tsv",
        edgelist= temp('alignment/all_pairs_for_alignment.tsv')
    params:
        treshold=config['pre_cluster_treshold'],
        min_identity=config['pre_cluster_min_identity'],
        N=config['subset_size_alignments'],
    resources:
        mem=config['mem']['large']

    run:

        import pandas as pd

        Q= gd.load_quality(input.quality)
        Q['quality_score']= Q.eval(config['quality_score'])

        Msh= gd.load_mash(input.dists)
        Msh_high= Msh.query(f"Identity>={params.treshold}")
        G= gd.to_graph(Msh_high)

        #map genomes to clusters by single linkage
        Clustering= {}
        for i,cc in enumerate(gd.nx.connected_components(G)):
            for g in cc:
                  Clustering[g]=i
        Clustering =gd.best_genome_from_table(pd.Series(Clustering),Q.quality_score)

        # don't foget  Genomes qith no connection
        lonly_genomes= Q.index.difference(Clustering.index)
        Clustering.loc[lonly_genomes]=lonly_genomes

        #save clustering
        Clustering.name='Precluster_representative'
        Clustering.index.name='Genome'
        Clustering.to_csv(output.preclustering,sep='\t')

        #Subset_pairwise distances
        Cluster_representatives= Clustering.unique()

        Msh_sub= Msh.loc[ (
        Msh.index.levels[0].intersection(Cluster_representatives),
        Msh.index.levels[1].intersection(Cluster_representatives)
         ),:].query(f"Identity>={params.min_identity}")

        print(f"From {Clustering.shape[0]} genomes {Cluster_representatives.shape[0]} "
              f"({Cluster_representatives.shape[0]/Clustering.shape[0]*100:.2f}%) are representatives.\n"
              f"This decreases the number of interaction to {Msh_sub.shape[0]:d} ({Msh_sub.shape[0]/Msh.shape[0]:.2g})"
             )



        # save list of all comparison to perform alignment on them
        G= gd.to_graph(Msh_sub)
        if hasattr(G,'selfloop_edges'):
            G.remove_edges_from(G.selfloop_edges())


        gd.nx.write_edgelist(G,output.edgelist,delimiter='\t',data=False,comments=None)
