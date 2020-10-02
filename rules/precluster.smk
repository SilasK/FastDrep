
def precluster_sketch_max_dist(wildcards):

    if wildcards.genomeset=="all_genomes":

        treshold= config['pre_cluster_treshold']

        assert treshold >= 0.95, "'pre_cluster_treshold' should be high >= 0.95"
        assert treshold < 1, "'pre_cluster_treshold' should be a float lower than 1"

    elif wildcards.genomeset=="precluster_representatives":

        treshold= config['pre_cluster_min_identity']

        assert treshold < 0.95, "'pre_cluster_treshold' should be lower than 0.95"

    else:
        raise Exception("genome set '{wildcards.genomeset}' not understood.")

    return 1-treshold


localrules: list_genomes
rule list_genomes:
    input:
        genome_folder
    output:
        temp("sketches/all_genomes.list")
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
        "sketches/{genomeset}.list"
    output:
        "sketches/{genomeset}.msh"
    params:
        k= config['sketch_k'],
        s= config['sketch_size'],
        out_name = lambda wc, output: os.path.splitext(output[0])[0]
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/sketch/mash_sketch_{genomeset}.log"
    benchmark:
        "logs/benchmark/sketch/mash_sketch_{genomeset}.txt"
    shell:
        "mash sketch -o {params.out_name} -p {threads} -s {params.s} -k {params.k} -l {input[0]} 2> {log}"


rule mash_calculate_dist_precluster:
    input:
        genomes="sketches/{genomeset}.msh"
    output:
        temp("precluster/mash_dists_{genomeset}.tsv")
    params:
        d= precluster_sketch_max_dist,
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/precluster/mash_dist_{genomeset}.log"
    benchmark:
        "logs/benchmark/sketch/mash_dists_{genomeset}.tsv"
    shell:
        "mash dist -p {threads} -d {params.d} "
        "{input.genomes} {input.genomes} > {output[0]} 2> {log}"



rule mash_calculate_dist:
    input:
        genomes="sketches/all_Genomes.msh"
    output:
        "tables/mash_dists.tsv"
    params:
        d= config.get('sketch_max_dist',0.2),
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




rule bindash_sketch:
    input:
        "sketches/{genomeset}.list"
    output:
        temp(expand("sketches/{{genomeset}}_K{{k}}.{ext}",ext=['bdsh','bdsh.dat','bdsh.txt'])),
    params:
        sketchsize64= int(config['sketch_size'])//64,
        extra=config.get('bindash_extra',"")
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/sketch/bindash_sketch_{genomeset}_K{k}.log"
    benchmark:
        "logs/benchmark/sketch/bindash_sketch_{genomeset}_K{k}.txt"
    shell:
        "bindash sketch "
        "--outfname={output[0]} "
        "--nthreads={threads} "
        "--sketchsize64={params.sketchsize64} "
        "--kmerlen={wildcards.k} "
        "{params.extra} "
        "--listfname={input[0]} 2> {log}"

rule bindash_dist_precluster:
    input:
        expand("sketches/{{genomeset}}_K{k}.{ext}",k=config[sketch_k],ext=['bdsh','bdsh.dat','bdsh.txt'])
    output:
        temp("precluster/bindash_dists_{genomeset}.tsv")
    params:
        d= precluster_sketch_max_dist
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/precluster/bindash_dists_{genomeset}.log"
    benchmark:
        "logs/benchmark/sketch/bindash_dist_{genomeset}.txt"
    shell:
        "bindash dist "
        "--nthreads={threads} "
        "--mthres={params.d} "
        "--outfname={output} {input[0]} 2> {log}"


rule bindash_dist:
    input:
        expand("sketches/{{genomeset}}_K{{k}}.{ext}",ext=['bdsh','bdsh.dat','bdsh.txt'])
    output:
        "tables/bindash_dists_K{k}.tsv"
    params:
        d= config.get('sketch_max_dist',0.2)
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/dist_K{k}.log"
    benchmark:
        "logs/benchmark/bindash_dist_K{k}.txt"
    shell:
        "bindash dist "
        "--nthreads={threads} "
        "--mthres={params.d} "
        "--outfname={output} {input[0]} 2> {log}"



rule precluster:
    input:
        dists=f"precluster/{config['sketcher']}_dists_all_genomes.tsv",
        quality="tables/Genome_quality.tsv",
    output:
        preclustering="precluster/preclustering.tsv",
        cluster_list="sketches/precluster_representatives.list"
    params:
        treshold=config['pre_cluster_treshold'],
        N=config['subset_size_alignments'],
        genome_folder=genome_folder,
    # resources:
    #     mem=config['mem']['large'],
    #     time=config['runtime']['precluster']
    benchmark:
        "logs/benchmark/precluster.txt"
    log:
        'logs/precluster/log.txt',
        stats='logs/precluster/stats.txt'
    script:
        "../scripts/precluster.py"



checkpoint get_alignment_subsets:
    input:
        f"precluster/{config['sketcher']}_dists_precluster_representatives.tsv",
        stats="logs/precluster/stats.txt"
    output:
        dir=temp(directory('alignment/subsets'))
    params:
        N=config['subset_size_alignments'],
        min_identity=config['pre_cluster_min_identity'],
    resources:
        mem=config['mem']['large'],
        time=config['runtime']['precluster']
    benchmark:
        "logs/benchmark/precluster/get_alignment_subsets.txt"
    run:

        import networkx as nx

        Msh= gd.load_mash(input[0])
        Msh= Msh.query(f"Identity>={params.min_identity}")
        G= gd.to_graph(Msh)

        G.remove_edges_from(nx.selfloop_edges(G))


        os.makedirs(output.dir)

        with open(input.stats,'a') as logfile:
            logfile.write(f"{len(G.edges())} pairwise comparisons have to be calculated for preclustering.\n"
             )


        # save pairs ordered alphabetically in multiple files
        fout=None
        for i,pair in enumerate(G.edges()):
            if (i % int(params.N)) ==0:
                n_file= int(i // params.N )+1
                if fout is not None: fout.close()
                fout= open(f"{output.dir}/subset_{n_file}.txt",'w')

            fout.write("\t".join(sorted(pair))+'\n')


        fout.close()
