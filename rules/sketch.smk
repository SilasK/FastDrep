localrules: list_genomes
rule list_genomes:
    input:
        genome_folder
    output:
        temp("sketches/genome_list_{sketcher}.txt")
    run:
        from glob import glob
        with open(output[0],'w') as f:
            f.write('\n'.join(glob(f"{input[0]}/*{config['fasta_extension']}"))+'\n' )


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

rule tsv2parquet:
    input:
        "tables/{tool}_dists.tsv"
    output:
        "tables/{tool}_dists.parquet"
    resources:
        mem_mb=config['mem']['large'] *1000
    threads:
        1
    run:

        if wildcards.tool == "mummer":
            open_function= gd.load_mummer
        elif wildcards.tool == "minimap":
            open_function= gd.load_minimap
        elif wildcards.tool == "bindash":
            open_function= gd.load_bindash
        elif wildcards.tool == "fastani":
            open_function= gd.load_fastani
        elif wildcards.tool == "mash":
            open_function= gd.load_mash
        else:
            raise Exception(
                f"Don't know how to load table from tool : {wildcards.tool}"
            )


        M= open_function(input[0]).drop(['Identity'],axis=1)
        M.to_parquet(output[0],engine="pyarrow")



checkpoint filter_sketch:
    input:
        f"tables/{config['sketcher']}_dists.tsv"
    output:
        temp(directory('{aligner}/subsets'))
    params:
        threshold=config['pre_cluster_threshold'],
        N=config['subset_size_alignments']
    benchmark:
        "logs/benchmark/filter_sketch_{aligner}.txt"
    resources:
        mem_mb=config['mem']['large']*1000
    run:

        F= gd.load_mash(input[0])
        G= gd.to_graph(F.query(f"Distance<={params.threshold}"))
        if hasattr(G,'selfloop_edges'):
            G.remove_edges_from(G.selfloop_edges())

        os.makedirs(output[0])

        fout=None
        for i,e in enumerate(G.edges()):
            if (i % int(params.N)) ==0:
                n= int(i // params.N )+1
                if fout is not None: fout.close()
                fout= open(f"{output[0]}/subset_{n}.txt",'w')

            fout.write("\t".join(sorted(e))+'\n')

def get_alignment_subsets(**wildcards):
    subset_dir= checkpoints.filter_sketch.get(**wildcards).output[0]

    subsets= glob_wildcards(os.path.join(subset_dir,'{subset}.txt')).subset

    return subsets
