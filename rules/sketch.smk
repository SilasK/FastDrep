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

    log:
        'logs/precluster/log.txt',
        stats='logs/precluster/stats.txt'
    script:
        "../scripts/precluster.py"
